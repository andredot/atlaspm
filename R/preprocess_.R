#' Preprocess a mortality file down to avoidable causes of death
#'
#' Starting from a raw mortality table, this keeps only premature deaths
#' (age under \code{age_threshold}, default 75, following the OECD/Eurostat
#' definition of premature mortality), attaches the preventable-cause metadata
#' from the OECD/Eurostat crosswalk via longest-prefix ICD-10 matching, and
#' returns only the rows whose cause of death is on the preventable list.
#'
#' ICD-10 codes are matched hierarchically. Each code is first upper-cased and
#' stripped of any non-alphanumeric characters (so the dotless \code{"I679"}
#' and a dotted \code{"I67.9"} are treated identically). The first 4 characters
#' are matched against the 4-character keys in the lookup; if there is no
#' 4-character match, the first 3 characters are matched against the
#' 3-character keys. The longer match always wins, which is what distinguishes
#' carve-outs such as \code{A40.3} (preventable) from the rest of \code{A40}
#' (not preventable). Codes that match nothing, and missing/invalid codes
#' (e.g. \code{NA}, \code{"ZZZZ"}), are dropped by the final filter.
#'
#' The eight causes that the OECD/Eurostat list splits 50/50 between the
#' preventable and treatable categories (tuberculosis, cervical cancer,
#' diabetes, and five cardiovascular causes) carry \code{weight = 0.5}; all
#' other preventable causes carry \code{weight = 1}. To count deaths without
#' double-counting against a treatable analysis, weight by this column.
#'
#' @param mort_raw A data frame of individual death records (e.g. \code{mort_raw}).
#'   If your data is in a file, read it first (e.g. with
#'   \code{readr::read_csv()}) and pass the resulting data frame.
#' @param lookup The preventable crosswalk. Either a path to
#'   \code{preventable_lookup.csv} or a data frame containing at least the
#'   columns \code{key}, \code{group}, \code{cause}, \code{weight} and
#'   \code{mechanism}.
#' @param code_col Name of the ICD-10 cause-of-death column in \code{mort_raw}.
#'   Default \code{"code"}.
#' @param age_col Name of the age-in-years column in \code{mort_raw}.
#'   Default \code{"age"}. Records with a missing age are dropped by the
#'   \code{age < age_threshold} filter.
#' @param age_threshold Upper age bound, exclusive: deaths at or above this age
#'   are removed. Default \code{75}.
#'
#' @return A data frame of preventable premature deaths only: every original
#'   column of \code{mort_raw}, plus \code{match_key} (the ICD-10 prefix that
#'   matched) and the appended \code{group}, \code{cause}, \code{weight},
#'   \code{mechanism} (and \code{flag}, if present in the lookup).
#'
#' @examples
#' \dontrun{
#' preventable <- filter_preventable(mort_raw, "preventable_lookup.csv")
#'
#' # Weighted death counts by cause (50/50 causes count as 0.5):
#' preventable |>
#'   dplyr::count(cause, wt = weight, sort = TRUE)
#' }
#'
#' @importFrom dplyr filter mutate left_join select case_when any_of all_of |>
#' @importFrom stringr str_remove_all str_to_upper str_sub str_length
#' @importFrom readr read_csv
#' @importFrom rlang .data
#' @export
preprocess_mortality <- function(mort_raw,
                                 lookup,
                                 code_col      = "code",
                                 age_col       = "age",
                                 age_threshold = 75) {

  # Allow a CSV path or a ready-made data frame for the lookup
  if (is.character(lookup)) {
    lookup <- readr::read_csv(lookup, show_col_types = FALSE)
  }
  lookup <- dplyr::select(
    lookup,
    dplyr::any_of(c("key", "group", "cause", "type", "weight", "mechanism",
                    "flag"))
  )
  validate_lookup(lookup)

  # Split the lookup by key length so the 4-char match can take priority
  key4_set <- lookup$key[stringr::str_length(lookup$key) == 4]
  key3_set <- lookup$key[stringr::str_length(lookup$key) == 3]

  mort_raw |>
    # 1. premature deaths only
    dplyr::filter(.data[[age_col]] < age_threshold) |>
    # 2. normalise the code and derive the 4- and 3-char prefixes
    dplyr::mutate(
      .code_norm = stringr::str_to_upper(
        stringr::str_remove_all(.data[[code_col]], "[^A-Za-z0-9]")
      ),
      .key4 = stringr::str_sub(.data[[".code_norm"]], 1L, 4L),
      .key3 = stringr::str_sub(.data[[".code_norm"]], 1L, 3L)
    ) |>
    # 3. longest-prefix resolution: try 4-char key, then 3-char key
    dplyr::mutate(
      match_key = dplyr::case_when(
        .data[[".key4"]] %in% key4_set ~ .data[[".key4"]],
        .data[[".key3"]] %in% key3_set ~ .data[[".key3"]],
        TRUE                           ~ NA_character_
      )
    ) |>
    # 4. append the avoidable-cause metadata
    dplyr::left_join(lookup, by = c("match_key" = "key")) |>
    # 5. keep only the preventable causes (the matched rows)
    dplyr::filter(!is.na(.data[["cause"]])) |>
    dplyr::select(-dplyr::all_of(c(".code_norm", ".key4", ".key3")))
}
#' Crude mortality rate table by area, wide over cause / group / mechanism
#'
#' Builds a wide crude-mortality-rate table with one row per area (comune by
#' default, but any spatial unit - e.g. Milan's NILs - via
#' \code{mort_col}/\code{pop_col}) and one column per category of each
#' classification, namely \code{cause} (columns prefixed \code{"C_"}),
#' \code{group} (prefixed \code{"G_"}) and \code{mechanism} (prefixed
#' \code{"M_"}). Every cell is the crude rate \code{deaths / population * per}
#' (deaths per \code{per} residents, 100 000 by default).
#'
#' The area key is kept as a plain column (\code{group_var}) rather than baked
#' into the rates, so the same area-level table can later be re-aggregated to
#' ASST, distretto, etc. by joining a crosswalk and re-deriving rates.
#'
#' Deaths are aggregated as the \strong{sum of the \code{weight} column}, so the
#' causes the OECD/Eurostat list splits 50/50 each count as 0.5 of a death. Set
#' \code{weight_col = NULL} to count every record as a whole death.
#'
#' The denominator is read from a semicolon-separated population CSV. Only the
#' year \code{pop_year} (default 2023) is kept and \code{numero} is summed over
#' all ages and both sexes to give one population figure per area.
#'
#' \strong{Area key.} The deaths-side key is \code{mort_col} and the
#' population-side key is \code{pop_col}; the two must produce identical strings
#' for the join. By default the population key is the numeric ISTAT code
#' \code{Codice comune}, zero-padded to six digits (\code{pad_area = TRUE}) to
#' match the padded codes in the deaths file. For a mixed key containing
#' non-numeric values - e.g. an \code{area} column holding both \code{"015011"}
#' (comune) and \code{"015146_79"} (Milan comune + NIL) - set \code{pad_area =
#' FALSE} so the string passes through untouched (\code{as.integer("015146_79")}
#' would otherwise become \code{NA} and silently drop that area's denominator).
#'
#' Every area present in the population file is returned; areas with no deaths
#' get 0 in every rate column.
#'
#' @param mort_count Data frame of deaths, one row per death, containing the
#'   area column \code{mort_col}, the classification columns \code{cause},
#'   \code{group} and \code{mechanism}, and (unless \code{weight_col = NULL})
#'   the weight column.
#' @param population Path to the population CSV (or a pre-read data frame) with
#'   columns \code{pop_col}, \code{Eta}, \code{anno}, \code{sesso},
#'   \code{Comune}, \code{numero}.
#' @param group_var Name to give the area key in the output. Default
#'   \code{"comune"}.
#' @param mort_col Name of the area column in \code{mort_count} that matches the
#'   population key \code{pop_col}. Default \code{"comune_residenza"}.
#' @param pop_col Area column in \code{population} matching \code{mort_col}.
#'   Default \code{"Codice comune"}. Use \code{"area"} for a NIL-aware population
#'   table keyed on a mixed comune/NIL column.
#' @param pad_area Whether to zero-pad the population area key to six digits with
#'   \code{sprintf("\%06d", ...)}. Default \code{TRUE} (correct for a purely
#'   numeric ISTAT code). Set \code{FALSE} when \code{pop_col} holds non-numeric
#'   keys such as \code{"015146_79"}, which must not be coerced through
#'   \code{as.integer}.
#' @param class_vars Named character vector mapping classification columns in
#'   \code{mort_count} to their output column prefix. Default
#'   \code{c(cause = "C", group = "G", mechanism = "M")}.
#' @param pop_year Census/population year to keep from the CSV. Default
#'   \code{2023}.
#' @param weight_col Name of the weight column in \code{mort_count}; deaths are
#'   summed over it (50/50 causes contribute 0.5). \code{NULL} counts each
#'   record as one death. Default \code{"weight"}.
#' @param per Rate denominator multiplier. Default \code{100000}.
#'
#' @return A tibble with one row per area: \code{group_var}, \code{population},
#'   one crude-rate column per category named \code{<prefix>_<label>} with the
#'   label cleaned to snake_case (e.g. \code{C_lung_cancer}, \code{G_cancer},
#'   \code{M_immunisation_and_prophylaxis}), and a \code{total} column with the
#'   crude rate over all deaths in the area. Note the per-prefix columns do
#'   not sum to \code{total}: each classification (cause / group / mechanism)
#'   already partitions the same deaths, and split causes are counted under
#'   more than one mechanism.
#'
#' @examples
#' \dontrun{
#' # comune-level (defaults)
#' rates <- preprocess_cmr(mort_count, "popolazione.csv")
#'
#' # area-level: Milan split into NILs, mixed comune/NIL key, no padding
#' area_rates <- preprocess_cmr(mort_count_area, pop_area_table,
#'                              group_var = "area",
#'                              mort_col  = "area_residenza",
#'                              pop_col   = "area",
#'                              pad_area  = FALSE)
#'
#' # re-aggregate comuni to ASST later via a crosswalk:
#' rates |>
#'   dplyr::left_join(asst_crosswalk, by = "comune") |>
#'   dplyr::group_by(asst) # then recompute rates from counts at the new level
#' }
#'
#' @importFrom dplyr mutate group_by across all_of any_of summarise left_join
#'   coalesce rename filter select |>
#' @importFrom tidyr pivot_wider
#' @importFrom purrr imap_dfr
#' @importFrom readr read_delim
#' @importFrom rlang .data :=
#' @importFrom janitor make_clean_names
#' @export
preprocess_cmr <- function(mort_count,
                           population,
                           group_var  = "comune",
                           mort_col   = "comune_residenza",
                           pop_col    = "Codice comune",
                           pad_area   = TRUE,
                           class_vars = c(cause = "C", group = "G", mechanism = "M"),
                           pop_year   = 2023,
                           weight_col = "weight",
                           per        = 100000) {

  use_weight <- !is.null(weight_col) && weight_col %in% names(mort_count)

  # ---- 1. population denominator: one row per area, pop_year, all ages/sexes ----
  if (is.character(population)) {
    population <- readr::read_delim(population, delim = ";", show_col_types = FALSE)
  }

  pop <- population |>
    dplyr::filter(.data[["anno"]] == pop_year) |>
    dplyr::mutate(
      # match the deaths-side key: pad numeric comune codes to 6 digits, but let
      # mixed keys like "015146_79" pass through as-is (as.integer() -> NA).
      .area = if (pad_area) sprintf("%06d", as.integer(.data[[pop_col]]))
      else as.character(.data[[pop_col]])
    ) |>
    dplyr::group_by(.area) |>
    dplyr::summarise(population = sum(.data[["numero"]], na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::rename(!!group_var := ".area")

  # ---- 2. weighted deaths per area x (each classification, each level) ----
  m <- mort_count |>
    dplyr::mutate(
      .w    = if (use_weight) .data[[weight_col]] else 1,
      .area = as.character(.data[[mort_col]])
    )

  # long table of (area, column-name, deaths), looping over classifications
  counts <- purrr::imap_dfr(class_vars, function(prefix, var) {
    lvls   <- unique(as.character(m[[var]]))
    clean  <- stats::setNames(
      paste0(prefix, "_", janitor::make_clean_names(lvls)),
      lvls
    )
    m |>
      dplyr::group_by(.area, .lvl = as.character(.data[[var]])) |>
      dplyr::summarise(deaths = sum(.data[[".w"]], na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(col = clean[.lvl]) |>
      dplyr::select(.area, col, deaths)
  })

  # ---- 3. widen counts, join onto full population list, divide by pop ----
  # Wide counts for areas that actually have deaths (one column per category).
  wide_counts <- counts |>
    tidyr::pivot_wider(
      id_cols     = ".area",
      names_from  = "col",
      values_from = "deaths",
      values_fill = 0
    )

  # Total deaths per area (counted once, not summed across classifications).
  totals <- m |>
    dplyr::group_by(.area) |>
    dplyr::summarise(total = sum(.data[[".w"]], na.rm = TRUE), .groups = "drop")

  wide_counts <- dplyr::left_join(wide_counts, totals, by = ".area")

  rate_cols <- setdiff(names(wide_counts), ".area")

  # Left-join onto EVERY area in the population; areas with no deaths get 0.
  pop |>
    dplyr::rename(.area = dplyr::all_of(group_var)) |>
    dplyr::left_join(wide_counts, by = ".area") |>
    dplyr::mutate(dplyr::across(dplyr::all_of(rate_cols),
                                ~ dplyr::coalesce(.x, 0) / .data[["population"]] * per)) |>
    dplyr::rename(!!group_var := ".area")
}


#' Indirectly standardised mortality table by area, wide over cause / group / mechanism
#'
#' Companion to \code{\link{preprocess_cmr}} that returns indirectly
#' age-sex-standardised mortality instead of the crude rate. For every area
#' (comune by default, but any spatial unit - e.g. Milan's NILs - via
#' \code{mort_col}/\code{pop_col}) and every category of each classification -
#' \code{cause} (columns prefixed \code{"C_"}), \code{group} (\code{"G_"}) and
#' \code{mechanism} (\code{"M_"}) - it returns four columns: observed deaths
#' (\code{<prefix>_<label>_obs}), expected deaths (\code{<prefix>_<label>_exp}),
#' the SMR (\code{<prefix>_<label>_smr}, observed / expected) and the indirectly
#' standardised rate (\code{<prefix>_<label>_isr}, SMR x the standard crude rate,
#' per \code{per} residents).
#'
#' Indirect standardisation applies a standard schedule of age-sex specific
#' rates to each area's own age-sex population to obtain \emph{expected} deaths,
#' then forms SMR = observed / expected. The standard schedule is computed
#' \strong{internally} from the pooled data: for each single year of age x sex
#' (and category), the standard rate is the total deaths across all areas
#' divided by the total population across all areas. Strata are single years of
#' age (\code{eta} / \code{Eta}) crossed with sex (\code{sesso}).
#'
#' Deaths are weighted by \code{weight_col} (50/50 split causes contribute 0.5).
#' Observed and expected deaths are aggregated on the same weighting, so the SMR
#' is internally consistent. Every area present in the population file is
#' returned; areas with no deaths get SMR/ISR/observed of 0, but their
#' \emph{expected} stays positive (it must not be zeroed - a Poisson spatial
#' model uses \code{log(expected)} as its offset, and \code{log(0)} is
#' undefined).
#'
#' \strong{Area key.} The deaths-side key is \code{mort_col} and the
#' population-side key is \code{pop_col}; the two must produce identical strings
#' for a join to occur. By default the population key is the numeric ISTAT code
#' \code{Codice comune}, which is zero-padded to six digits (\code{pad_area =
#' TRUE}) to match the padded codes carried in the deaths file. For a mixed key
#' that contains non-numeric values - e.g. an \code{area} column holding both
#' \code{"015011"} (comune) and \code{"015146_79"} (Milan comune + NIL) - set
#' \code{pad_area = FALSE} so the string passes through untouched
#' (\code{as.integer("015146_79")} would otherwise become \code{NA} and silently
#' drop that area's denominator).
#'
#' @param mort_count Data frame of deaths, one row per death, containing the
#'   area column \code{mort_col}, the strata columns \code{age_col} and
#'   \code{sex_col}, the classification columns named in \code{class_vars}, and
#'   (unless \code{weight_col = NULL}) the weight column.
#' @param population Path to the semicolon-separated population CSV (or a
#'   pre-read data frame) with columns \code{pop_col}, \code{Eta}, \code{anno},
#'   \code{sesso}, \code{Comune}, \code{numero}. \code{Eta} is the single year of
#'   age and \code{numero} the population count.
#' @param group_var Name to give the area key in the output. Default
#'   \code{"comune"}.
#' @param mort_col Area column in \code{mort_count} matching the population key
#'   \code{pop_col}. Default \code{"comune_residenza"}.
#' @param pop_col Area column in \code{population} matching \code{mort_col}.
#'   Default \code{"Codice comune"}. Use \code{"area"} for a NIL-aware population
#'   table keyed on a mixed comune/NIL column.
#' @param pad_area Whether to zero-pad the population area key to six digits with
#'   \code{sprintf("\%06d", ...)}. Default \code{TRUE} (correct for a purely
#'   numeric ISTAT code). Set \code{FALSE} when \code{pop_col} holds non-numeric
#'   keys such as \code{"015146_79"}, which must not be coerced through
#'   \code{as.integer}.
#' @param age_col,sex_col Age (single year) and sex columns in \code{mort_count}.
#'   Defaults \code{"eta"} and \code{"sesso"}. They are matched to \code{Eta}
#'   and \code{sesso} in the population file.
#' @param class_vars Named character vector mapping classification columns in
#'   \code{mort_count} to their output column prefix. Default
#'   \code{c(cause = "C", group = "G", mechanism = "M")}.
#' @param pop_year Population year to keep from the CSV. Default \code{2023}.
#' @param weight_col Weight column in \code{mort_count}. \code{NULL} counts each
#'   record as one death. Default \code{"weight"}.
#' @param per Multiplier for the standardised rate. Default \code{100000}.
#'
#' @return A tibble with one row per area: \code{group_var}, \code{population},
#'   then for every category four columns \code{<prefix>_<label>_obs},
#'   \code{<prefix>_<label>_exp}, \code{<prefix>_<label>_smr} and
#'   \code{<prefix>_<label>_isr}, plus \code{total_obs} (observed deaths),
#'   \code{total_exp} (expected deaths under the standard schedule),
#'   \code{total_smr} and \code{total_isr} over all deaths. SMR is
#'   observed/expected; ISR is the SMR scaled by the standard population's crude
#'   rate for that category. The per-category \code{_obs}/\code{_exp} columns are
#'   the response and offset a Poisson spatial model (e.g. BYM2) needs to smooth
#'   each category separately.
#'
#' @examples
#' \dontrun{
#' # comune-level (defaults)
#' smr <- preprocess_smr(mort_count, "popolazione.csv")
#'
#' # area-level: Milan split into NILs, mixed comune/NIL key, no padding
#' area_smr <- preprocess_smr(mort_count_area, pop_area_table,
#'                            group_var = "area",
#'                            mort_col  = "area_residenza",
#'                            pop_col   = "area",
#'                            pad_area  = FALSE)
#' }
#'
#' @seealso \code{\link{preprocess_cmr}} for the crude-rate version.
#' @importFrom dplyr mutate group_by ungroup across all_of summarise left_join
#'   inner_join coalesce rename filter select distinct |>
#' @importFrom tidyr pivot_wider
#' @importFrom purrr imap_dfr
#' @importFrom readr read_delim
#' @importFrom rlang .data :=
#' @importFrom janitor make_clean_names
#' @export
preprocess_smr <- function(mort_count,
                           population,
                           group_var  = "comune",
                           mort_col   = "comune_residenza",
                           pop_col    = "Codice comune",
                           pad_area   = TRUE,
                           age_col    = "eta",
                           sex_col    = "sesso",
                           class_vars = c(cause = "C", group = "G", mechanism = "M"),
                           pop_year   = 2023,
                           weight_col = "weight",
                           per        = 100000) {

  use_weight <- !is.null(weight_col) && weight_col %in% names(mort_count)

  # ---- 1. population by area x age x sex (single year) ----
  if (is.character(population)) {
    population <- readr::read_delim(population, delim = ";", show_col_types = FALSE)
  }

  pop_strata <- population |>
    dplyr::filter(.data[["anno"]] == pop_year) |>
    dplyr::mutate(
      # match the deaths-side key: pad numeric comune codes to 6 digits, but let
      # mixed keys like "015146_79" pass through as-is (as.integer() -> NA).
      .area = if (pad_area) sprintf("%06d", as.integer(.data[[pop_col]]))
      else as.character(.data[[pop_col]]),
      .age  = as.integer(.data[["Eta"]]),
      .sex  = as.integer(.data[["sesso"]]),
      .pop  = as.numeric(.data[["numero"]])
    ) |>
    dplyr::group_by(.area, .age, .sex) |>
    dplyr::summarise(.pop = sum(.pop, na.rm = TRUE), .groups = "drop")

  # total resident population per area (all ages/sexes) for the output
  pop_total <- pop_strata |>
    dplyr::group_by(.area) |>
    dplyr::summarise(population = sum(.pop), .groups = "drop")

  # ---- 2. observed deaths per area x category, and per stratum x category ----
  m <- mort_count |>
    dplyr::mutate(
      .w    = if (use_weight) .data[[weight_col]] else 1,
      .area = as.character(.data[[mort_col]]),
      .age  = as.integer(.data[[age_col]]),
      .sex  = as.integer(.data[[sex_col]])
    )

  # standard schedule denominator: total population per stratum across ALL areas
  std_denom <- pop_strata |>
    dplyr::group_by(.age, .sex) |>
    dplyr::summarise(std_pop = sum(.pop), .groups = "drop")

  total_std_pop <- sum(pop_strata$.pop)

  # one tidy block per classification, then bind
  blocks <- purrr::imap_dfr(class_vars, function(prefix, var) {
    # clean DISTINCT levels once (avoids make_clean_names() uniquifying per row)
    lvls  <- unique(as.character(m[[var]]))
    clean <- stats::setNames(paste0(prefix, "_", janitor::make_clean_names(lvls)), lvls)

    mc <- m |>
      dplyr::mutate(.lvl = as.character(.data[[var]]))

    # standard age-sex-specific rate per (level, age, sex): pooled deaths / pooled pop
    std_rate <- mc |>
      dplyr::group_by(.lvl, .age, .sex) |>
      dplyr::summarise(std_deaths = sum(.data[[".w"]], na.rm = TRUE), .groups = "drop") |>
      dplyr::left_join(std_denom, by = c(".age", ".sex")) |>
      dplyr::mutate(std_rate = std_deaths / std_pop)

    # standard crude rate per level (for converting SMR -> ISR):
    #   total pooled deaths / total pooled population
    std_crude <- mc |>
      dplyr::group_by(.lvl) |>
      dplyr::summarise(std_crude = sum(.data[[".w"]], na.rm = TRUE) / total_std_pop,
                       .groups = "drop")

    # OBSERVED deaths per area x level
    observed <- mc |>
      dplyr::group_by(.area, .lvl) |>
      dplyr::summarise(observed = sum(.data[[".w"]], na.rm = TRUE), .groups = "drop")

    # EXPECTED deaths per area x level = sum over strata of
    #   (area pop in stratum) x (standard rate for that level/stratum)
    expected <- pop_strata |>
      dplyr::inner_join(std_rate, by = c(".age", ".sex"),
                        relationship = "many-to-many") |>
      dplyr::mutate(exp = .pop * std_rate) |>
      dplyr::group_by(.area, .lvl) |>
      dplyr::summarise(expected = sum(exp, na.rm = TRUE), .groups = "drop")

    # combine: SMR = observed / expected ; ISR = SMR x standard crude rate x per
    expected |>
      dplyr::left_join(observed, by = c(".area", ".lvl")) |>
      dplyr::left_join(std_crude, by = ".lvl") |>
      dplyr::mutate(
        observed = dplyr::coalesce(observed, 0),
        smr      = dplyr::if_else(expected > 0, observed / expected, NA_real_),
        isr      = smr * std_crude * per,
        col      = clean[.lvl]
      ) |>
      dplyr::select(.area, col, observed, expected, smr, isr)
  })

  # ---- 3. add an all-cause "total" block (deaths counted once) ----
  std_rate_tot <- m |>
    dplyr::group_by(.age, .sex) |>
    dplyr::summarise(std_deaths = sum(.data[[".w"]], na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(std_denom, by = c(".age", ".sex")) |>
    dplyr::mutate(std_rate = std_deaths / std_pop)
  std_crude_tot <- sum(m$.w, na.rm = TRUE) / total_std_pop

  observed_tot <- m |>
    dplyr::group_by(.area) |>
    dplyr::summarise(observed = sum(.data[[".w"]], na.rm = TRUE), .groups = "drop")
  expected_tot <- pop_strata |>
    dplyr::left_join(std_rate_tot, by = c(".age", ".sex")) |>
    dplyr::mutate(exp = .pop * std_rate) |>
    dplyr::group_by(.area) |>
    dplyr::summarise(expected = sum(exp, na.rm = TRUE), .groups = "drop")
  total_block <- expected_tot |>
    dplyr::left_join(observed_tot, by = ".area") |>
    dplyr::mutate(
      observed  = dplyr::coalesce(observed, 0),
      total_smr = dplyr::if_else(expected > 0, observed / expected, NA_real_),
      total_isr = total_smr * std_crude_tot * per
    ) |>
    # keep observed + expected counts: the Poisson BYM2 needs them as response
    # and offset (log(E)); they also let downstream code spot unstable areas.
    dplyr::select(.area, total_obs = "observed", total_exp = "expected",
                  total_smr, total_isr)

  # ---- 4. widen the per-category blocks (obs, exp, smr, isr), join totals + pop ----
  wide_obs <- blocks |>
    tidyr::pivot_wider(id_cols = ".area", names_from = "col",
                       values_from = "observed", names_glue = "{col}_obs",
                       values_fill = 0)          # no death in a cell -> 0 observed: correct
  wide_exp <- blocks |>
    tidyr::pivot_wider(id_cols = ".area", names_from = "col",
                       values_from = "expected", names_glue = "{col}_exp",
                       values_fill = NA_real_)   # missing expected is NOT 0 -> NA, never log(0)
  wide_smr <- blocks |>
    tidyr::pivot_wider(id_cols = ".area", names_from = "col",
                       values_from = "smr", names_glue = "{col}_smr",
                       values_fill = 0)
  wide_isr <- blocks |>
    tidyr::pivot_wider(id_cols = ".area", names_from = "col",
                       values_from = "isr", names_glue = "{col}_isr",
                       values_fill = 0)

  # expected columns must be shielded from the blanket coalesce-to-0 below
  exp_cols <- grep("_exp$", names(wide_exp), value = TRUE)

  pop_total |>
    dplyr::left_join(wide_smr,    by = ".area") |>
    dplyr::left_join(wide_isr,    by = ".area") |>
    dplyr::left_join(wide_obs,    by = ".area") |>
    dplyr::left_join(wide_exp,    by = ".area") |>
    dplyr::left_join(total_block, by = ".area") |>
    # areas with no deaths: SMR/ISR/observed are 0, but EXPECTED must stay as-is
    # (a death-free area still has positive expected; zeroing it breaks log(E)).
    dplyr::mutate(dplyr::across(
      -dplyr::all_of(c(".area", "population", exp_cols)),
      ~ dplyr::coalesce(.x, 0)
    )) |>
    dplyr::rename(!!group_var := ".area")
}

build_area_shp <- function(nil_shp, pop_shp, milan_code = "015146") {

  # 0. put NILs in pop_shp's CRS (pop_shp is the master)
  nil <- nil_shp |>
    sf::st_transform(sf::st_crs(pop_shp)) |>
    sf::st_make_valid()

  milan_geom <- pop_shp |>
    dplyr::filter(PRO_COM_T == milan_code) |>
    sf::st_make_valid() |>
    sf::st_geometry() |>
    sf::st_union()

  # 1. cookie-cut NILs to Milan's authoritative outline (removes overhang)
  nil_cut <- nil |>
    dplyr::select(ID_NIL) |>
    sf::st_intersection(milan_geom) |>
    sf::st_collection_extract("POLYGON") |>   # drop stray lines/points
    sf::st_make_valid()

  # 2. slivers of Milan left uncovered by any NIL (the boundary shortfall)
  gaps <- sf::st_difference(milan_geom, sf::st_union(nil_cut))

  if (length(gaps) > 0 && !all(sf::st_is_empty(gaps))) {
    gaps_sfc <- gaps |> sf::st_cast("MULTIPOLYGON") |> sf::st_cast("POLYGON")
    gaps_sf  <- sf::st_sf(ID_NIL = NA_integer_, geometry = gaps_sfc)
    gaps_sf$ID_NIL <- nil_cut$ID_NIL[sf::st_nearest_feature(gaps_sf, nil_cut)]
    nil_cut <- rbind(nil_cut, gaps_sf)
  }

  # 3. dissolve back to one polygon per NIL -> now tiles Milan exactly
  milan_nils <- nil_cut |>
    dplyr::group_by(ID_NIL) |>
    dplyr::summarise(.groups = "drop") |>
    dplyr::transmute(area = paste0(milan_code, "_", ID_NIL))

  # 4. everything else keeps its pop_shp geometry, keyed the same way
  others <- pop_shp |>
    dplyr::filter(PRO_COM_T != milan_code) |>
    dplyr::transmute(area = PRO_COM_T)

  # 5. single layer, `area` == mort_count$area_residenza
  rbind(others, milan_nils) |> sf::st_make_valid()
}

#' Build a binary spatial adjacency matrix from comune geometries
#'
#' Thin wrapper around \code{geostan::shape2mat()} that returns the binary
#' contiguity matrix \code{C} used as the spatial weights in a BYM2/ICAR model.
#' The row order of \code{C} matches the row order of \code{geo}, so the same
#' object must be passed to \code{\link{fit_bym2}} as the model data.
#'
#' @param geo An \code{sf} of comune polygons, e.g.
#'   \code{add_geo(preprocess_smr(...), comuni)}.
#' @return A sparse binary adjacency matrix (\code{"B"} style), one row/column
#'   per comune, in the row order of \code{geo}.
#' @examples
#' \dontrun{
#' geo <- add_geo(mort_smr, pop_shp, data_key = "comune")
#' C   <- build_adjacency(geo)
#' }
#' @importFrom geostan shape2mat
#' @export
build_adjacency <- function(geo) {
  geostan::shape2mat(geo, style = "B")
}



#' Population-weighted quintiles
#'
#' Assigns each element to one of \code{n} groups so that the groups hold
#' (approximately) equal shares of total weight \code{w} when ordered by
#' \code{x}. This reproduces the "quintili di popolazione" used by the Italian
#' Deprivation Index: cut points fall at equal shares of \emph{population}, not
#' at equal numbers of areas.
#'
#' @param x Numeric vector to rank.
#' @param w Non-negative weights (e.g. municipal population).
#' @param n Number of groups. Default 5.
#' @return Integer vector of group indices (1 = lowest \code{x}).
#' @export
wtd_quantile_group <- function(x, w, n = 5) {
  o    <- order(x)
  cumw <- cumsum(w[o]) / sum(w[o])
  g    <- cut(cumw, breaks = c(-Inf, seq_len(n - 1) / n, Inf), labels = FALSE)
  out  <- integer(length(x))
  out[o] <- g
  out
}


#' Build the four-indicator deprivation proxy by municipality
#'
#' Aggregates census-section counts to municipality (\code{PROCOM}), forms the
#' four proxy indicators, standardises each against the \strong{national}
#' distribution, sums them into a single deprivation score, derives national
#' population quintiles, and only then filters to the municipalities present in
#' the mortality data.
#'
#' This is the permanent-census \emph{proxy} reformulation of the Italian
#' Deprivation Index, \strong{not} the validated five-indicator 2011 index. Two
#' original indicators (non-home-ownership and single-parent family) are absent
#' from the permanent-census release and are replaced by foreign-population
#' share; overcrowding is approximated by occupants per dwelling and
#' unemployment by non-employment. Treat it as a contemporary sensitivity
#' comparator, not as a newer edition of the same construct.
#'
#' Standardisation is national by design: z-scores use the mean and SD across
#' \emph{every} Italian municipality in \code{census}, so the study-area filter
#' is applied last and never influences the score or the quintile cut points.
#' All four indicators are oriented so that a higher value means more
#' disadvantage, so a plain sum of z-scores is the index.
#'
#' @param census Section-level table from \code{\link{import_census_2023}}.
#' @param mort_raw Mortality table whose \code{mort_col} lists the study
#'   municipalities to retain.
#' @param mort_col Municipality-code column in \code{mort_raw}. Default
#'   \code{"comune_residenza"}. Matched after zero-padding both sides.
#'
#' @return A tibble, one row per retained municipality: \code{comune} (6-digit
#'   key), \code{population}, the four raw indicators, the continuous
#'   \code{di_score}, and national \code{di_quintile}.
#' @seealso \code{\link{import_census_2023}}, \code{\link{wtd_quantile_group}}.
#' @export
build_deprivation_proxy <- function(census, mort_raw,
                                    mort_col = "comune_residenza") {

  # 1. aggregate section counts to municipality ------------------------------
  agg <- census |>
    dplyr::rename(comune = PROCOM) |>
    dplyr::group_by(comune) |>
    dplyr::summarise(dplyr::across(dplyr::everything(),
                                   ~ sum(.x, na.rm = TRUE)),
                     .groups = "drop")

  # 2. four proxy indicators (higher = more disadvantage) --------------------
  pop_1564 <- rowSums(
    dplyr::select(agg, dplyr::all_of(paste0("P", 17:26))),
    na.rm = TRUE
  )

  ind <- agg |>
    dplyr::mutate(
      edu_low = (P86 + P87 + P88) / P83,   # low education, pop 9+
      nonemp  = 1 - P101 / pop_1564,        # non-employment, 15-64
      foreign = ST1 / P1,                   # foreign-resident share
      crowd   = P1 / A2                     # occupants per dwelling
    )

  # 3. national standardisation (across all municipalities) ------------------
  z <- function(v) (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)

  ind <- ind |>
    dplyr::mutate(
      di_score = z(edu_low) + z(nonemp) + z(foreign) + z(crowd),
      # 4. national population quintiles (before any filtering)
      di_quintile = wtd_quantile_group(di_score, P1)
    )

  # 5. filter to the mortality study area ------------------------------------
  keep <- pad(as.character(unique(mort_raw[[mort_col]])))

  matched <- sum(ind$comune %in% keep)
  message(matched, " of ", length(keep),
          " study municipalities matched in the census data.")

  ind |>
    dplyr::filter(.data[["comune"]] %in% keep) |>
    dplyr::select("comune", population = "P1",
                  "edu_low", "nonemp", "foreign", "crowd",
                  "di_score", "di_quintile")
}

build_pop_area_table <- function(pop_finale, pop_nil, mort_count,
                                 milan_code = "015146") {

  years <- sort(unique(mort_count$anno))          # analysis years, e.g. 2022:2024

  # --- comuni: keep pop_finale as-is, drop Milan's aggregate, key on `area` ---
  finale_area <- pop_finale |>
    dplyr::mutate(area = sprintf("%06d", as.integer(`Codice comune`)),
                  Eta  = sprintf("%03d", as.integer(Eta)),
                  anno = as.character(anno),
                  sesso = as.integer(sesso)) |>
    dplyr::filter(area != milan_code) |>          # Milan is replaced by its NILs
    dplyr::select(area, Eta, anno, sesso, Comune, numero)

  # --- NILs: reshape to the same schema, replicate across years ---
  nil_area <- pop_nil |>
    rlang::set_names(c("nil_label", "genere", "eta_label", "n")) |>
    dplyr::mutate(
      nil_num = as.integer(stringr::str_extract(nil_label, "^[0-9]+")),
      area    = paste0(milan_code, "_", nil_num),
      sesso   = dplyr::recode(genere, "Maschi" = 1L, "Femmine" = 2L),
      Eta     = sprintf("%03d", as.integer(stringr::str_extract(eta_label, "[0-9]+")))
    ) |>
    dplyr::group_by(area, Eta, sesso) |>
    dplyr::summarise(numero = sum(n), .groups = "drop") |>
    tidyr::crossing(anno = as.character(years)) |>   # snapshot -> every analysis year
    dplyr::mutate(Comune = "Milano") |>
    dplyr::select(area, Eta, anno, sesso, Comune, numero)

  dplyr::bind_rows(finale_area, nil_area)
}


#' Validate the avoidable-cause lookup before it is joined
#'
#' The 50/50 causes are represented as \strong{two rows sharing one ICD key},
#' one \code{Preventable} and one \code{Treatable}, each carrying
#' \code{weight = 0.5} and a different \code{mechanism}. The fan-out that
#' produces the duplicated death records in \code{mort_count} is therefore a
#' property of the lookup, not an explicit duplication step in the code. If a
#' split cause is ever collapsed to a single row, the second mechanism
#' disappears silently and every weighted count is wrong by half that cause.
#'
#' This function fails loudly on the three ways that can break:
#' \enumerate{
#'   \item a key whose weights do not sum to exactly 1;
#'   \item a duplicated key whose rows share a \code{mechanism} (the join would
#'     fan out into indistinguishable rows);
#'   \item a weight that is neither 1 nor 0.5.
#' }
#'
#' @param lookup The lookup table, after column selection.
#' @return \code{lookup}, invisibly, if every check passes.
#' @importFrom dplyr group_by summarise filter n_distinct n |>
#' @export
validate_lookup <- function(lookup) {

  stopifnot(all(c("key", "weight", "mechanism") %in% names(lookup)))

  bad_w <- setdiff(unique(lookup$weight), c(1, 0.5))
  if (length(bad_w)) {
    stop("Unexpected weight(s) in lookup: ",
         paste(bad_w, collapse = ", "), call. = FALSE)
  }

  by_key <- lookup |>
    dplyr::group_by(.data[["key"]]) |>
    dplyr::summarise(
      total_w  = sum(.data[["weight"]]),
      n_rows   = dplyr::n(),
      n_mech   = dplyr::n_distinct(.data[["mechanism"]]),
      .groups  = "drop"
    )

  unbalanced <- dplyr::filter(by_key, abs(.data[["total_w"]] - 1) > 1e-9)
  if (nrow(unbalanced)) {
    stop("Lookup keys whose weights do not sum to 1: ",
         paste(unbalanced$key, collapse = ", "), call. = FALSE)
  }

  collided <- dplyr::filter(by_key, .data[["n_rows"]] > .data[["n_mech"]])
  if (nrow(collided)) {
    stop("Lookup keys duplicated within the same mechanism: ",
         paste(collided$key, collapse = ", "), call. = FALSE)
  }

  invisible(lookup)
}

#' One row per decedent, with preventability and mechanism as columns
#'
#' \code{mort_count} has one row per \emph{death x avoidability arm}: the 50/50
#' causes occupy two rows apiece, one preventable and one treatable, each with
#' \code{weight = 0.5}. Its row count therefore overstates the number of people
#' who died, and de-duplicating it with \code{distinct()} silently discards one
#' of the two arms.
#'
#' This function pivots instead of de-duplicating. Every death appears exactly
#' once, and the information that was spread across rows is moved into columns:
#' whether the death was preventable, whether it was treatable, and the
#' mechanism assigned within each arm. Nothing is lost.
#'
#' Use the result for anything describing \strong{decedents}: Table 1, the age
#' and sex distribution, the total N in the abstract. For anything describing
#' \strong{deaths attributable to a layer}, keep using \code{mort_count} and sum
#' \code{weight}, since a 50/50 death contributes half a death to each of two
#' layers and that is the quantity a layer analysis needs.
#'
#' @section Columns added:
#' \describe{
#'   \item{\code{preventable}}{Factor, \code{Yes}/\code{No}.}
#'   \item{\code{treatable}}{Factor, \code{Yes}/\code{No}.}
#'   \item{\code{avoidability}}{Factor with three levels: \code{Preventable
#'     only}, \code{Treatable only}, \code{Preventable and treatable}. The last
#'     is exactly the set of 50/50 causes.}
#'   \item{\code{mechanism_preventable}}{Factor. The preventive service function
#'     the death is attributed to, or \code{NA} (or \code{na_label}) if the
#'     death is not preventable.}
#'   \item{\code{mechanism_treatable}}{As above, for the treatable arm.}
#'   \item{\code{mechanism_any}}{Factor. A single mechanism column for
#'     convenience: the preventable mechanism where one exists, otherwise the
#'     treatable one. Use this when you want one mechanism per death and are
#'     content to let the preventable arm take precedence, which is the
#'     convention the OECD/Eurostat list itself uses.}
#'   \item{\code{weight_preventable}, \code{weight_treatable}}{The weights, kept
#'     so that layer totals can be recovered from this table without returning
#'     to \code{mort_count}.}
#' }
#' Column labels are attached with \code{attr(x, "label")}, so
#' \code{gtsummary::tbl_summary()} picks them up without a \code{label}
#' argument.
#'
#' @param mort_count Output of \code{\link{preprocess_mortality}}, with a
#'   \code{death_id} assigned before the 50/50 fan-out.
#' @param type_col Name of the column distinguishing the preventable from the
#'   treatable arm. Default \code{"type"}. Values are matched case-insensitively
#'   against \code{"preventable"} and \code{"treatable"}.
#' @param na_label Optional string used in place of \code{NA} in the two
#'   mechanism columns, e.g. \code{"Not applicable"}. Supplying it stops
#'   \code{gtsummary} reporting those cells as \code{Unknown}, which is
#'   misleading: the mechanism is not missing, it does not exist for that death.
#'   Default \code{NULL} keeps \code{NA}.
#'
#' @return A tibble with one row per \code{death_id}.
#'
#' @examples
#' \dontrun{
#' deaths <- build_deaths(mort_count, na_label = "Not applicable")
#'
#' # Table 1 over decedents. Labels come from the attributes.
#' deaths |>
#'   dplyr::select(sesso, eta, group, cause, avoidability,
#'                 mechanism_preventable) |>
#'   dplyr::mutate(cause = forcats::fct_lump_n(cause, 15)) |>
#'   gtsummary::tbl_summary(
#'     by = avoidability,
#'     missing = "no"
#'   ) |>
#'   gtsummary::add_overall() |>
#'   gtsummary::bold_labels()
#'
#' # Mechanism distribution among preventable deaths only
#' deaths |>
#'   dplyr::filter(preventable == "Yes") |>
#'   dplyr::select(mechanism_preventable, sesso, eta) |>
#'   gtsummary::tbl_summary(by = mechanism_preventable)
#'
#' # Layer totals still come from the weights, and agree with mort_count
#' deaths |>
#'   dplyr::count(mechanism_preventable, wt = weight_preventable)
#' }
#'
#' @importFrom dplyr arrange distinct select any_of all_of mutate across if_else
#' @importFrom dplyr left_join n_distinct group_by summarise filter pull |>
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom rlang .data
#' @export
build_deaths <- function(mort_count,
                         type_col = "type",
                         na_label = NULL) {

  # ---- validation ---------------------------------------------------------

  if (!"death_id" %in% names(mort_count)) {
    stop("`mort_count` has no `death_id`. Assign it in preprocess_mortality() ",
         "BEFORE the 50/50 fan-out, otherwise the two arms of a split cause ",
         "cannot be recognised as the same death.", call. = FALSE)
  }
  if (!type_col %in% names(mort_count)) {
    stop("`mort_count` has no `", type_col, "` column, so the preventable and ",
         "treatable arms cannot be told apart. Check that the lookup carries ",
         "a `type` column and that preprocess_mortality() retains it.",
         call. = FALSE)
  }
  if (!"mechanism" %in% names(mort_count)) {
    stop("`mort_count` has no `mechanism` column.", call. = FALSE)
  }
  if (!"weight" %in% names(mort_count)) {
    stop("`mort_count` has no `weight` column.", call. = FALSE)
  }

  n_deaths <- dplyr::n_distinct(mort_count$death_id)
  w_total  <- sum(mort_count$weight)
  if (abs(w_total - n_deaths) > 1e-6) {
    stop("sum(weight) = ", w_total, " but there are ", n_deaths,
         " distinct death_id. A split cause has lost one of its two rows, ",
         "or a death_id was assigned after the lookup join.", call. = FALSE)
  }

  arm <- tolower(trimws(as.character(mort_count[[type_col]])))
  bad <- setdiff(unique(arm), c("preventable", "treatable"))
  if (length(bad)) {
    stop("Unexpected value(s) in `", type_col, "`: ",
         paste(shQuote(bad), collapse = ", "),
         ". Expected only 'preventable' and 'treatable'.", call. = FALSE)
  }

  dup <- mort_count |>
    dplyr::mutate(.arm = arm) |>
    dplyr::count(.data[["death_id"]], .data[[".arm"]]) |>
    dplyr::filter(.data[["n"]] > 1L)
  if (nrow(dup)) {
    stop(nrow(dup), " death_id x arm combination(s) appear more than once, ",
         "e.g. death_id ", dup$death_id[1], ". Each death should have at most ",
         "one preventable row and one treatable row.", call. = FALSE)
  }

  # Columns that vary within a death (the fanned-out ones) versus those that
  # identify the decedent. Anything not in `arm_cols` must be constant within
  # death_id; if it is not, distinct() below would pick arbitrarily.
  arm_cols  <- c("mechanism", type_col, "weight", "flag")
  base_cols <- setdiff(names(mort_count), arm_cols)

  inconstant <- mort_count |>
    dplyr::group_by(.data[["death_id"]]) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(setdiff(base_cols, "death_id")),
                    ~ dplyr::n_distinct(.x, na.rm = FALSE) > 1L),
      .groups = "drop"
    ) |>
    dplyr::select(-dplyr::all_of("death_id"))
  offenders <- names(inconstant)[vapply(inconstant, any, logical(1))]
  if (length(offenders)) {
    stop("These columns vary within a single death_id and would be resolved ",
         "arbitrarily: ", paste(offenders, collapse = ", "),
         ". Add them to the fan-out handling or drop them before calling ",
         "build_deaths().", call. = FALSE)
  }

  # ---- one row per decedent ----------------------------------------------

  base <- mort_count |>
    dplyr::arrange(.data[["death_id"]]) |>
    dplyr::distinct(.data[["death_id"]], .keep_all = TRUE) |>
    dplyr::select(dplyr::all_of(base_cols))

  # ---- pivot the arm-varying columns to wide ------------------------------

  wide <- mort_count |>
    dplyr::mutate(.arm = arm) |>
    dplyr::select(dplyr::all_of(c("death_id", ".arm", "mechanism", "weight"))) |>
    tidyr::pivot_wider(
      id_cols     = dplyr::all_of("death_id"),
      names_from  = ".arm",
      values_from = dplyr::all_of(c("mechanism", "weight")),
      names_glue  = "{.value}_{.arm}"
    )

  # pivot_wider omits a column entirely if no death has that arm
  for (nm in c("mechanism_preventable", "mechanism_treatable")) {
    if (!nm %in% names(wide)) wide[[nm]] <- NA_character_
  }
  for (nm in c("weight_preventable", "weight_treatable")) {
    if (!nm %in% names(wide)) wide[[nm]] <- NA_real_
  }

  mech_levels <- sort(unique(stats::na.omit(as.character(mort_count$mechanism))))

  out <- base |>
    dplyr::left_join(wide, by = "death_id") |>
    dplyr::mutate(
      preventable = factor(
        dplyr::if_else(is.na(.data[["mechanism_preventable"]]), "No", "Yes"),
        levels = c("No", "Yes")
      ),
      treatable = factor(
        dplyr::if_else(is.na(.data[["mechanism_treatable"]]), "No", "Yes"),
        levels = c("No", "Yes")
      ),
      avoidability = factor(
        dplyr::case_when(
          !is.na(.data[["mechanism_preventable"]]) &
            !is.na(.data[["mechanism_treatable"]]) ~ "Preventable and treatable",
          !is.na(.data[["mechanism_preventable"]]) ~ "Preventable only",
          !is.na(.data[["mechanism_treatable"]])   ~ "Treatable only"
        ),
        levels = c("Preventable only", "Treatable only",
                   "Preventable and treatable")
      ),
      mechanism_any = factor(
        dplyr::coalesce(.data[["mechanism_preventable"]],
                        .data[["mechanism_treatable"]]),
        levels = mech_levels
      ),
      mechanism_preventable = factor(.data[["mechanism_preventable"]],
                                     levels = mech_levels),
      mechanism_treatable   = factor(.data[["mechanism_treatable"]],
                                     levels = mech_levels),
      weight_preventable = tidyr::replace_na(.data[["weight_preventable"]], 0),
      weight_treatable   = tidyr::replace_na(.data[["weight_treatable"]], 0)
    )

  # Explicit level instead of NA, so gtsummary does not report "Unknown" for
  # cells where the mechanism does not exist rather than being missing.
  if (!is.null(na_label)) {
    out <- out |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(c("mechanism_preventable", "mechanism_treatable")),
          ~ forcats::fct_na_value_to_level(.x, level = na_label)
        )
      )
  }

  # ---- post-conditions ----------------------------------------------------

  stopifnot(nrow(out) == n_deaths)
  stopifnot(!any(is.na(out$avoidability)))
  stopifnot(
    abs(sum(out$weight_preventable) + sum(out$weight_treatable) - n_deaths) < 1e-6
  )

  # ---- labels for gtsummary ----------------------------------------------

  labs <- c(
    preventable           = "Preventable",
    treatable             = "Treatable",
    avoidability          = "Avoidability category",
    mechanism_preventable = "Preventive service function",
    mechanism_treatable   = "Treatment service function",
    mechanism_any         = "Service function"
  )
  for (nm in names(labs)) {
    if (nm %in% names(out)) attr(out[[nm]], "label") <- unname(labs[nm])
  }

  out
}
