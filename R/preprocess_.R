#' Preprocess a mortality file down to preventable causes of death
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
    dplyr::any_of(c("key", "group", "cause", "weight", "mechanism"))
  )

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
    # 4. append the preventable-cause metadata
    dplyr::left_join(lookup, by = c("match_key" = "key")) |>
    # 5. keep only the preventable causes (the matched rows)
    dplyr::filter(!is.na(.data[["cause"]])) |>
    dplyr::select(-dplyr::all_of(c(".code_norm", ".key4", ".key3")))
}

#' Crude mortality rate table by comune, wide over cause / group / mechanism
#'
#' Builds a wide crude-mortality-rate table with one row per area (comune by
#' default) and one column per category of each classification, namely
#' \code{cause} (columns prefixed \code{"C_"}), \code{group} (prefixed
#' \code{"G_"}) and \code{mechanism} (prefixed \code{"M_"}). Every cell is the
#' crude rate \code{deaths / population * per} (deaths per \code{per} residents,
#' 100 000 by default).
#'
#' The area key is kept as a plain column (\code{group_var}) rather than baked
#' into the rates, so the same comune-level table can later be re-aggregated to
#' ASST, distretto, etc. by joining a crosswalk and re-deriving rates.
#'
#' Deaths are aggregated as the \strong{sum of the \code{weight} column}, so the
#' causes the OECD/Eurostat list splits 50/50 each count as 0.5 of a death. Set
#' \code{weight_col = NULL} to count every record as a whole death.
#'
#' The denominator is read from a semicolon-separated population CSV with the
#' columns \code{Codice comune;Eta;anno;sesso;Comune;numero}. Only the year
#' \code{pop_year} (default 2023) is kept and \code{numero} is summed over all
#' ages and both sexes to give one population figure per comune. The population
#' \code{Codice comune} is matched to \code{mort_col} in \code{mort_count}
#' (default \code{"comune_residenza_2023"}).
#'
#' Every comune present in the population file is returned; comuni with no
#' deaths get 0 in every rate column.
#'
#' @param mort_count Data frame of deaths, one row per death, containing the
#'   area column \code{mort_col}, the classification columns \code{cause},
#'   \code{group} and \code{mechanism}, and (unless \code{weight_col = NULL})
#'   the weight column.
#' @param population Path to the population CSV (or a pre-read data frame) with
#'   columns \code{Codice comune}, \code{Eta}, \code{anno}, \code{sesso},
#'   \code{Comune}, \code{numero}.
#' @param group_var Name to give the area key in the output. Default
#'   \code{"comune"}.
#' @param mort_col Name of the comune column in \code{mort_count} that matches
#'   the population \code{Codice comune}. Default \code{"comune_residenza_2023"}.
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
#' @return A tibble with one row per comune: \code{group_var}, \code{population},
#'   one crude-rate column per category named \code{<prefix>_<label>} with the
#'   label cleaned to snake_case (e.g. \code{C_lung_cancer}, \code{G_cancer},
#'   \code{M_immunisation_and_prophylaxis}), and a \code{total} column with the
#'   crude rate over all deaths in the comune. Note the per-prefix columns do
#'   not sum to \code{total}: each classification (cause / group / mechanism)
#'   already partitions the same deaths, and split causes are counted under
#'   more than one mechanism.
#'
#' @examples
#' \dontrun{
#' rates <- preprocess_cmr(mort_count, "popolazione.csv")
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
                           mort_col   = "comune_residenza_2023",
                           class_vars = c(cause = "C", group = "G", mechanism = "M"),
                           pop_year   = 2023,
                           weight_col = "weight",
                           per        = 100000) {

  use_weight <- !is.null(weight_col) && weight_col %in% names(mort_count)

  # ---- 1. population denominator: one row per comune, 2023, all ages/sexes ----
  if (is.character(population)) {
    population <- readr::read_delim(population, delim = ";", show_col_types = FALSE)
  }

  pop <- population |>
    dplyr::filter(.data[["anno"]] == pop_year) |>
    dplyr::mutate(`Codice comune` = sprintf("%06d", as.integer(`Codice comune`))) |>
    dplyr::group_by(.data[["Codice comune"]]) |>
    dplyr::summarise(population = sum(.data[["numero"]], na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::rename(!!group_var := "Codice comune")

  # ---- 2. weighted deaths per comune x (each classification, each level) ----
  m <- mort_count |>
    dplyr::mutate(
      .w    = if (use_weight) .data[[weight_col]] else 1,
      .area = as.character(.data[[mort_col]])
    )

  # long table of (comune, column-name, deaths), looping over classifications
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
  # Wide counts for comuni that actually have deaths (one column per category).
  wide_counts <- counts |>
    tidyr::pivot_wider(
      id_cols     = ".area",
      names_from  = "col",
      values_from = "deaths",
      values_fill = 0
    )

  # Total deaths per comune (counted once, not summed across classifications).
  totals <- m |>
    dplyr::group_by(.area) |>
    dplyr::summarise(total = sum(.data[[".w"]], na.rm = TRUE), .groups = "drop")

  wide_counts <- dplyr::left_join(wide_counts, totals, by = ".area")

  rate_cols <- setdiff(names(wide_counts), ".area")

  # Left-join onto EVERY comune in the population; comuni with no deaths get 0.
  pop |>
    dplyr::rename(.area = dplyr::all_of(group_var)) |>
    dplyr::left_join(wide_counts, by = ".area") |>
    dplyr::mutate(dplyr::across(dplyr::all_of(rate_cols),
                                ~ dplyr::coalesce(.x, 0) / .data[["population"]] * per)) |>
    dplyr::rename(!!group_var := ".area")
}
#' Indirectly standardised mortality table by comune, wide over cause / group / mechanism
#'
#' Companion to \code{\link{preprocess_cmr}} that returns indirectly
#' age-sex-standardised mortality instead of the crude rate. For every area
#' (comune by default) and every category of each classification - \code{cause}
#' (columns prefixed \code{"C_"}), \code{group} (\code{"G_"}) and
#' \code{mechanism} (\code{"M_"}) - it returns two columns: the SMR
#' (\code{<prefix>_<label>_smr}, observed / expected) and the indirectly
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
#' is internally consistent. Every comune present in the population file is
#' returned; comuni with no deaths get SMR/ISR of 0.
#'
#' @param mort_count Data frame of deaths, one row per death, containing the
#'   area column \code{mort_col}, the strata columns \code{age_col} and
#'   \code{sex_col}, the classification columns named in \code{class_vars}, and
#'   (unless \code{weight_col = NULL}) the weight column.
#' @param population Path to the semicolon-separated population CSV (or a
#'   pre-read data frame) with columns \code{Codice comune}, \code{Eta},
#'   \code{anno}, \code{sesso}, \code{Comune}, \code{numero}. \code{Eta} is the
#'   single year of age and \code{numero} the population count.
#' @param group_var Name to give the area key in the output. Default
#'   \code{"comune"}.
#' @param mort_col Comune column in \code{mort_count} matching the population
#'   \code{Codice comune}. Default \code{"comune_residenza_2023"}.
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
#' @return A tibble with one row per comune: \code{group_var}, \code{population},
#'   then for every category two columns \code{<prefix>_<label>_smr} and
#'   \code{<prefix>_<label>_isr}, plus \code{total_obs} (observed deaths),
#'   \code{total_exp} (expected deaths under the standard schedule),
#'   \code{total_smr} and \code{total_isr} over all deaths. SMR is
#'   observed/expected; ISR is the SMR scaled by the standard population's crude
#'   rate for that category. \code{total_obs}/\code{total_exp} are the response
#'   and offset a Poisson spatial model (e.g. BYM2) needs.
#'
#' @examples
#' \dontrun{
#' smr <- preprocess_smr(mort_count, "popolazione.csv")
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
                           mort_col   = "comune_residenza_2023",
                           age_col    = "eta",
                           sex_col    = "sesso",
                           class_vars = c(cause = "C", group = "G", mechanism = "M"),
                           pop_year   = 2023,
                           weight_col = "weight",
                           per        = 100000) {

  use_weight <- !is.null(weight_col) && weight_col %in% names(mort_count)

  # ---- 1. population by comune x age x sex (single year), padded code ----
  if (is.character(population)) {
    population <- readr::read_delim(population, delim = ";", show_col_types = FALSE)
  }

  pop_strata <- population |>
    dplyr::filter(.data[["anno"]] == pop_year) |>
    dplyr::mutate(
      .area = sprintf("%06d", as.integer(`Codice comune`)),  # match the 6-digit deaths code
      .age  = as.integer(`Eta`),
      .sex  = as.integer(.data[["sesso"]]),
      .pop  = as.numeric(.data[["numero"]])
    ) |>
    dplyr::group_by(.area, .age, .sex) |>
    dplyr::summarise(.pop = sum(.pop, na.rm = TRUE), .groups = "drop")

  # total resident population per comune (all ages/sexes) for the output
  pop_total <- pop_strata |>
    dplyr::group_by(.area) |>
    dplyr::summarise(population = sum(.pop), .groups = "drop")

  # ---- 2. observed deaths per comune x category, and per stratum x category ----
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

    # OBSERVED deaths per comune x level
    observed <- mc |>
      dplyr::group_by(.area, .lvl) |>
      dplyr::summarise(observed = sum(.data[[".w"]], na.rm = TRUE), .groups = "drop")

    # EXPECTED deaths per comune x level = sum over strata of
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
      dplyr::select(.area, col, smr, isr)
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

  # ---- 4. widen the per-category blocks (smr and isr), join totals + pop ----
  wide_smr <- blocks |>
    tidyr::pivot_wider(id_cols = ".area", names_from = "col",
                       values_from = "smr", names_glue = "{col}_smr",
                       values_fill = 0)
  wide_isr <- blocks |>
    tidyr::pivot_wider(id_cols = ".area", names_from = "col",
                       values_from = "isr", names_glue = "{col}_isr",
                       values_fill = 0)

  pop_total |>
    dplyr::left_join(wide_smr,    by = ".area") |>
    dplyr::left_join(wide_isr,    by = ".area") |>
    dplyr::left_join(total_block, by = ".area") |>
    # comuni with no deaths: every standardised value is 0
    dplyr::mutate(dplyr::across(-dplyr::all_of(c(".area", "population")),
                                ~ dplyr::coalesce(.x, 0))) |>
    dplyr::rename(!!group_var := ".area")
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
#'   \code{"comune_residenza_2023"}. Matched after zero-padding both sides.
#'
#' @return A tibble, one row per retained municipality: \code{comune} (6-digit
#'   key), \code{population}, the four raw indicators, the continuous
#'   \code{di_score}, and national \code{di_quintile}.
#' @seealso \code{\link{import_census_2023}}, \code{\link{wtd_quantile_group}}.
#' @export
build_deprivation_proxy <- function(census, mort_raw,
                                    mort_col = "comune_residenza_2023") {

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
