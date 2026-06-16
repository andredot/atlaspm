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

#' Crude preventable mortality rate by area of residence
#'
#' Takes the preventable-deaths table (the output of
#' \code{filter_preventable()}) together with a population denominator and
#' returns the crude preventable mortality rate for each area of residence.
#'
#' Deaths are aggregated as the \strong{sum of the \code{weight} column}, so the
#' eight causes the OECD/Eurostat list splits 50/50 between the preventable and
#' treatable categories (tuberculosis, cervical cancer, diabetes and five
#' cardiovascular causes) each count as 0.5 of a death. Set
#' \code{weight_col = NULL} to count every record as a whole death instead.
#'
#' The rate is \emph{crude} (not age-standardised): the input is assumed to be a
#' single premature-age band (0-74) already produced by
#' \code{filter_preventable()}. The crude rate is
#' \code{deaths / population * per}, i.e. per \code{per} residents
#' (100 000 by default). Every area present in \code{population} is returned,
#' including areas with no preventable deaths (rate 0).
#'
#' @param preventable Data frame of preventable deaths, one row per death, as
#'   returned by \code{filter_preventable()}. Must contain the grouping column
#'   \code{group_var} and (unless \code{weight_col = NULL}) the weight column.
#' @param population Population lookup table keyed by area of residence: either
#'   a path to a CSV or a data frame, with one row per area. Must contain
#'   \code{group_var} and \code{pop_col}.
#' @param group_var Name of the area-of-residence column, present in both
#'   \code{preventable} and \code{population}. Default \code{"area_residenza"}.
#' @param pop_col Name of the population column in \code{population}.
#'   Default \code{"population"}.
#' @param weight_col Name of the weight column in \code{preventable}; deaths are
#'   summed over it (50/50 causes contribute 0.5). Set to \code{NULL} to count
#'   each record as one death. Default \code{"weight"}.
#' @param per Rate denominator multiplier. Default \code{100000} gives the rate
#'   per 100 000 residents.
#'
#' @return A tibble with one row per area, sorted by descending rate, with the
#'   columns: \code{group_var}, \code{deaths} (weighted death count),
#'   \code{population}, and \code{crude_rate} (per \code{per} residents).
#'
#' @examples
#' \dontrun{
#' preventable <- filter_preventable(mort_raw, "preventable_lookup.csv")
#'
#' # population denominator, one row per area:
#' pop <- tibble::tibble(
#'   area_residenza = c("Lazio", "Lombardia"),
#'   population     = c(5700000, 9900000)
#' )
#'
#' preproces_cmr(preventable, pop)
#' # rate per 10,000 instead, with a differently named area column:
#' preproces_cmr(preventable, pop, group_var = "asl", per = 10000)
#' }
#'
#' @seealso \code{\link{filter_preventable}}
#' @importFrom dplyr mutate group_by across all_of summarise left_join
#'   coalesce arrange desc |>
#' @importFrom readr read_csv
#' @importFrom rlang .data
#' @export
preprocess_cmr <- function(preventable,
                          population,
                          group_var  = "area_residenza",
                          pop_col    = "population",
                          weight_col = "weight",
                          per        = 100000) {

  use_weight <- !is.null(weight_col) && weight_col %in% names(preventable)

  # population is a lookup table: accept a CSV path or a data frame
  if (is.character(population)) {
    population <- readr::read_csv(population, show_col_types = FALSE)
  }

  # 1. weighted death count per area (50/50 causes contribute 0.5)
  deaths <- preventable |>
    dplyr::mutate(.w = if (use_weight) .data[[weight_col]] else 1) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_var))) |>
    dplyr::summarise(deaths = sum(.data[[".w"]], na.rm = TRUE),
                     .groups = "drop")

  # 2. tidy the denominator and standardise its column name to "population"
  pop2 <- dplyr::select(population, dplyr::all_of(c(group_var, pop_col)))
  names(pop2)[names(pop2) == pop_col] <- "population"

  # 3. join denominator to deaths and compute the crude rate
  pop2 |>
    dplyr::left_join(deaths, by = group_var) |>
    dplyr::mutate(
      deaths     = dplyr::coalesce(.data[["deaths"]], 0),
      crude_rate = .data[["deaths"]] / .data[["population"]] * per
    ) |>
    dplyr::arrange(dplyr::desc(.data[["crude_rate"]]))
}
