
# Here below put your small tiny supporting functions -------------
view_in_excel <- function(.data) {
  if (interactive()) {
    tmp <- fs::file_temp("excel", ext = "csv")
    readr::write_excel_csv(.data, tmp)
    fs::file_show(tmp)
  }
  invisible(.data)
}


extract_fct_names <- function(path) {
  readr::read_lines(path) |>
    stringr::str_extract_all("^.*(?=`? ?<- ?function)") |>
    unlist() |>
    purrr::compact() |>
    stringr::str_remove_all("[\\s`]+")
}



get_input_data_path <- function(x = "") {
  file.path(
    Sys.getenv("PRJ_SHARED_PATH"),
    Sys.getenv("INPUT_DATA_FOLDER"),
    x
  ) |>
    normalizePath()
}

get_output_data_path <- function(x = "") {
  file.path(
    Sys.getenv("PRJ_SHARED_PATH"),
    Sys.getenv("OUTPUT_DATA_FOLDER"),
    x
  ) |>
    normalizePath(mustWork = FALSE)
}


share_objects <- function(obj_list) {
  now <- lubridate::now() |>
    stringr::str_remove_all("\\W+") |>
    stringr::str_sub(1, 12)

  file_name_now <- stringr::str_c(
    names(obj_list), "-", now, ".rds"
  )

  file_name_latest <- stringr::str_c(
    names(obj_list), "-", "latest", ".rds"
  )

  obj_paths_now <- get_output_data_path(file_name_now) |>
    normalizePath(mustWork = FALSE) |>
    purrr::set_names(names(obj_list))

  obj_paths_latest <- get_output_data_path(file_name_latest) |>
    normalizePath(mustWork = FALSE) |>
    purrr::set_names(names(obj_list))

  # Those must be RDS
  obj_list |>
    purrr::walk2(obj_paths_now, readr::write_rds)
  obj_list |>
    purrr::walk2(obj_paths_latest, readr::write_rds)

  obj_paths_latest
}

pad  <- function(x) sprintf(paste0("%0", 6, "d"), as.integer(x))


#' Round half away from zero
#'
#' Base R's [round()] implements IEEE 754 round-half-to-even ("banker's
#' rounding"), so `round(0.5)` is `0` and `round(2.5)` is `2`. The OECD/Eurostat
#' fractional allocation described in the methods rounds half **up**, and the
#' difference is not academic: half-weight counts are exactly the values at
#' which the two rules disagree, and they are common in this data (every one of
#' the nine split causes produces them).
#'
#' @param x Numeric vector.
#' @param digits Number of decimal places. Default `0`.
#'
#' @return Numeric vector, rounded half away from zero.
#' @examples
#' round_half_up(c(0.5, 1.5, 2.5, -0.5))   # 1 2 3 -1
#' round(c(0.5, 1.5, 2.5, -0.5))           # 0 2 2  0
#' @export
round_half_up <- function(x, digits = 0) {
  m <- 10^digits
  sign(x) * trunc(abs(x) * m + 0.5) / m
}


#' Build the wide-table column name for a classification level
#'
#' [preprocess_cmr()] and [preprocess_smr()] name their output columns
#' `<prefix>_<janitor::make_clean_names(label)>`. Downstream code that needs one
#' specific column - the cerebrovascular tracer, the all-cancer negative
#' control - must derive the name the same way rather than hardcoding a guess,
#' because `make_clean_names()` transliterates in ways that are hard to predict
#' ("Lifestyle and NCDs" does not become `lifestyle_and_ncds`).
#'
#' @param prefix Classification prefix: `"C"` (cause), `"G"` (group) or `"M"`
#'   (mechanism).
#' @param label The level label exactly as it appears in the lookup, e.g.
#'   `"Cerebrovascular diseases"`.
#' @param suffix Optional suffix, e.g. `"_obs"` or `"_exp"`.
#'
#' @return A length-1 character column name.
#' @examples
#' smr_col("C", "Cerebrovascular diseases", "_obs")
#' @importFrom janitor make_clean_names
#' @export
smr_col <- function(prefix, label, suffix = "") {
  paste0(prefix, "_", janitor::make_clean_names(label), suffix)
}


#' Assert that a set of columns exists, failing with the near misses
#'
#' A plain `stopifnot()` on column presence tells you what is missing but not
#' what is there instead, which is the information you actually need when a
#' name has been transliterated or a join has silently renamed something.
#'
#' @param data A data frame or `sf`.
#' @param cols Character vector of required column names.
#' @param what Short description of `data`, used in the error message.
#'
#' @return `data`, invisibly.
#' @export
require_cols <- function(data, cols, what = "data") {
  missing <- setdiff(cols, names(data))
  if (length(missing)) {
    near <- names(data)[
      vapply(names(data),
             function(n) any(utils::adist(n, missing) <= 4L),
             logical(1))
    ]
    stop(
      "`", what, "` is missing: ", paste(missing, collapse = ", "), ".\n",
      if (length(near)) paste0("  Closest available: ",
                               paste(near, collapse = ", "), "\n") else "",
      "  All available: ", paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(data)
}


#' Map register sex codes onto the ISTAT numeric convention
#'
#' The ReNCaM extract codes sex as `"M"`/`"F"`; `pop_finale.csv` and
#' `pop_nil.csv` use the ISTAT convention `1 = Maschi`, `2 = Femmine`. Indirect
#' standardisation joins deaths to population on this code, so the two sides
#' must agree. Anything unrecognised is a hard error rather than an `NA`: a
#' silently dropped sex would bias every expected count without warning.
#'
#' @param x Character, factor or numeric vector of sex codes. Accepts
#'   `"M"`/`"F"`, `"MASCHIO"`/`"FEMMINA"`, `"1"`/`"2"` and `1`/`2`, in any case
#'   and with surrounding whitespace.
#' @param allow_na Whether missing values are permitted. Default `FALSE`.
#'
#' @return Integer vector of `1L` (male) and `2L` (female).
#' @examples
#' recode_sex(c("M", "F", "m"))
#' @export
recode_sex <- function(x, allow_na = FALSE) {
  s <- toupper(trimws(as.character(x)))

  out <- rep(NA_integer_, length(s))
  out[s %in% c("M", "1", "MASCHIO", "MASCHI", "MALE")]     <- 1L
  out[s %in% c("F", "2", "FEMMINA", "FEMMINE", "FEMALE")]  <- 2L

  bad <- is.na(out) & !(is.na(s) & allow_na)
  if (any(bad)) {
    stop("Unrecognised sex code(s): ",
         paste(unique(shQuote(s[bad])), collapse = ", "),
         ". Expected M/F or 1/2 (ISTAT: 1 = male, 2 = female).",
         call. = FALSE)
  }
  out
}


#' Check that every requested population year is present
#'
#' Silently returning a smaller denominator because one year of the population
#' file is absent would inflate every rate in the study without any visible
#' symptom. Fail instead.
#'
#' @param population Population table with an `anno` column.
#' @param pop_year Requested year(s).
#'
#' @return `population`, with `anno` coerced to the same type as `pop_year`.
#' @export
check_pop_years <- function(population, pop_year) {
  require_cols(population, "anno", "population")

  # `anno` arrives as character from read_delim() but pop_year is usually
  # numeric; compare on a common type rather than relying on coercion inside
  # %in%, which would match nothing and filter the table to zero rows.
  population[["anno"]] <- as.character(population[["anno"]])
  pop_year_chr <- as.character(pop_year)

  missing <- setdiff(pop_year_chr, unique(population[["anno"]]))
  if (length(missing)) {
    stop("Population file has no rows for year(s): ",
         paste(missing, collapse = ", "), ".\n",
         "  Available: ", paste(sort(unique(population[["anno"]])),
                                collapse = ", "),
         call. = FALSE)
  }
  population
}


#' Assert that deaths and population use the same sex coding
#'
#' Indirect standardisation joins on `(age, sex)`. If the two sides disagree
#' about what `1` means, the join still succeeds and the expected counts are
#' quietly wrong for every area. This checks the value sets match and that
#' neither side is a subset of the other.
#'
#' @param deaths_sex,pop_sex Integer vectors of sex codes.
#'
#' @return `TRUE`, invisibly.
#' @export
assert_sex_alignment <- function(deaths_sex, pop_sex) {
  d <- sort(unique(stats::na.omit(deaths_sex)))
  p <- sort(unique(stats::na.omit(pop_sex)))

  if (!identical(as.integer(d), as.integer(p))) {
    stop("Sex codes differ between deaths (", paste(d, collapse = "/"),
         ") and population (", paste(p, collapse = "/"), ").\n",
         "  Both must use the ISTAT convention 1 = male, 2 = female. ",
         "See recode_sex().", call. = FALSE)
  }
  if (any(is.na(deaths_sex))) {
    stop(sum(is.na(deaths_sex)), " death record(s) have a missing sex code; ",
         "they would be dropped from the standard schedule.", call. = FALSE)
  }
  invisible(TRUE)
}
