#' Import and clean the mortality register
#'
#' Reads the semicolon-delimited RENCAM mortality export, keeps the fields
#' relevant to the preventable-mortality analysis, and applies the cleaning
#' steps needed before record linkage and spatial modelling.
#'
#' @details
#' All columns are read as character to preserve identifier codes (e.g. the
#' leading zeros in ISTAT municipality codes such as `015146`). Only the
#' fields used downstream are retained; `istatres`, `b61prov` and
#' `icd2`–`icd7` are dropped.
#'
#' The function applies four transformations:
#' \itemize{
#'   \item \code{sesso} is converted to a factor.
#'   \item \code{eta} is coerced to integer.
#'   \item \code{icd1} may contain several space-separated codes; only the
#'     first is kept.
#'   \item \code{nil_2023} equal to \code{"999"} marks a resident living
#'     outside the municipality of Milan and is recoded to \code{NA}.
#' }
#'
#' It also derives \code{area_residenza}, a single residence-geography
#' identifier: for residents of Milan (comune \code{015146}) it is the
#' municipality code joined to the NIL (e.g. \code{"015146_35"}); for all
#' other residents it is the municipality code alone. The NIL is prefixed
#' with the comune code so the identifier is globally unique and cannot
#' collide with a municipality code.
#'
#' @param file_path Character scalar. Path to the RENCAM CSV export
#'   (semicolon-delimited). Typically supplied by an upstream `targets`
#'   target that resolves the file location.
#'
#' @return A [tibble][tibble::tibble] with one row per death and the columns:
#' \describe{
#'   \item{sesso}{Sex, as a factor.}
#'   \item{icd1}{First ICD descriptor code of the cause of death (character).}
#'   \item{comune_residenza_2023}{ISTAT municipality of residence code (character).}
#'   \item{asst_2023}{ASST code — the higher-level health authority (character).}
#'   \item{nil_2023}{Milan NIL (local area) code, or `NA` if resident outside Milan (character).}
#'   \item{distretto_2023}{District code — the most proximal local health authority (character).}
#'   \item{eta}{Age at death (integer).}
#'   \item{causa}{Cause-of-death code (character).}
#'   \item{area_residenza}{Combined residence geography: NIL within Milan, otherwise the municipality (character).}
#' }
#'
#' @examples
#' \dontrun{
#'   import_mortality("data/rencam_2023.csv")
#' }
#'
#' @export
import_mortality <- function(file_path) {
  readr::read_delim(
    file_path,
    delim     = ";",
    col_select = c(
      sesso, icd1,
      comune_residenza_2023, asst_2023, nil_2023, distretto_2023,
      eta, causa
    ),
    col_types = readr::cols(.default = readr::col_character()),
    na        = c("", "NA"),
    trim_ws   = TRUE
  ) |>
    dplyr::mutate(
      # sex as a categorical variable
      sesso = factor(sesso),

      # age at death
      eta = as.integer(eta),

      # icd1 may hold several space-separated codes; keep only the first
      icd1 = stringr::str_extract(icd1, "^\\S+"),

      # nil = 999 is the placeholder for "resident outside Milan" -> NA
      nil_2023 = dplyr::na_if(nil_2023, "999"),

      # single residence geography:
      #   Milan (015146) -> the specific NIL; elsewhere -> the municipality
      area_residenza = dplyr::if_else(
        comune_residenza_2023 == "015146" & !is.na(nil_2023),
        paste0(comune_residenza_2023, "_", nil_2023),
        comune_residenza_2023
      )
    )
}

#' Import the ISTAT social and material vulnerability index (IVSM)
#'
#' Reads a municipal IVSM table (ISTAT \emph{8milaCensus} / Ministero
#' dell'Interno release) and returns one row per municipality, keyed on a
#' zero-padded 6-digit \code{comune} code so the result joins directly against
#' the geometry and rate tables produced elsewhere in the package (the same
#' key convention used by \code{\link{add_geo}} and \code{\link{import_population}}).
#'
#' The source file is the ISTAT-distributed IVSM, which summarises municipal
#' vulnerability through seven elementary indicators spanning the "material"
#' and "social" dimensions, with the synthetic index expressed relative to a
#' national average of 100. Only the municipality code and the index value are
#' retained by default; pass \code{keep_indicators = TRUE} to carry the seven
#' component indicators through as well.
#'
#' The municipality code is read as character and left-padded with zeros to 6
#' characters via \code{\link{pad}}, so a numeric \code{15002} in the source
#' becomes \code{"015002"} and matches the \code{PRO_COM_T}-style codes used by
#' the ISTAT boundary layer. This mirrors the padding done inside
#' \code{add_geo}, so an IVSM table imported here can be passed straight to
#' \code{add_geo(ivsm, comuni)} with the default \code{data_key = "comune"}.
#'
#' @param path Path to the IVSM file. ISTAT distributes it as a
#'   semicolon-separated CSV with a decimal comma; both are handled by the
#'   default \code{\link[readr]{locale}} (\code{"it"}).
#' @param code_col Name of the municipality-code column in the source file.
#'   Defaults to \code{"Codice comune"}, matching the population file. Inspect
#'   the header with \code{readr::read_lines(path, n_max = 1)} if unsure.
#' @param ivsm_col Name of the IVSM value column in the source file. Defaults
#'   to \code{"IVSM"}.
#' @param keep_indicators Logical; if \code{TRUE}, the seven elementary
#'   indicators are retained alongside the synthetic index. Defaults to
#'   \code{FALSE} (index only).
#'
#' @return A \code{tibble} with one row per municipality: the key column
#'   \code{comune} (character, 6-digit zero-padded) and \code{ivsm} (numeric).
#'   When \code{keep_indicators = TRUE}, the component indicator columns are
#'   appended with cleaned snake_case names.
#'
#' @export
import_ivsm <- function(path,
                        code_col = "PROCOM",
                        ivsm_col = "IVSM",
                        keep_indicators = FALSE) {
  pad  <- function(x) sprintf(paste0("%0", 6, "d"), as.integer(x))
  it_locale <- readr::locale(decimal_mark = ",", grouping_mark = ".")

  raw <- readr::read_delim(
    path,
    delim = ",",
    quote = "\"",
    locale = it_locale,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )

  if (!code_col %in% names(raw)) {
    stop("Code column '", code_col, "' not found. Columns are: ",
         paste(names(raw), collapse = ", "), call. = FALSE)
  }
  if (!ivsm_col %in% names(raw)) {
    stop("IVSM column '", ivsm_col, "' not found. Columns are: ",
         paste(names(raw), collapse = ", "), call. = FALSE)
  }

  out <- raw |>
    dplyr::transmute(
      comune = pad(.data[[code_col]]),
      ivsm   = readr::parse_number(
        .data[[ivsm_col]],
        locale = it_locale
      )
    )

  if (keep_indicators) {
    indicators <- raw |>
      dplyr::select(-dplyr::all_of(c(code_col, ivsm_col))) |>
      janitor::clean_names() |>
      dplyr::mutate(
        dplyr::across(
          dplyr::everything(),
          ~ readr::parse_number(.x, locale = it_locale)
        )
      )
    out <- dplyr::bind_cols(out, indicators)
  }

  out
}

#' Import 2023 census section data (regional xlsx files)
#'
#' Reads the per-region ISTAT 2023 permanent-census section workbooks
#' (\code{R01_..._2023_sezioni.xlsx} ... \code{R20_..._2023_sezioni.xlsx}) from a
#' directory and row-binds them into a single national section-level table. Only
#' the columns needed by \code{\link{build_deprivation_proxy}} are kept, so the
#' result stays small even though the inputs span ~350,000 census sections.
#'
#' The record layout (TRACCIATO) is shared across regions, so one column
#' selection applies to every file. The TRACCIATO workbook and any other
#' non-"_sezioni" file in the directory are skipped by the default \code{pattern}.
#' The municipality national code \code{PROCOM} is read as character and
#' zero-padded to 6 digits via \code{\link{pad}}, matching the key convention
#' used elsewhere in the package.
#'
#' @param dir Directory containing the regional xlsx files.
#' @param pattern Regex selecting the regional data files. Default
#'   \code{"_2023_sezioni\\\\.xlsx$"}, which excludes the TRACCIATO workbook.
#' @param sheet Worksheet to read from each file (name or index). Default 1.
#' @param cols Columns to keep. Defaults to the key (\code{PROCOM}) plus the
#'   counts required to build the proxy: total and 9+ population, low-education
#'   counts, employed 15-64, the ten 15-64 age bands, foreign residents, and
#'   occupied dwellings.
#'
#' @return A tibble, one row per census section: \code{PROCOM} (character,
#'   6-digit) and the selected count columns as numeric.
#' @export
import_census_2023 <- function(dir,
                               pattern = "_2023_sezioni\\.xlsx$",
                               sheet   = 1,
                               cols    = c("PROCOM",
                                           "P1", "P83",
                                           "P86", "P87", "P88",
                                           "P101",
                                           paste0("P", 17:26),  # pop 15-64 bands
                                           "ST1", "A2")) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    stop("No files matching '", pattern, "' in ", dir, call. = FALSE)
  }

  read_one <- function(f) {
    d <- readxl::read_excel(f, sheet = sheet)
    missing <- setdiff(cols, names(d))
    if (length(missing) > 0L) {
      stop("In ", basename(f), " these columns are absent: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    d <- dplyr::select(d, dplyr::all_of(cols))
    count_cols <- setdiff(cols, "PROCOM")
    dplyr::mutate(
      d,
      PROCOM = pad(as.character(PROCOM)),
      dplyr::across(dplyr::all_of(count_cols),
                    ~ suppressWarnings(as.numeric(as.character(.x))))
    )
  }

  purrr::map(files, read_one) |>
    purrr::list_rbind()
}

