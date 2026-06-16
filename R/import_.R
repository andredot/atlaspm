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
