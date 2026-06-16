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

#' Import the comune geometries (with population) as an sf for mapping
#'
#' Reads the ISTAT comune boundaries with \code{sf}, attaches the resident
#' population for a single year (all ages and both sexes) from the population
#' CSV, and returns an \code{sf} keyed by comune code. The result is meant to be
#' handed to \code{\link{preprocess_cmr}} / \code{\link{preprocess_smr}} via
#' their \code{geometry} argument, which left-joins it onto the rate table so
#' the output can be mapped with \code{ggplot2::geom_sf()}.
#'
#' The join is an \strong{inner} join on the comune code, so only comuni present
#' in both the population file and the shapefile are kept. Because this analysis
#' is limited to the metropolitan area of Milan, comuni outside it are dropped.
#'
#' The ISTAT code is matched between the population \code{code_col} and the
#' shapefile \code{shp_key} (\code{"PRO_COM_T"} by default, the zero-padded text
#' code such as \code{"015002"}); both are coerced to a zero-padded
#' \code{code_width}-character string, and exposed as a \code{comune} column so
#' it matches the rate functions' \code{group_var}.
#'
#' @param population Path to the semicolon-separated population CSV (or a
#'   pre-read data frame) with columns \code{Codice comune}, \code{Eta},
#'   \code{anno}, \code{sesso}, \code{Comune}, \code{numero}.
#' @param shp Path to the ISTAT comuni shapefile
#' @param pop_year Population year to keep from the CSV. Default \code{2023}.
#' @param code_col Name of the comune-code column in the population file.
#'   Default \code{"Codice comune"}.
#' @param shp_key Name of the comune-code column in the shapefile. Default
#'   \code{"PRO_COM_T"}.
#' @param group_var Name to give the comune-code column in the output, matching
#'   the rate functions' \code{group_var}. Default \code{"comune"}.
#' @param code_width Width to zero-pad both codes to before joining. Default
#'   \code{6}.
#'
#' @return An \code{sf} object, one row per matched comune, with the geometry,
#'   the shapefile attributes, a \code{comune} key column and a
#'   \code{population} column (total residents in \code{pop_year}).
#'
#' @examples
#' \dontrun{
#' geom <- import_population("popolazione.csv")
#'
#' # map population directly
#' library(ggplot2)
#' ggplot(geom) + geom_sf(aes(fill = population))
#'
#' # attach geometry to a rate table for mapping
#' crude <- preprocess_cmr(mort_count, "popolazione.csv", geometry = geom)
#' ggplot(crude) + geom_sf(aes(fill = total))
#' }
#'
#' @seealso \code{\link{preprocess_cmr}}, \code{\link{preprocess_smr}}
#' @importFrom dplyr filter mutate group_by summarise inner_join select |>
#' @importFrom readr read_delim
#' @importFrom sf st_read st_make_valid
#' @importFrom rlang .data :=
#' @export
import_population <- function(population,
                              shp,
                              pop_year   = 2023,
                              code_col   = "Codice comune",
                              shp_key    = "PRO_COM_T",
                              group_var  = "comune",
                              code_width = 6) {

  pad <- function(x) sprintf(paste0("%0", code_width, "d"), as.integer(x))

  if (is.character(population)) {
    population <- readr::read_delim(population, delim = ";", show_col_types = FALSE)
  }

  pop <- population |>
    dplyr::filter(.data[["anno"]] == pop_year) |>
    dplyr::mutate(.code = pad(.data[[code_col]])) |>
    dplyr::group_by(.code) |>
    dplyr::summarise(population = sum(.data[["numero"]], na.rm = TRUE),
                     .groups = "drop")

  comuni <- sf::st_read(shp, quiet = TRUE) |>
    sf::st_make_valid() |>
    dplyr::mutate(.code = pad(.data[[shp_key]]))

  comuni |>
    dplyr::inner_join(pop, by = ".code") |>
    dplyr::rename(!!group_var := ".code")
}
