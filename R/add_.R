#' Attach comune geometry to a tabular layer for mapping
#'
#' Joins any comune-level table - the population denominator, a crude-rate table
#' from \code{\link{preprocess_cmr}}, a standardised table from
#' \code{\link{preprocess_smr}}, etc. - onto the ISTAT comune boundaries and
#' returns an \code{sf} ready for \code{ggplot2::geom_sf()}. Geometry is joined
#' on the fly, so the rate functions stay purely tabular.
#'
#' The join keys are matched as zero-padded \code{code_width}-character strings,
#' so a numeric shapefile code (\code{15002}) and a text table code
#' (\code{"015002"}) still line up. By default it is a left join from the table,
#' keeping every row of \code{data} (set \code{keep = "geometry"} to instead
#' keep every comune in the shapefile).
#'
#' @param data A comune-level data frame (e.g. output of \code{preprocess_cmr()}
#'   or \code{preprocess_smr()}, or the population table).
#' @param shp The comune geometries as an already-read \code{sf} object
#'   (e.g. \code{sf::st_read("geodata/Com01012025_g/Com01012025_g_WGS84.shp")}).
#' @param data_key Name of the comune-code column in \code{data}. Default
#'   \code{"comune"}.
#' @param shp_key Name of the comune-code column in the shapefile. Default
#'   \code{"PRO_COM_T"}.
#' @param keep Which side drives the rows: \code{"data"} (default, left join
#'   keeping every row of \code{data}) or \code{"geometry"} (keep every comune
#'   in the shapefile, e.g. to show the whole map with NA where there is no
#'   data).
#' @param code_width Width to zero-pad both keys to before joining. Default
#'   \code{6}.
#'
#' @return An \code{sf} object: the columns of \code{data} plus the shapefile
#'   geometry (and its attribute columns), one row per matched comune.
#'
#' @examples
#' \dontrun{
#' comuni <- sf::st_read("geodata/Com01012025_g/Com01012025_g_WGS84.shp")
#' crude  <- preprocess_cmr(mort_count, "popolazione.csv")
#'
#' crude_map <- add_geo(crude, comuni)               # comune <-> PRO_COM_T
#' library(ggplot2)
#' ggplot(crude_map) + geom_sf(aes(fill = total))
#'
#' add_geo(smr, comuni)                              # same for the standardised table
#' add_geo(pop, comuni, data_key = "Codice comune")  # or the population table
#' }
#'
#' @seealso \code{\link{preprocess_cmr}}, \code{\link{preprocess_smr}}
#' @importFrom dplyr mutate left_join right_join select all_of |>
#' @importFrom sf st_as_sf
#' @importFrom rlang .data
#' @export
add_geo <- function(data,
                    shp,
                    data_key   = "Codice comune",
                    shp_key    = "PRO_COM_T",
                    keep       = c("data", "geometry"),
                    code_width = 6) {

  keep <- match.arg(keep)
  pad  <- function(x) sprintf(paste0("%0", code_width, "d"), as.integer(x))

  # normalise both keys to a common padded string under one shared name
  data <- dplyr::mutate(data, .key = pad(.data[[data_key]]))
  shp  <- dplyr::mutate(shp,  .key = pad(.data[[shp_key]]))

  joined <- if (keep == "data") {
    dplyr::left_join(data, shp, by = ".key")    # every row of data kept
  } else {
    dplyr::right_join(data, shp, by = ".key")   # every comune in the shapefile kept
  }

  joined |>
    dplyr::select(-".key") |>
    sf::st_as_sf()
}
