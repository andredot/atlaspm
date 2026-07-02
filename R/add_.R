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
                    data_key   = "comune",
                    shp_key    = "PRO_COM_T",
                    keep       = c("data", "geometry"),
                    code_width = 6) {

  keep <- match.arg(keep)
  pad  <- function(x) sprintf(paste0("%0", code_width, "d"), as.integer(x))

  # `targets` storage (qs/rds) can strip the sf class on read, leaving a plain
  # tibble with a bare sfc geometry column - which makes dplyr::mutate() below
  # fail with "all columns must be vectors / geometry is an sfc object".
  # Re-promote defensively so the function works regardless of how shp arrived.
  shp <- sf::st_as_sf(shp)

  # normalise both keys to a common padded string under one shared name
  data <- dplyr::mutate(data, .key = pad(.data[[data_key]]))
  shp  <- dplyr::mutate(shp,  .key = pad(.data[[shp_key]]))

  # Keep the sf (shp) on the LEFT of the join. dplyr only preserves the sf class
  # when the sf is the left table; joining a plain tibble on the left returns a
  # demoted tibble carrying a bare sfc, which later breaks dplyr::mutate() with
  # "all columns must be vectors / geometry is an sfc object".
  joined <- if (keep == "data") {
    # keep every row of `data`  -> inner-style on data's keys: right_join(shp, data)
    dplyr::right_join(shp, data, by = ".key")
  } else {
    # keep every comune in `shp`
    dplyr::left_join(shp, data, by = ".key")
  }

  joined |>
    dplyr::select(-".key") |>
    sf::st_as_sf()
}

#' Attach a standardised covariate to the comune sf, preserving row order
#'
#' Left-joins a one-row-per-comune covariate table onto \code{geo}, optionally
#' imputes any missing values from the mean of contiguous neighbours, and adds a
#' z-scored copy of the covariate. The left join keeps every row of \code{geo}
#' in its original order, so the result stays aligned with an adjacency matrix
#' built from the same \code{geo} (essential for BYM2/ICAR, where C is matched
#' to the data by position). Imputation keeps all rows, so the adjacency matrix
#' built from \code{geo} remains valid without rebuilding.
#'
#' @param geo The comune \code{sf} (e.g. \code{smr_geo}), defining row order.
#' @param cov A data frame with the join key and the covariate, one row per
#'   comune (e.g. \code{ivsm_raw}: \code{comune}, \code{ivsm}).
#' @param var Name of the covariate column in \code{cov}. Default \code{"ivsm"}.
#' @param by Join key present in both. Default \code{"comune"}.
#' @param z_name Name for the standardised column. Default \code{paste0(var, "_z")}.
#' @param impute_missing If \code{TRUE} (default), comuni with a missing
#'   covariate have it filled with the mean of their contiguous (queen)
#'   neighbours' values, computed on \code{geo}'s own geometry. Keeps all rows,
#'   so the adjacency matrix stays valid. If \code{FALSE}, missing values are
#'   left as NA (which will break a Stan fit).
#' @return \code{geo} with the (possibly imputed) raw covariate and its z-scored
#'   version added, one row per comune in the original order.
#' @importFrom dplyr select all_of left_join mutate |>
#' @importFrom sf st_as_sf
#' @importFrom spdep poly2nb
#' @importFrom rlang .data :=
#' @export
add_covariate <- function(geo, cov, var = "ivsm", by = "comune",
                          z_name = paste0(var, "_z"),
                          impute_missing = TRUE) {

  # `targets` storage can strip the sf class on read, leaving a plain tibble
  # with a bare sfc geometry column - which makes dplyr::mutate() below fail with
  # "all columns must be vectors / geometry is an sfc object". Re-promote first.
  geo <- sf::st_as_sf(geo)

  cov1 <- dplyr::select(cov, dplyr::all_of(c(by, var)))
  out  <- dplyr::left_join(geo, cov1, by = by)

  x       <- out[[var]]
  missing <- which(is.na(x))

  if (length(missing) > 0) {
    if (impute_missing) {
      # contiguity neighbours from the polygons (same definition as the model's
      # adjacency); fill each NA with the mean of its non-missing neighbours.
      nb <- spdep::poly2nb(out)
      filled <- x
      for (i in missing) {
        nb_vals <- x[nb[[i]]]
        filled[i] <- mean(nb_vals, na.rm = TRUE)
      }
      still_na <- is.na(filled)
      if (any(still_na)) {
        warning(sprintf(
          "%d comuni still have no '%s' after neighbour imputation (no neighbour had a value); filled with the overall mean instead.",
          sum(still_na), var), call. = FALSE)
        filled[still_na] <- mean(x, na.rm = TRUE)
      }
      warning(sprintf(
        "%d comuni had a missing '%s'; imputed from the mean of their contiguous neighbours.",
        length(missing), var), call. = FALSE)
      out[[var]] <- filled
    } else {
      warning(sprintf(
        "%d comuni have no '%s' value (left as NA); the Stan fit will error unless handled.",
        length(missing), var), call. = FALSE)
    }
  }

  # z-score AFTER imputation, so all modelled comuni contribute to the scaling
  dplyr::mutate(out, !!z_name := as.numeric(scale(.data[[var]])))
}
