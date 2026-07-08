#' Attach area geometry to a tabular layer for mapping
#'
#' Joins any area-level table - the population denominator, a crude-rate table
#' from \code{\link{preprocess_cmr}}, a standardised table from
#' \code{\link{preprocess_smr}}, etc. - onto the boundary geometries and returns
#' an \code{sf} ready for \code{ggplot2::geom_sf()}. Geometry is joined on the
#' fly, so the rate functions stay purely tabular.
#'
#' By default the join keys are matched as zero-padded \code{code_width}-character
#' strings (\code{pad_keys = TRUE}), so a numeric shapefile code (\code{15002})
#' and a text table code (\code{"015002"}) still line up. This assumes a purely
#' \strong{numeric} ISTAT comune code. For a mixed key that contains non-numeric
#' values - e.g. an \code{area} column holding both \code{"015011"} (comune) and
#' \code{"015146_79"} (Milan comune + NIL) - set \code{pad_keys = FALSE}: the
#' keys are then compared as-is (both sides are already padded strings), avoiding
#' the \code{as.integer("015146_79") -> NA} collapse that would merge every NIL
#' into a single \code{"NA"} key.
#'
#' By default it is a left join from the table, keeping every row of \code{data}
#' (set \code{keep = "geometry"} to instead keep every area in the shapefile).
#'
#' @param data An area-level data frame (e.g. output of \code{preprocess_cmr()}
#'   or \code{preprocess_smr()}, or the population table).
#' @param shp The area geometries as an already-read \code{sf} object. For the
#'   comune analysis this is \code{pop_shp} (key \code{"PRO_COM_T"}); for the
#'   NIL-aware analysis it is \code{area_shp}, whose key column is \code{"area"}.
#' @param data_key Name of the area-code column in \code{data}. Default
#'   \code{"comune"}.
#' @param shp_key Name of the area-code column in the shapefile. Default
#'   \code{"PRO_COM_T"}. Use \code{"area"} for \code{area_shp}.
#' @param keep Which side drives the rows: \code{"data"} (default, left join
#'   keeping every row of \code{data}) or \code{"geometry"} (keep every area
#'   in the shapefile, e.g. to show the whole map with NA where there is no
#'   data).
#' @param pad_keys Whether to zero-pad both keys to \code{code_width} digits
#'   before joining. Default \code{TRUE} (correct for a purely numeric ISTAT
#'   code). Set \code{FALSE} when either key holds non-numeric values such as
#'   \code{"015146_79"}, which must not be coerced through \code{as.integer}.
#' @param code_width Width to zero-pad both keys to before joining, when
#'   \code{pad_keys = TRUE}. Default \code{6}.
#'
#' @return An \code{sf} object: the columns of \code{data} plus the shapefile
#'   geometry (and its attribute columns), one row per matched area.
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
#'
#' # NIL-aware layer: mixed comune/NIL keys on both sides, no padding
#' add_geo(area_smr, area_shp,
#'         data_key = "area", shp_key = "area", pad_keys = FALSE)
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
                    pad_keys   = TRUE,
                    code_width = 6) {
  keep <- match.arg(keep)
  # numeric ISTAT codes get zero-padded to a common width; mixed comune/NIL keys
  # (e.g. "015146_79") are compared as-is - as.integer() would turn them to NA.
  pad  <- if (pad_keys) {
    function(x) sprintf(paste0("%0", code_width, "d"), as.integer(x))
  } else {
    function(x) as.character(x)
  }
  # `targets` storage (qs/rds) can strip the sf class on read, leaving a plain
  # tibble with a bare sfc geometry column - which makes dplyr::mutate() below
  # fail with "all columns must be vectors / geometry is an sfc object".
  # Re-promote defensively so the function works regardless of how shp arrived.
  shp <- sf::st_as_sf(shp)
  # normalise both keys to a common string under one shared name
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
    # keep every area in `shp`
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

#' Expand a comune-level covariate to area (comune/NIL) keys
#'
#' Re-keys a one-row-per-comune covariate table (e.g. \code{ivsm_raw},
#' \code{deprivation}) onto the mixed comune/NIL \code{area} keys used by the
#' NIL-aware analysis, so it can be joined by \code{\link{add_covariate}} with
#' \code{by = "area"}. Each area key encodes its parent comune as the substring
#' before the first underscore: \code{"015146_79"} (a Milan NIL) maps to comune
#' \code{"015146"}, while a plain comune key such as \code{"015011"} maps to
#' itself. Every area therefore inherits its parent comune's covariate value -
#' so all Milan NILs receive the single Milan-comune value, which is the correct
#' behaviour when the covariate is only defined at comune level (there is no NIL
#' definition of IVSM).
#'
#' This is a deliberate alternative to letting \code{add_covariate()} impute the
#' NILs from neighbours: neighbour imputation would fill each NIL from adjacent
#' NILs, which are themselves all missing, cascading to surrounding comuni and
#' producing a muddy value. Assigning the true comune value is cleaner and
#' honest about the covariate's resolution.
#'
#' @param cov A comune-level covariate data frame with the join key \code{by}
#'   and the covariate column(s), one row per comune.
#' @param area_keys Character vector of all area keys to produce a row for,
#'   typically \code{area_shp$area} or \code{area_geo$area}.
#' @param by Name of the comune-code column in \code{cov}. Default
#'   \code{"comune"}.
#'
#' @return A tibble with one row per element of \code{area_keys}: an \code{area}
#'   column plus every non-key column of \code{cov}, with each area carrying its
#'   parent comune's values. Areas whose parent comune is absent from \code{cov}
#'   come through as \code{NA} (so a downstream \code{add_covariate()} still has
#'   the chance to impute or warn).
#'
#' @examples
#' \dontrun{
#' ivsm_area <- expand_cov_to_area(ivsm_raw, area_shp$area)   # comune -> area
#' add_covariate(area_geo, ivsm_area, var = "ivsm", by = "area")
#' }
#'
#' @seealso \code{\link{add_covariate}}
#' @importFrom dplyr mutate left_join select all_of |>
#' @importFrom tibble tibble
#' @importFrom rlang :=
#' @export
expand_cov_to_area <- function(cov, area_keys, by = "comune") {
  tibble::tibble(area = as.character(area_keys)) |>
    # parent comune = substring before the first "_"  ("015146_79" -> "015146",
    # "015011" -> "015011"); a plain comune key is left unchanged.
    dplyr::mutate(!!by := sub("_.*$", "", .data[["area"]])) |>
    dplyr::left_join(cov, by = by) |>
    dplyr::select(-dplyr::all_of(by))
}
