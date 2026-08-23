#' Build a manifest of the downloaded air-quality rasters
#'
#' Scans a directory for the EEA interpolated concentration GeoTIFFs and
#' returns one row per raster, so the pipeline carries no hard-coded filenames.
#' The default \code{pattern} matches the local naming convention
#' \code{<pollutant>_avg_<yy>.tif} (e.g. \code{pm25_avg_23.tif},
#' \code{no2_avg_13.tif}).
#'
#' Only \code{.tif}/\code{.tiff} files are considered, so the browsable
#' \code{.png} previews that ship alongside each raster are ignored without
#' needing to be deleted.
#'
#' \strong{Two-digit years are expanded as \code{2000 + yy}}, which is correct
#' for this archive (the EEA PM2.5 series begins well after 2000) but would
#' silently mis-date a 1990s file. The expansion is checked against
#' \code{year_range} and errors rather than guessing.
#'
#' @param dir Directory holding the rasters. Searched recursively.
#' @param pattern Case-insensitive regex with two capture groups: pollutant and
#'   two-digit year. Override if the files are renamed.
#' @param year_range Plausible bounds for the expanded year, as a sanity check.
#'   Default \code{c(2000L, 2035L)}.
#'
#' @return A tibble with one row per raster and the columns \code{pollutant}
#'   (normalised to \code{"pm25"} / \code{"no2"}), \code{year} (integer, four
#'   digits), \code{file} (basename) and \code{path}.
#'
#' @examples
#' \dontrun{
#' manifest <- discover_aq_files(get_input_data_path("eea_aq"))
#' }
#'
#' @seealso \code{\link{build_pollution_area}}
#' @importFrom dplyr tibble arrange |>
#' @importFrom stringr str_match str_to_lower
#' @importFrom rlang .data
#' @export
discover_aq_files <- function(dir,
                              pattern    = "^(pm25|no2)_avg_(\\d{2})\\.tiff?$",
                              year_range = c(2000L, 2035L)) {

  if (!dir.exists(dir)) {
    stop("Raster directory not found: ", dir, call. = FALSE)
  }

  files <- list.files(dir, pattern = "\\.tiff?$", full.names = TRUE,
                      ignore.case = TRUE, recursive = TRUE)
  if (length(files) == 0L) {
    stop("No GeoTIFFs found under ", dir, call. = FALSE)
  }

  m    <- stringr::str_match(stringr::str_to_lower(basename(files)), pattern)
  keep <- !is.na(m[, 1L])

  if (!any(keep)) {
    stop("No filename under ", dir, " matched the expected pattern. ",
         "Supply `pattern`, or pass a manifest directly to ",
         "build_pollution_area().", call. = FALSE)
  }
  if (any(!keep)) {
    message("discover_aq_files(): ignoring ", sum(!keep),
            " unmatched raster(s): ",
            paste(basename(files[!keep]), collapse = ", "))
  }

  year <- 2000L + as.integer(m[keep, 3L])
  if (any(year < year_range[1L] | year > year_range[2L])) {
    stop("Two-digit year expansion produced an implausible year: ",
         paste(sort(unique(year)), collapse = ", "),
         ". Rename the files with four-digit years or widen `year_range`.",
         call. = FALSE)
  }

  dplyr::tibble(
    pollutant = m[keep, 2L],
    year      = year,
    file      = basename(files[keep]),
    path      = files[keep]
  ) |>
    dplyr::arrange(.data[["pollutant"]], .data[["year"]])
}


#' Read one concentration raster, cropped to the study area
#'
#' Two operations here affect the extracted values, not merely the runtime.
#'
#' \strong{Cropping first.} The source grids cover the whole of Europe at 1 km;
#' extracting ~260 small polygons from the full grid reads far more than is
#' needed. The crop is padded by \code{buffer_m} so that cells only partially
#' overlapping a boundary unit survive and can still contribute their share.
#'
#' \strong{Negative values become \code{NA}.} The EEA GeoTIFFs carry a large
#' negative sentinel outside the mapped domain, and a concentration cannot be
#' negative, so any negative cell is nodata by definition. This matters because
#' the failure is silent: a single sentinel cell clipped by a boundary polygon
#' would drag that unit's areal mean to an absurd value while still returning a
#' perfectly finite number, and the unit would only surface much later as an
#' influential point in \code{pareto_k_table()}.
#'
#' @param path Path to a single GeoTIFF.
#' @param geo An \code{sf} of the study area, used only for the crop extent.
#'   Any CRS; it is transformed internally.
#' @param buffer_m Padding in metres around the study-area bounding box.
#'   Default 5000.
#'
#' @return A single-layer \code{SpatRaster} in the raster's native CRS
#'   (EPSG:3035 for the EEA products), cropped and with nodata set.
#'
#' @seealso \code{\link{extract_aq}}
#' @importFrom terra rast nlyr crs crop ext classify
#' @importFrom sf st_as_sf st_transform st_bbox
#' @export
import_aq_raster <- function(path, geo, buffer_m = 5000) {

  geo <- sf::st_as_sf(geo)

  r <- terra::rast(path)
  if (terra::nlyr(r) > 1L) r <- r[[1L]]

  box <- geo |>
    sf::st_transform(terra::crs(r)) |>
    sf::st_bbox()

  box[["xmin"]] <- box[["xmin"]] - buffer_m
  box[["ymin"]] <- box[["ymin"]] - buffer_m
  box[["xmax"]] <- box[["xmax"]] + buffer_m
  box[["ymax"]] <- box[["ymax"]] + buffer_m

  r <- terra::crop(
    r,
    terra::ext(box[["xmin"]], box[["xmax"]], box[["ymin"]], box[["ymax"]])
  )

  # negative sentinel -> NA (right = FALSE keeps exact zeros)
  terra::classify(r, cbind(-Inf, 0, NA), right = FALSE)
}


#' Areal means of every raster in a manifest
#'
#' Computes one areal mean per areal unit per raster, with partial-cell
#' weighting via \code{exactextractr::exact_extract}. Partial-cell weighting is
#' not cosmetic at NIL scale: on a 1 km grid several central NILs are smaller
#' than a single cell, so a centroid rule would hand identical values to
#' adjacent NILs and a majority-cell rule would discard most of the overlap
#' information.
#'
#' Each concentration column is accompanied by a \code{coverage_*} column
#' giving the fraction of the unit's area that fell on a non-\code{NA} cell.
#' A unit with low coverage has a mean computed from part of itself only, which
#' is worth knowing before it turns up as an outlier.
#'
#' Row order follows \code{geo}, which is what keeps the result aligned with an
#' adjacency matrix built from the same \code{geo} - the BYM2/ESF models match
#' \code{C} to the data by position, not by key.
#'
#' @param geo An \code{sf} of areal units carrying the key column \code{area}
#'   (typically \code{area_shp}).
#' @param manifest A tibble from \code{\link{discover_aq_files}}, or any tibble
#'   with \code{pollutant}, \code{year} and \code{path}.
#'
#' @return A tibble with one row per unit of \code{geo}, in the same order:
#'   \code{area}, then \code{<pollutant>_<year>} and
#'   \code{coverage_<pollutant>_<year>} for each raster.
#'
#' @seealso \code{\link{build_pollution_area}}, \code{\link{import_aq_raster}}
#' @importFrom sf st_as_sf st_transform
#' @importFrom terra crs
#' @importFrom exactextractr exact_extract
#' @importFrom dplyr tibble
#' @export
extract_aq <- function(geo, manifest) {

  geo <- sf::st_as_sf(geo)

  if (!"area" %in% names(geo)) {
    stop("`geo` must carry an 'area' key column.", call. = FALSE)
  }
  if (nrow(manifest) == 0L) {
    stop("`manifest` is empty - nothing to extract.", call. = FALSE)
  }

  out <- dplyr::tibble(area = as.character(geo[["area"]]))

  for (i in seq_len(nrow(manifest))) {

    r     <- import_aq_raster(manifest[["path"]][i], geo)
    geo_p <- sf::st_transform(geo, terra::crs(r))
    nm    <- paste0(manifest[["pollutant"]][i], "_", manifest[["year"]][i])

    res <- exactextractr::exact_extract(
      r, geo_p,
      fun = function(df) {
        ok  <- !is.na(df[["value"]])
        w   <- df[["coverage_fraction"]]
        data.frame(
          mean = if (any(ok)) sum(df[["value"]][ok] * w[ok]) / sum(w[ok])
                 else NA_real_,
          coverage = if (sum(w) > 0) sum(w[ok]) / sum(w) else NA_real_
        )
      },
      summarize_df = TRUE,
      progress     = FALSE
    )

    out[[nm]]                      <- res[["mean"]]
    out[[paste0("coverage_", nm)]] <- res[["coverage"]]
  }

  out
}


#' Build the area-level pollution table
#'
#' Top-level entry point for the pollution layer: discovers the rasters,
#' extracts them onto the areal units, and checks the result. With the 2013 and
#' 2023 maps in place the returned columns are \code{pm25_2013},
#' \code{pm25_2023}, \code{no2_2013} and \code{no2_2023}, plus their
#' \code{coverage_*} companions.
#'
#' \strong{The two years are not two exposure measurements.} 2023 is the
#' exposure, sitting at the midpoint of the 2022-2024 mortality window; 2013
#' exists solely so that \code{\link{check_aq_ranks}} can test whether the
#' spatial ordering of areas is persistent. Do not average the two years
#' together and do not fit a model on 2013.
#'
#' @param area_shp An \code{sf} of areal units with the key column \code{area}.
#' @param dir Directory of downloaded rasters. Ignored when \code{manifest} is
#'   supplied.
#' @param manifest Optional pre-built manifest, bypassing
#'   \code{\link{discover_aq_files}}.
#' @param min_coverage Units whose non-\code{NA} raster coverage falls below
#'   this trigger a warning. Default 0.95.
#'
#' @return A tibble keyed on \code{area}, in \code{area_shp} row order.
#'
#' @examples
#' \dontrun{
#' pollution_area <- build_pollution_area(
#'   area_shp,
#'   dir = get_input_data_path("eea_aq")
#' )
#' }
#'
#' @seealso \code{\link{check_aq_ranks}}, \code{\link{add_pollution}}
#' @importFrom rlang .data
#' @export
build_pollution_area <- function(area_shp, dir = NULL, manifest = NULL,
                                 min_coverage = 0.95) {

  if (is.null(manifest)) {
    if (is.null(dir)) {
      stop("Supply either `dir` or `manifest`.", call. = FALSE)
    }
    manifest <- discover_aq_files(dir)
  }

  out <- extract_aq(area_shp, manifest)

  cov_cols <- grep("^coverage_", names(out), value = TRUE)
  val_cols <- setdiff(names(out), c("area", cov_cols))

  low <- vapply(out[cov_cols],
                function(x) sum(x < min_coverage, na.rm = TRUE), integer(1L))
  if (any(low > 0L)) {
    warning("Areal units below ", min_coverage, " raster coverage: ",
            paste(names(low)[low > 0L], low[low > 0L],
                  sep = "=", collapse = ", "),
            ". Their means are computed from part of the unit only.",
            call. = FALSE)
  }

  n_na <- vapply(out[val_cols], function(x) sum(is.na(x)), integer(1L))
  if (any(n_na > 0L)) {
    warning("Missing concentrations after extraction: ",
            paste(names(n_na)[n_na > 0L], n_na[n_na > 0L],
                  sep = "=", collapse = ", "),
            ". Check that the cropped extent covers the whole ATS territory.",
            call. = FALSE)
  }

  out
}


#' Attach the pollution covariates to a modelling \code{sf}
#'
#' Thin wrapper that pushes each pollutant column through
#' \code{\link{add_covariate}}, so the pollution layer enters the models by
#' exactly the same route as IVSM and the deprivation index: left join on
#' \code{area} preserving row order, neighbour imputation for any missing unit,
#' and a \code{_z} column standardised \emph{after} imputation.
#'
#' Only the exposure year is attached by default. The 2013 columns are
#' diagnostic and have no business in a design matrix.
#'
#' @param geo The modelling \code{sf} (e.g. \code{smr_geo}), defining row order.
#' @param pollution_area Output of \code{\link{build_pollution_area}}.
#' @param vars Character vector of pollutant columns to attach. Default
#'   \code{c("pm25_2023", "no2_2023")}.
#' @param impute_missing Passed to \code{\link{add_covariate}}. Default
#'   \code{TRUE}.
#'
#' @return \code{geo} with each requested column and its \code{_z} version
#'   added, in the original row order.
#'
#' @examples
#' \dontrun{
#' smr_geo_poll <- add_pollution(smr_geo, pollution_area)
#' # -> pm25_2023, pm25_2023_z, no2_2023, no2_2023_z
#' }
#'
#' @seealso \code{\link{add_covariate}}
#' @importFrom sf st_as_sf
#' @export
add_pollution <- function(geo, pollution_area,
                          vars = c("pm25_2023", "no2_2023"),
                          impute_missing = TRUE) {

  geo <- sf::st_as_sf(geo)

  missing_vars <- setdiff(vars, names(pollution_area))
  if (length(missing_vars) > 0L) {
    stop("Not present in `pollution_area`: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  for (v in vars) {
    geo <- add_covariate(geo, pollution_area, var = v, by = "area",
                         impute_missing = impute_missing)
  }

  geo
}


#' Rank stability of the spatial pattern between two years
#'
#' The Section 8.3 selection rule prefers the pollutant whose spatial pattern is more
#' persistent, on the reasoning that chronic-disease mortality responds to
#' long-run exposure, so a field that reshuffles between years is a poor proxy
#' for the exposure that actually generated the deaths.
#'
#' Spearman rather than Pearson, because what must be stable is the
#' \strong{ordering} of areas, not the level. Concentrations fell substantially
#' across the decade for both pollutants; in a cross-sectional model that
#' common shift is absorbed by the intercept and is irrelevant to
#' identification.
#'
#' \strong{Read the coefficient asymmetrically.} The two maps were produced a
#' decade apart under successive revisions of the EEA mapping methodology and
#' with different underlying model inputs, so part of any observed rank
#' movement is processing rather than the field. A high correlation is
#' therefore trustworthy evidence of persistence; a low one is ambiguous
#' between genuine instability and method drift, and should be reported as
#' ambiguous rather than treated as evidence against the pollutant.
#'
#' @param pollution_area Output of \code{\link{build_pollution_area}}.
#' @param years Length-2 integer vector, older year first. Default
#'   \code{c(2013L, 2023L)}.
#'
#' @return A tibble with one row per pollutant: \code{pollutant}, the two
#'   column names, the mean and SD in each year, and Spearman \code{rho}.
#'
#' @examples
#' \dontrun{
#' check_aq_ranks(pollution_area)
#' }
#'
#' @importFrom dplyr tibble bind_rows
#' @importFrom stats cor sd
#' @export
check_aq_ranks <- function(pollution_area, years = c(2013L, 2023L)) {

  if (length(years) != 2L) {
    stop("`years` must be length 2, older first.", call. = FALSE)
  }

  value_cols <- grep("^[a-z0-9]+_\\d{4}$",
                     setdiff(names(pollution_area), "area"), value = TRUE)
  value_cols <- value_cols[!grepl("^coverage_", value_cols)]
  stems      <- unique(sub("_\\d{4}$", "", value_cols))

  rows <- lapply(stems, function(p) {

    a <- paste0(p, "_", years[1L])
    b <- paste0(p, "_", years[2L])
    if (!all(c(a, b) %in% names(pollution_area))) return(NULL)

    dplyr::tibble(
      pollutant = p,
      col_old   = a,
      col_new   = b,
      mean_old  = mean(pollution_area[[a]], na.rm = TRUE),
      mean_new  = mean(pollution_area[[b]], na.rm = TRUE),
      sd_old    = stats::sd(pollution_area[[a]], na.rm = TRUE),
      sd_new    = stats::sd(pollution_area[[b]], na.rm = TRUE),
      rho       = stats::cor(pollution_area[[a]], pollution_area[[b]],
                             method = "spearman", use = "complete.obs")
    )
  })

  dplyr::bind_rows(rows)
}


#' Pre-specified pollutant selection diagnostic (Section 8.3)
#'
#' Every criterion in this table is computed on the \strong{exposure alone}.
#' That is the entire point of the rule: selecting the pollutant that shows the
#' stronger association with mortality is outcome-dependent selection and would
#' forfeit the pre-specification defence for the whole analysis. Run this,
#' record the output with the date, and only then fit an outcome model.
#'
#' Preferred direction by column:
#' \describe{
#'   \item{\code{sd}, \code{range}, \code{cv}}{Higher - more contrast to
#'     estimate from.}
#'   \item{\code{sd_milan}}{Higher - tests whether the 1 km grid buys anything
#'     across the NILs. Near zero means the disaggregation gained nothing for
#'     this pollutant.}
#'   \item{\code{cor_z_low}}{\strong{Lower} - a covariate strongly correlated
#'     with the least-shrunk canonical regressor is precisely the configuration
#'     that produces spatial confounding (Section 10.3).}
#'   \item{\code{cor_di}}{Lower - less mutual attenuation against deprivation.}
#'   \item{\code{rho}}{Higher - a more persistent long-run proxy.}
#' }
#'
#' The prior expectation is attached to the result as the
#' \code{"prior_expectation"} attribute so that it travels with the numbers and
#' cannot be quietly rewritten after the fact.
#'
#' @param dat A data frame with \code{area}, the pollutant columns for
#'   \code{year}, and optionally the deprivation column.
#' @param year The exposure year the models will use. Default \code{2023L}.
#' @param di_col Name of the deprivation column. Default \code{"di_score"}.
#' @param z_low Optional numeric vector: the least-shrunk canonical regressor
#'   from the Section 10.3 diagnostic, in \code{dat} row order. When \code{NULL},
#'   \code{cor_z_low} is \code{NA} and the decision cannot yet be closed - this
#'   is usually the criterion that separates the two pollutants.
#' @param ranks Optional output of \code{\link{check_aq_ranks}}, merged in.
#'
#' @return A tibble, one row per pollutant, with a \code{"prior_expectation"}
#'   attribute.
#'
#' @examples
#' \dontrun{
#' pollutant_selection_table(
#'   dplyr::left_join(pollution_area, deprivation_area, by = "area"),
#'   ranks = aq_ranks
#' )
#' }
#'
#' @seealso \code{\link{check_aq_ranks}}
#' @importFrom dplyr tibble bind_rows left_join all_of select
#' @importFrom stats cor sd
#' @export
pollutant_selection_table <- function(dat, year = 2023L,
                                      di_col = "di_score",
                                      z_low  = NULL,
                                      ranks  = NULL) {

  is_nil <- grepl("^015146_", dat[["area"]])
  cols   <- grep(paste0("^[a-z0-9]+_", year, "$"),
                 setdiff(names(dat), "area"), value = TRUE)
  cols   <- cols[!grepl("^coverage_", cols)]

  if (length(cols) == 0L) {
    stop("No pollutant columns found for year ", year, ".", call. = FALSE)
  }

  out <- dplyr::bind_rows(lapply(cols, function(cl) {

    v <- dat[[cl]]

    dplyr::tibble(
      pollutant = sub("_\\d{4}$", "", cl),
      column    = cl,
      n_missing = sum(is.na(v)),
      mean      = mean(v, na.rm = TRUE),
      sd        = stats::sd(v, na.rm = TRUE),
      min       = min(v, na.rm = TRUE),
      max       = max(v, na.rm = TRUE),
      range     = max(v, na.rm = TRUE) - min(v, na.rm = TRUE),
      cv        = stats::sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE),
      sd_milan  = if (sum(is_nil) > 1L) stats::sd(v[is_nil], na.rm = TRUE)
                  else NA_real_,
      cor_z_low = if (!is.null(z_low))
                    stats::cor(v, z_low, use = "complete.obs")
                  else NA_real_,
      cor_di    = if (di_col %in% names(dat))
                    stats::cor(v, dat[[di_col]], use = "complete.obs")
                  else NA_real_
    )
  }))

  if (!is.null(ranks) && nrow(ranks) > 0L) {
    out <- dplyr::left_join(
      out,
      dplyr::select(ranks, dplyr::all_of(c("pollutant", "rho"))),
      by = "pollutant"
    )
  }

  attr(out, "prior_expectation") <- paste(
    "NO2 is expected to be the better-identified exposure: traffic-driven,",
    "with steep gradients, loading onto higher-frequency eigenvectors that the",
    "ICAR shrinks hard. PM2.5 in the Po basin is largely secondary aerosol",
    "forming at basin scale, so it varies at low spatial frequency and",
    "competes with the least-shrunk canonical regressors. Caveat: in the EEA",
    "product NO2's fine structure is partly generated by a road-data",
    "regressor, so a low cor_z_low for NO2 is weaker evidence than it appears.",
    "If the diagnostic contradicts this expectation, follow the diagnostic and",
    "say so."
  )

  out
}
