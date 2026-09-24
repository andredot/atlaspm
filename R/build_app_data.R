# The Shiny app's data bundle, as a pipeline product -----------------------------
#
# This is a target, not a script. Everything it needs arrives as an argument, so
# `targets` can see the dependency edges and rebuild the bundle when a fit,
# a covariate or the geography changes. Nothing here calls tar_read(): a target
# that reads other targets behind the pipeline's back is invisible to the
# dependency graph, and the failure mode is an app that renders yesterday's
# results perfectly.
#
# WHY A BUNDLE RATHER THAN LETTING THE APP READ THE STORE
#
# The app needs quantities that live inside the stanfit objects. Reading those
# at app start would pull hundreds of megabytes of posterior draws into every
# session to recompute something that changes only when the pipeline reruns.
# The bundle is written once per `tar_make()` and read in about a second.


#' Stop with a legible message when expected columns are absent
#'
#' The default error from a downstream `$` on a missing column names neither the
#' object nor the fix, and this function is assembling fifteen of them.
#'
#' @param x A data frame or `sf`.
#' @param cols Required column names.
#' @param what Name to use for `x` in the message.
#' @return `x`, invisibly.
#' @noRd
.need_cols <- function(x, cols, what) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) {
    stop("`", what, "` is missing: ", paste(miss, collapse = ", "), ".\n",
         "  Columns present: ",
         paste(utils::head(names(x), 25), collapse = ", "),
         if (length(names(x)) > 25) ", ..." else "", call. = FALSE)
  }
  invisible(x)
}


#' Human-readable, searchable names for the modelling keys
#'
#' The modelling key is `"015146_79"` for a Milan NIL and a bare ISTAT code
#' elsewhere. Nobody types that into a search box, so both sides are resolved to
#' names. A failed lookup degrades to the key rather than to `NA`: an ugly entry
#' in the search list is a nuisance, a blank one is a hole the user cannot find
#' an area through.
#'
#' @param area Character vector of modelling keys.
#' @param nil_shp,pop_shp The NIL and comune shapefiles.
#' @return A character vector of labels, unique and non-empty.
#' @noRd
.area_labels <- function(area, nil_shp, pop_shp) {

  nil_tab <- sf::st_drop_geometry(nil_shp)
  nil_name_col <- intersect(c("NIL", "nil", "NIL_NAME", "DENOM_NIL",
                              "Nome", "NOME", "nome"), names(nil_tab))[1]
  if (is.na(nil_name_col)) {
    warning("No NIL name column in `nil_shp`; Milan entries will show ID ",
            "numbers in the app's search box. Columns present: ",
            paste(names(nil_tab), collapse = ", "), call. = FALSE)
    nil_lookup <- data.frame(id = nil_tab[["ID_NIL"]],
                             nm = paste("NIL", nil_tab[["ID_NIL"]]))
  } else {
    nil_lookup <- data.frame(id = nil_tab[["ID_NIL"]],
                             nm = as.character(nil_tab[[nil_name_col]]))
  }
  nil_lookup <- nil_lookup[!duplicated(nil_lookup[["id"]]), , drop = FALSE]

  com_tab <- sf::st_drop_geometry(pop_shp)
  com_name_col <- intersect(c("COMUNE", "Comune", "NOME_COM", "COMUNE_A",
                              "DEN_COM"), names(com_tab))[1]
  com_lookup <- data.frame(
    code = as.character(com_tab[["PRO_COM_T"]]),
    nm   = if (is.na(com_name_col)) as.character(com_tab[["PRO_COM_T"]])
           else as.character(com_tab[[com_name_col]])
  )
  com_lookup <- com_lookup[!duplicated(com_lookup[["code"]]), , drop = FALSE]

  area   <- as.character(area)
  is_nil <- grepl("_", area, fixed = TRUE)
  out    <- character(length(area))

  id <- suppressWarnings(as.integer(sub("^.*_", "", area[is_nil])))
  nm <- nil_lookup[["nm"]][match(id, nil_lookup[["id"]])]
  out[is_nil]  <- paste0(ifelse(is.na(nm), paste("NIL", id), nm), " (Milan)")
  out[!is_nil] <- com_lookup[["nm"]][match(area[!is_nil], com_lookup[["code"]])]

  blank <- is.na(out) | !nzchar(out)
  out[blank] <- area[blank]

  # Two comuni can share a name across provinces; the key disambiguates.
  dup <- unique(out[duplicated(out)])
  if (length(dup)) {
    hit <- out %in% dup
    out[hit] <- paste0(out[hit], " [", area[hit], "]")
  }
  out
}


#' Simplify polygons for the browser, in a unit that means something
#'
#' `st_simplify()`'s `dTolerance` is in the units of the CRS. Simplifying after
#' a transform to EPSG:4326 applies a tolerance in DEGREES, and a tolerance of
#' 10 there would erase Lombardy - without emptying the geometries, so an
#' `st_is_empty()` check would not catch it. Simplify in the projected CRS,
#' then transform; if the input is already geographic, convert the tolerance
#' instead of applying the wrong one.
#'
#' @param geo An `sf`.
#' @param metres Display tolerance in metres.
#' @return `geo` simplified where safe, in EPSG:4326.
#' @noRd
.simplify_for_web <- function(geo, metres = 15) {

  if (sf::st_is_longlat(geo)) {
    tol <- metres / 111320          # ~111 km per degree of latitude
  } else {
    u <- tryCatch(as.character(sf::st_crs(geo)$units),
                  error = function(e) NA_character_)
    if (is.na(u) || !identical(u, "m")) {
      warning("Projected CRS is not in metres (units: ",
              if (is.na(u)) "unknown" else u, "); skipping simplification ",
              "rather than guessing a tolerance.", call. = FALSE)
      return(sf::st_transform(geo, 4326))
    }
    tol <- metres
  }

  small <- tryCatch(
    sf::st_simplify(geo, dTolerance = tol, preserveTopology = TRUE),
    error = function(e) NULL
  )
  if (is.null(small) || nrow(small) != nrow(geo) ||
      any(sf::st_is_empty(small))) {
    return(sf::st_transform(geo, 4326))
  }
  sf::st_transform(sf::st_make_valid(small), 4326)
}


#' Assemble the data bundle the exploratory Shiny app reads
#'
#' Collects the pipeline outputs the app needs into one compact file. Called
#' from a `format = "file"` target so that `targets` tracks the written bundle
#' and rebuilds it whenever an input changes.
#'
#' Every input is an argument rather than a `tar_read()` call, which is what
#' makes the bundle a real member of the dependency graph instead of a
#' side-effect that quietly goes stale.
#'
#' @param path Where to write the bundle. The directory is created if needed.
#' @param smr_geo_full The modelling frame. Row order must match the fits.
#' @param nil_shp,pop_shp Shapefiles, used only to resolve area names.
#' @param gp_density_area Primary-care indicator, for the unweighted variant
#'   that [add_covariate()] does not carry into the modelling frame.
#' @param crude_geo Crude rates by area.
#' @param aug_total [augment_bym2()] output for the all-causes model (M3).
#' @param aug_models Named list of [augment_bym2()] outputs for the adjusted
#'   models, e.g. `list(M5 = aug_M5, M6 = aug_M6)`. The names become the app's
#'   model selector.
#' @param mech_bym2 Per-mechanism smoothed surfaces (`smr_geo_mech_bym2`).
#' @param excess Named list of [residual_excess()] outputs, one per model, with
#'   the same names as `aug_models`.
#' @param attribution Named list of [attribution_grid()] tibbles, same names.
#' @param smr_geo_allage The all-age frame carrying the I63 outcome.
#' @param aug_i63,aug_i63_adj Unadjusted and adjusted I63 surfaces.
#' @param stroke_centres Geocoded stroke centres.
#' @param n_years,rr_threshold,rr_threshold_stroke,prob_cutoff Study constants,
#'   passed through to the app so its captions cannot drift from the pipeline.
#' @param aq_year Air-quality year, selecting the NO2 and PM2.5 columns.
#' @param simplify_m Display simplification tolerance, in metres.
#' @param compress Passed to [saveRDS()]. `TRUE` (gzip) is the default because
#'   this now runs on every `tar_make()`; `"xz"` is smaller but much slower.
#'
#' @return `path`, invisibly, so the target can be `format = "file"`.
#'
#' @examples
#' \dontrun{
#' # from _targets.R, with every dependency declared
#' tar_target(app_data, build_app_data("inst/shiny/atlas/app_data.rds",
#'                                     smr_geo_full = smr_geo_full, ...),
#'            format = "file")
#' }
#' @seealso [run_atlas()]
#' @export
build_app_data <- function(path,
                           smr_geo_full,
                           nil_shp,
                           pop_shp,
                           gp_density_area,
                           crude_geo,
                           aug_total,
                           aug_models,
                           mech_bym2,
                           excess,
                           attribution,
                           smr_geo_allage,
                           aug_i63,
                           aug_i63_adj,
                           stroke_centres,
                           n_years             = 3,
                           rr_threshold        = 1.10,
                           rr_threshold_stroke = 1.20,
                           prob_cutoff         = 0.80,
                           aq_year             = 2023L,
                           simplify_m          = 15,
                           compress            = TRUE) {

  no2  <- sprintf("no2_%d",  aq_year)
  pm25 <- sprintf("pm25_%d", aq_year)

  models <- names(aug_models)
  if (is.null(models) || !length(models) || anyDuplicated(models)) {
    stop("`aug_models` must be a named list with unique names; the names ",
         "become the app's model selector.", call. = FALSE)
  }
  # A mismatch here would put one model's excess under another's label, which
  # is the kind of error that looks entirely plausible on a map.
  for (nm in c("excess", "attribution")) {
    got <- names(get(nm))
    if (!setequal(got, models)) {
      stop("`", nm, "` is named ", paste(got, collapse = ", "),
           " but `aug_models` is named ", paste(models, collapse = ", "),
           ". They must describe the same models.", call. = FALSE)
    }
  }

  # --- 1. geography and labels --------------------------------------------------

  geo <- sf::st_as_sf(smr_geo_full)
  .need_cols(geo, c("area", "total_obs", "total_exp", "total_smr",
                    "di_score_z", "gp_density_z", paste0(no2, "_z"),
                    "person_years"), "smr_geo_full")
  n_geo <- nrow(geo)

  geo[["label"]] <- .area_labels(geo[["area"]], nil_shp, pop_shp)

  # --- 2. the covariate alternative the modelling frame does not carry ----------
  #
  # add_covariate() takes one variable at a time, so the unweighted primary-care
  # indicator - the predecessor of the one in the model - never reaches
  # smr_geo_full. The app's first tab exists to compare the two.

  gp <- dplyr::select(
    gp_density_area,
    dplyr::any_of(c("area", "gp_density_unweighted", "gp_caseload",
                    "gp_weight_ratio"))
  )
  geo <- dplyr::left_join(geo, gp, by = "area")
  if (nrow(geo) != n_geo) {
    stop("Joining `gp_density_area` changed the row count from ", n_geo,
         " to ", nrow(geo), ": the indicator has duplicate `area` keys.",
         call. = FALSE)
  }
  if ("gp_density_unweighted" %in% names(geo)) {
    geo[["gp_density_unweighted_z"]] <-
      as.numeric(scale(geo[["gp_density_unweighted"]]))
  } else {
    warning("`gp_density_area` has no `gp_density_unweighted` column; the ",
            "right-hand map of the primary-care row will say so.",
            call. = FALSE)
  }

  # --- 3. crude rates and smoothed surfaces --------------------------------------

  crude <- tibble::as_tibble(sf::st_drop_geometry(crude_geo))
  .need_cols(crude, c("area", "total"), "crude_geo")

  aug_pair <- function(aug, prefix, what) {
    a <- sf::st_drop_geometry(aug)
    .need_cols(a, c("area", "bym2_rr", "bym2_exceed"), what)
    out <- tibble::tibble(area = as.character(a[["area"]]))
    out[[paste0(prefix, "_rr")]]     <- a[["bym2_rr"]]
    out[[paste0(prefix, "_exceed")]] <- a[["bym2_exceed"]]
    out
  }

  surfaces <- aug_pair(aug_total, "M3", "aug_total")
  for (m in models) {
    surfaces <- dplyr::left_join(
      surfaces, aug_pair(aug_models[[m]], m, paste0("aug_models$", m)),
      by = "area"
    )
  }

  mech_surface <- tibble::as_tibble(sf::st_drop_geometry(mech_bym2))
  mech_surface <- dplyr::select(
    mech_surface,
    dplyr::all_of("area"), dplyr::matches("^M_.*_bym2(_exc)?$")
  )
  mech_surface[["area"]] <- as.character(mech_surface[["area"]])

  # --- 4. excess and attribution, stacked with a model column --------------------

  excess_tbl <- dplyr::bind_rows(lapply(models, function(m) {
    pa <- excess[[m]][["per_area"]]
    .need_cols(pa, c("area", "observed", "expected", "predicted", "excess",
                     "excess_low", "excess_high", "rr_resid", "p_exceed",
                     "flagged"), paste0("excess$", m, "$per_area"))
    tibble::tibble(
      area = as.character(pa[["area"]]), model = m,
      observed = pa[["observed"]], expected = pa[["expected"]],
      predicted = pa[["predicted"]], excess = pa[["excess"]],
      excess_low = pa[["excess_low"]], excess_high = pa[["excess_high"]],
      rr_resid = pa[["rr_resid"]], p_exceed = pa[["p_exceed"]],
      flagged = pa[["flagged"]]
    )
  }))

  attribution_tbl <- dplyr::bind_rows(lapply(models, function(m) {
    at <- attribution[[m]]
    .need_cols(at, c("area", "component", "label", "deaths_yr", "pct",
                     "dominant"), paste0("attribution$", m))
    at[["area"]]  <- as.character(at[["area"]])
    at[["model"]] <- m
    at
  }))

  # --- 5. stroke sub-model --------------------------------------------------------

  allage <- sf::st_as_sf(smr_geo_allage)
  .need_cols(allage, c("area", "i63_obs", "i63_exp"), "smr_geo_allage")

  unadj <- sf::st_drop_geometry(aug_i63)
  adj   <- sf::st_drop_geometry(aug_i63_adj)
  .need_cols(unadj, c("bym2_rr", "bym2_exceed"), "aug_i63")
  .need_cols(adj,   c("bym2_rr", "bym2_exceed"), "aug_i63_adj")

  stroke <- tibble::tibble(
    area         = as.character(allage[["area"]]),
    i63_obs      = as.numeric(allage[["i63_obs"]]),
    i63_exp      = as.numeric(allage[["i63_exp"]]),
    unadj_rr     = unadj[["bym2_rr"]],
    unadj_exceed = unadj[["bym2_exceed"]],
    adj_rr       = adj[["bym2_rr"]],
    adj_exceed   = adj[["bym2_exceed"]]
  )

  # Points, so take the coordinates directly. st_centroid() on a POINT layer is
  # a no-op that warns about assuming planar geometry.
  centres <- sf::st_transform(sf::st_as_sf(stroke_centres), 4326)
  xy <- sf::st_coordinates(sf::st_geometry(centres))
  centre_tab <- tibble::tibble(
    facility = as.character(centres[["facility"]]),
    is_hub   = if ("is_hub" %in% names(centres))
      as.logical(centres[["is_hub"]]) else TRUE,
    lon = xy[, 1], lat = xy[, 2]
  )

  # --- 6. classification menus ----------------------------------------------------
  #
  # Built from the columns that exist, not a hard-coded list, so a change to the
  # lookup shows up in the app rather than producing a menu entry that points at
  # a column that is not there.

  tab <- sf::st_drop_geometry(geo)
  menu_for <- function(prefix) {
    stems <- sub("_smr$", "",
                 grep(sprintf("^%s_.*_smr$", prefix), names(tab), value = TRUE))
    if (!length(stems)) {
      return(tibble::tibble(stem = character(), label = character(),
                            deaths = numeric()))
    }
    d <- vapply(paste0(stems, "_obs"), function(v) {
      if (v %in% names(tab)) sum(tab[[v]], na.rm = TRUE) else NA_real_
    }, numeric(1), USE.NAMES = FALSE)
    m <- tibble::tibble(stem = stems, label = stratum_label(stems), deaths = d)
    m[order(-m[["deaths"]]), , drop = FALSE]
  }
  menus <- list(group = menu_for("G"), cause = menu_for("C"),
                mechanism = menu_for("M"))

  # --- 7. trim, simplify, write ----------------------------------------------------

  keep <- unique(c(
    "area", "label", "population", "person_years",
    "total_obs", "total_exp", "total_smr", "total_isr",
    "di_score", "di_score_z", "ivsm", "ivsm_z",
    "gp_density", "gp_density_z", "gp_caseload",
    "gp_density_unweighted", "gp_density_unweighted_z",
    paste0(no2, c("", "_z")), paste0(pm25, c("", "_z")),
    "t_hub_mean", "t_hub_mean_z", "t_centre_mean", "t_centre_mean_z",
    grep("^[CGM]_.*_(obs|exp|smr)$", names(geo), value = TRUE)
  ))
  keep <- intersect(keep, names(geo))

  app_data <- list(
    geo          = .simplify_for_web(geo[, keep], metres = simplify_m),
    crude        = crude,
    surfaces     = surfaces,
    mech_surface = mech_surface,
    excess       = excess_tbl,
    attribution  = attribution_tbl,
    stroke       = stroke,
    centres      = centre_tab,
    menus        = menus,
    meta = list(
      n_years          = n_years,
      rr_threshold     = rr_threshold,
      rr_threshold_str = rr_threshold_stroke,
      prob_cutoff      = prob_cutoff,
      models           = models,
      no2 = no2, pm25 = pm25,
      built_at = Sys.time()
    )
  )

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(app_data, path, compress = compress)
  invisible(path)
}
