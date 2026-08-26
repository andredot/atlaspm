# Stroke sub-model reporting helpers -------------------------------------------
#
# Everything the stroke report needs that is not already produced by
# R/stroke.R. Kept separate because these are presentation and diagnosis, not
# construction: they do not build the exposure, they interrogate it.
#
# The question they exist to answer is not "is there an association" but "could
# this design have detected one". A null result from an exposure with no
# contrast is not evidence of no effect, and the two are indistinguishable in a
# coefficient table.


#' Stroke centre table for visual verification
#'
#' Every centre with its coordinates, how it was geocoded, and a link that
#' opens the exact point on a map. Geocoding is the step in this pipeline most
#' likely to be silently wrong - a facility resolved to the wrong side of a
#' city changes the catchment of every area around it - and the only reliable
#' check is a human looking at the pin.
#'
#' The DGR node list gives a comune, not a street address, so the coordinate
#' and the map link are the address for verification purposes.
#'
#' @param centres `sf` from [geocode_stroke_centres()].
#' @param link Whether to include a map URL column. Default `TRUE`.
#' @param zoom Zoom level for the link.
#'
#' @return A tibble ordered hubs first, then by comune.
#' @examples
#' \dontrun{
#' stroke_centre_table(tar_read(stroke_centres))
#' }
#' @export
stroke_centre_table <- function(centres, link = TRUE, zoom = 17) {

  xy <- sf::st_coordinates(sf::st_as_sf(centres))
  tab <- sf::st_drop_geometry(centres)

  out <- tibble::tibble(
    centre_id = tab[["centre_id"]],
    facility  = tab[["facility"]],
    level     = ifelse(tab[["is_hub"]], "Hub (SU II, thrombectomy)",
                       "Spoke (SU I)"),
    ente      = if ("ente" %in% names(tab)) tab[["ente"]] else NA_character_,
    comune    = tab[["comune"]],
    provincia = tab[["provincia"]],
    lat       = round(xy[, 2], 5),
    lon       = round(xy[, 1], 5),
    source    = if ("geocode_source" %in% names(tab))
      tab[["geocode_source"]] else NA_character_
  )

  if (link) {
    out[["map"]] <- sprintf(
      "[check](https://www.openstreetmap.org/?mlat=%.5f&mlon=%.5f#map=%d/%.5f/%.5f)",
      out[["lat"]], out[["lon"]], zoom, out[["lat"]], out[["lon"]]
    )
  }

  out[order(!grepl("^Hub", out[["level"]]), out[["comune"]]), , drop = FALSE]
}


#' Does the exposure have enough contrast to detect anything?
#'
#' The first question to ask of a null result, and the one a coefficient table
#' cannot answer. If every area in the study sits within the therapeutic window,
#' there is no exposure gradient, and no amount of modelling will find one.
#'
#' Reports the spread of each accessibility measure, how many areas and how much
#' population fall beyond clinically meaningful thresholds, and how spatially
#' smooth the exposure is - the last because a smooth exposure competes with the
#' BYM2 random effect for the same variation, which is a separate reason a
#' coefficient can be driven toward zero.
#'
#' @param geo Modelling `sf` carrying the accessibility columns.
#' @param C Adjacency matrix, for the smoothness diagnostic. Optional.
#' @param thresholds Minute thresholds to count areas beyond.
#' @param vars Accessibility columns to summarise.
#'
#' @return A tibble, one row per measure.
#' @seealso [mde_per_unit()]
#' @export
exposure_contrast <- function(geo,
                              C = NULL,
                              thresholds = c(30, 45, 60),
                              vars = c("t_centre_mean", "t_hub_mean",
                                       "t_hub_p90")) {

  tab  <- sf::st_drop_geometry(geo)
  vars <- vars[vars %in% names(tab)]
  if (!length(vars)) {
    stop("No accessibility columns found in `geo`.", call. = FALSE)
  }
  pop <- if ("population" %in% names(tab)) tab[["population"]] else
    rep(1, nrow(tab))

  rows <- lapply(vars, function(v) {
    x <- as.numeric(tab[[v]])
    q <- stats::quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE,
                         names = FALSE)

    thr <- vapply(thresholds, function(t) sum(x > t, na.rm = TRUE), numeric(1))
    pthr <- vapply(thresholds, function(t) {
      100 * sum(pop[x > t], na.rm = TRUE) / sum(pop, na.rm = TRUE)
    }, numeric(1))

    moran <- if (!is.null(C)) {
      tryCatch(moran_test_raw(x, C, nsim = 999)$statistic,
               error = function(e) NA_real_)
    } else NA_real_

    out <- tibble::tibble(
      measure = v,
      min     = min(x, na.rm = TRUE),
      p05     = q[1], p25 = q[2], median = q[3], p75 = q[4], p95 = q[5],
      max     = max(x, na.rm = TRUE),
      iqr     = q[4] - q[2],
      sd      = stats::sd(x, na.rm = TRUE),
      # The ratio that matters: an SD of 3 minutes across a 20-minute range
      # means the "per SD" coefficient is describing a clinically trivial
      # contrast, however tight its interval looks.
      range_min = max(x, na.rm = TRUE) - min(x, na.rm = TRUE),
      moran_i   = moran
    )
    for (i in seq_along(thresholds)) {
      out[[paste0("n_over_", thresholds[i])]]   <- thr[i]
      out[[paste0("pct_pop_over_", thresholds[i])]] <- pthr[i]
    }
    out
  })

  do.call(rbind, rows)
}


#' Re-express a standardised coefficient in clinical units
#'
#' A coefficient "per standard deviation" is uninterpretable when the reader
#' needs to know what a ten-minute delay does. Worse, it disguises a
#' range-restriction problem: an SD of two minutes makes a null look like
#' strong evidence of no effect, when it is really evidence about a contrast
#' nobody cares about.
#'
#' Converts the posterior to a relative risk per `per` units of the raw
#' covariate, and reports the largest effect the data can rule out - the upper
#' credible bound - which is the honest way to present a null.
#'
#' @param fit A fitted model.
#' @param geo The `sf` it was fitted to.
#' @param var Raw (unstandardised) covariate column, e.g. `"t_hub_mean"`.
#' @param z_var The standardised column actually in the model. Derived by
#'   appending `_z` when `NULL`.
#' @param per Units of `var` to express the effect per. Default `10`.
#' @param probs Credible-interval bounds.
#'
#' @return A one-row tibble.
#' @examples
#' \dontrun{
#' mde_per_unit(control_fits$tracer, smr_geo_tracer, "t_hub_mean", per = 10)
#' }
#' @export
mde_per_unit <- function(fit, geo, var, z_var = NULL, per = 10,
                         probs = c(0.025, 0.975)) {

  if (is.null(z_var)) z_var <- paste0(var, "_z")
  require_cols(geo, c(var, z_var), "geo")

  sd_raw <- stats::sd(sf::st_drop_geometry(geo)[[var]], na.rm = TRUE)
  draws  <- as.numeric(as.matrix(fit$stanfit, pars = .beta_pars(fit, z_var)))

  # beta is per 1 SD of the raw variable; convert to per `per` raw units.
  b <- draws * (per / sd_raw)
  q <- stats::quantile(b, probs = probs, names = FALSE)

  tibble::tibble(
    variable   = var,
    per        = per,
    sd_raw     = sd_raw,
    rr         = exp(mean(b)),
    ci_low     = exp(q[1]),
    ci_high    = exp(q[2]),
    # The largest effect compatible with the data, in either direction.
    rules_out_above = exp(max(abs(q))),
    p_direction = max(mean(b > 0), mean(b < 0))
  )
}


#' Map the stroke network over the study area
#'
#' @param area_shp Modelling geography.
#' @param centres `sf` of stroke centres.
#' @param caption Passed to `labs()`.
#' @return A ggplot.
#' @export
plot_stroke_centres <- function(area_shp, centres, caption = NULL) {

  a <- sf::st_as_sf(area_shp)
  p <- sf::st_as_sf(centres) |> sf::st_transform(sf::st_crs(a))
  p[["Level"]] <- ifelse(p[["is_hub"]], "Hub (thrombectomy)", "Spoke")

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = a, fill = "grey94", colour = "white",
                     linewidth = 0.2) +
    ggplot2::geom_sf(data = p,
                     ggplot2::aes(shape = .data[["Level"]],
                                  fill = .data[["Level"]]),
                     size = 2.6, colour = "grey15", stroke = 0.4) +
    ggplot2::scale_shape_manual(values = c("Hub (thrombectomy)" = 24,
                                           "Spoke" = 21)) +
    ggplot2::scale_fill_manual(values = c("Hub (thrombectomy)" = "#8c2d04",
                                          "Spoke" = "white")) +
    ggplot2::labs(caption = caption) +
    ggplot2::coord_sf(
      xlim   = c(sf::st_bbox(a)[["xmin"]], sf::st_bbox(a)[["xmax"]]),
      ylim   = c(sf::st_bbox(a)[["ymin"]], sf::st_bbox(a)[["ymax"]]),
      expand = TRUE
    ) +
    theme_atlas(map = TRUE)
}


#' Map a travel-time surface
#'
#' @param geo Modelling `sf`.
#' @param value Column to map.
#' @param centres Optional `sf` of centres to overlay.
#' @param legend_title,caption Labels.
#' @return A ggplot.
#' @export
plot_travel_time <- function(geo, value = "t_hub_mean", centres = NULL,
                             legend_title = "Minutes", caption = NULL) {

  g <- sf::st_as_sf(geo)
  require_cols(g, value, "geo")

  p <- ggplot2::ggplot(g) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[value]]), colour = NA) +
    ggplot2::scale_fill_viridis_c(option = "cividis", direction = -1,
                                  name = legend_title)

  if (!is.null(centres)) {
    cen <- sf::st_as_sf(centres) |> sf::st_transform(sf::st_crs(g))
    p <- p +
      ggplot2::geom_sf(data = cen[cen[["is_hub"]], ], shape = 24,
                       fill = "#d94801", colour = "grey10", size = 2.2)
  }

  p + ggplot2::labs(caption = caption) + theme_atlas(map = TRUE)
}


#' Scatter of the tracer outcome against the exposure, before any model
#'
#' The plot to look at before believing a coefficient. If there is no
#' relationship visible here, a null coefficient is unsurprising; if there is
#' one here but not in the model, the spatial term has absorbed it, and that is
#' a different problem with a different remedy.
#'
#' @param geo Modelling `sf` with the observed/expected pair attached.
#' @param x Exposure column.
#' @param obs_col,exp_col Observed and expected counts. Default to the
#'   cerebrovascular 0-74 tracer; pass `"i63_obs"` / `"i63_exp"` for the all-age
#'   cerebral-infarction outcome.
#' @param ylab Y-axis label. Derived from `obs_col` when `NULL`.
#' @param caption Passed to `labs()`.
#' @return A ggplot.
#' @export
plot_tracer_raw <- function(geo, x = "t_hub_mean",
                            obs_col = "cvd_obs", exp_col = "cvd_exp",
                            ylab = NULL, caption = NULL) {

  tab <- sf::st_drop_geometry(geo)
  require_cols(tab, c(x, obs_col, exp_col), "geo")

  if (is.null(ylab)) {
    ylab <- switch(obs_col,
                   cvd_obs    = "Cerebrovascular SMR (observed / expected)",
                   i63_obs    = "I63 SMR, all ages (observed / expected)",
                   haem_obs   = "Haemorrhagic stroke SMR (observed / expected)",
                   cvdall_obs = "All-cerebrovascular SMR, all ages",
                   "SMR (observed / expected)")
  }

  d <- data.frame(
    x   = as.numeric(tab[[x]]),
    smr = tab[[obs_col]] / tab[[exp_col]],
    w   = tab[[exp_col]]
  )
  rho <- stats::cor(d$x, d$smr, method = "spearman", use = "complete.obs")

  ggplot2::ggplot(d, ggplot2::aes(x = .data[["x"]], y = .data[["smr"]])) +
    ggplot2::geom_hline(yintercept = 1, linetype = "22", colour = "grey45") +
    ggplot2::geom_point(ggplot2::aes(size = .data[["w"]]), alpha = 0.5,
                        colour = "#8c2d04") +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                         colour = "grey25", linewidth = 0.6) +
    ggplot2::scale_size_area(name = "Expected deaths", max_size = 5) +
    ggplot2::labs(
      x = "Population-weighted travel time (minutes)",
      y = ylab,
      subtitle = sprintf("Spearman rank correlation = %.3f", rho),
      caption = caption
    ) +
    theme_atlas()
}


# All-age outcomes from the raw register ---------------------------------------
#
# The avoidable-mortality frame stops at 75 because that is what the
# OECD/Eurostat definition specifies. That restriction is right for the study's
# primary outcome and wrong for a tracer: an ischaemic stroke at 90 is as
# treatable as one at 60, and excluding it discards most of the events the
# exposure could plausibly act on. These functions go back to the register.


#' Build an indirectly standardised outcome for an arbitrary ICD set
#'
#' Bypasses the avoidable-cause lookup entirely and standardises a chosen set of
#' ICD-10 prefixes over any age range, against the same population denominator
#' the rest of the study uses.
#'
#' The standard is internal, as elsewhere: age-sex specific rates are computed
#' across the whole study territory and applied to each area's own age-sex
#' structure. The expected counts are therefore on the same footing as
#' `total_exp`, and the resulting SMRs are comparable with the main analysis
#' even though the age range differs.
#'
#' @param mort_raw Output of [import_mortality()] - the unfiltered register,
#'   before [preprocess_mortality()] applies the age and cause restrictions.
#' @param pop_area_table Population by area, age, sex and year.
#' @param areas Character vector of study areas to retain.
#' @param prefixes ICD-10 prefixes to match, e.g. `"I63"` or `c("I60","I61")`.
#'   Matched against the leading characters of the normalised cause code, so
#'   `"I63"` captures I63.0 through I63.9.
#' @param age_min,age_max Age range, inclusive. Defaults to all ages.
#' @param pop_year Years to sum for person-years.
#' @param label Stem for the output columns.
#'
#' @return A tibble: `area`, `<label>_obs`, `<label>_exp`, `<label>_smr`.
#'   Attributes record the definition used.
#'
#' @examples
#' \dontrun{
#' i63 <- build_icd_outcome(mort_raw, pop_area_table, area_shp$area,
#'                          prefixes = "I63", label = "i63")
#' }
#' @seealso [outcome_feasibility()], [preprocess_smr()]
#' @export
build_icd_outcome <- function(mort_raw, pop_area_table, areas,
                              prefixes,
                              age_min  = 0,
                              age_max  = Inf,
                              pop_year = 2022:2024,
                              label    = "outcome") {

  require_cols(mort_raw, c("causa", "eta", "sesso", "area_residenza"),
               "mort_raw")
  require_cols(pop_area_table, c("area", "Eta", "sesso", "anno", "numero"),
               "pop_area_table")

  # Same normalisation as preprocess_mortality(): strip punctuation, upcase,
  # then match on the leading characters. Doing it differently here would let
  # the tracer and the main outcome disagree about what an ICD code is.
  code <- toupper(gsub("[^A-Za-z0-9]", "", as.character(mort_raw[["causa"]])))
  hit  <- Reduce(`|`, lapply(prefixes, function(p) startsWith(code, p)))

  d <- mort_raw[hit, , drop = FALSE]
  d <- d[!is.na(d[["eta"]]) &
           d[["eta"]] >= age_min & d[["eta"]] <= age_max, , drop = FALSE]
  d <- d[d[["area_residenza"]] %in% areas, , drop = FALSE]

  if (!nrow(d)) {
    stop("No deaths matched prefix(es) ", paste(prefixes, collapse = ", "),
         " in ages ", age_min, "-", age_max, ".", call. = FALSE)
  }

  deaths <- d |>
    dplyr::mutate(.age = as.integer(.data[["eta"]]),
                  .sex = as.integer(.data[["sesso"]]),
                  .area = as.character(.data[["area_residenza"]])) |>
    dplyr::count(.data[[".area"]], .data[[".age"]], .data[[".sex"]],
                 name = "obs")
  names(deaths)[1:3] <- c(".area", ".age", ".sex")

  pop <- check_pop_years(pop_area_table, pop_year) |>
    dplyr::filter(.data[["anno"]] %in% as.character(pop_year),
                  .data[["area"]] %in% areas) |>
    dplyr::mutate(.age = as.integer(.data[["Eta"]]),
                  .sex = as.integer(.data[["sesso"]]),
                  .area = as.character(.data[["area"]]),
                  .pop = as.numeric(.data[["numero"]])) |>
    dplyr::filter(.data[[".age"]] >= age_min, .data[[".age"]] <= age_max) |>
    dplyr::group_by(.data[[".area"]], .data[[".age"]], .data[[".sex"]]) |>
    dplyr::summarise(.pop = sum(.data[[".pop"]], na.rm = TRUE),
                     .groups = "drop")
  names(pop)[1:3] <- c(".area", ".age", ".sex")

  assert_sex_alignment(deaths[[".sex"]], pop[[".sex"]])

  # Internal standard: age-sex specific rates across the whole territory.
  std <- deaths |>
    dplyr::group_by(.data[[".age"]], .data[[".sex"]]) |>
    dplyr::summarise(std_deaths = sum(.data[["obs"]]), .groups = "drop") |>
    dplyr::full_join(
      pop |>
        dplyr::group_by(.data[[".age"]], .data[[".sex"]]) |>
        dplyr::summarise(std_pop = sum(.data[[".pop"]]), .groups = "drop"),
      by = c(".age", ".sex")
    ) |>
    dplyr::mutate(
      std_deaths = dplyr::coalesce(.data[["std_deaths"]], 0),
      std_rate   = dplyr::if_else(.data[["std_pop"]] > 0,
                                  .data[["std_deaths"]] / .data[["std_pop"]], 0)
    )

  expected <- pop |>
    dplyr::left_join(std, by = c(".age", ".sex")) |>
    dplyr::group_by(.data[[".area"]]) |>
    dplyr::summarise(exp = sum(.data[[".pop"]] * .data[["std_rate"]]),
                     .groups = "drop")

  observed <- deaths |>
    dplyr::group_by(.data[[".area"]]) |>
    dplyr::summarise(obs = sum(.data[["obs"]]), .groups = "drop")

  out <- tibble::tibble(area = areas) |>
    dplyr::left_join(observed, by = c("area" = ".area")) |>
    dplyr::left_join(expected, by = c("area" = ".area")) |>
    dplyr::mutate(
      obs = dplyr::coalesce(.data[["obs"]], 0),
      exp = dplyr::coalesce(.data[["exp"]], NA_real_),
      smr = .data[["obs"]] / .data[["exp"]]
    )

  names(out) <- c("area", paste0(label, c("_obs", "_exp", "_smr")))

  attr(out, "definition") <- list(
    prefixes = prefixes, age_min = age_min, age_max = age_max,
    pop_year = pop_year, n_deaths = sum(observed[["obs"]]),
    n_areas = length(areas)
  )
  message(sprintf(
    "%s: %d deaths, ICD %s, ages %d-%s, across %d areas.",
    label, sum(observed[["obs"]]), paste(prefixes, collapse = "/"),
    age_min, if (is.finite(age_max)) as.character(age_max) else "max",
    length(areas)))

  out
}


#' Is an outcome dense enough to model at this resolution?
#'
#' Answers the question that should precede fitting, not follow it. A BYM2 on a
#' surface that is zero across most units estimates the prior rather than the
#' data, and produces a smooth, plausible, uninformative map.
#'
#' The thresholds are conventional rather than exact: an expected count below
#' about 5 is where a Poisson SMR stops behaving, and a majority of zero
#' observations is where a spatial smoother has nothing left to smooth.
#'
#' @param outcome Output of [build_icd_outcome()].
#' @param label The column stem used.
#'
#' @return A one-row tibble with a `verdict` column.
#' @export
outcome_feasibility <- function(outcome, label) {

  obs <- outcome[[paste0(label, "_obs")]]
  exp <- outcome[[paste0(label, "_exp")]]
  n   <- length(obs)

  pct_zero <- 100 * sum(obs == 0, na.rm = TRUE) / n
  pct_exp5 <- 100 * sum(exp < 5, na.rm = TRUE) / n

  verdict <- if (pct_zero > 50) {
    "NOT MODELLABLE - most areas have no events; a spatial model would return the prior"
  } else if (pct_zero > 30 || pct_exp5 > 75) {
    "MARGINAL - fit it, but expect wide intervals and read the exceedance map with caution"
  } else {
    "MODELLABLE"
  }

  tibble::tibble(
    outcome        = label,
    deaths         = sum(obs, na.rm = TRUE),
    areas          = n,
    mean_per_area  = mean(obs, na.rm = TRUE),
    median_per_area = stats::median(obs, na.rm = TRUE),
    max_per_area   = max(obs, na.rm = TRUE),
    n_zero         = sum(obs == 0, na.rm = TRUE),
    pct_zero       = pct_zero,
    median_expected = stats::median(exp, na.rm = TRUE),
    pct_expected_under5 = pct_exp5,
    verdict        = verdict
  )
}


#' Overlay the stroke network on an existing map
#'
#' Adds centres to a map from [plot_smr_map()] or [plot_exceedance_map()],
#' with a legend, without disturbing the fill scale.
#'
#' Two things this handles that a bare `geom_sf()` does not.
#'
#' \strong{Extent.} Routing deliberately includes centres well outside the ATS
#' Milano boundary, because residents near the edge may be closer to them. Drawn
#' naively, those points expand the plot to cover Varese, Bergamo and Pavia and
#' the study area shrinks into a corner. The map is cropped back to the
#' choropleth's own bounding box, so out-of-area centres inform the routing but
#' do not dictate the frame.
#'
#' \strong{Legend.} Marker type is mapped to the `shape` aesthetic rather than
#' set as a fixed parameter, which is what produces a legend entry. `fill` stays
#' a fixed parameter because the panel's discrete SMR or exceedance palette
#' already owns that scale.
#'
#' @param p A ggplot from one of the map functions.
#' @param centres `sf` of stroke centres.
#' @param hubs_only Show only level-II hubs. Default `FALSE`, so the spoke
#'   network is visible too - the distinction between the two carries the
#'   argument throughout this analysis, and a map that hides it invites the
#'   reader to assume every marker offers thrombectomy.
#' @param crop Clip the map to the underlying layer's bounding box. Default
#'   `TRUE`.
#' @param size,fill,colour Marker appearance. White with a dark outline reads
#'   against every tier of both palettes.
#' @param legend_title Title for the shape legend. `NULL` for none.
#'
#' @return The ggplot, with the centre layer and its legend added.
#' @examples
#' \dontrun{
#' plot_smr_map(aug_cvd, value = "bym2_rr") |> overlay_hubs(stroke_centres)
#' }
#' @export
overlay_hubs <- function(p, centres, hubs_only = FALSE, crop = TRUE,
                         size = 2.4, fill = "white", colour = "grey10",
                         legend_title = NULL) {

  cen <- sf::st_as_sf(centres)
  if (hubs_only && "is_hub" %in% names(cen)) cen <- cen[cen[["is_hub"]], ]

  cen[[".marker"]] <- if ("is_hub" %in% names(cen)) {
    factor(ifelse(cen[["is_hub"]],
                  "Thrombectomy hub (SU II)", "Stroke unit (SU I)"),
           levels = c("Thrombectomy hub (SU II)", "Stroke unit (SU I)"))
  } else {
    factor("Stroke centre")
  }

  out <- p + ggplot2::geom_sf(
    data = cen, inherit.aes = FALSE,
    ggplot2::aes(shape = .data[[".marker"]]),
    size = size, fill = fill, colour = colour, stroke = 0.5
  ) +
    ggplot2::scale_shape_manual(
      values = c("Thrombectomy hub (SU II)" = 24,   # filled triangle
                 "Stroke unit (SU I)"       = 21,   # filled circle
                 "Stroke centre"            = 21),
      name   = legend_title,
      drop   = TRUE
    ) +
    ggplot2::guides(shape = ggplot2::guide_legend(order = 2))

  # Crop to the choropleth, not to the union of choropleth and centres.
  if (crop && !is.null(p$data)) {
    bb <- tryCatch(sf::st_bbox(sf::st_as_sf(p$data)), error = function(e) NULL)
    if (!is.null(bb)) {
      out <- out + ggplot2::coord_sf(
        xlim   = c(bb[["xmin"]], bb[["xmax"]]),
        ylim   = c(bb[["ymin"]], bb[["ymax"]]),
        expand = TRUE
      )
    }
  }
  out
}


#' Two maps side by side, sharing one legend
#'
#' For raw-versus-smoothed comparisons, where the whole point is that the two
#' surfaces are on the same scale. Collecting the guide is not cosmetic: two
#' separate legends invite the reader to compare colours that happen to look
#' alike but belong to different scales, which is the specific mistake this
#' pairing exists to prevent.
#'
#' @param p1,p2 ggplots, normally both from [plot_smr_map()].
#' @param titles Optional length-2 character vector of panel titles.
#' @param position Where to put the collected legend.
#'
#' @return A patchwork object, or `p1` with a message if patchwork is absent.
#' @examples
#' \dontrun{
#' map_pair(plot_smr_map(g, "cvd_smr"), plot_smr_map(aug_cvd, "bym2_rr"),
#'          titles = c("Raw", "Smoothed"))
#' }
#' @export
map_pair <- function(p1, p2, titles = NULL, position = "bottom") {

  if (!is.null(titles)) {
    p1 <- p1 + ggplot2::ggtitle(titles[1])
    p2 <- p2 + ggplot2::ggtitle(titles[2])
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("patchwork is not installed; returning the first panel only. ",
            "install.packages('patchwork') to see the pair.")
    return(p1)
  }

  styled <- ggplot2::theme(
    plot.title      = ggplot2::element_text(face = "bold", hjust = 0.5,
                                            size = ggplot2::rel(0.95)),
    legend.position = position
  )

  patchwork::wrap_plots(p1, p2, nrow = 1) +
    patchwork::plot_layout(guides = "collect") &
    styled
}


#' Deaths attributable to an exposure, under the fitted model
#'
#' Answers "how many of these deaths does the model attribute to travel time?"
#' rather than "is the coefficient different from zero". The two are not the
#' same question, and for a policy audience the first is the useful one.
#'
#' The counterfactual holds everything else fixed - including the spatial
#' random effect - and moves only the exposure, to `reference`. The difference
#' between fitted and counterfactual counts, summed over areas and annualised,
#' is the attributable excess.
#'
#' \strong{Read the interval, not the point estimate.} When the coefficient is
#' null the point estimate will be near zero, and that is not the finding. The
#' finding is the upper bound: the largest number of deaths per year the data
#' allow to be attributed to access. A tight bound around zero is a
#' substantive, quantified negative result; a wide one means the design cannot
#' say, and should be reported that way.
#'
#' \strong{This is model-based attribution, not a causal effect.} It inherits
#' every identification problem of the fit it is given - most importantly, that
#' hub placement is endogenous to population need, so an observational contrast
#' between near and far areas is confounded by the reasons hubs were sited
#' where they are.
#'
#' @param fit A fitted model containing `var_z`.
#' @param geo The `sf` it was fitted to.
#' @param var Raw exposure column, e.g. `"t_hub_mean"`.
#' @param z_var Standardised column in the model. Derived when `NULL`.
#' @param reference Counterfactual exposure value. `"min"` uses the best-served
#'   area (what would be avoided if everyone had that access); a number uses
#'   that value, e.g. `45` for the regional centralisation threshold.
#' @param obs_col Observed-count column, for the percentage denominator.
#' @param n_years Years the counts span.
#' @param probs Credible-interval bounds.
#'
#' @return A one-row tibble: attributable deaths per year with interval, as a
#'   count and a share of the outcome.
#' @examples
#' \dontrun{
#' attributable_excess(control_fits$tracer, smr_geo_tracer, "t_hub_mean",
#'                     obs_col = "cvd_obs")
#' }
#' @seealso [mde_per_unit()], [residual_excess()]
#' @export
attributable_excess <- function(fit, geo, var, z_var = NULL,
                                reference = "min",
                                obs_col   = "cvd_obs",
                                n_years   = 3,
                                probs     = c(0.025, 0.975)) {

  if (is.null(z_var)) z_var <- paste0(var, "_z")
  tab <- sf::st_drop_geometry(geo)
  require_cols(tab, c(var, z_var, obs_col), "geo")

  x      <- as.numeric(tab[[var]])
  sd_raw <- stats::sd(x, na.rm = TRUE)
  ref    <- if (identical(reference, "min")) min(x, na.rm = TRUE) else
    as.numeric(reference)

  # Fitted counts, and draws of the coefficient.
  mu   <- rstan::extract(fit$stanfit, pars = "fitted")$fitted   # [draws, areas]
  beta <- as.numeric(as.matrix(fit$stanfit, pars = .beta_pars(fit, z_var)))

  # Counterfactual: same random effect, same intercept, exposure moved to the
  # reference. On the log scale that is a shift of beta * (x_ref - x_i) / sd.
  shift  <- outer(beta, (ref - x) / sd_raw)     # [draws, areas]
  mu_cf  <- mu * exp(shift)

  attr_draws <- rowSums(mu - mu_cf) / n_years
  q <- stats::quantile(attr_draws, probs = probs, names = FALSE)

  total_yr <- sum(round_half_up(tab[[obs_col]])) / n_years

  tibble::tibble(
    variable       = var,
    reference      = ref,
    exposure_max   = max(x, na.rm = TRUE),
    deaths_per_year = mean(attr_draws),
    ci_low         = q[1],
    ci_high        = q[2],
    pct_of_outcome = 100 * mean(attr_draws) / total_yr,
    pct_ci_high    = 100 * q[2] / total_yr,
    upper_bound    = max(abs(q)),
    total_per_year = total_yr
  )
}


# Reporting summaries for the stroke sub-model ---------------------------------

#' What the BYM2 smoothing did to an outcome surface
#'
#' The raw ratio and the smoothed relative risk are the same quantity seen
#' before and after shrinkage, and the gap between their ranges is the honest
#' statement of how much of the raw map was noise. Reporting the two ranges
#' side by side with `rho` stops a reader taking a dramatic raw map at face
#' value: with a median of two or three expected deaths per area, one extra
#' death moves a small unit's ratio by a large fraction.
#'
#' @param geo The `sf` carrying the observed and expected counts.
#' @param aug Output of [augment_bym2()] for the same outcome.
#' @param diag Output of [check_bym2_fit()] for the same fit.
#' @param obs_col,exp_col Observed and expected count columns in `geo`.
#' @param rr_col Smoothed relative-risk column in `aug`.
#'
#' @return A two-column tibble, `Quantity` and `Value`, ready for `kable()`.
#' @examples
#' \dontrun{
#' smoothing_summary(smr_geo_allage, aug_i63, diag_i63_bym2,
#'                   obs_col = "i63_obs", exp_col = "i63_exp")
#' }
#' @seealso [exceedance_count()]
#' @export
smoothing_summary <- function(geo, aug, diag,
                              obs_col = "cvd_obs", exp_col = "cvd_exp",
                              rr_col  = "bym2_rr") {

  tab <- sf::st_drop_geometry(geo)
  require_cols(tab, c(obs_col, exp_col), "geo")
  aug_tab <- sf::st_drop_geometry(aug)
  require_cols(aug_tab, rr_col, "aug")

  raw <- tab[[obs_col]] / tab[[exp_col]]
  rr  <- aug_tab[[rr_col]]

  rng2 <- function(x) sprintf("%.2f to %.2f",
                              min(x, na.rm = TRUE), max(x, na.rm = TRUE))

  tibble::tibble(
    Quantity = c(
      "Median expected deaths per area",
      "Mixing parameter rho (share of residual variation that is spatially structured)",
      "Range of raw SMR",
      "Range of smoothed relative risk",
      "Max R-hat",
      "Min effective sample size"),
    Value = c(
      sprintf("%.1f", stats::median(tab[[exp_col]], na.rm = TRUE)),
      sprintf("%.3f", diag[["rho_mean"]]),
      rng2(raw),
      rng2(rr),
      sprintf("%.4f", diag[["rhat_max"]]),
      format(round(diag[["neff_min"]]), big.mark = ","))
  )
}


#' Count areas above posterior exceedance cut-offs
#'
#' @param aug Output of [augment_bym2()].
#' @param col Exceedance-probability column.
#' @param cutoffs Posterior probabilities to count above.
#'
#' @return A one-row tibble: `n_areas`, then one `p80`-style column per cut-off,
#'   with the relative-risk threshold carried in the `"threshold"` attribute
#'   when `aug` records one.
#' @examples
#' \dontrun{
#' exceedance_count(aug_i63)
#' }
#' @seealso [smoothing_summary()]
#' @export
exceedance_count <- function(aug, col = "bym2_exceed",
                             cutoffs = c(0.80, 0.95)) {

  tab <- sf::st_drop_geometry(aug)
  require_cols(tab, col, "aug")
  x <- tab[[col]]

  out <- tibble::tibble(n_areas = sum(!is.na(x)))
  for (cut in cutoffs) {
    out[[paste0("p", round(100 * cut))]] <- sum(x > cut, na.rm = TRUE)
  }
  attr(out, "threshold") <- attr(aug, "threshold")
  out
}


#' Cerebrovascular deaths by subtype and the care channel each belongs to
#'
#' The point of the table is not epidemiological description but outcome
#' choice. Thrombectomy acts on cerebral infarction alone; neurosurgical and
#' neurocritical care act on the haemorrhagic subtypes. An outcome that pools
#' them dilutes any access effect toward zero by construction, and the size of
#' that dilution is what this table makes visible.
#'
#' @param mort_raw Output of [import_mortality()], unfiltered.
#' @param age_max Upper age, inclusive. `NULL` (default) keeps all ages.
#'
#' @return A tibble: `subtype`, `deaths`, `pct`, ordered by `deaths`.
#' @examples
#' \dontrun{
#' stroke_subtype_table(mort_raw)
#' }
#' @export
stroke_subtype_table <- function(mort_raw, age_max = NULL) {

  require_cols(mort_raw, c("causa", "eta"), "mort_raw")

  d <- mort_raw
  d[["key"]] <- toupper(gsub("[^A-Za-z0-9]", "", as.character(d[["causa"]])))
  d <- d[startsWith(d[["key"]], "I6"), , drop = FALSE]
  if (!is.null(age_max)) {
    d <- d[!is.na(d[["eta"]]) & d[["eta"]] <= age_max, , drop = FALSE]
  }

  k3 <- substr(d[["key"]], 1, 3)
  subtype <- dplyr::case_when(
    k3 == "I63"                     ~ "I63 Cerebral infarction (thrombectomy)",
    k3 %in% c("I60", "I61", "I62")  ~ "I60-I62 Haemorrhagic (neurosurgery)",
    k3 == "I64"                     ~ "I64 Stroke, type not specified",
    k3 == "I69"                     ~ "I69 Sequelae",
    TRUE                            ~ paste(k3, "Other cerebrovascular")
  )

  tibble::tibble(subtype = subtype) |>
    dplyr::count(subtype, name = "deaths") |>
    dplyr::mutate(pct = 100 * .data[["deaths"]] / sum(.data[["deaths"]])) |>
    dplyr::arrange(dplyr::desc(.data[["deaths"]]))
}
