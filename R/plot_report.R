# Figures for the thesis -------------------------------------------------------
#
# Every function here returns a ggplot with no title baked in. Titles live in
# the report's figure captions, because a title rendered inside the image cannot
# be styled, cannot be cross-referenced, and gets duplicated the moment someone
# adds a caption in Word. Subtitles and captions carrying methodological
# information (threshold, period, source) do stay in the plot, since those must
# travel with the image if it is ever pulled out of the document.


#' Shared theme for the report figures
#'
#' One theme, applied everywhere, so the figures look like they belong to the
#' same document. Base size 11 to survive being scaled down to a single column.
#'
#' @param base_size Base font size.
#' @param map Whether this is a choropleth (drops axes and grid entirely).
#' @return A ggplot2 theme.
#' @export
theme_atlas <- function(base_size = 11, map = FALSE) {
  base <- if (map) {
    ggplot2::theme_void(base_size = base_size)
  } else {
    ggplot2::theme_minimal(base_size = base_size)
  }

  base + ggplot2::theme(
    plot.subtitle   = ggplot2::element_text(colour = "grey30",
                                            margin = ggplot2::margin(b = 8)),
    plot.caption    = ggplot2::element_text(colour = "grey40", hjust = 0,
                                            size = ggplot2::rel(0.8),
                                            margin = ggplot2::margin(t = 8)),
    legend.position = if (map) "right" else "bottom",
    legend.title    = ggplot2::element_text(size = ggplot2::rel(0.9)),
    strip.text      = ggplot2::element_text(face = "bold",
                                            size = ggplot2::rel(0.9)),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin      = ggplot2::margin(8, 8, 8, 8)
  )
}


#' Forest plot of relative risks
#'
#' The figure Results "Covariate association" is missing, and the one the
#' control panel needs. Log x-axis, because a relative risk of 0.5 and one of
#' 2.0 are the same distance from the null and a linear axis would hide that.
#'
#' @param coefs A tibble from [collect_coefficients()] or [collect_controls()].
#' @param facet Optional faceting column, e.g. `"model"` or `"role"`.
#' @param y Column to use for the row labels. Default `"label"`.
#' @param xlab Axis label.
#' @param subtitle,caption Passed through to `labs()`.
#' @param null_line Where to draw the reference line. Default `1`.
#'
#' @return A ggplot.
#' @examples
#' \dontrun{
#' plot_forest(collect_coefficients(fits), facet = "model")
#' }
#' @export
plot_forest <- function(coefs,
                        facet     = NULL,
                        y         = "label",
                        xlab      = "Relative risk per standard deviation (95% credible interval)",
                        subtitle  = NULL,
                        caption   = NULL,
                        null_line = 1) {

  require_cols(coefs, c(y, "estimate", "ci_low", "ci_high"), "coefs")

  # Validate at BUILD time. A ggplot defers evaluation until it is printed, so
  # an empty or facet-less input produces a target that stores "successfully"
  # and only fails later, inside the report, as "Faceting variables must have
  # at least one value" - a message about ggplot rather than about the empty
  # coefficient table that caused it.
  if (!nrow(coefs)) {
    stop("`coefs` has no rows, so there is nothing to plot. The coefficient ",
         "table upstream came back empty; fix that rather than this figure.",
         call. = FALSE)
  }
  if (!is.null(facet)) {
    require_cols(coefs, facet, "coefs")
    if (!length(stats::na.omit(unique(coefs[[facet]])))) {
      stop("Faceting column `", facet, "` is present but entirely NA.",
           call. = FALSE)
    }
  }

  coefs[[".row"]] <- factor(coefs[[y]], levels = rev(unique(coefs[[y]])))
  coefs[[".sig"]] <- !coefs[["crosses_null"]] %in% TRUE

  p <- ggplot2::ggplot(
    coefs,
    ggplot2::aes(x = .data[["estimate"]], y = .data[[".row"]])
  ) +
    ggplot2::geom_vline(xintercept = null_line, linetype = "22",
                        colour = "grey45") +
    # geom_linerange(), not geom_errorbarh(): the latter is deprecated as of
    # ggplot2 3.5 and the caps were suppressed with height = 0 anyway.
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data[["ci_low"]], xmax = .data[["ci_high"]]),
      linewidth = 0.6, colour = "grey25"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data[[".sig"]]),
      shape = 21, size = 2.8, colour = "grey15", stroke = 0.5
    ) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#8c2d04", `FALSE` = "white"),
      guide  = "none"
    ) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::labs(x = xlab, y = NULL, subtitle = subtitle, caption = caption) +
    theme_atlas()

  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)),
                                 ncol = 1, scales = "free_y")
  }
  p
}


#' Choropleth of the crude mortality rate
#'
#' Results compares the crude and standardised surfaces; the pipeline mapped
#' only the standardised one. Deliberately a continuous scale rather than the
#' binned SMR scale: the point of showing the crude map is that it is unstable
#' and driven by small denominators, and binning would tidy that away.
#'
#' @param geo An `sf` carrying the rate column.
#' @param value Column to map. Default `"total"` (the all-cause crude rate from
#'   [preprocess_cmr()]).
#' @param subtitle,caption Passed to `labs()`.
#' @param legend_title Legend title.
#'
#' @return A ggplot.
#' @export
plot_cmr_map <- function(geo,
                         value        = "total",
                         subtitle     = NULL,
                         caption      = NULL,
                         legend_title = "Deaths per 100,000\nperson-years") {

  require_cols(geo, value, "geo")

  ggplot2::ggplot(sf::st_as_sf(geo)) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[value]]), colour = NA) +
    ggplot2::scale_fill_viridis_c(option = "magma", direction = -1,
                                  name = legend_title) +
    ggplot2::labs(subtitle = subtitle, caption = caption) +
    theme_atlas(map = TRUE)
}


#' Map of residual excess deaths
#'
#' The primary output for practice: not the relative risk, but the count of
#' deaths beyond what the covariates predict. Counts rather than ratios because
#' a doubled risk in an area of 4,000 people and a 10% excess in an area of
#' 90,000 mean very different things to a service planner, and the ratio map
#' hides that.
#'
#' A diverging scale centred on zero, because the sign matters here in a way it
#' does not on a relative-risk map.
#'
#' @param excess Output of [residual_excess()].
#' @param geo The `sf` the model was fitted to.
#' @param flagged_only Shade only the areas above the probability cutoff,
#'   leaving the rest grey. Default `FALSE`.
#' @param subtitle,caption Passed to `labs()`.
#'
#' @return A ggplot.
#' @export
plot_excess_map <- function(excess, geo, flagged_only = FALSE,
                            subtitle = NULL, caption = NULL) {

  g <- sf::st_as_sf(geo)
  g[["excess"]]  <- excess$per_area[["excess"]]
  g[["flagged"]] <- excess$per_area[["flagged"]]

  if (flagged_only) g[["excess"]][!g[["flagged"]]] <- NA_real_

  lim <- max(abs(g[["excess"]]), na.rm = TRUE)

  if (is.null(caption)) {
    caption <- sprintf(
      paste0("Excess = posterior mean of (fitted \u2212 covariate-only ",
             "predicted) deaths per year.%s"),
      if (flagged_only)
        sprintf(" Only areas with P(residual RR > %.2f) > %.2f are shaded.",
                excess$threshold, excess$prob_cutoff) else ""
    )
  }

  ggplot2::ggplot(g) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[["excess"]]), colour = NA) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "#f7f7f7", high = "#8c2d04",
      midpoint = 0, limits = c(-lim, lim),
      na.value = "grey88",
      name = "Excess deaths\nper year"
    ) +
    ggplot2::labs(subtitle = subtitle, caption = caption) +
    theme_atlas(map = TRUE)
}


#' Compare two risk surfaces
#'
#' Results "Sensitivity to specifications" needs a BYM2-versus-ESF comparison.
#' A scatter answers the question the section actually asks - do the two
#' specifications rank areas the same way - more directly than two maps side by
#' side, where the eye cannot do the comparison.
#'
#' The Spearman correlation is printed in the subtitle, since that is the number
#' the sentence will quote.
#'
#' @param geo_a,geo_b Augmented `sf` objects from [augment_bym2()].
#' @param labels Axis labels for the two specifications.
#' @param value Column to compare. Default `"bym2_rr"`.
#' @param caption Passed to `labs()`.
#'
#' @return A ggplot.
#' @export
plot_rr_compare <- function(geo_a, geo_b,
                            labels = c("BYM2", "ESF"),
                            value  = "bym2_rr",
                            caption = NULL) {

  a <- sf::st_drop_geometry(geo_a)[[value]]
  b <- sf::st_drop_geometry(geo_b)[[value]]

  if (length(a) != length(b)) {
    stop("The two surfaces have ", length(a), " and ", length(b),
         " areas. They must be the same geography in the same order.",
         call. = FALSE)
  }

  rho <- stats::cor(a, b, method = "spearman")
  d   <- data.frame(a = a, b = b)

  ggplot2::ggplot(d, ggplot2::aes(x = .data[["a"]], y = .data[["b"]])) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22",
                         colour = "grey45") +
    ggplot2::geom_point(alpha = 0.55, size = 1.8, colour = "#8c2d04") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = paste(labels[1], "smoothed relative risk"),
      y = paste(labels[2], "smoothed relative risk"),
      subtitle = sprintf("Spearman rank correlation = %.3f", rho),
      caption  = caption
    ) +
    theme_atlas()
}


#' Rank-concordance heatmap across the mechanism strata
#'
#' Results asks which clusters persist across strata. [rank_concordance()]
#' produces the matrix; this draws it.
#'
#' @param concordance Output of [rank_concordance()].
#' @param caption Passed to `labs()`.
#'
#' @return A ggplot.
#' @export
plot_concordance <- function(concordance, caption = NULL) {

  m  <- concordance$matrix
  df <- expand.grid(a = rownames(m), b = colnames(m),
                    stringsAsFactors = FALSE)
  df[["rho"]] <- as.vector(m)

  ggplot2::ggplot(df, ggplot2::aes(x = .data[["a"]], y = .data[["b"]],
                                   fill = .data[["rho"]])) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data[["rho"]])),
                       size = 3, colour = "grey15") +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac", mid = "#f7f7f7", high = "#8c2d04",
      midpoint = 0, limits = c(-1, 1),
      name = "Spearman\ncorrelation"
    ) +
    ggplot2::labs(x = NULL, y = NULL, caption = caption) +
    theme_atlas() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      panel.grid  = ggplot2::element_blank()
    )
}


#' Save a figure at publication resolution
#'
#' Consistent dimensions and 300 dpi, so the figures do not have to be
#' re-exported one at a time when the thesis is assembled.
#'
#' @param plot A ggplot.
#' @param name File stem, without extension.
#' @param dir Output directory.
#' @param width Width in millimetres. Default `190` (full text width).
#' @param height Height in millimetres.
#' @param dpi Resolution.
#'
#' @return The file path, invisibly.
#' @export
save_figure <- function(plot, name, dir = "output/figures",
                        width = 190, height = 150, dpi = 300) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".png"))

  # ragg does proper font fallback through systemfonts, so characters the
  # default sans lacks - the >= in the exceedance legends, the en dash in the
  # probability bands, the minus sign in the excess caption - resolve to some
  # installed font instead of rendering as tofu. grDevices::png() on Windows
  # goes through GDI and does not fall back at all.
  args <- list(filename = path, plot = plot, width = width, height = height,
               units = "mm", dpi = dpi, bg = "white")
  if (requireNamespace("ragg", quietly = TRUE)) {
    args[["device"]] <- ragg::agg_png
  }
  do.call(ggplot2::ggsave, args)

  invisible(path)
}


#' NO2 and PM2.5 side by side
#'
#' The figure behind the pollutant selection argument. Both surfaces are in the
#' same units, so two presentations are possible and they answer different
#' questions.
#'
#' \strong{Free scales} (the default) show each pollutant's own spatial pattern:
#' where in the territory each is high. Use this to see whether the two
#' pollutants pick out the same places.
#'
#' \strong{Shared scale} puts both on one common ug/m3 legend. This is the
#' honest version of the selection argument, because it shows directly what the
#' methods claim - that PM2.5 varies far less across the study area than NO2
#' does, and therefore carries less contrast for a model to use. On free scales
#' that difference is invisible, since each panel is stretched to fill its own
#' range.
#'
#' @param geo Modelling `sf` carrying the pollutant columns. Column names are
#'   detected with [exposure_columns()], so the exposure year does not need to
#'   be hardcoded.
#' @param shared_scale Put both panels on one legend. Default `FALSE`.
#' @param caption Passed to `labs()`; a default noting the source is used when
#'   `NULL`.
#' @param cor_note Append the rank correlation between the two surfaces to the
#'   caption. Default `TRUE`.
#'
#' @return A patchwork object when `shared_scale = FALSE` and \pkg{patchwork} is
#'   available, otherwise a single faceted ggplot.
#' @examples
#' \dontrun{
#' plot_pollution_pair(smr_geo_full)
#' plot_pollution_pair(smr_geo_full, shared_scale = TRUE)
#' }
#' @seealso [exposure_columns()], [plot_travel_time()]
#' @export
plot_pollution_pair <- function(geo,
                                shared_scale = FALSE,
                                caption      = NULL,
                                cor_note     = TRUE) {

  g   <- sf::st_as_sf(geo)
  tab <- sf::st_drop_geometry(g)
  cols <- exposure_columns(tab)

  if (length(cols) < 2L) {
    stop("Need both a NO2 and a PM2.5 column. Found: ",
         if (length(cols)) paste(cols, collapse = ", ") else "none",
         ". Was add_pollution() run?", call. = FALSE)
  }

  # Subscripts as plotmath rather than U+2082/U+2085. R draws a plotmath
  # subscript by positioning a normal digit, so it cannot fall back to tofu the
  # way a literal subscript glyph does when the device font lacks it - which is
  # what the default png() device on Windows does with these two codepoints.
  # Unmatched columns are quoted so that every label is parseable and
  # label_parsed below can never be handed something it cannot evaluate.
  labs_pretty <- c("NO2 (ug/m3)" = "NO[2]", "PM2.5 (ug/m3)" = "PM[2.5]")
  nm <- ifelse(names(cols) %in% names(labs_pretty),
               labs_pretty[names(cols)],
               sprintf('"%s"', names(cols)))

  if (is.null(caption)) {
    caption <- paste0(
      "Annual mean, ", sub(".*_", "", cols[1]),
      ". EEA interpolated concentration maps at 1 km, area-weighted to the ",
      "modelling units."
    )
  }
  if (cor_note) {
    rho <- stats::cor(tab[[cols[1]]], tab[[cols[2]]], method = "spearman",
                      use = "complete.obs")
    caption <- paste0(caption, "\nSpearman rank correlation between the two ",
                      "surfaces: ", sprintf("%.3f", rho), ".")
  }

  # ---- one shared legend: the selection argument, made visible -------------
  if (shared_scale) {
    long <- do.call(rbind, lapply(seq_along(cols), function(i) {
      d <- g[, cols[i], drop = FALSE]
      names(d)[1] <- "value"
      d[["pollutant"]] <- factor(nm[i], levels = nm)
      d
    }))

    return(
      ggplot2::ggplot(long) +
        ggplot2::geom_sf(ggplot2::aes(fill = .data[["value"]]), colour = NA) +
        ggplot2::facet_wrap(~ .data[["pollutant"]], nrow = 1,
                            labeller = ggplot2::label_parsed) +
        ggplot2::scale_fill_viridis_c(
          option = "inferno", direction = -1,
          name = expression(paste(mu, "g/m"^3))
        ) +
        ggplot2::labs(caption = caption) +
        theme_atlas(map = TRUE) +
        ggplot2::theme(legend.position = "right")
    )
  }

  # ---- free scales: each pollutant's own pattern ---------------------------
  panel <- function(col, title) {
    ggplot2::ggplot(g) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data[[col]]), colour = NA) +
      ggplot2::scale_fill_viridis_c(
        option = "inferno", direction = -1,
        name = expression(paste(mu, "g/m"^3))
      ) +
      ggplot2::labs(title = parse(text = title)[[1]]) +
      theme_atlas(map = TRUE) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.key.width = ggplot2::unit(14, "mm"),
        legend.key.height = ggplot2::unit(3, "mm"),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )
  }

  p1 <- panel(unname(cols[1]), nm[1])
  p2 <- panel(unname(cols[2]), nm[2])

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("patchwork is not installed; falling back to a shared scale. ",
            "Install patchwork for independent legends per panel.")
    return(plot_pollution_pair(geo, shared_scale = TRUE, caption = caption,
                               cor_note = FALSE))
  }

  patchwork::wrap_plots(p1, p2, nrow = 1) +
    patchwork::plot_annotation(caption = caption,
                               theme = theme_atlas())
}
