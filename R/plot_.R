#' Map indirectly standardised preventable mortality by area
#'
#' Draws a choropleth of the all-cause standardised mortality ratio
#' (\code{total_smr}) from a geometry-attached \code{smr} table (the output of
#' \code{preprocess_smr()} passed through \code{add_geo()}). The SMR is binned
#' into 5 classes anchored on 1.0 (= deaths matching the age-sex expectation) and
#' shown with an earthy diverging palette, no area borders and a clean
#' \code{theme_void}.
#'
#' @param smr An \code{sf} object with one row per area and a numeric SMR
#'   column, e.g. \code{add_geo(preprocess_smr(...), comuni)}.
#' @param value Name of the SMR column to map. Default \code{"total_smr"}.
#' @param breaks Numeric cut points defining the 5 classes. Default
#'   \code{c(-Inf, 0.9, 0.95, 1.05, 1.1, Inf)}, straddling 1.0 so the middle class
#'   is the "as expected" band.
#' @param title,subtitle,caption Plot annotations. Sensible defaults describe an
#'   indirectly age-sex standardised preventable-mortality map.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' smr_geo <- add_geo(preprocess_smr(mort_count, "popolazione.csv"), comuni)
#' plot_smr_map(smr_geo)
#' plot_smr_map(smr_geo, value = "C_lung_cancer_smr",
#'              title = "Lung-cancer mortality vs expectation")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_sf scale_fill_manual labs theme_void
#'   theme element_text margin unit
#' @importFrom dplyr mutate
#' @importFrom sf st_as_sf
#' @importFrom rlang .data
#' @export
plot_smr_map <- function(smr,
                         value    = "total_smr",
                         breaks   = c(-Inf, 0.9, 0.95, 1.05, 1.1, Inf),
                         title    = "Preventable mortality by area",
                         subtitle = "Indirectly standardised mortality ratio (SMR) for preventable causes",
                         caption  = "SMR = 1 means deaths match the area-wide age-sex expectation.") {

  labels <- c("< 0.90", "0.90 \u2013 0.95", "0.95 \u2013 1.05", "1.05 \u2013 1.10", "> 1.10")
  pal <- c(
    "< 0.90"        = "#5a8a7d",  # deep sage  (well below expected)
    "0.90 \u2013 0.95" = "#a3c4b5",  # soft green
    "0.95 \u2013 1.05" = "#f2e8d5",  # pale sand  (\u2248 expected)
    "1.05 \u2013 1.10" = "#dca678",  # warm clay
    "> 1.10"        = "#b5651d"   # terracotta (well above expected)
  )

  # Ensure a properly registered sf: a table whose geometry is a bare sfc
  # list-column (e.g. demoted by an upstream join) would make dplyr::mutate()
  # below fail with "all columns must be vectors / geometry is an sfc object".
  smr <- sf::st_as_sf(smr)

  smr <- dplyr::mutate(
    smr,
    .smr_class = cut(.data[[value]], breaks = breaks, labels = labels, right = FALSE)
  )

  ggplot2::ggplot(smr) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[".smr_class"]]), colour = NA) +
    ggplot2::scale_fill_manual(
      values   = pal,
      na.value = "grey85",
      drop     = FALSE,
      name     = "SMR\n(observed / expected)"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption) +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 15,
                                              margin = ggplot2::margin(b = 4)),
      plot.subtitle   = ggplot2::element_text(colour = "grey30",
                                              margin = ggplot2::margin(b = 10)),
      plot.caption    = ggplot2::element_text(colour = "grey45", size = 9, hjust = 0),
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.9, "lines"),
      plot.margin     = ggplot2::margin(12, 12, 12, 12)
    )
}


#' Faceted map of standardised mortality by category
#'
#' Draws a small-multiples choropleth of the standardised mortality ratio across
#' several categories at once - by default the seven mechanism columns
#' (\code{M_*_smr}) produced by \code{preprocess_smr()} and passed through
#' \code{add_geo()}. Each facet is one category; all share the same SMR bins
#' (anchored on 1.0), the same earthy diverging palette, no area borders and a
#' single shared legend, so the panels are directly comparable.
#'
#' The columns to map are selected with a \code{tidyselect} expression
#' (\code{cols}); their facet titles are the column names with the prefix and the
#' \code{"_smr"} suffix stripped and underscores turned to spaces.
#'
#' @param smr An \code{sf} object with one row per area and several SMR
#'   columns, e.g. \code{add_geo(preprocess_smr(...), comuni)}.
#' @param cols A \code{tidyselect} expression choosing the SMR columns to map.
#'   Default \code{dplyr::matches("^M_.*_smr$")} (the seven mechanism SMRs;
#'   note a bare \code{starts_with("M_")} would also catch the \code{_isr}
#'   columns).
#' @param breaks Numeric cut points defining the 5 classes. Default
#'   \code{c(-Inf, 0.9, 0.95, 1.05, 1.1, Inf)}, straddling 1.0.
#' @param strip_prefix,strip_suffix Regular expressions removed from each column
#'   name to make the facet label. Defaults strip a leading \code{"X_"} prefix
#'   and a trailing \code{"_smr"}.
#' @param ncol Number of facet columns. Default \code{4} (so 7 maps fill two
#'   rows).
#' @param title,subtitle,caption Plot annotations.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' smr_geo <- add_geo(preprocess_smr(mort_count, "popolazione.csv"), comuni)
#' plot_smr_facets(smr_geo)                                       # 7 mechanism maps
#' plot_smr_facets(smr_geo, cols = dplyr::matches("^G_.*_smr$"))  # by cause group
#' }
#'
#' @seealso \code{\link{plot_smr_map}}, \code{\link{preprocess_smr}}
#' @importFrom dplyr select all_of mutate matches row_number left_join |>
#' @importFrom tidyr pivot_longer
#' @importFrom tibble tibble
#' @importFrom sf st_as_sf st_drop_geometry st_geometry
#' @importFrom stringr str_remove str_replace_all
#' @importFrom ggplot2 ggplot aes geom_sf scale_fill_manual labs facet_wrap
#'   theme_void theme element_text margin unit
#' @importFrom rlang .data
#' @export
plot_smr_facets <- function(smr,
                            cols         = dplyr::matches("^M_.*_smr$"),
                            breaks       = c(-Inf, 0.9, 0.95, 1.05, 1.1, Inf),
                            strip_prefix = "^[A-Z]_",
                            strip_suffix = "_smr$",
                            ncol         = 4,
                            title    = "Standardised preventable mortality by mechanism, by area",
                            subtitle = "Indirectly age-sex standardised mortality ratio (SMR); 1 = deaths match the age-sex expectation",
                            caption  = "Each panel standardised on the same age-sex schedule; bins shared across panels.") {

  labels <- c("< 0.90", "0.90 \u2013 0.95", "0.95 \u2013 1.05", "1.05 \u2013 1.10", "> 1.10")
  pal <- c(
    "< 0.90"        = "#5a8a7d",  # deep sage  (well below expected)
    "0.90 \u2013 0.95" = "#a3c4b5",  # soft green
    "0.95 \u2013 1.05" = "#f2e8d5",  # pale sand  (\u2248 expected)
    "1.05 \u2013 1.10" = "#dca678",  # warm clay
    "> 1.10"        = "#b5651d"   # terracotta (well above expected)
  )

  smr <- sf::st_as_sf(smr)   # defend against a demoted (non-sf) input

  # The geometry (one polygon per area) is the same across every category, so
  # pivot the PLAIN attribute table to long, then attach a single geometry copy
  # per area back by row index. Reshaping an sf directly trips
  # "column `geometry` is an sfc object", so we keep geometry out of the pivot.
  geom_lookup <- tibble::tibble(
    .row     = seq_len(nrow(smr)),
    geometry = sf::st_geometry(smr)
  )

  long <- smr |>
    sf::st_drop_geometry() |>
    dplyr::mutate(.row = dplyr::row_number()) |>
    dplyr::select(".row", {{ cols }}) |>
    tidyr::pivot_longer(-".row", names_to = "category", values_to = "smr") |>
    dplyr::mutate(
      smr_class = cut(.data[["smr"]], breaks = breaks, labels = labels, right = FALSE),
      category  = stringr::str_remove(.data[["category"]], strip_prefix),
      category  = stringr::str_remove(.data[["category"]], strip_suffix),
      category  = stringr::str_replace_all(.data[["category"]], "_", " ")
    ) |>
    dplyr::left_join(geom_lookup, by = ".row") |>
    sf::st_as_sf(sf_column_name = "geometry")

  ggplot2::ggplot(long) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[["smr_class"]]), colour = NA) +
    ggplot2::facet_wrap(~ category, ncol = ncol) +
    ggplot2::scale_fill_manual(
      values   = pal,
      na.value = "grey85",
      drop     = FALSE,
      name     = "SMR\n(observed / expected)"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 15,
                                              margin = ggplot2::margin(b = 4)),
      plot.subtitle   = ggplot2::element_text(colour = "grey30",
                                              margin = ggplot2::margin(b = 10)),
      plot.caption    = ggplot2::element_text(colour = "grey45", size = 9, hjust = 0),
      strip.text      = ggplot2::element_text(face = "bold", size = 11,
                                              margin = ggplot2::margin(4, 0, 4, 0)),
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.9, "lines"),
      plot.margin     = ggplot2::margin(12, 12, 12, 12)
    )
}

#' Crude vs standardised mortality scatter (overall)
#'
#' One point per area: crude mortality rate on the x-axis, indirectly
#' standardised rate (ISR) on the y-axis, coloured by the SMR. The dashed
#' diagonal marks where the standardised rate equals the crude rate, so points
#' above it are pulled up by standardisation (younger-than-standard population)
#' and points below pulled down.
#'
#' @param cmr Crude-rate table from \code{preprocess_cmr()} (one row per area).
#' @param smr Standardised table from \code{preprocess_smr()} (one row per area).
#' @param group_var area key present in both. Default \code{"area"}.
#' @param crude,isr,smr_col Column names to use. Defaults \code{"total"},
#'   \code{"total_isr"}, \code{"total_smr"}.
#' @param title,subtitle Plot annotations.
#' @return A \code{ggplot} object.
#' @importFrom dplyr select all_of inner_join |>
#' @importFrom ggplot2 ggplot aes geom_abline geom_point scale_colour_gradient2
#'   labs theme_minimal theme element_text
#' @importFrom rlang .data
#' @export
plot_cmr_isr <- function(cmr, smr,
                         group_var = "area",
                         crude     = "total",
                         isr       = "total_isr",
                         smr_col   = "total_smr",
                         title     = "Crude vs standardised preventable mortality, by area",
                         subtitle  = "Each point a area; colour = SMR (observed / expected)") {

  d <- dplyr::inner_join(
    dplyr::select(cmr, dplyr::all_of(c(group_var, crude)))   |> stats::setNames(c(group_var, "crude")),
    dplyr::select(smr, dplyr::all_of(c(group_var, isr, smr_col))) |> stats::setNames(c(group_var, "isr", "smr")),
    by = group_var
  )

  ggplot2::ggplot(d, ggplot2::aes(x = .data[["crude"]], y = .data[["isr"]], colour = .data[["smr"]])) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_point(size = 2.4, alpha = 0.8) +
    ggplot2::geom_point(
      ggplot2::aes(text = .data[[group_var]]),
      size = 2.4, alpha = 0.8
    ) +
    ggplot2::scale_colour_gradient2(
      low = "#5a8a7d", mid = "#f2e8d5", high = "#b5651d", midpoint = 1,
      name = "SMR"
    ) +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = "Crude mortality rate (per 100,000)",
      y = "Indirectly standardised rate (per 100,000)"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey30")
    )
}


#' Crude vs standardised mortality scatter, faceted by category
#'
#' As \code{plot_cmr_isr()} but one panel per category (the seven mechanisms by
#' default): crude rate vs ISR, coloured by SMR, with free scales so each
#' category is readable on its own magnitude. Crude columns come from \code{cmr}
#' (e.g. \code{M_lifestyle_and_ncds}); the matching ISR/SMR columns come from
#' \code{smr} (\code{..._isr}, \code{..._smr}).
#'
#' @param cmr,smr Crude and standardised tables (one row per area).
#' @param group_var area key present in both. Default \code{"area"}.
#' @param prefix Column prefix selecting the category block. Default \code{"M_"}.
#' @param ncol Number of facet columns. Default \code{4}.
#' @param title,subtitle Plot annotations.
#' @return A \code{ggplot} object.
#' @importFrom dplyr select all_of matches mutate inner_join |>
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_remove str_replace_all
#' @importFrom ggplot2 ggplot aes geom_abline geom_point scale_colour_gradient2
#'   labs facet_wrap theme_minimal theme element_text
#' @importFrom rlang .data
#' @export
plot_cmr_isr_facets <- function(cmr, smr,
                                group_var = "area",
                                prefix    = "M_",
                                ncol      = 4,
                                title     = "Crude vs standardised preventable mortality by mechanism",
                                subtitle  = "Each point a area; colour = SMR. Free scales per panel.") {

  crude_long <- cmr |>
    dplyr::select(dplyr::all_of(group_var), dplyr::matches(paste0("^", prefix))) |>
    tidyr::pivot_longer(-dplyr::all_of(group_var), names_to = "category", values_to = "crude")

  isr_long <- smr |>
    dplyr::select(dplyr::all_of(group_var), dplyr::matches(paste0("^", prefix, ".*_isr$"))) |>
    tidyr::pivot_longer(-dplyr::all_of(group_var), names_to = "category", values_to = "isr") |>
    dplyr::mutate(category = stringr::str_remove(.data[["category"]], "_isr$"))

  smr_long <- smr |>
    dplyr::select(dplyr::all_of(group_var), dplyr::matches(paste0("^", prefix, ".*_smr$"))) |>
    tidyr::pivot_longer(-dplyr::all_of(group_var), names_to = "category", values_to = "smr") |>
    dplyr::mutate(category = stringr::str_remove(.data[["category"]], "_smr$"))

  d <- crude_long |>
    dplyr::inner_join(isr_long, by = c(group_var, "category")) |>
    dplyr::inner_join(smr_long, by = c(group_var, "category")) |>
    dplyr::mutate(
      category = stringr::str_remove(.data[["category"]], "^[A-Z]_"),
      category = stringr::str_replace_all(.data[["category"]], "_", " ")
    )

  ggplot2::ggplot(d, ggplot2::aes(x = .data[["crude"]], y = .data[["isr"]], colour = .data[["smr"]])) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_point(size = 1.8, alpha = 0.7) +
    ggplot2::facet_wrap(~ category, ncol = ncol, scales = "free") +
    ggplot2::scale_colour_gradient2(
      low = "#5a8a7d", mid = "#f2e8d5", high = "#b5651d", midpoint = 1,
      name = "SMR"
    ) +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = "Crude mortality rate (per 100,000)",
      y = "Indirectly standardised rate (per 100,000)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey30"),
      strip.text    = ggplot2::element_text(face = "bold")
    )
}

#' Standardised mortality rate vs a vulnerability/deprivation index scatter
#'
#' One point per area: a contextual index (IVSM, Deprivation Index, ...) on
#' the x-axis and the overall indirectly standardised mortality rate (ISR,
#' \emph{not} the SMR ratio) on the y-axis. A positive slope indicates that more
#' deprived/vulnerable comuni also carry higher standardised mortality. The ISR
#' comes from \code{\link{preprocess_smr}} (the \code{total_isr} column by
#' default); the index comes from \code{\link{import_ivsm}} or
#' \code{\link{build_deprivation_proxy}}. The two are joined on the shared
#' \code{area} key.
#'
#' The function is index-agnostic: point the \code{index_col} argument at the
#' relevant column and set \code{ref_line} to the index's natural anchor — 100
#' for IVSM (national average), 0 for the Deprivation Index (a sum of national
#' z-scores). Leave \code{ref_line = NULL} to omit the line.
#'
#' @param smr Standardised table from \code{preprocess_smr()} (one row per
#'   area), supplying the standardised rate.
#' @param index Index table (one row per area) from \code{import_ivsm()} or
#'   \code{build_deprivation_proxy()}.
#' @param group_var area key present in both. Default \code{"area"}.
#' @param isr Column in \code{smr} holding the overall standardised rate.
#'   Default \code{"total_isr"}.
#' @param index_col Column in \code{index} holding the index value. Default
#'   \code{"ivsm"}.
#' @param ref_line Numeric x-position for a dashed vertical reference line, or
#'   \code{NULL} for none. Default \code{NULL}. Use 100 for IVSM, 0 for the DI.
#' @param xlab X-axis label. Default a generic index label.
#' @param smooth Logical; add a linear trend line. Default \code{TRUE}.
#' @param title,subtitle Plot annotations.
#' @return A \code{ggplot} object.
#' @importFrom dplyr select all_of inner_join
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_vline labs
#'   theme_minimal theme element_text
#' @importFrom rlang .data
#' @export
plot_scatter_smr_index <- function(smr, index,
                                   group_var = "area",
                                   isr       = "total_isr",
                                   index_col = "ivsm",
                                   ref_line  = NULL,
                                   xlab      = "Vulnerability / deprivation index",
                                   smooth    = TRUE,
                                   title     = "Standardised mortality vs index, by area",
                                   subtitle  = "Each point a area; x = index, y = indirectly standardised rate") {

  d <- dplyr::inner_join(
    dplyr::select(smr,   dplyr::all_of(c(group_var, isr)))       |> stats::setNames(c(group_var, "isr")),
    dplyr::select(index, dplyr::all_of(c(group_var, index_col))) |> stats::setNames(c(group_var, "index")),
    by = group_var
  )

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[["index"]], y = .data[["isr"]]))

  if (!is.null(ref_line)) {
    p <- p + ggplot2::geom_vline(xintercept = ref_line,
                                 linetype = "dashed", colour = "grey70")
  }

  p <- p + ggplot2::geom_point(size = 2.4, alpha = 0.8, colour = "#3b528b")

  if (smooth) {
    p <- p + ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                                  se = TRUE, colour = "#b5651d", fill = "#f2e8d5")
  }

  p +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = xlab,
      y = "Indirectly standardised mortality rate (per 100,000)"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey30")
    )
}

#' Strategic prioritisation map from BYM2 exceedance probabilities
#'
#' Maps the exceedance-probability column from \code{\link{augment_bym2}} as four
#' tiers using a calibrated-likelihood lexicon (almost certain / likely /
#' possible / unlikely), with an earthy palette. The tiers express how strong the
#' evidence is that an area's relative risk exceeds the threshold - which is what
#' an exceedance probability measures - rather than the point estimate alone.
#'
#' The \code{threshold} here only sets the *label* - the underlying probabilities
#' were fixed when \code{augment_bym2()} ran. To keep label and data in sync, the
#' function reads the threshold stored on the \code{sf} by \code{augment_bym2()}
#' when \code{threshold} is left \code{NULL}, falling back to \code{1.10} with a
#' message if no stored value is found.
#'
#' @param geo An \code{sf} from \code{augment_bym2()} with an exceedance column.
#' @param value Name of the exceedance-probability column. Default
#'   \code{"bym2_exceed"}.
#' @param threshold Relative-risk threshold used for the label. \code{NULL}
#'   (default) reads the value recorded by \code{augment_bym2()}.
#' @param title,subtitle Plot annotations. \code{subtitle = NULL} is built from
#'   the threshold.
#' @param threshold_label Optional manual override of the threshold text.
#' @param legend_title Legend heading. Default
#'   \code{"Evidence of excess mortality"}.
#' @return A \code{ggplot} object.
#' @importFrom dplyr mutate case_when
#' @importFrom sf st_as_sf
#' @importFrom stats setNames
#' @importFrom ggplot2 ggplot aes geom_sf scale_fill_manual labs theme_void
#'   theme element_text margin unit
#' @importFrom rlang .data
#' @export
plot_exceedance_map <- function(geo,
                                value           = "bym2_exceed",
                                threshold       = NULL,
                                title           = "Strategic prioritisation map",
                                subtitle        = NULL,
                                threshold_label = NULL,
                                legend_title    = "Evidence of excess mortality") {

  # resolve the threshold: explicit arg > attribute from augment_bym2() > 1.10
  if (is.null(threshold)) {
    threshold <- attr(geo, "bym2_exceed_threshold")
    if (is.null(threshold)) {
      threshold <- 1.10
      message("No `threshold` supplied or stored on `geo`; labelling as RR > 1.10. ",
              "Pass `threshold` to match the value used in augment_bym2().")
    }
  }
  if (is.null(threshold_label)) {
    pct <- round((threshold - 1) * 100)
    threshold_label <- sprintf("+%d%% (RR > %.2f)", pct, threshold)
  }
  if (is.null(subtitle)) {
    subtitle <- paste0("Probability that preventable mortality exceeds ", threshold_label)
  }

  # labels defined ONCE, worst-to-best; palette derives from them so the two
  # can never drift out of sync (the bug that drops a tier to grey)
  tier_labels <- c(
    "Almost certain (\u2265 95%)",
    "Likely (80\u201395%)",
    "Possible (50\u201380%)",
    "Unlikely (< 50%)"
  )
  pal <- stats::setNames(
    c("#8c2d04",   # burnt umber
      "#d94801",   # terracotta
      "#fdb863",   # ochre
      "#eaeaea"),  # light stone
    tier_labels
  )

  geo <- sf::st_as_sf(geo) |>
    dplyr::mutate(
      .tier = factor(
        dplyr::case_when(
          .data[[value]] >= 0.95 ~ tier_labels[1],
          .data[[value]] >= 0.80 ~ tier_labels[2],
          .data[[value]] >= 0.50 ~ tier_labels[3],
          TRUE                   ~ tier_labels[4]
        ),
        levels = tier_labels
      )
    )

  ggplot2::ggplot(geo) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[".tier"]]), colour = NA) +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE, name = legend_title) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 15,
                                              margin = ggplot2::margin(b = 4)),
      plot.subtitle   = ggplot2::element_text(colour = "grey30",
                                              margin = ggplot2::margin(b = 10)),
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.9, "lines"),
      plot.margin     = ggplot2::margin(12, 12, 12, 12)
    )
}

#' Faceted strategic-prioritisation map from per-mechanism exceedance probabilities
#'
#' The small-multiples counterpart to \code{\link{plot_exceedance_map}}: instead
#' of one panel for a single \code{bym2_exceed} column, it facets over the
#' per-mechanism exceedance columns produced by
#' \code{\link{augment_bym2_mechanisms}} (named \code{<stem>_bym2_exc}). Each
#' panel bins its comuni into the same four calibrated-likelihood tiers
#' (almost certain / likely / possible / unlikely) on a shared palette and a
#' single legend, so the mechanisms are directly comparable by how strong the
#' evidence of excess mortality is.
#'
#' Columns are selected with a \code{tidyselect} expression (\code{cols}); facet
#' titles are the column names with the prefix and the \code{_bym2_exc} suffix
#' stripped and underscores turned to spaces.
#'
#' @param geo An \code{sf} from \code{\link{augment_bym2_mechanisms}} carrying
#'   the \code{<stem>_bym2_exc} columns.
#' @param cols A \code{tidyselect} expression choosing the exceedance columns.
#'   Default \code{dplyr::matches("^M_.*_bym2_exc$")} (the seven mechanisms).
#' @param threshold Relative-risk threshold the probabilities refer to. Default
#'   \code{NULL} reads the value stored by \code{augment_bym2_mechanisms()},
#'   falling back to \code{1.10}.
#' @param strip_prefix,strip_suffix Regexes removed from each column name to make
#'   the facet label. Defaults strip a leading \code{"M_"} and a trailing
#'   \code{"_bym2_exc"}.
#' @param ncol Number of facet columns. Default \code{4} (seven panels -> two
#'   rows).
#' @param title,subtitle,caption Plot annotations. \code{subtitle = NULL} is
#'   built from the threshold.
#' @param threshold_label Optional manual override of the threshold text.
#' @param legend_title Legend heading. Default
#'   \code{"Evidence of excess mortality"}.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' plot_exceedance_facets(smr_geo_mech_bym2)
#' }
#' @importFrom dplyr select all_of mutate matches row_number left_join |>
#' @importFrom tidyr pivot_longer
#' @importFrom tibble tibble
#' @importFrom sf st_as_sf st_drop_geometry st_geometry
#' @importFrom stringr str_remove str_replace_all
#' @importFrom stats setNames
#' @importFrom ggplot2 ggplot aes geom_sf scale_fill_manual labs facet_wrap
#'   theme_void theme element_text margin unit
#' @importFrom rlang .data
#' @export
plot_exceedance_facets <- function(geo,
                                   cols            = dplyr::matches("^M_.*_bym2_exc$"),
                                   threshold       = NULL,
                                   strip_prefix    = "^M_",
                                   strip_suffix    = "_bym2_exc$",
                                   ncol            = 4,
                                   title    = "Strategic prioritisation map by mechanism, by area",
                                   subtitle = NULL,
                                   caption  = "Per-mechanism BYM2 exceedance probabilities; tiers shared across panels.",
                                   threshold_label = NULL,
                                   legend_title    = "Evidence of excess mortality") {

  # resolve the threshold: explicit arg > attribute > 1.10
  if (is.null(threshold)) {
    threshold <- attr(geo, "bym2_exceed_threshold")
    if (is.null(threshold)) {
      threshold <- 1.10
      message("No `threshold` supplied or stored on `geo`; labelling as RR > 1.10.")
    }
  }
  if (is.null(threshold_label)) {
    pct <- round((threshold - 1) * 100)
    threshold_label <- sprintf("+%d%% (RR > %.2f)", pct, threshold)
  }
  if (is.null(subtitle)) {
    subtitle <- paste0("Probability that preventable mortality exceeds ", threshold_label)
  }

  # tiers defined ONCE, worst-to-best; palette derives from them
  tier_labels <- c(
    "Almost certain (\u2265 95%)",
    "Likely (80\u201395%)",
    "Possible (50\u201380%)",
    "Unlikely (< 50%)"
  )
  pal <- stats::setNames(
    c("#8c2d04", "#d94801", "#fdb863", "#eaeaea"),
    tier_labels
  )

  geo <- sf::st_as_sf(geo)   # defend against a demoted (non-sf) input

  # Reshaping an sf directly trips "column `geometry` is an sfc object", so pivot
  # the PLAIN attribute table long, then re-attach one geometry copy per area
  # by row index (geometry is identical across every mechanism panel).
  geom_col <- attr(geo, "sf_column")
  geo_idx  <- dplyr::mutate(geo, .row = dplyr::row_number())

  attr_long <- sf::st_drop_geometry(geo_idx) |>
    dplyr::select(.row, !!cols) |>
    tidyr::pivot_longer(
      cols      = -.row,
      names_to  = "mechanism",
      values_to = ".exceed"
    ) |>
    dplyr::mutate(
      mechanism = stringr::str_replace_all(
        stringr::str_remove(stringr::str_remove(mechanism, strip_prefix),
                            strip_suffix),
        "_", " "),
      .tier = factor(
        dplyr::case_when(
          .exceed >= 0.95 ~ tier_labels[1],
          .exceed >= 0.80 ~ tier_labels[2],
          .exceed >= 0.50 ~ tier_labels[3],
          TRUE            ~ tier_labels[4]
        ),
        levels = tier_labels
      )
    )

  # one geometry per area, keyed by the same row index
  geo_lookup <- tibble::tibble(.row = geo_idx$.row)
  geo_lookup[[geom_col]] <- sf::st_geometry(geo)

  plot_df <- attr_long |>
    dplyr::left_join(geo_lookup, by = ".row") |>
    sf::st_as_sf(sf_column_name = geom_col)

  ggplot2::ggplot(plot_df) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[".tier"]]), colour = NA) +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE, name = legend_title) +
    ggplot2::facet_wrap(~ mechanism, ncol = ncol) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption) +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 15,
                                              margin = ggplot2::margin(b = 4)),
      plot.subtitle   = ggplot2::element_text(colour = "grey30",
                                              margin = ggplot2::margin(b = 10)),
      plot.caption    = ggplot2::element_text(colour = "grey45", size = 9, hjust = 0),
      strip.text      = ggplot2::element_text(face = "bold", size = 11,
                                              margin = ggplot2::margin(b = 4)),
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.9, "lines"),
      plot.margin     = ggplot2::margin(12, 12, 12, 12)
    )
}
