#' Map indirectly standardised preventable mortality by comune
#'
#' Draws a choropleth of the all-cause standardised mortality ratio
#' (\code{total_smr}) from a geometry-attached \code{smr} table (the output of
#' \code{preprocess_smr()} passed through \code{add_geo()}). The SMR is binned
#' into 5 classes anchored on 1.0 (= deaths matching the age-sex expectation) and
#' shown with an earthy diverging palette, no comune borders and a clean
#' \code{theme_void}.
#'
#' @param smr An \code{sf} object with one row per comune and a numeric SMR
#'   column, e.g. \code{add_geo(preprocess_smr(...), comuni)}.
#' @param value Name of the SMR column to map. Default \code{"total_smr"}.
#' @param breaks Numeric cut points defining the 5 classes. Default
#'   \code{c(-Inf, 0.5, 0.8, 1.25, 2, Inf)}, straddling 1.0 so the middle class
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
                         breaks   = c(-Inf, 0.5, 0.8, 1.25, 2, Inf),
                         title    = "Preventable mortality by Comune",
                         subtitle = "Indirectly standardised mortality ratio (SMR) for preventable causes",
                         caption  = "SMR = 1 means deaths match the area-wide age-sex expectation.") {

  labels <- c("< 0.50", "0.50 \u2013 0.80", "0.80 \u2013 1.25", "1.25 \u2013 2.00", "> 2.00")
  pal <- c(
    "< 0.50"        = "#5a8a7d",  # deep sage  (well below expected)
    "0.50 \u2013 0.80" = "#a3c4b5",  # soft green
    "0.80 \u2013 1.25" = "#f2e8d5",  # pale sand  (\u2248 expected)
    "1.25 \u2013 2.00" = "#dca678",  # warm clay
    "> 2.00"        = "#b5651d"   # terracotta (well above expected)
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
#' (anchored on 1.0), the same earthy diverging palette, no comune borders and a
#' single shared legend, so the panels are directly comparable.
#'
#' The columns to map are selected with a \code{tidyselect} expression
#' (\code{cols}); their facet titles are the column names with the prefix and the
#' \code{"_smr"} suffix stripped and underscores turned to spaces.
#'
#' @param smr An \code{sf} object with one row per comune and several SMR
#'   columns, e.g. \code{add_geo(preprocess_smr(...), comuni)}.
#' @param cols A \code{tidyselect} expression choosing the SMR columns to map.
#'   Default \code{dplyr::matches("^M_.*_smr$")} (the seven mechanism SMRs;
#'   note a bare \code{starts_with("M_")} would also catch the \code{_isr}
#'   columns).
#' @param breaks Numeric cut points defining the 5 classes. Default
#'   \code{c(-Inf, 0.5, 0.8, 1.25, 2, Inf)}, straddling 1.0.
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
                            breaks       = c(-Inf, 0.5, 0.8, 1.25, 2, Inf),
                            strip_prefix = "^[A-Z]_",
                            strip_suffix = "_smr$",
                            ncol         = 4,
                            title    = "Standardised preventable mortality by mechanism, by comune",
                            subtitle = "Indirectly age-sex standardised mortality ratio (SMR); 1 = deaths match the age-sex expectation",
                            caption  = "Each panel standardised on the same age-sex schedule; bins shared across panels.") {

  labels <- c("< 0.50", "0.50 \u2013 0.80", "0.80 \u2013 1.25", "1.25 \u2013 2.00", "> 2.00")
  pal <- c(
    "< 0.50"        = "#5a8a7d",
    "0.50 \u2013 0.80" = "#a3c4b5",
    "0.80 \u2013 1.25" = "#f2e8d5",
    "1.25 \u2013 2.00" = "#dca678",
    "> 2.00"        = "#b5651d"
  )

  smr <- sf::st_as_sf(smr)   # defend against a demoted (non-sf) input

  # The geometry (one polygon per comune) is the same across every category, so
  # pivot the PLAIN attribute table to long, then attach a single geometry copy
  # per comune back by row index. Reshaping an sf directly trips
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
#' One point per comune: crude mortality rate on the x-axis, indirectly
#' standardised rate (ISR) on the y-axis, coloured by the SMR. The dashed
#' diagonal marks where the standardised rate equals the crude rate, so points
#' above it are pulled up by standardisation (younger-than-standard population)
#' and points below pulled down.
#'
#' @param cmr Crude-rate table from \code{preprocess_cmr()} (one row per comune).
#' @param smr Standardised table from \code{preprocess_smr()} (one row per comune).
#' @param group_var Comune key present in both. Default \code{"comune"}.
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
                         group_var = "comune",
                         crude     = "total",
                         isr       = "total_isr",
                         smr_col   = "total_smr",
                         title     = "Crude vs standardised preventable mortality, by comune",
                         subtitle  = "Each point a comune; colour = SMR (observed / expected)") {

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
#' @param cmr,smr Crude and standardised tables (one row per comune).
#' @param group_var Comune key present in both. Default \code{"comune"}.
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
                                group_var = "comune",
                                prefix    = "M_",
                                ncol      = 4,
                                title     = "Crude vs standardised preventable mortality by mechanism",
                                subtitle  = "Each point a comune; colour = SMR. Free scales per panel.") {

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
