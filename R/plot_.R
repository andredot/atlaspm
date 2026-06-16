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
