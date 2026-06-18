#' Fit a BYM2 (ICAR) spatial model to standardised mortality
#'
#' Fits a Poisson \code{geostan::stan_icar} model with \code{type = "bym2"} to
#' smooth the comune-level mortality relative risk. The response is the observed
#' death count and the offset is \code{log(expected)}, so the latent quantity is
#' the relative risk RR = observed / expected on the smoothed scale (i.e. a
#' model-based SMR). This base model has no covariates - just the spatial random
#' effect; covariates (IVSM, etc.) can be added later via \code{formula}.
#'
#' \code{geo} and \code{C} must be in the same row order (both built from the
#' same \code{sf}); \code{\link{build_adjacency}} guarantees this.
#'
#' @param geo An \code{sf} with one row per comune containing the observed-count
#'   and expected-count columns (\code{obs_col}, \code{exp_col}), in the same
#'   order as \code{C}.
#' @param C Binary adjacency matrix from \code{\link{build_adjacency}}.
#' @param formula Model formula. Default
#'   \code{total_obs ~ offset(log(total_exp))} - intercept + spatial effect only.
#' @param obs_col,exp_col Names of the observed- and expected-count columns,
#'   used only to guard against zero/NA expected. Defaults \code{"total_obs"},
#'   \code{"total_exp"}.
#' @param chains,iter Sampler settings passed to \code{stan_icar}. Defaults
#'   \code{4} and \code{4000}.
#' @param ... Further arguments forwarded to \code{geostan::stan_icar()}
#'   (e.g. \code{prior}, \code{control}, \code{refresh}).
#' @return A fitted \code{geostan} model object (\code{stanfit} wrapper).
#' @examples
#' \dontrun{
#' fit <- fit_bym2(geo, C)
#' }
#' @importFrom geostan stan_icar
#' @importFrom stats poisson as.formula offset
#' @export
fit_bym2 <- function(geo, C,
                     formula = total_obs ~ offset(log(total_exp)),
                     obs_col = "total_obs",
                     exp_col = "total_exp",
                     chains  = 4,
                     iter    = 4000,
                     ...) {

  # expected deaths must be positive for the log-offset; a comune with E == 0
  # contributes no information and would produce log(0) = -Inf.
  if (any(geo[[exp_col]] <= 0 | is.na(geo[[exp_col]]))) {
    warning("Some comuni have expected deaths <= 0 or NA; their log-offset is undefined. ",
            "Consider filtering them before fitting.", call. = FALSE)
  }

  geostan::stan_icar(
    formula = formula,
    data    = geo,
    family  = stats::poisson(),
    C       = C,
    type    = "bym2",
    chains  = chains,
    iter    = iter,
    ...
  )
}


#' Attach BYM2 smoothed estimates back onto the comune sf
#'
#' Pulls the posterior fitted values from a \code{\link{fit_bym2}} model and adds
#' three columns to the geometry table:
#' \itemize{
#'   \item \code{bym2_rr} - posterior mean smoothed relative risk (model SMR);
#'   \item \code{bym2_isr} - smoothed indirectly standardised rate, i.e. the
#'     posterior mean fitted deaths divided by population, per \code{per};
#'   \item \code{bym2_exceed} - exceedance probability P(RR > \code{threshold}),
#'     the share of posterior draws above the threshold (1 by default).
#' }
#' Row order of \code{fit} and \code{geo} must match (both came from the same
#' \code{sf}); no join is done, columns are bound by position.
#'
#' @param geo The \code{sf} passed to \code{\link{fit_bym2}} (same row order).
#' @param fit The fitted model from \code{\link{fit_bym2}}.
#' @param exp_col,pop_col Names of the expected-count and population columns in
#'   \code{geo}. Defaults \code{"total_exp"}, \code{"population"}.
#' @param per Rate multiplier for \code{bym2_isr}. Default \code{100000}.
#' @param threshold Relative-risk threshold for the exceedance probability.
#'   Default \code{1} (more deaths than expected).
#' @return \code{geo} with the three \code{bym2_*} columns added.
#' @examples
#' \dontrun{
#' geo2 <- augment_bym2(geo, fit)
#' plot_smr_map(geo2, value = "bym2_isr",
#'              title = "BYM2-smoothed preventable mortality")
#' }
#' @importFrom rstan extract
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
augment_bym2 <- function(geo, fit,
                         exp_col   = "total_exp",
                         pop_col   = "population",
                         per       = 100000,
                         threshold = 1) {

  # geostan stores the posterior expected COUNTS under the Stan parameter
  # "fitted"; pull them straight from the underlying stanfit. (geostan's S3
  # fitted() method is registered but not exported, so geostan::fitted() errors,
  # and posterior_epred() has no geostan method - hence the direct extract.)
  mu <- rstan::extract(fit$stanfit, pars = "fitted")$fitted   # [draws, n_areas]
  E  <- geo[[exp_col]]

  # relative-risk draws = fitted count / expected count
  rr <- sweep(mu, MARGIN = 2, STATS = E, FUN = "/")           # [draws, n_areas]

  geo |>
    dplyr::mutate(
      bym2_rr     = colMeans(rr),
      bym2_isr    = colMeans(mu) / .data[[pop_col]] * per,
      bym2_exceed = colMeans(rr > threshold)
    )
}


#' Strategic prioritisation map from BYM2 exceedance probabilities
#'
#' Maps the exceedance probability column from \code{\link{augment_bym2}} as
#' four action tiers (urgent / high / watchlist / normal), using an earthy
#' palette. The tiers express how *certain* the model is that a comune's
#' relative risk exceeds the threshold, which is more decision-relevant than the
#' point estimate alone.
#'
#' @param geo An \code{sf} from \code{augment_bym2()} with an exceedance column.
#' @param value Name of the exceedance-probability column. Default
#'   \code{"bym2_exceed"}.
#' @param threshold_label Text describing the risk threshold, used in the
#'   subtitle. Default \code{"+20% (RR > 1.10)"} - set to match the
#'   \code{threshold} used in \code{augment_bym2()}.
#' @param title,subtitle Plot annotations.
#' @return A \code{ggplot} object.
#' @importFrom dplyr mutate case_when
#' @importFrom sf st_as_sf
#' @importFrom ggplot2 ggplot aes geom_sf scale_fill_manual labs theme_void
#'   theme element_text margin unit
#' @importFrom rlang .data
#' @export
plot_exceedance_map <- function(geo,
                                value           = "bym2_exceed",
                                threshold_label = "+10% (RR > 1.10)",
                                title    = "Strategic prioritisation map",
                                subtitle = NULL) {

  if (is.null(subtitle)) {
    subtitle <- paste0("Probability that preventable mortality exceeds ", threshold_label)
  }

  tiers <- c("Urgent (\u2265 95%)", "High (80\u201395%)", "Watchlist (50\u201380%)", "Low (< 50%)")
  pal <- c(
    "Urgent (\u2265 95%)"     = "#8c2d04",  # burnt umber
    "High (80\u201395%)"      = "#d94801",  # terracotta
    "Watchlist (50\u201380%)" = "#fdb863",  # ochre
    "Low (< 50%)"          = "#eaeaea"   # light stone
  )

  geo <- sf::st_as_sf(geo) |>
    dplyr::mutate(
      .tier = factor(
        dplyr::case_when(
          .data[[value]] >= 0.95 ~ tiers[1],
          .data[[value]] >= 0.80 ~ tiers[2],
          .data[[value]] >= 0.50 ~ tiers[3],
          TRUE                   ~ tiers[4]
        ),
        levels = tiers
      )
    )

  ggplot2::ggplot(geo) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[".tier"]]), colour = NA) +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE, name = "Intervention priority") +
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
