#' Compute the BYM2 scaling factor for a connected component
#'
#' Internal helper implementing the Riebler et al. (2016) / Morris et al. (2019)
#' scaling: the geometric mean of the marginal variances of the (sum-to-zero
#' constrained) ICAR precision matrix built from \code{C}.
#'
#' @param C Binary adjacency matrix for a *single* connected component.
#' @return A scalar scaling factor.
#' @keywords internal
#' @noRd
scale_c <- function(C) {
  geomean <- function(x) exp(mean(log(x)))
  Cs <- methods::as(C, "CsparseMatrix")
  N  <- nrow(Cs)
  Q  <- Matrix::Diagonal(N, Matrix::rowSums(Cs)) - Cs
  Q_pert <- Q + Matrix::Diagonal(N) * max(Matrix::diag(Q)) * sqrt(.Machine$double.eps)
  Sigma  <- Matrix::solve(Q_pert)                       # dense N x N inverse
  SA     <- Matrix::rowSums(Sigma)                      # Sigma %*% 1
  diag_c <- Matrix::diag(Sigma) - SA^2 / sum(SA)        # sum-to-zero correction
  geomean(as.numeric(diag_c))
}


#' Prepare the BYM2 scale factor from an adjacency matrix
#'
#' Produces the \code{k}-length \code{scale_factor} vector (one entry per
#' connected component of the graph) that \code{geostan::stan_icar()} needs for a
#' *properly scaled* BYM2 model. Without this, \code{stan_icar} silently sets the
#' scale factor to a vector of ones, which turns BYM2 back into an unscaled BYM
#' and breaks the interpretation of \code{rho} (the proportion of variance that
#' is spatially structured) and the portability of the priors. Singletons /
#' islands receive a scale factor of 1, following Freni-Sterrantino et al. (2018)
#' and Donegan & Morris (2021).
#'
#' @param C Binary adjacency matrix from \code{\link{build_adjacency}}, in the
#'   same row order as the \code{sf} you will fit.
#' @return Numeric vector of length \code{k} (number of connected components),
#'   suitable for the \code{scale_factor} argument of \code{geostan::stan_icar()}.
#' @examples
#' \dontrun{
#' S <- compute_scale_factor(C)
#' fit <- fit_bym2(geo, C, scale_factor = S)
#' }
#' @importFrom geostan prep_icar_data
#' @export
compute_scale_factor <- function(C) {
  icar_data <- geostan::prep_icar_data(C)
  k  <- icar_data$k
  sf <- numeric(k)
  for (j in seq_len(k)) {
    idx <- which(icar_data$comp_id == j)
    if (length(idx) == 1L) {       # singleton / island: no neighbours to scale
      sf[j] <- 1
      next
    }
    sf[j] <- scale_c(C[idx, idx, drop = FALSE])
  }
  sf
}


#' Fit a (properly scaled) BYM2 spatial model to standardised mortality
#'
#' Fits a Poisson \code{geostan::stan_icar} model with \code{type = "bym2"} to
#' smooth comune-level mortality relative risk. The response is the observed
#' death count and the offset is \code{log(expected)}, so the latent quantity is
#' the relative risk RR = observed / expected on the smoothed scale (a model-based
#' SMR). The base model has only the spatial random effect; covariates (IVSM,
#' etc.) can be added via \code{formula}.
#'
#' Differences from a bare \code{stan_icar} call, all of which affect the output:
#' \itemize{
#'   \item \strong{Scaling.} If \code{scale_factor} is not supplied it is computed
#'     with \code{\link{compute_scale_factor}}, so this is a *true*
#'     BYM2 in which \code{rho} is the share of spatially-structured variance and
#'     priors are graph-portable. Pass \code{scale_factor = NULL} explicitly only
#'     if you understand you are then fitting an unscaled model.
#'   \item \strong{Measurement error.} \code{ME} is forwarded to \code{stan_icar}.
#'     When a covariate is modelled with error, \code{stan_icar} treats the
#'     observed covariate as a noisy reading of a latent true value and estimates
#'     \code{beta} on that latent value, so the covariate's uncertainty
#'     propagates into \code{beta} (and hence into the adjusted RR) instead of
#'     being ignored. Build \code{ME} with \code{geostan::prep_me_data()}; the
#'     column names in its \code{se} data frame must match the covariate names in
#'     \code{formula}.
#'   \item \strong{Sampler tuning.} \code{control} defaults to a higher
#'     \code{adapt_delta} and \code{max_treedepth} because BYM2's funnel geometry
#'     (the \code{spatial_scale}-by-\code{phi} interaction) routinely produces
#'     divergent transitions at Stan's defaults. Divergences bias the posterior
#'     and especially the tail-sensitive exceedance probabilities, so always run
#'     \code{\link{check_bym2_fit}} afterwards.
#' }
#'
#' \code{geo} and \code{C} must be in the same row order (both built from the
#' same \code{sf}); \code{\link{build_adjacency}} guarantees this.
#'
#' @param geo An \code{sf} with one row per comune containing the observed- and
#'   expected-count columns (\code{obs_col}, \code{exp_col}), in the same order
#'   as \code{C}.
#' @param C Binary adjacency matrix from \code{\link{build_adjacency}}.
#' @param formula Model formula. Default
#'   \code{total_obs ~ offset(log(total_exp))} - intercept + spatial effect only.
#' @param ME Optional measurement-error specification from
#'   \code{geostan::prep_me_data()}; \code{NULL} (default) for none.
#' @param scale_factor BYM2 scaling vector. If \code{NULL} (the explicit default
#'   marker) it is computed from \code{C} via \code{\link{compute_scale_factor}}.
#' @param obs_col,exp_col Names of the observed- and expected-count columns,
#'   used only to guard against zero/NA expected. Defaults \code{"total_obs"},
#'   \code{"total_exp"}.
#' @param chains,iter Sampler settings passed to \code{stan_icar}. Defaults
#'   \code{4} and \code{4000}.
#' @param control List of Stan control parameters. Default
#'   \code{list(adapt_delta = 0.97, max_treedepth = 12)}.
#' @param ... Further arguments forwarded to \code{geostan::stan_icar()}
#'   (e.g. \code{prior}, \code{slx}, \code{re}, \code{censor_point},
#'   \code{refresh}). Do NOT pass \code{slim = TRUE} or \code{drop} the
#'   \code{"fitted"} quantities - \code{\link{augment_bym2}} needs them.
#' @return A fitted \code{geostan} model object (\code{stanfit} wrapper).
#' @examples
#' \dontrun{
#' # base spatial-only model (scale factor computed automatically)
#' fit <- fit_bym2(geo, C)
#'
#' # with IVSM measured with error, so its uncertainty enters beta
#' ME  <- geostan::prep_me_data(se = data.frame(IVSM = geo$IVSM_se))
#' fit <- fit_bym2(geo, C,
#'                 formula = total_obs ~ offset(log(total_exp)) + IVSM,
#'                 ME = ME)
#' check_bym2_fit(fit)
#' }
#' @importFrom geostan stan_icar
#' @importFrom stats poisson as.formula offset
#' @export
fit_bym2 <- function(geo, C,
                     formula = total_obs ~ offset(log(total_exp)),
                     ME      = NULL,
                     scale_factor = NULL,
                     obs_col = "total_obs",
                     exp_col = "total_exp",
                     chains  = 4,
                     iter    = 4000,
                     control = list(adapt_delta = 0.97, max_treedepth = 12),
                     ...) {

  # expected deaths must be positive for the log-offset; a comune with E == 0
  # contributes no information and would produce log(0) = -Inf.
  if (any(geo[[exp_col]] <= 0 | is.na(geo[[exp_col]]))) {
    warning("Some comuni have expected deaths <= 0 or NA; their log-offset is undefined. ",
            "Consider filtering them before fitting.", call. = FALSE)
  }

  # Proper BYM2 scaling: compute from the graph unless the caller supplied one.
  if (is.null(scale_factor)) {
    scale_factor <- compute_scale_factor(C)
  }

  geostan::stan_icar(
    formula      = formula,
    data         = geo,
    family       = stats::poisson(),
    C            = C,
    type         = "bym2",
    scale_factor = scale_factor,
    ME           = ME,
    chains       = chains,
    iter         = iter,
    control      = control,
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
#'     the share of posterior draws above the threshold.
#' }
#' The \code{threshold} used is also stored as an attribute
#' (\code{"bym2_exceed_threshold"}) on the returned object so that
#' \code{\link{plot_exceedance_map}} can label the map consistently without you
#' restating it.
#'
#' Row order of \code{fit} and \code{geo} must match (both came from the same
#' \code{sf}); no join is done, columns are bound by position.
#'
#' @param geo The \code{sf} passed to \code{\link{fit_bym2}} (same row order).
#' @param fit The fitted model from \code{\link{fit_bym2}}.
#' @param exp_col,pop_col Names of the expected-count and population columns in
#'   \code{geo}. Defaults \code{"total_exp"}, \code{"population"}.
#' @param per Rate multiplier for \code{bym2_isr}. Default \code{100000}. Affects
#'   \code{bym2_isr} only.
#' @param threshold Relative-risk threshold for the exceedance probability.
#'   Default \code{1.10} (a +10\% excess in mortality). Affects \code{bym2_exceed}
#'   only.
#' @return \code{geo} with the three \code{bym2_*} columns added and the
#'   threshold recorded as an attribute.
#' @examples
#' \dontrun{
#' geo2 <- augment_bym2(geo, fit)                 # P(RR > 1.10)
#' plot_exceedance_map(geo2, value = "bym2_isr")
#' }
#' @importFrom rstan extract
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @export
augment_bym2 <- function(geo, fit,
                         exp_col   = "total_exp",
                         pop_col   = "population",
                         per       = 100000,
                         threshold = 1.10) {

  # geostan stores the posterior expected COUNTS under the Stan parameter
  # "fitted"; pull them straight from the underlying stanfit. (geostan's S3
  # fitted() method is registered but not exported, so geostan::fitted() errors,
  # and posterior_epred() has no geostan method - hence the direct extract.)
  mu <- rstan::extract(fit$stanfit, pars = "fitted")$fitted   # [draws, n_areas]
  E  <- geo[[exp_col]]

  # relative-risk draws = fitted count / expected count
  rr <- sweep(mu, MARGIN = 2, STATS = E, FUN = "/")           # [draws, n_areas]

  out <- geo |>
    dplyr::mutate(
      bym2_rr     = colMeans(rr),
      bym2_isr    = colMeans(mu) / .data[[pop_col]] * per,
      bym2_exceed = colMeans(rr > threshold)
    )

  # record the threshold so the map can label itself consistently
  attr(out, "bym2_exceed_threshold") <- threshold
  out
}


#' Diagnostic report for a fitted BYM2 model
#'
#' One call that gathers the checks you actually need before trusting a BYM2 fit:
#' the raw \code{rstan::check_hmc_diagnostics()} messages, plus quantified
#' divergences, max-treedepth saturation, per-chain E-BFMI (energy), global
#' convergence (worst R-hat and effective sample size across *all* monitored
#' quantities), a focal-parameter table (intercept, \code{spatial_scale},
#' \code{rho}, plus any you name), the BYM2 mixing interpretation of \code{rho},
#' and WAIC where available. Each check prints with an \code{OK}/\code{WARN}
#' verdict and everything is also returned invisibly for programmatic use.
#'
#' Why these in particular for BYM2:
#' \itemize{
#'   \item \strong{Divergences / treedepth / E-BFMI} catch the funnel pathology
#'     that BYM2 is prone to; non-zero divergences or low E-BFMI (< ~0.3) mean the
#'     posterior - and especially the tail-based \code{bym2_exceed} - is suspect.
#'   \item \strong{spatial_scale} is the global smoothing parameter; check its
#'     posterior is well identified (good R-hat/ESS) and look at its magnitude.
#'   \item \strong{rho} is the share of random-effect variance that is spatially
#'     structured; its posterior mean is reported as a percentage.
#' }
#'
#' @param fit A fitted model from \code{\link{fit_bym2}}.
#' @param pars Focal parameters to tabulate. Default
#'   \code{c("intercept", "spatial_scale", "rho")}; names not present are
#'   skipped (e.g. add \code{"beta[1]"} when you have covariates).
#' @param rhat_warn,ess_warn,ebfmi_warn Thresholds that trigger a WARN verdict.
#'   Defaults \code{1.01}, \code{400}, \code{0.30}.
#' @param probs Quantiles for the focal-parameter table. Default
#'   \code{c(0.025, 0.5, 0.975)}.
#' @param digits Rounding for the focal-parameter table. Default \code{3}.
#' @param print If \code{TRUE} (default) print the report.
#' @return Invisibly, a list with the divergence/treedepth/E-BFMI counts, global
#'   \code{rhat_max}/\code{neff_min} and how many parameters breach the
#'   thresholds, \code{rho_mean}, \code{waic}, the focal-parameter table
#'   (\code{focal_summary}) and the full \code{rstan::summary} matrix
#'   (\code{full_summary}).
#' @examples
#' \dontrun{
#' d <- check_bym2_fit(fit)
#' d$rhat_max
#' d$focal_summary
#' }
#' @importFrom rstan check_hmc_diagnostics get_sampler_params
#' @export
check_bym2_fit <- function(fit,
                           pars       = c("intercept", "spatial_scale", "rho"),
                           rhat_warn  = 1.01,
                           ess_warn   = 400,
                           ebfmi_warn = 0.30,
                           probs      = c(0.025, 0.5, 0.975),
                           digits     = 3,
                           print      = TRUE) {

  sf <- fit$stanfit

  ## --- raw rstan diagnostic messages (side effect) -------------------------
  if (print) {
    cat("rstan::check_hmc_diagnostics():\n")
    rstan::check_hmc_diagnostics(sf)
    cat("\n")
  }

  ## --- sampler-level diagnostics -------------------------------------------
  sp       <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  n_chains <- length(sp)
  n_perchn <- nrow(sp[[1]])
  n_draws  <- n_chains * n_perchn

  n_div   <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  pct_div <- 100 * n_div / n_draws

  max_td  <- tryCatch(sf@stan_args[[1]]$control$max_treedepth,
                      error = function(e) NULL)
  if (is.null(max_td)) max_td <- 10L
  n_td    <- sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= max_td),
                        numeric(1)))
  pct_td  <- 100 * n_td / n_draws

  # E-BFMI per chain (standard estimator from the energy diagnostic)
  ebfmi <- vapply(sp, function(x) {
    e <- x[, "energy__"]
    sum(diff(e)^2) / sum((e - mean(e))^2)
  }, numeric(1))
  low_ebfmi <- which(ebfmi < ebfmi_warn)

  ## --- convergence over ALL monitored quantities ---------------------------
  full     <- rstan::summary(sf, probs = probs)$summary
  rhat_all <- full[, "Rhat"]
  neff_all <- full[, "n_eff"]
  rhat_max <- max(rhat_all, na.rm = TRUE)
  neff_min <- min(neff_all, na.rm = TRUE)
  n_bad_rhat <- sum(rhat_all > rhat_warn, na.rm = TRUE)
  n_bad_ess  <- sum(neff_all < ess_warn,  na.rm = TRUE)
  worst_rhat_par <- rownames(full)[which.max(rhat_all)]
  worst_ess_par  <- rownames(full)[which.min(neff_all)]

  ## --- focal-parameter table -----------------------------------------------
  key     <- pars[pars %in% rownames(full)]
  key_tab <- round(full[key, , drop = FALSE], digits)

  ## --- BYM2 mixing interpretation ------------------------------------------
  rho_mean <- if ("rho" %in% rownames(full)) full["rho", "mean"] else NA_real_

  ## --- WAIC (guarded; optional) --------------------------------------------
  waic_val <- tryCatch({
    w <- geostan::waic(fit)
    if (is.numeric(w) && "WAIC" %in% names(w))      as.numeric(w[["WAIC"]])
    else if (is.list(w) && !is.null(w$WAIC))        as.numeric(w$WAIC)
    else NA_real_
  }, error = function(e) NA_real_)

  verdict <- function(flag) if (isTRUE(flag)) "OK  " else "WARN"

  if (print) {
    line <- strrep("-", 66)
    cat(line, "\n", "BYM2 model diagnostics\n", line, "\n", sep = "")
    cat(sprintf("Chains: %d   Post-warmup draws/chain: %d   Total: %d\n\n",
                n_chains, n_perchn, n_draws))

    cat("HMC diagnostics\n")
    cat(sprintf("  [%s] Divergences   : %d (%.2f%%)\n",
                verdict(n_div == 0), n_div, pct_div))
    cat(sprintf("  [%s] Treedepth     : %d hit limit %d (%.2f%%)\n",
                verdict(n_td == 0), n_td, max_td, pct_td))
    cat(sprintf("  [%s] E-BFMI (min)  : %.3f%s\n",
                verdict(length(low_ebfmi) == 0), min(ebfmi),
                if (length(low_ebfmi))
                  sprintf("   [chains %s < %.2f]",
                          paste(low_ebfmi, collapse = ","), ebfmi_warn) else ""))
    cat("\n")

    cat("Convergence (all monitored quantities)\n")
    cat(sprintf("  [%s] Max R-hat     : %.4f  (%s)\n",
                verdict(rhat_max <= rhat_warn), rhat_max, worst_rhat_par))
    cat(sprintf("  [%s] Min ESS       : %.0f  (%s)\n",
                verdict(neff_min >= ess_warn), neff_min, worst_ess_par))
    cat(sprintf("       %d param(s) R-hat > %.2f ; %d param(s) ESS < %d\n\n",
                n_bad_rhat, rhat_warn, n_bad_ess, ess_warn))

    cat("Focal parameters\n")
    print(key_tab)
    cat("\n")

    if (!is.na(rho_mean)) {
      cat(sprintf("BYM2 mixing: posterior mean rho = %.3f\n", rho_mean))
      cat(sprintf("  -> ~%.0f%% of the random-effect variance is spatially structured\n",
                  100 * rho_mean))
      cat("     (rho -> 1 neighbour-driven smoothing; rho -> 0 unstructured)\n\n")
    }

    if (!is.na(waic_val)) cat(sprintf("WAIC: %.1f\n\n", waic_val))
    cat(line, "\n", sep = "")
  }

  invisible(list(
    n_chains       = n_chains,
    n_draws        = n_draws,
    divergences    = n_div,
    pct_divergent  = pct_div,
    max_treedepth  = max_td,
    treedepth_hits = n_td,
    ebfmi          = ebfmi,
    rhat_max       = rhat_max,
    neff_min       = neff_min,
    n_rhat_above   = n_bad_rhat,
    n_ess_below    = n_bad_ess,
    rho_mean       = rho_mean,
    waic           = waic_val,
    focal_summary  = key_tab,
    full_summary   = full
  ))
}


#' Side-by-side comparison of two or more fitted BYM2 models
#'
#' Pulls the focal parameter posteriors (intercept, any \code{beta}, the
#' \code{spatial_scale} and the BYM2 mixing parameter \code{rho}) from a set of
#' \code{\link{fit_bym2}} fits into one tidy table, and runs a proper LOO
#' comparison so you get elpd differences \emph{with} their standard errors -
#' which is what you need to judge whether one model actually predicts better,
#' rather than eyeballing raw WAIC numbers.
#'
#' The pointwise log-likelihood is reconstructed from each model's posterior
#' fitted counts under the Poisson likelihood
#' (\eqn{y_i \sim \mathrm{Poisson}(\mu_i)}, with \eqn{\mu_i} = the \code{fitted}
#' parameter, i.e. expected deaths including the offset), so this does not rely
#' on geostan exposing a \code{log_lik} object. Chains are kept separate
#' (\code{permuted = FALSE}) so \code{loo} can compute \code{r_eff} correctly.
#'
#' Two cautions on interpretation. (1) The fits are separate models, so the
#' difference in a posterior mean (e.g. \code{spatial_scale} base vs IVSM) has no
#' clean joint test - read each estimate with its own sd, don't difference them.
#' (2) The decision-relevant quantity is \code{elpd_diff / se_diff};
#' \code{loo_compare} orders the models best-first and reports every other
#' model's difference from the best. A common rule of thumb treats |ratio| < 2
#' as indistinguishable and > 4 as clear.
#'
#' @param fits A \emph{named} list of two or more \code{\link{fit_bym2}}
#'   objects, e.g. \code{list(base = m0, ivsm = m1, di = m2)}. The names are
#'   used as labels.
#' @param data A list the SAME length as \code{fits} and in the SAME order, each
#'   element the \code{sf}/data frame its fit was built from, containing
#'   \code{obs_col}.
#' @param obs_col Name of the observed-count column. Default \code{"total_obs"}.
#' @param pars Focal parameters to tabulate; names absent from a fit are skipped
#'   (so \code{"beta[1]"} appears only for covariate models). Default
#'   \code{c("intercept", "beta[1]", "spatial_scale", "rho")}.
#' @param param_labels Optional named character vector to relabel raw Stan names
#'   in the table, e.g. \code{c("beta[1]" = "ivsm_z")}. Default \code{NULL}.
#' @param probs Quantiles for the parameter table. Default
#'   \code{c(0.025, 0.5, 0.975)}.
#' @return An object of class \code{"bym2_comparison"}: a list with \code{params}
#'   (tidy parameter table), \code{loo} (the per-model \code{loo} objects) and
#'   \code{loo_compare} (the \code{loo::loo_compare} matrix). Has a print method.
#' @note Memory: the log-likelihood array is draws x areas in size; for several
#'   thousand comuni this is a few hundred MB per model. Fine on the machine that
#'   ran the 4-chain fits, but it scales with the number of models compared.
#' @examples
#' \dontrun{
#' cmp <- compare_bym2(
#'   fits = list(base = model_base, ivsm = model_ivsm, di = model_di),
#'   data = list(base = smr_geo,    ivsm = smr_geo_ivsm, di = smr_geo_di),
#'   param_labels = c("beta[1]" = "index_z")
#' )
#' cmp                      # prints the side-by-side report
#' cmp$loo_compare          # the raw elpd-difference matrix
#' }
#' @importFrom rstan summary extract
#' @importFrom stats dpois
#' @importFrom loo relative_eff loo loo_compare
#' @export
compare_bym2 <- function(fits,
                         data,
                         obs_col      = "total_obs",
                         pars         = c("intercept", "beta[1]",
                                          "spatial_scale", "rho"),
                         param_labels = NULL,
                         probs        = c(0.025, 0.5, 0.975)) {

  stopifnot(length(fits) >= 2L, length(fits) == length(data))
  labels <- names(fits)
  if (is.null(labels)) labels <- paste0("model_", seq_along(fits))

  ## ---- focal-parameter table ----------------------------------------------
  one_param_tab <- function(fit, label) {
    avail <- rownames(rstan::summary(fit$stanfit)$summary)
    keep  <- pars[pars %in% avail]
    s <- rstan::summary(fit$stanfit, pars = keep, probs = probs)$summary
    nm <- rownames(s)
    if (!is.null(param_labels)) {
      hit <- nm %in% names(param_labels)
      nm[hit] <- param_labels[nm[hit]]
    }
    data.frame(
      model     = label,
      parameter = nm,
      mean      = s[, "mean"],
      sd        = s[, "sd"],
      q2.5      = s[, "2.5%"],
      q97.5     = s[, "97.5%"],
      rhat      = s[, "Rhat"],
      n_eff     = s[, "n_eff"],
      row.names = NULL,
      check.names = FALSE
    )
  }
  par_tab <- do.call(rbind, Map(one_param_tab, fits, labels))

  ## ---- pointwise log-likelihood -> LOO ------------------------------------
  loglik_array <- function(fit, y) {
    mu <- rstan::extract(fit$stanfit, pars = "fitted", permuted = FALSE)
    d  <- dim(mu); n_it <- d[1]; n_ch <- d[2]; N <- d[3]
    ll <- array(0, dim = c(n_it, n_ch, N))
    for (i in seq_len(N)) {
      ll[, , i] <- matrix(stats::dpois(y[i], as.vector(mu[, , i]), log = TRUE),
                          nrow = n_it, ncol = n_ch)
    }
    ll
  }

  loo_one <- function(fit, dat) {
    y  <- as.numeric(dat[[obs_col]])
    ll <- loglik_array(fit, y)
    r  <- loo::relative_eff(exp(ll))
    loo::loo(ll, r_eff = r)
  }
  loo_list <- Map(loo_one, fits, data)
  names(loo_list) <- labels

  comp <- loo::loo_compare(loo_list)

  out <- list(
    labels      = labels,
    params      = par_tab,
    loo         = loo_list,
    loo_compare = comp
  )
  class(out) <- "bym2_comparison"
  out
}


#' @rdname compare_bym2
#' @param x A \code{"bym2_comparison"} object.
#' @param digits Rounding for the printed tables. Default \code{3}.
#' @param ... Unused.
#' @export
print.bym2_comparison <- function(x, digits = 3, ...) {
  line <- strrep("-", 64)
  cat(line, "\n", "BYM2 comparison: ", paste(x$labels, collapse = " vs "),
      "\n", line, "\n", sep = "")

  cat("Focal parameters\n")
  p <- x$params
  num <- vapply(p, is.numeric, logical(1))
  p[num] <- lapply(p[num], function(z) round(z, digits))
  print(p, row.names = FALSE)

  cat("\nLOO comparison (elpd_loo: higher is better; looic ~ WAIC scale)\n")
  print(round(x$loo_compare, digits))

  # loo_compare is ordered best-first; row 1 is the reference, rows 2..n hold
  # each other model's elpd difference from the best, with its se.
  verdict_for <- function(ratio) {
    if (is.na(ratio)) "se = 0"
    else if (ratio < 2) "indistinguishable (|diff| < 2 SE)"
    else if (ratio < 4) "weak (2-4 SE)"
    else "clear (> 4 SE)"
  }

  best   <- rownames(x$loo_compare)[1]
  others <- rownames(x$loo_compare)[-1]
  cat(sprintf("\nBest: %s. Each other model vs best:\n", best))
  for (m in others) {
    ed    <- x$loo_compare[m, "elpd_diff"]
    se    <- x$loo_compare[m, "se_diff"]
    ratio <- if (se > 0) abs(ed) / se else NA_real_
    cat(sprintf("  %-14s elpd_diff = %8.2f (se = %6.2f) | |diff|/se = %4s -> %s\n",
                m, ed, se,
                if (is.na(ratio)) "NA" else sprintf("%.1f", ratio),
                verdict_for(ratio)))
  }
  invisible(x)
}

#' Fit one BYM2 model per mechanism (or per any set of obs/exp column pairs)
#'
#' Loops over the per-mechanism observed/expected column pairs on \code{geo} and
#' fits a separate \code{\link{fit_bym2}} model to each, reusing the shared
#' adjacency matrix \code{C} and \code{scale_factor} (both depend only on the
#' graph, so they are computed once and passed in). Each mechanism is smoothed on
#' its own indirectly standardised expected counts, so a smoothed RR of 1 means
#' "matches the region-wide age-sex expectation for that mechanism".
#' geostan's Poisson response must be integer; split-cause weights (0.5) make
#' M_*_obs fractional, so it's rounded to the nearest integer and
#' shift is reported.
#'
#' Mechanism "stems" are discovered from the \code{_obs} columns matching
#' \code{obs_pattern}; for each stem \code{S} the model is
#' \code{S_obs ~ offset(log(S_exp))}. A stem whose \code{_exp} column contains any
#' non-positive or \code{NA} value is fatal (the log-offset would be undefined) -
#' this guards against an upstream \code{values_fill = 0} silently injecting
#' \code{log(0)} into the offset for death-free comuni.
#'
#' @param geo An \code{sf}, one row per comune, in the SAME order as \code{C},
#'   carrying \code{<stem>_obs} and \code{<stem>_exp} columns.
#' @param C Binary adjacency matrix from \code{\link{build_adjacency}}.
#' @param scale_factor BYM2 scaling vector from \code{\link{compute_scale_factor}}
#'   (graph-only, shared across all mechanisms).
#' @param obs_pattern Regex selecting the observed-count columns. Default
#'   \code{"^M_.*_obs$"} (the seven mechanisms).
#' @param chains,iter,control Passed to \code{\link{fit_bym2}}.
#' @param ... Further arguments forwarded to \code{\link{fit_bym2}}.
#' @return A named list of fitted models, one per mechanism stem.
#' @examples
#' \dontrun{
#' fits <- fit_bym2_mechanisms(smr_geo, C_matrix, scale_factor)
#' check_bym2_fit(fits[["M_screening"]])
#' }
#' @importFrom stats as.formula
#' @export
fit_bym2_mechanisms <- function(geo, C, scale_factor,
                                obs_pattern = "^M_.*_obs$",
                                chains  = 4,
                                iter    = 4000,
                                control = list(adapt_delta = 0.97,
                                               max_treedepth = 12),
                                ...) {

  obs_cols <- grep(obs_pattern, names(geo), value = TRUE)
  if (!length(obs_cols)) {
    stop("No columns matched obs_pattern = '", obs_pattern, "'. ",
         "Check that preprocess_smr() emits per-mechanism `_obs`/`_exp` columns.",
         call. = FALSE)
  }
  stems <- sub("_obs$", "", obs_cols)

  fits <- vector("list", length(stems))
  names(fits) <- stems

  for (stem in stems) {
    oc <- paste0(stem, "_obs")
    ec <- paste0(stem, "_exp")
    if (!ec %in% names(geo)) {
      warning("No expected column '", ec, "' to match '", oc, "'; skipping ",
              stem, ".", call. = FALSE)
      next
    }

    # the offset trap: expected must be strictly positive for every comune
    E <- geo[[ec]]
    bad <- is.na(E) | E <= 0
    if (any(bad)) {
      stop(sprintf(
        "Mechanism '%s': %d comuni have expected <= 0 or NA, so log(expected) is undefined.\n  Either this mechanism is too sparse to smooth, or `%s` was filled with 0 upstream\n  (a death-free comune still has POSITIVE expected - it must not be zero-filled).",
        stem, sum(bad), ec), call. = FALSE)
    }

    # geostan's Poisson response must be integer; split-cause weights (0.5) make
    # M_*_obs fractional. Round to the nearest integer and report the shift.
    y_raw <- geo[[oc]]
    y_int <- round(y_raw)
    shift <- sum(abs(y_raw - y_int))
    if (shift > 0) {
      message(sprintf("  %s: rounded observed to integer (total absolute shift %.1f deaths over %d comuni).",
                      stem, shift, sum(y_raw != y_int)))
    }
    geo[[oc]] <- y_int

    f <- stats::as.formula(sprintf("%s ~ offset(log(%s))", oc, ec))
    message("Fitting BYM2 for ", stem, " ...")
    fits[[stem]] <- fit_bym2(
      geo, C,
      formula      = f,
      scale_factor = scale_factor,
      obs_col      = oc,
      exp_col      = ec,
      chains       = chains,
      iter         = iter,
      control      = control,
      ...
    )
  }

  fits[!vapply(fits, is.null, logical(1))]
}


#' Collect per-mechanism smoothed RR (and exceedance) into one wide sf
#'
#' Augments \code{geo} with one smoothed-relative-risk column per mechanism,
#' named \code{<stem><out_suffix>} (default \code{<stem>_bym2}), plus a matching
#' exceedance column \code{<stem><out_suffix>_exc} = P(RR > \code{threshold}).
#' The result is shaped exactly like the input to \code{\link{plot_smr_facets}},
#' so the faceted small-multiples map needs no new plotting code - just point its
#' \code{cols} at the \code{_bym2} columns.
#'
#' @param geo The same \code{sf} passed to \code{\link{fit_bym2_mechanisms}}.
#' @param fits The named list returned by \code{\link{fit_bym2_mechanisms}}.
#' @param threshold Relative-risk threshold for the exceedance columns. Default
#'   \code{1.10}.
#' @param out_suffix Suffix for the smoothed-RR columns. Default \code{"_bym2"}.
#' @return \code{geo} with \code{<stem>_bym2} and \code{<stem>_bym2_exc} columns
#'   added, and the threshold recorded as an attribute.
#' @examples
#' \dontrun{
#' smr_geo_mech_bym2 <- augment_bym2_mechanisms(smr_geo, fits)
#' plot_smr_facets(smr_geo_mech_bym2,
#'                 cols   = dplyr::matches("^M_.*_bym2$"),
#'                 breaks = c(-Inf, 0.90, 0.95, 1.05, 1.10, Inf),
#'                 strip_suffix = "_bym2$")
#' }
#' @export
augment_bym2_mechanisms <- function(geo, fits,
                                    threshold  = 1.10,
                                    out_suffix = "_bym2") {
  out <- geo
  for (stem in names(fits)) {
    ec  <- paste0(stem, "_exp")
    aug <- augment_bym2(geo, fits[[stem]], exp_col = ec, threshold = threshold)
    out[[paste0(stem, out_suffix)]]          <- aug[["bym2_rr"]]
    out[[paste0(stem, out_suffix, "_exc")]]  <- aug[["bym2_exceed"]]
  }
  attr(out, "bym2_exceed_threshold") <- threshold
  out
}
