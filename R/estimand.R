# The estimand and the spatial diagnostics -------------------------------------
#
# Everything in Results "Residual excess", plus the two Moran's I values the
# methods promise (one before fitting, one on residuals).
#
# THE QUANTITY
#
# The methods define the primary estimand as "the number of deaths occurring
# beyond what the area's age and sex structure, deprivation, air pollution and
# primary care capacity jointly predict". Age and sex are in the offset;
# deprivation, pollution and primary care are the covariates. So the prediction
# to compare against is
#
#     pred_i = E_i * exp(alpha + x_i' beta)
#
# WITHOUT the random effects. That is the whole point: the random effects are
# the excess, so including them would make the estimand identically zero.
#
# Two versions are returned. The model-smoothed excess, mu_i - pred_i, is the
# coherent Bayesian quantity and is what should be reported: it borrows strength
# across neighbours in the same way the risk surface does. The raw excess,
# y_i - pred_i, is the arithmetic one; it sums to something very close to the
# smoothed version but is noisier per area. Reporting both, and noting that they
# agree, is cheap reassurance.


#' Parameter names of a stanfit, cheaply
#'
#' `rstan::summary()` computes R-hat and effective sample size for every
#' monitored quantity just to hand back its row names. With 279 areas that is
#' several thousand parameters - `fitted[]`, `log_lik[]`, `phi[]`, `theta[]` -
#' recomputed once per model per call. The flat names are already stored on the
#' object.
#'
#' @param fit A fitted geostan model.
#' @return Character vector of parameter names.
#' @noRd
.stan_par_names <- function(fit) {
  nm <- tryCatch(fit$stanfit@sim$fnames_oi, error = function(e) NULL)
  if (is.null(nm) || !length(nm)) {
    nm <- rownames(rstan::summary(fit$stanfit)$summary)
  }
  nm
}


#' Resolve the names geostan gave the coefficient parameters
#'
#' geostan does not name the covariate coefficients the same way across
#' versions and engines: they may appear in the stanfit as `beta[1]`, `beta[2]`
#' ... or under the covariate's own column name. Assuming one and grepping for
#' it produces an empty coefficient table rather than an error, which then
#' surfaces a long way downstream as a ggplot faceting failure.
#'
#' Resolve by looking at what is actually there, and fail with the available
#' names if neither convention matches.
#'
#' @param fit A fitted geostan model.
#' @param covariates Character vector of covariate names, in the order the
#'   model was fitted.
#'
#' @return Character vector of stanfit parameter names, same length and order
#'   as `covariates`. `character(0)` when the model has no covariates.
#' @noRd
.beta_pars <- function(fit, covariates) {

  if (!length(covariates)) return(character(0))

  avail <- .stan_par_names(fit)
  k     <- length(covariates)

  # 1. named after the covariates themselves
  if (all(covariates %in% avail)) return(covariates)

  # 2. the indexed vector, in order
  idx <- paste0("beta[", seq_len(k), "]")
  if (all(idx %in% avail)) return(idx)

  # 3. any beta[j] block of the right length, whatever it is indexed by
  hit <- grep("^beta\\[[0-9]+\\]$", avail, value = TRUE)
  if (length(hit) == k) return(hit)

  # 4. some geostan versions prefix the covariate name
  pref <- paste0("beta_", covariates)
  if (all(pref %in% avail)) return(pref)

  # 5. gamma, used by a few of the geostan model variants
  gam <- paste0("gamma[", seq_len(k), "]")
  if (all(gam %in% avail)) return(gam)

  # Nothing matched. This is almost always a sign that the model was fitted
  # WITHOUT the covariates - an intercept-only fit has no coefficient block to
  # find - so report the formula alongside the parameter names.
  spec <- attr(fit, "atlaspm_spec")
  interesting <- grep("^(fitted|log_lik|phi|theta|alpha_re|yrep)\\[",
                      avail, value = TRUE, invert = TRUE)

  stop("Cannot locate coefficient parameters for: ",
       paste(covariates, collapse = ", "), "\n",
       "  Tried: their own names, beta[1..", k, "], any beta[j] block, ",
       "beta_<name>, gamma[1..", k, "].\n",
       "  Model parameters present: ",
       paste(utils::head(interesting, 30), collapse = ", "),
       if (length(interesting) > 30) ", ..." else "", "\n",
       if (!is.null(spec))
         paste0("  Formula recorded on the fit: ", spec$formula, "\n") else "",
       "  If no coefficient block appears above, the model was fitted ",
       "intercept-only and the covariates never reached geostan.",
       call. = FALSE)
}


#' Decompose a fitted model into its fixed and random parts
#'
#' Internal workhorse for [residual_excess()], [variance_decomposition()] and
#' [moran_test()]. Returns posterior draws of the fitted counts, of the
#' covariate-only prediction, and of the implied area random effect.
#'
#' @param fit A model from [fit_model()] or [fit_bym2()].
#' @param geo The `sf` the model was fitted to, same row order.
#' @param covariates Character vector of covariate columns in `geo`. Taken from
#'   the fit's specification when `NULL`.
#' @param exp_col Expected-count column.
#'
#' @return A list of `[draws x areas]` matrices: `mu` (fitted counts),
#'   `pred` (covariate-only prediction), `re` (log(mu) - log(pred)), plus the
#'   observed vector `y` and the expected vector `E`.
#' @noRd
.model_parts <- function(fit, geo, covariates = NULL, exp_col = NULL) {

  sp <- tryCatch(model_spec(fit), error = function(e) NULL)
  if (is.null(covariates)) covariates <- if (is.null(sp)) character(0) else sp$covariates
  if (is.null(exp_col))    exp_col    <- if (is.null(sp)) "total_exp" else sp$exp_col
  obs_col <- if (is.null(sp)) "total_obs" else sp$obs_col

  require_cols(geo, c(exp_col, obs_col, covariates), "geo")

  mu <- rstan::extract(fit$stanfit, pars = "fitted")$fitted   # [draws, areas]
  if (ncol(mu) != nrow(geo)) {
    stop("The fit has ", ncol(mu), " areas but `geo` has ", nrow(geo),
         " rows. These must be the same object, in the same order.",
         call. = FALSE)
  }

  E <- as.numeric(geo[[exp_col]])
  y <- round_half_up(as.numeric(geo[[obs_col]]))

  alpha <- as.numeric(rstan::extract(fit$stanfit, pars = "intercept")$intercept)

  # log(pred) = alpha + X beta, broadcast over draws
  eta <- matrix(alpha, nrow = length(alpha), ncol = nrow(geo))

  if (length(covariates)) {
    X <- as.matrix(sf::st_drop_geometry(geo)[, covariates, drop = FALSE])
    if (anyNA(X)) {
      stop("Covariate matrix contains NA for ", sum(!stats::complete.cases(X)),
           " area(s): ", paste(covariates, collapse = ", "),
           ". The model cannot have been fitted on these rows.", call. = FALSE)
    }
    beta <- as.matrix(fit$stanfit, pars = .beta_pars(fit, covariates))
    if (is.null(dim(beta))) beta <- matrix(beta, ncol = 1L)
    if (ncol(beta) != ncol(X)) {
      stop("The fit has ", ncol(beta), " coefficient(s) but ", ncol(X),
           " covariate column(s) were supplied (",
           paste(covariates, collapse = ", "),
           "). Pass `covariates` in the order the model was fitted.",
           call. = FALSE)
    }
    eta <- eta + beta %*% t(X)
  }

  pred <- sweep(exp(eta), MARGIN = 2, STATS = E, FUN = "*")

  list(mu = mu, pred = pred, re = log(mu) - log(pred), y = y, E = E,
       covariates = covariates)
}


#' Residual excess of avoidable deaths
#'
#' The primary estimand. For each area, the number of deaths occurring beyond
#' what its age-sex structure (in the offset) and the model's covariates jointly
#' predict, with a posterior credible interval; and the total across the study
#' area, annualised.
#'
#' Also returns the residual relative risk - the exponentiated area random
#' effect - and the posterior probability that it exceeds `threshold`. This is
#' the exceedance quantity Results reports at a probability cutoff of 0.80,
#' and it differs from [augment_bym2()]'s `bym2_exceed`: that one is the
#' probability that *total* relative risk exceeds the threshold, this one is the
#' probability that risk exceeds the threshold *after* the covariates have had
#' their say.
#'
#' @param fit A fitted model, normally M5.
#' @param geo The `sf` it was fitted to.
#' @param covariates Covariate columns. Taken from the fit when `NULL`.
#' @param exp_col Expected-count column. Taken from the fit when `NULL`.
#' @param n_years Number of years the death counts span, used to annualise.
#'   Default `3`.
#' @param threshold Residual relative-risk threshold for the exceedance
#'   probability. Default `1.10`.
#' @param prob_cutoff Probability above which an area is flagged. Default
#'   `0.80` (Richardson et al. 2004).
#' @param probs Credible-interval bounds. Default `c(0.025, 0.975)`.
#'
#' @return A list:
#'   \describe{
#'     \item{`total`}{One-row tibble: annualised excess deaths with credible
#'       interval, as a count and as a percentage of all avoidable deaths.}
#'     \item{`per_area`}{One row per area: excess, credible interval, residual
#'       RR, exceedance probability and the flag.}
#'     \item{`n_flagged`}{Number of areas above `prob_cutoff`.}
#'     \item{`raw_total`}{The arithmetic version, `sum(y) - sum(pred)`, for
#'       comparison.}
#'   }
#'
#' @examples
#' \dontrun{
#' ex <- residual_excess(fit_M5, smr_geo_full, n_years = 3)
#' ex$total
#' subset(ex$per_area, flagged)
#' }
#' @seealso [variance_decomposition()], [augment_bym2()]
#' @importFrom rstan extract
#' @importFrom rlang .data
#' @export
residual_excess <- function(fit, geo,
                            covariates  = NULL,
                            exp_col     = NULL,
                            n_years     = 3,
                            threshold   = 1.10,
                            prob_cutoff = 0.80,
                            probs       = c(0.025, 0.975)) {

  p <- .model_parts(fit, geo, covariates, exp_col)

  excess_draws <- (p$mu - p$pred) / n_years          # [draws, areas]
  total_draws  <- rowSums(excess_draws)              # [draws]
  rr_draws     <- exp(p$re)

  total_deaths_yr <- sum(p$y) / n_years
  q <- function(x) stats::quantile(x, probs = probs, names = FALSE)

  total <- tibble::tibble(
    excess_per_year = mean(total_draws),
    ci_low          = q(total_draws)[1],
    ci_high         = q(total_draws)[2],
    total_per_year  = total_deaths_yr,
    pct_of_total    = 100 * mean(total_draws) / total_deaths_yr,
    pct_ci_low      = 100 * q(total_draws)[1] / total_deaths_yr,
    pct_ci_high     = 100 * q(total_draws)[2] / total_deaths_yr,
    n_years         = n_years
  )

  area_key <- if ("area" %in% names(geo)) geo[["area"]] else seq_len(nrow(geo))

  per_area <- tibble::tibble(
    area        = area_key,
    observed    = p$y,
    expected    = p$E,
    predicted   = colMeans(p$pred),
    excess      = colMeans(excess_draws),
    excess_low  = apply(excess_draws, 2, function(x) q(x)[1]),
    excess_high = apply(excess_draws, 2, function(x) q(x)[2]),
    rr_resid    = colMeans(rr_draws),
    rr_low      = apply(rr_draws, 2, function(x) q(x)[1]),
    rr_high     = apply(rr_draws, 2, function(x) q(x)[2]),
    p_exceed    = colMeans(rr_draws > threshold)
  )
  per_area[["flagged"]] <- per_area[["p_exceed"]] > prob_cutoff

  list(
    total       = total,
    per_area    = per_area,
    n_flagged   = sum(per_area[["flagged"]]),
    raw_total   = (sum(p$y) - sum(colMeans(p$pred))) / n_years,
    threshold   = threshold,
    prob_cutoff = prob_cutoff
  )
}


#' How much between-area variation the covariates absorb
#'
#' Results reports that "adding the covariates to the unadjusted model reduced
#' between-area variance by \[N\]%". This computes that number.
#'
#' The comparison is between the *residual* variance in each model: the variance
#' across areas of the log area random effect, evaluated per posterior draw and
#' then summarised. Doing it per draw rather than on posterior means matters -
#' the variance of the means is not the mean of the variances, and the former
#' understates the spread by an amount that grows with how sparse the counts
#' are, which is exactly the regime this study is in.
#'
#' @param fit_null The unadjusted model (M3).
#' @param fit_adj The adjusted model (M5).
#' @param geo The `sf` both were fitted to.
#' @param labels Names for the two models in the output.
#' @param probs Credible-interval bounds.
#'
#' @return A list with `by_model` (one row per model: residual SD and variance
#'   with credible intervals, plus the BYM2 mixing parameter where available)
#'   and `reduction` (percentage reduction in residual variance, with a
#'   credible interval derived per draw).
#'
#' @examples
#' \dontrun{
#' vd <- variance_decomposition(fit_M3, fit_M5, smr_geo_full)
#' vd$reduction
#' }
#' @export
variance_decomposition <- function(fit_null, fit_adj, geo,
                                   labels = c("M3", "M5"),
                                   probs  = c(0.025, 0.975)) {

  p0 <- .model_parts(fit_null, geo)
  p1 <- .model_parts(fit_adj,  geo)

  v0 <- apply(p0$re, 1, stats::var)
  v1 <- apply(p1$re, 1, stats::var)

  n <- min(length(v0), length(v1))
  if (length(v0) != length(v1)) {
    warning("The two fits have different draw counts (", length(v0), " vs ",
            length(v1), "); the reduction is computed on the first ", n,
            " of each. The draws are not paired, so read its interval as ",
            "indicative.", call. = FALSE)
  }
  red <- 100 * (1 - v1[seq_len(n)] / v0[seq_len(n)])

  q <- function(x) stats::quantile(x, probs = probs, names = FALSE)
  rho_of <- function(f) {
    s <- rstan::summary(f$stanfit)$summary
    if ("rho" %in% rownames(s)) s["rho", "mean"] else NA_real_
  }

  by_model <- tibble::tibble(
    model     = labels,
    resid_var = c(mean(v0), mean(v1)),
    var_low   = c(q(v0)[1], q(v1)[1]),
    var_high  = c(q(v0)[2], q(v1)[2]),
    resid_sd  = c(mean(sqrt(v0)), mean(sqrt(v1))),
    sd_low    = c(q(sqrt(v0))[1], q(sqrt(v1))[1]),
    sd_high   = c(q(sqrt(v0))[2], q(sqrt(v1))[2]),
    rho       = c(rho_of(fit_null), rho_of(fit_adj))
  )

  list(
    by_model  = by_model,
    reduction = tibble::tibble(
      pct_reduction = mean(red),
      ci_low        = q(red)[1],
      ci_high       = q(red)[2]
    )
  )
}


#' Coerce an adjacency matrix into a spatial weights list
#'
#' `geostan::shape2mat(style = "B")` returns a **pattern** sparse matrix
#' (`ngCMatrix`), which records only where the non-zeros are and has no `@x`
#' slot. `spdep::mat2listw()` reaches for that slot and fails with
#' "no slot of name \"x\" for this object of class \"ngCMatrix\"".
#'
#' Densifying is free at this size - 279 areas is a 279x279 matrix - and avoids
#' depending on which Matrix coercion paths are current.
#'
#' @noRd
.listw <- function(C, style = "W") {
  if (!is.matrix(C)) C <- as.matrix(C)
  storage.mode(C) <- "double"
  spdep::mat2listw(C, style = style, zero.policy = TRUE)
}


#' Moran's I on a plain numeric vector
#'
#' The pre-model version the methods call for: computed on the crude
#' standardised ratios, before any spatial term has had a chance to absorb the
#' structure it is meant to detect.
#'
#' @param x Numeric vector, one value per area, in the row order of `C`.
#' @param C Binary adjacency matrix from [build_adjacency()].
#' @param nsim Permutations for the Monte Carlo p-value. Default `9999`.
#' @param style Spatial weights style passed to `spdep::mat2listw()`. Default
#'   `"W"` (row-standardised).
#'
#' @return A one-row tibble: `statistic`, `expectation`, `p_value`, `n_sim`.
#' @examples
#' \dontrun{
#' moran_test_raw(smr_geo$total_smr, C_matrix)
#' }
#' @importFrom spdep mat2listw moran.mc
#' @export
moran_test_raw <- function(x, C, nsim = 9999, style = "W") {

  ok <- is.finite(x)
  if (!all(ok)) {
    stop(sum(!ok), " non-finite value(s) in `x`. Moran's I is undefined on ",
         "them; decide explicitly whether those areas belong in the test.",
         call. = FALSE)
  }

  lw  <- .listw(C, style)
  res <- spdep::moran.mc(x, lw, nsim = nsim, zero.policy = TRUE)

  tibble::tibble(
    statistic   = unname(res$statistic),
    expectation = -1 / (length(x) - 1),
    p_value     = res$p.value,
    n_sim       = nsim
  )
}


#' Moran's I on model residuals
#'
#' The post-model version. Confirms that the spatial term has absorbed the
#' structure that [moran_test_raw()] found in the crude ratios.
#'
#' Two things this does that the ad-hoc version in the old `reports/explore.R`
#' did not. It works on **Pearson residuals**, \eqn{(y - \mu)/\sqrt{\mu}}, not
#' on the smoothed relative risks: smoothed RRs are spatially autocorrelated by
#' construction, so testing them tells you only that the smoother smoothed. And
#' it propagates posterior uncertainty, returning the distribution of I across
#' draws alongside the permutation test on the posterior-mean residual.
#'
#' @param fit A fitted model.
#' @param geo The `sf` it was fitted to.
#' @param C Adjacency matrix.
#' @param nsim Permutations for the p-value on the mean residual.
#' @param n_draws Draws to use for the posterior distribution of I. Default
#'   `500`; the full posterior is rarely worth the time here.
#' @param style Weights style.
#'
#' @return A one-row tibble: `statistic` (on the posterior-mean residual),
#'   `p_value`, and `post_mean` / `post_low` / `post_high` summarising I across
#'   draws.
#' @seealso [moran_test_raw()]
#' @export
moran_test <- function(fit, geo, C, nsim = 9999, n_draws = 500, style = "W") {

  p  <- .model_parts(fit, geo)
  lw <- .listw(C, style)

  pearson <- sweep(-p$mu, MARGIN = 2, STATS = p$y, FUN = "+") / sqrt(p$mu)

  r_mean <- colMeans(pearson)
  res    <- spdep::moran.mc(r_mean, lw, nsim = nsim, zero.policy = TRUE)

  idx <- sample.int(nrow(pearson), min(n_draws, nrow(pearson)))
  I_draws <- vapply(
    idx,
    function(i) unname(spdep::moran(pearson[i, ], lw,
                                    n = length(r_mean),
                                    S0 = spdep::Szero(lw),
                                    zero.policy = TRUE)$I),
    numeric(1)
  )

  tibble::tibble(
    statistic   = unname(res$statistic),
    expectation = -1 / (length(r_mean) - 1),
    p_value     = res$p.value,
    n_sim       = nsim,
    post_mean   = mean(I_draws),
    post_low    = stats::quantile(I_draws, 0.025, names = FALSE),
    post_high   = stats::quantile(I_draws, 0.975, names = FALSE)
  )
}


#' Test a conditional-independence implication on model residuals
#'
#' The DAG implies that primary care capacity is conditionally independent of
#' the outcome given deprivation and location. A testable version: after fitting
#' a model that contains deprivation and the spatial term but *not* primary
#' care, the residuals should be uncorrelated with primary care capacity.
#'
#' This is a *testable implication*, not a test of the DAG. Failing it says the
#' graph is wrong somewhere; passing it says only that this particular
#' implication survived, which is weaker than it sounds and should be reported
#' as such.
#'
#' @param fit The model omitting `var` (normally M4).
#' @param geo The `sf` it was fitted to.
#' @param var Column in `geo` to correlate the residuals against.
#' @param method Correlation method. Default `"spearman"`, which does not
#'   assume the residuals are normal.
#'
#' @return A one-row tibble: `variable`, `estimate`, `p_value`, `method`, `n`.
#' @export
test_conditional_independence <- function(fit, geo, var, method = "spearman") {
  require_cols(geo, var, "geo")

  p       <- .model_parts(fit, geo)
  pearson <- sweep(-p$mu, MARGIN = 2, STATS = p$y, FUN = "+") / sqrt(p$mu)
  r_mean  <- colMeans(pearson)
  x       <- as.numeric(sf::st_drop_geometry(geo)[[var]])

  ok <- is.finite(r_mean) & is.finite(x)
  ct <- suppressWarnings(stats::cor.test(r_mean[ok], x[ok], method = method))

  tibble::tibble(
    variable = var,
    estimate = unname(ct$estimate),
    p_value  = ct$p.value,
    method   = method,
    n        = sum(ok)
  )
}
