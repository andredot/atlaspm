# Attribution: decomposing the excess into its named parts ---------------------
#
# Section 4 says how much excess there is and where. This file says what it is
# made of. Three questions, three entry points:
#
#   attribution_local()   For one area (normally the flagged one): what share of
#                         its excess is deprivation, pollution, primary care,
#                         and the spatial random effect.
#   attribution_grid()    The same decomposition for every area, posterior means
#                         only, so the dominant driver can be mapped.
#   stratum_shares()      What share of avoidable deaths sits in each OECD
#                         cause group and each service function.
#
# THE REFERENCE, AND WHY IT IS NOT THE ONE IN residual_excess()
#
# residual_excess() compares the fitted count against a prediction that already
# contains the covariates, so its estimand is the random effect alone. That is
# the right estimand for "what can the model not explain", and the wrong one for
# "what is this area's excess made of", because it has already spent the
# covariates before the question is asked.
#
# The decomposition here uses a reference area instead: same age and sex
# structure, study-average deprivation, study-average NO2, study-average primary
# care density, and no residual spatial effect. Since the covariates enter the
# model as z-scores, "study average" is exactly x = 0, so that reference is
#
#     B_i = E_i * exp(alpha)
#
# and every component is a departure from it. The two estimands reconcile
# exactly: the ICAR and IID parts of a *sequential* decomposition that adds the
# covariates first sum to residual_excess()'s excess. That sequential figure is
# returned alongside the main one so the reconciliation can be checked rather
# than asserted.
#
# WHY SHAPLEY AND NOT beta * x
#
# The obvious decomposition - contribution of k is beta_k * x_ik - is a
# decomposition of the log relative risk, not of deaths, and the two are not the
# same shape: on the death scale the components interact, because a 20% excess
# from deprivation and a 20% excess from pollution together produce 44%, not
# 40%. Splitting deaths by beta * x therefore does not add up, and splitting
# them sequentially makes the answer depend on the order the components are
# named in.
#
# Shapley values are the unique attribution that is exact (the parts sum to the
# whole), symmetric (two components with identical contributions get identical
# credit) and order-independent. With five components that is 32 subsets, which
# costs nothing. The log-scale share is returned too, in `pct_log`, because it
# is what a reader will reach for to check the arithmetic by hand.
#
# A WARNING ABOUT SIGNS
#
# A component whose local value is better than the study average contributes
# NEGATIVE deaths: it is holding the excess down. Its share is then negative,
# and the positive shares sum to more than 100%. That is the correct reading and
# it must survive into the figure. Clamping contributions at zero before
# plotting them as filled proportions - which is what the feasibility notebook
# did with pmax(0, .) and position = "fill" - deletes the protective factors and
# rescales what is left, so a determinant that was reducing deaths disappears
# and the remaining bars grow to cover for it.


#' Drop geometry if there is any
#'
#' The tabular functions here work on the modelling frame, which is an `sf`, but
#' nothing they do needs the geometry. Accepting a plain data frame too means
#' they can be pointed at an intermediate table - and tested - without an `sf`
#' round trip.
#'
#' @param x An `sf` or a data frame.
#' @return A data frame.
#' @noRd
.plain <- function(x) {
  if (inherits(x, "sf")) sf::st_drop_geometry(x) else as.data.frame(x)
}


#' Posterior draws in one consistent ordering
#'
#' `rstan::extract()` permutes draws and merges chains; `as.matrix()` does not.
#' Both are fine on their own, and mixing them is not: a quantity built from an
#' `extract()`ed `fitted` and an `as.matrix()`ed `beta` pairs each draw of one
#' with an unrelated draw of the other. Marginal summaries survive that - the
#' mean of a permuted sample is the same mean - but anything joint does not, and
#' every quantity in this file is joint.
#'
#' So: one call, one ordering, no exceptions.
#'
#' @param fit A fitted geostan model.
#' @param pars Character vector of parameter names or stems.
#' @return A `[draws x parameters]` matrix.
#' @noRd
.fit_draws <- function(fit, pars) {
  as.matrix(fit$stanfit, pars = pars)
}


#' Pull an n-length parameter vector out of a draws matrix
#'
#' @param D Draws matrix from [.fit_draws()].
#' @param stem Parameter stem, e.g. `"phi"`.
#' @param n Expected length.
#' @return A `[draws x n]` matrix, or `NULL` if the stem is absent.
#' @noRd
.draw_block <- function(D, stem, n) {
  cols <- grep(sprintf("^%s\\[[0-9]+\\]$", stem), colnames(D))
  if (!length(cols)) return(NULL)
  if (length(cols) != n) {
    stop("Parameter `", stem, "` has ", length(cols), " elements but `geo` has ",
         n, " rows. These must be the same model.", call. = FALSE)
  }
  # grep returns columns in name order, which is also index order for [1]..[n]
  # only up to 9 areas; sort on the parsed index instead.
  idx <- as.integer(sub(".*\\[([0-9]+)\\]$", "\\1", colnames(D)[cols]))
  D[, cols[order(idx)], drop = FALSE]
}


#' Split a BYM2 fit into the log contributions of its named parts
#'
#' Returns, per posterior draw and per area, the additive log-scale contribution
#' of each covariate and of each random-effect component, plus the baseline
#' count they are departures from.
#'
#' The BYM2 convolution is `phi + theta` where `phi` is the spatially structured
#' (ICAR) term already scaled by `spatial_scale * sqrt(rho) / sqrt(scale_factor)`
#' and `theta` is the unstructured term scaled by `spatial_scale * sqrt(1 - rho)`
#' (geostan's `convolve_bym2()`). Both are monitored, so the two can be reported
#' separately rather than as one lump. They are not the same object for policy:
#' a spatially structured excess is shared with the neighbours and points at
#' something territorial, an unstructured one is that area alone.
#'
#' The identity `mu = B * exp(sum of contributions)` is checked rather than
#' assumed. Where it does not hold - a disconnected graph puts singleton areas
#' on a different convolution branch, and a measurement-error model estimates
#' beta against a latent covariate rather than the observed column - the gap is
#' returned as an explicit `unexplained` component instead of being silently
#' distributed over the others.
#'
#' @param fit A BYM2 fit from [fit_model()] or [fit_bym2()].
#' @param geo The `sf` it was fitted to, same row order.
#' @param covariates Covariate columns. Taken from the fit when `NULL`.
#' @param exp_col Expected-count column. Taken from the fit when `NULL`.
#' @param tol Tolerance, on the log scale, before the reconstruction gap is
#'   promoted to its own component. Default `1e-6`.
#'
#' @return A list: `d` (named list of `[draws x areas]` log contributions),
#'   `base` (`[draws x areas]` baseline counts `E * exp(alpha)`), `mu`
#'   (`[draws x areas]` fitted counts), `E`, `y`, `covariates`, and
#'   `max_gap` (the largest absolute reconstruction error on the log scale).
#' @noRd
#' @param spatial `"split"` keeps the ICAR and IID parts of a BYM2 convolution
#'   apart; `"combined"` reports one `spatial` component. `"auto"` (default)
#'   splits where `phi` and `theta` exist and combines otherwise, which is what
#'   an ESF or plain-GLM fit needs.
#' @noRd
.bym2_components <- function(fit, geo, covariates = NULL, exp_col = NULL,
                             spatial = c("auto", "split", "combined"),
                             tol = 1e-6) {

  spatial <- match.arg(spatial)
  sp <- tryCatch(model_spec(fit), error = function(e) NULL)
  if (is.null(covariates)) {
    covariates <- if (is.null(sp)) character(0) else sp$covariates
  }
  if (is.null(exp_col)) exp_col <- if (is.null(sp)) "total_exp" else sp$exp_col
  obs_col <- if (is.null(sp)) "total_obs" else sp$obs_col

  require_cols(geo, c(exp_col, covariates), "geo")
  n <- nrow(geo)

  avail     <- .stan_par_names(fit)
  has_phi   <- any(grepl("^phi\\[", avail))
  has_theta <- any(grepl("^theta\\[", avail))

  if (spatial == "split" && !has_phi) {
    stop("This fit has no `phi` parameter, so its spatial term cannot be ",
         "split into structured and unstructured parts. An ESF fit spreads ",
         "its spatial variation over eigenvectors and a GLM has none at all; ",
         "use spatial = \"combined\" to report one lumped term instead.",
         call. = FALSE)
  }
  split <- spatial != "combined" && has_phi

  beta_pars <- if (length(covariates)) .beta_pars(fit, covariates) else character(0)
  D <- .fit_draws(fit, c("intercept", beta_pars,
                         if (split) "phi", if (split && has_theta) "theta",
                         "fitted"))

  alpha <- D[, "intercept"]
  mu    <- .draw_block(D, "fitted", n)

  E    <- as.numeric(geo[[exp_col]])
  base <- outer(exp(alpha), E)                    # [draws x areas]

  d <- list()
  if (length(covariates)) {
    X <- as.matrix(.plain(geo)[, covariates, drop = FALSE])
    if (anyNA(X)) {
      stop("Covariate matrix contains NA for ",
           sum(!stats::complete.cases(X)), " area(s). The model cannot have ",
           "been fitted on these rows.", call. = FALSE)
    }
    for (k in seq_along(covariates)) {
      d[[covariates[k]]] <- outer(D[, beta_pars[k]], X[, k])
    }
  }

  if (split) {
    d[["icar"]] <- .draw_block(D, "phi", n)
    if (has_theta) d[["iid"]] <- .draw_block(D, "theta", n)
  }

  # The identity, checked. Do not trust it because the Stan code says so; the
  # Stan code has branches (islands, measurement error, centred covariates)
  # that this reconstruction does not.
  gap     <- log(mu) - log(base) - Reduce(`+`, d)
  max_gap <- max(abs(gap))

  if (!split) {
    # Everything the intercept and the covariates do not account for is the
    # spatial term, whatever its parameterisation. Exact by construction, which
    # is why no tolerance check applies on this branch.
    d[["spatial"]] <- gap
    max_gap <- 0
  } else if (max_gap > tol) {
    d[["unexplained"]] <- gap
    warning("The fitted counts are not exactly reproduced by intercept + ",
            "covariates + phi + theta (largest log-scale gap ",
            format(max_gap, digits = 3), "). The gap is reported as an ",
            "`unexplained` component. Usual causes: a disconnected adjacency ",
            "graph (singleton areas take a different convolution branch), or ",
            "a measurement-error model, whose beta multiplies a latent ",
            "covariate rather than the observed column.", call. = FALSE)
  }

  list(d = d, base = base, mu = mu, E = E,
       y = if (obs_col %in% names(geo)) round_half_up(as.numeric(geo[[obs_col]])) else NULL,
       covariates = covariates, max_gap = max_gap)
}


#' Shapley decomposition of a multiplicative excess
#'
#' Attributes `base * (exp(sum(d)) - 1)` - the excess over `base` produced by a
#' set of additive log-scale contributions - to the individual contributions.
#'
#' The value of a coalition `S` is `v(S) = base * (exp(sum_{k in S} d_k) - 1)`,
#' the excess that would exist if only those components were switched on. The
#' Shapley value of component `k` averages its marginal contribution
#' `v(S + k) - v(S)` over every coalition it could join, weighted so that every
#' ordering of the components counts equally. The results are exact - they sum
#' to `v(full set)` - and do not depend on the order the components are listed.
#'
#' The marginal contribution simplifies: `v(S + k) - v(S) = base * e_S * (e_k - 1)`
#' where `e_S` is the product of `exp(d_j)` over `S`. So the whole thing costs
#' `K` exponentials and `2^K` products, not `2^K` exponentials.
#'
#' @param d Named list of numeric matrices (or vectors), all the same shape:
#'   the additive log-scale contributions.
#' @param base Numeric of the same shape as each element of `d`: the counts the
#'   excess is measured against.
#' @param max_components Refuse to run beyond this many components, since the
#'   cost is `2^K`. Default `12`.
#'
#' @return A named list of the same shape as `d`, holding attributed counts.
#'   `Reduce("+", result)` equals `base * (exp(Reduce("+", d)) - 1)` to
#'   floating-point precision.
#' @examples
#' d <- list(a = log(1.2), b = log(1.5))
#' s <- shapley_excess(d, base = 100)
#' # exact: the parts sum to the whole
#' all.equal(s$a + s$b, 100 * (1.2 * 1.5 - 1))
#' @export
shapley_excess <- function(d, base, max_components = 12) {

  K <- length(d)
  if (!K) return(list())
  if (K > max_components) {
    stop("Shapley attribution over ", K, " components needs 2^", K,
         " coalitions. Group them first, or raise `max_components` ",
         "deliberately.", call. = FALSE)
  }
  if (is.null(names(d)) || anyDuplicated(names(d))) {
    stop("`d` must be a list with unique names; the names become the ",
         "component labels.", call. = FALSE)
  }
  .shape <- function(x) if (is.null(dim(x))) length(x) else dim(x)
  shapes <- lapply(d, .shape)
  if (!all(vapply(shapes, identical, logical(1), shapes[[1]]))) {
    stop("All elements of `d` must have the same shape.", call. = FALSE)
  }
  if (!identical(.shape(base), shapes[[1]])) {
    stop("`base` must have the same shape as the elements of `d` (it is a ",
         "draw-by-area matrix too, because the intercept varies by draw).",
         call. = FALSE)
  }

  # Shapley weights by coalition size: s! (K - s - 1)! / K!
  w  <- vapply(0:(K - 1),
               function(s) exp(lfactorial(s) + lfactorial(K - s - 1) -
                                 lfactorial(K)),
               numeric(1))
  ed  <- lapply(d, exp)
  out <- lapply(d, function(x) x * 0)
  bit <- 2^(seq_len(K) - 1)

  for (mask in seq_len(2^K) - 1L) {
    members <- which(bitwAnd(mask, bit) > 0)
    s       <- length(members)
    if (s == K) next                       # the full coalition joins nothing
    eS      <- if (!s) 1 else Reduce(`*`, ed[members])
    common  <- w[s + 1L] * base * eS
    for (k in setdiff(seq_len(K), members)) {
      out[[k]] <- out[[k]] + common * (ed[[k]] - 1)
    }
  }
  out
}


#' Sequential decomposition, for reconciliation
#'
#' Adds the components one at a time in the given order and records what each
#' one moved. Unlike the Shapley values this depends on the order, which is
#' exactly why it is here: with the covariates added first and the random
#' effects last, the random-effect entries sum to `mu - pred`, which is the
#' estimand [residual_excess()] reports. Two different decompositions of the
#' same fit should be seen to agree where they are supposed to.
#'
#' @param d Named list of log contributions, in the order to add them.
#' @param base Baseline counts, same shape.
#' @return Named list of attributed counts, same shape.
#' @noRd
.sequential_excess <- function(d, base) {
  out  <- vector("list", length(d))
  names(out) <- names(d)
  cum  <- base
  csum <- 0
  for (k in names(d)) {
    csum     <- csum + d[[k]]
    nxt      <- base * exp(csum)
    out[[k]] <- nxt - cum
    cum      <- nxt
  }
  out
}


#' Readable labels for the decomposition components
#'
#' @param x Component names as they come out of [.bym2_components()].
#' @param extra Optional named character vector of overrides.
#' @return A character vector of labels.
#' @examples
#' component_label(c("di_score_z", "icar", "gp_density_z"))
#' @export
component_label <- function(x, extra = NULL) {
  known <- c(
    di_score_z      = "Deprivation",
    ivsm_z          = "Deprivation (IVSM)",
    gp_density_z    = "Primary care density",
    t_hub_mean_z    = "Travel time to hub",
    t_centre_mean_z = "Travel time to stroke centre",
    icar            = "Spatial random effect (ICAR)",
    iid             = "Unstructured random effect (IID)",
    spatial         = "Spatial term (combined)",
    unexplained     = "Unexplained (reconstruction gap)"
  )
  # Pollutant columns carry their year, so match them by pattern rather than
  # enumerating every vintage.
  out <- unname(known[x])
  out[is.na(out) & grepl("^no2_",  x)] <- "Air pollution (NO2)"
  out[is.na(out) & grepl("^pm25_", x)] <- "Air pollution (PM2.5)"
  if (!is.null(extra)) {
    hit <- match(x, names(extra))
    out[!is.na(hit)] <- unname(extra[hit[!is.na(hit)]])
  }
  ifelse(is.na(out), gsub("_", " ", sub("_z$", "", x)), out)
}


#' Decompose one area's excess into its determinants
#'
#' The local diagnostic grid, for the areas that matter. For each named area,
#' how many of its excess avoidable deaths per year are attributable to
#' deprivation, to air pollution, to primary care capacity, and to the spatially
#' structured and unstructured random effects - with credible intervals, and as
#' a percentage of that area's excess.
#'
#' Excess here is measured against a reference area of the same age and sex
#' structure sitting at the study-area average on every covariate with no
#' residual spatial effect, i.e. `E * exp(alpha)`. See the note at the top of
#' this file for why that reference and not the one in [residual_excess()], and
#' for how the two reconcile.
#'
#' @section Reading the percentages:
#' They are signed and they sum to 100. A determinant on which the area does
#' better than the study average contributes a negative number of deaths and a
#' negative share: it is holding the excess down, and the other shares exceed
#' 100 to compensate. Report it that way; do not clamp it.
#'
#' The interval on a percentage is a ratio of two posterior quantities and is
#' only meaningful while the denominator stays away from zero. Where the area's
#' total excess has a credible interval spanning zero the percentage intervals
#' are returned as `NA`, because a share of approximately nothing is not a
#' number worth printing.
#'
#' @param fit A BYM2 fit, normally M5.
#' @param geo The `sf` it was fitted to.
#' @param area Area key(s) to decompose. When `NULL`, taken from `excess`.
#' @param excess Output of [residual_excess()]; used to pick the flagged areas
#'   when `area` is `NULL`, and to check the reconciliation.
#' @param covariates,exp_col Taken from the fit when `NULL`.
#' @param n_years Years the death counts span, used to annualise. Default `3`.
#' @param probs Credible-interval bounds.
#' @param spatial `"split"` keeps the ICAR and IID components apart, `"combined"`
#'   reports one lumped spatial term, and `"auto"` (the default) splits a BYM2
#'   and lumps anything else. An ESF fit has no ICAR/IID split to make.
#' @param labels Optional named character vector overriding the component
#'   labels, passed to [component_label()].
#'
#' @return A list:
#'   \describe{
#'     \item{`by_component`}{One row per area x component: log contribution,
#'       relative-risk multiplier, attributable deaths per year, and share of
#'       the area's excess, each with an interval; plus `pct_log` (the
#'       log-scale share) and `deaths_seq` (the sequential attribution).}
#'     \item{`totals`}{One row per area: observed, expected, baseline, fitted
#'       and total excess per year, and the total relative risk.}
#'     \item{`reconciliation`}{One row per area: the random-effect part of the
#'       sequential decomposition against [residual_excess()]'s excess, which
#'       are the same quantity and should match.}
#'   }
#'
#' @examples
#' \dontrun{
#' # the flagged area, picked up automatically
#' att <- attribution_local(fit_M5, smr_geo_full, excess = excess, n_years = 3)
#' att$by_component
#' att$reconciliation
#' }
#' @seealso [attribution_grid()], [plot_attribution()], [residual_excess()]
#' @importFrom rlang .data
#' @export
attribution_local <- function(fit, geo,
                              area       = NULL,
                              excess     = NULL,
                              covariates = NULL,
                              exp_col    = NULL,
                              n_years    = 3,
                              probs      = c(0.025, 0.975),
                              spatial    = c("auto", "split", "combined"),
                              labels     = NULL) {

  spatial <- match.arg(spatial)

  keys <- if (!is.null(geo[["area"]])) as.character(geo[["area"]]) else
    as.character(seq_len(nrow(geo)))

  if (is.null(area)) {
    if (is.null(excess)) {
      stop("Name the area(s) to decompose, or pass `excess` so the flagged ",
           "ones can be taken from it.", call. = FALSE)
    }
    area <- as.character(excess$per_area$area[excess$per_area$flagged])
    if (!length(area)) {
      stop("`excess` flags no areas at its cutoff, so there is nothing to ",
           "decompose. Name an area explicitly to look at it anyway.",
           call. = FALSE)
    }
  }
  area <- as.character(area)
  idx  <- match(area, keys)
  if (anyNA(idx)) {
    stop("Area(s) not found in `geo`: ", paste(area[is.na(idx)], collapse = ", "),
         ". Available keys look like: ",
         paste(utils::head(keys, 3), collapse = ", "), ".", call. = FALSE)
  }

  p <- .bym2_components(fit, geo, covariates, exp_col, spatial = spatial)

  # Subset to the areas of interest BEFORE the 2^K loop.
  sub  <- function(m) m[, idx, drop = FALSE]
  d    <- lapply(p$d, sub)
  base <- sub(p$base)
  mu   <- sub(p$mu)

  shap <- shapley_excess(d, base)
  seqn <- .sequential_excess(d, base)

  total <- (mu - base) / n_years                       # [draws x areas]
  q     <- function(x) stats::quantile(x, probs = probs, names = FALSE)

  # Share intervals are only defined while the denominator stays off zero.
  tq        <- apply(total, 2, q)
  denom_ok  <- tq[1, ] * tq[2, ] > 0

  rows <- lapply(names(d), function(k) {
    a_yr  <- shap[[k]] / n_years
    s_yr  <- seqn[[k]] / n_years
    pct_d <- 100 * a_yr / total
    tibble::tibble(
      area        = area,
      component   = k,
      label       = component_label(k, extra = labels),
      log_mean    = colMeans(d[[k]]),
      log_low     = apply(d[[k]], 2, function(x) q(x)[1]),
      log_high    = apply(d[[k]], 2, function(x) q(x)[2]),
      rr          = colMeans(exp(d[[k]])),
      rr_low      = apply(exp(d[[k]]), 2, function(x) q(x)[1]),
      rr_high     = apply(exp(d[[k]]), 2, function(x) q(x)[2]),
      deaths_yr   = colMeans(a_yr),
      deaths_low  = apply(a_yr, 2, function(x) q(x)[1]),
      deaths_high = apply(a_yr, 2, function(x) q(x)[2]),
      # Point estimate as a ratio of posterior means, not the mean of the
      # ratios: the latter is dominated by draws where the denominator is small.
      pct         = 100 * colMeans(a_yr) / colMeans(total),
      pct_low     = ifelse(denom_ok, apply(pct_d, 2, function(x) q(x)[1]), NA_real_),
      pct_high    = ifelse(denom_ok, apply(pct_d, 2, function(x) q(x)[2]), NA_real_),
      deaths_seq  = colMeans(s_yr)
    )
  })
  by_component <- do.call(rbind, rows)

  # Log-scale share, for the reader who wants to check beta * x by hand.
  log_tot <- colMeans(Reduce(`+`, d))
  by_component[["pct_log"]] <-
    100 * by_component[["log_mean"]] / log_tot[match(by_component[["area"]], area)]

  by_component <- by_component[order(match(by_component[["area"]], area),
                                     -abs(by_component[["deaths_yr"]])), ,
                               drop = FALSE]

  totals <- tibble::tibble(
    area           = area,
    observed       = if (is.null(p$y)) NA_real_ else p$y[idx],
    expected       = p$E[idx],
    baseline_yr    = colMeans(base) / n_years,
    fitted_yr      = colMeans(mu) / n_years,
    excess_yr      = colMeans(total),
    excess_low     = tq[1, ],
    excess_high    = tq[2, ],
    rr_total       = colMeans(mu) / p$E[idx],
    pct_defined    = denom_ok,
    n_years        = n_years
  )

  # The sequential random-effect part is residual_excess()'s estimand.
  re_names <- intersect(c("icar", "iid", "spatial"), names(d))
  re_seq   <- if (!length(re_names)) rep(NA_real_, length(area)) else
    colSums(do.call(rbind, lapply(re_names, function(k)
      colMeans(seqn[[k]]) / n_years)), na.rm = FALSE)
  recon <- tibble::tibble(
    area              = area,
    excess_sequential = as.numeric(re_seq),
    excess_residual   = if (is.null(excess)) NA_real_ else
      excess$per_area$excess[match(area, as.character(excess$per_area$area))]
  )
  recon[["difference"]] <- recon[["excess_sequential"]] - recon[["excess_residual"]]

  list(by_component = by_component, totals = totals, reconciliation = recon,
       max_gap = p$max_gap)
}


#' The same decomposition, for every area
#'
#' Posterior means only, and no intervals. The Shapley loop over 2^K coalitions
#' is run on the full posterior - the means are means of the attributions, not
#' attributions of the means, which are not the same number - but keeping draws
#' for 279 areas x 5 components would cost several hundred megabytes to produce
#' a table nobody reads a credible interval out of. Use [attribution_local()]
#' for the areas that get named in the text.
#'
#' @inheritParams attribution_local
#' @return A tibble, one row per area x component: `area`, `component`,
#'   `label`, `deaths_yr`, `log_mean`, `pct`, and `dominant` (TRUE for the
#'   component with the largest positive attribution in that area).
#' @examples
#' \dontrun{
#' grid <- attribution_grid(fit_M5, smr_geo_full, n_years = 3)
#' plot_dominant_driver(grid, smr_geo_full)
#' }
#' @seealso [attribution_local()], [plot_dominant_driver()]
#' @importFrom rlang .data
#' @export
attribution_grid <- function(fit, geo,
                             covariates = NULL,
                             exp_col    = NULL,
                             n_years    = 3,
                             spatial    = c("auto", "split", "combined"),
                             labels     = NULL) {

  spatial <- match.arg(spatial)

  keys <- if (!is.null(geo[["area"]])) as.character(geo[["area"]]) else
    as.character(seq_len(nrow(geo)))

  p     <- .bym2_components(fit, geo, covariates, exp_col, spatial = spatial)
  shap  <- shapley_excess(p$d, p$base)
  total <- colMeans(p$mu - p$base) / n_years

  out <- do.call(rbind, lapply(names(p$d), function(k) {
    tibble::tibble(
      area      = keys,
      component = k,
      label     = component_label(k, extra = labels),
      deaths_yr = colMeans(shap[[k]]) / n_years,
      log_mean  = colMeans(p$d[[k]]),
      # A share of an excess that rounds to nothing is not a number. Areas
      # sitting on the reference get NA rather than a ratio of two roundings.
      pct       = ifelse(abs(total) < 1e-9, NA_real_,
                         100 * (colMeans(shap[[k]]) / n_years) / total)
    )
  }))

  # "Dominant" means largest POSITIVE contribution. An area whose every
  # component is protective has no dominant driver and gets none, rather than
  # being labelled with whichever factor is least protective.
  out <- out |>
    dplyr::group_by(.data[["area"]]) |>
    dplyr::mutate(
      dominant = .data[["deaths_yr"]] > 0 &
        .data[["deaths_yr"]] == max(.data[["deaths_yr"]])
    ) |>
    dplyr::ungroup()

  out
}


#' Bar chart of one area's attributed excess
#'
#' Deaths per year, signed, with credible intervals - not a stacked proportion
#' bar. Proportions cannot show a protective determinant, and the interval is
#' most of the information: in a study this size several components will have
#' attributions whose sign is not established, and a bar without a whisker
#' invites the reader to believe otherwise.
#'
#' @param att Output of [attribution_local()].
#' @param subtitle,caption Passed to `labs()`.
#' @return A ggplot, faceted by area when more than one was decomposed.
#' @examples
#' \dontrun{
#' plot_attribution(attribution_local(fit_M5, smr_geo_full, excess = excess))
#' }
#' @export
plot_attribution <- function(att, subtitle = NULL, caption = NULL) {

  d <- att$by_component
  if (!nrow(d)) {
    stop("`att$by_component` is empty; there is nothing to plot.", call. = FALSE)
  }
  d[[".row"]] <- factor(d[["label"]],
                        levels = rev(unique(d[["label"]][order(d[["deaths_yr"]])])))

  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = .data[["deaths_yr"]], y = .data[[".row"]])
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "22", colour = "grey45") +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data[["deaths_low"]], xmax = .data[["deaths_high"]]),
      linewidth = 0.6, colour = "grey25"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data[["deaths_low"]] * .data[["deaths_high"]] > 0),
      shape = 21, size = 2.8, colour = "grey15", stroke = 0.5
    ) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#8c2d04", `FALSE` = "white"),
                               guide = "none") +
    ggplot2::labs(
      x = "Attributable avoidable deaths per year (95% credible interval)",
      y = NULL, subtitle = subtitle, caption = caption
    ) +
    theme_atlas()

  if (dplyr::n_distinct(d[["area"]]) > 1) {
    p <- p + ggplot2::facet_wrap(~ .data[["area"]], ncol = 1, scales = "free_y")
  }
  p
}


#' Map of which determinant dominates in each area
#'
#' The decision-support artefact the proposal's scenario matrix asks for: not
#' how much excess an area has, but which lever is the one to pull there.
#'
#' Read it with the exceedance map beside it. This map answers "which driver"
#' for every area including the ones whose excess is not distinguishable from
#' zero, where the answer is arbitrary. Colour is a categorical claim and it
#' will look equally confident everywhere; the confidence lives in the other
#' map.
#'
#' @param grid Output of [attribution_grid()].
#' @param geo The `sf` carrying the geometries and an `area` key.
#' @param min_excess Areas whose total attributed excess per year is below this
#'   are drawn as "no material excess" rather than given a driver colour.
#'   Default `0`.
#' @param title,subtitle,caption Passed to `labs()`.
#' @return A ggplot.
#' @examples
#' \dontrun{
#' plot_dominant_driver(attribution_grid(fit_M5, smr_geo_full), smr_geo_full)
#' }
#' @export
plot_dominant_driver <- function(grid, geo, min_excess = 0,
                                 title    = NULL,
                                 subtitle = "Largest positive contribution to the area's excess",
                                 caption  = NULL) {

  require_cols(grid, c("area", "label", "deaths_yr", "dominant"), "grid")
  if (is.null(geo[["area"]])) {
    stop("`geo` needs an `area` key to join the grid onto.", call. = FALSE)
  }

  tot <- grid |>
    dplyr::group_by(.data[["area"]]) |>
    dplyr::summarise(excess_yr = sum(.data[["deaths_yr"]]), .groups = "drop")

  # An exact tie between two components leaves both flagged; take the first
  # rather than dropping the area, and accept that a tie is arbitrary.
  dom <- grid[grid[["dominant"]] %in% TRUE,
              c("area", "label"), drop = FALSE]
  dom <- dom[!duplicated(dom[["area"]]), , drop = FALSE]
  dom <- dom |>
    dplyr::left_join(tot, by = "area") |>
    dplyr::mutate(driver = dplyr::if_else(
      .data[["excess_yr"]] > min_excess, .data[["label"]],
      "No material excess"
    ))

  g <- geo |>
    dplyr::mutate(area = as.character(.data[["area"]])) |>
    dplyr::left_join(dom, by = "area") |>
    dplyr::mutate(driver = factor(dplyr::coalesce(.data[["driver"]],
                                                  "No material excess")))

  ggplot2::ggplot(g) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[["driver"]]),
                     colour = "white", linewidth = 0.1) +
    ggplot2::scale_fill_brewer(palette = "Set2", na.value = "grey90",
                               name = "Dominant driver") +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption) +
    theme_atlas(map = TRUE)
}


# Cause groups and service functions -------------------------------------------


#' Readable label for a wide-table column stem
#'
#' [mechanism_label()] knows the seven service functions. This adds the OECD
#' cause groups (the `G_` columns) and falls through to a de-snaked version for
#' anything else, so an unexpected stem is legible rather than absent.
#'
#' @param stem Column stem, e.g. `"G_diseases_of_the_circulatory_system"`.
#' @return A display label.
#' @examples
#' stratum_label(c("M_screening", "G_cancer"))
#' @export
stratum_label <- function(stem) {
  groups <- c(
    "Infectious diseases", "Cancer", "Endocrine and metabolic diseases",
    "Diseases of the nervous system", "Diseases of the circulatory system",
    "Diseases of the respiratory system", "Diseases of the digestive system",
    "Diseases of the genitourinary system",
    "Pregnancy, childbirth and perinatal period", "Congenital malformations",
    "Adverse effects of medical and surgical care", "Injuries",
    "Alcohol-related and drug-related deaths",
    "Provisional assignment of new diseases"
  )
  key <- stats::setNames(groups, smr_col("G", groups))

  out <- unname(key[stem])
  is_mech <- is.na(out) & grepl("^M_", stem)
  out[is_mech] <- mechanism_label(stem[is_mech])
  ifelse(is.na(out),
         sub("^(.)", "\\U\\1", gsub("_", " ", sub("^[CGM]_", "", stem)),
             perl = TRUE),
         out)
}


#' Share of avoidable deaths, and of the excess, in each stratum
#'
#' Answers both halves of "what proportion of avoidable deaths is attributable
#' to each cause group / each service function": the composition of the deaths,
#' and the composition of the departure from the age-sex expectation. The two
#' are different questions and they usually give different answers - a stratum
#' can be a large share of deaths and none of the excess, which is the whole
#' point of standardising.
#'
#' No model is fitted. The `<stem>_obs` and `<stem>_exp` columns that
#' [preprocess_smr()] emits are additive across strata by construction - the
#' per-stratum expected counts are built from stratum-specific standard rates
#' applied to the same age-sex denominators, so they sum to the all-cause
#' expected - which means the shares are exact and are available for every
#' stratum, including the ones far too sparse to smooth. Fitting a BYM2 per
#' cause group would give intervals, and would give them for perhaps half the
#' groups: a stratum carrying forty deaths across 279 areas estimates its prior.
#' A share is a descriptive quantity and is better computed descriptively.
#'
#' @section Weights:
#' A death from one of the seven causes the OECD list splits 50/50 contributes
#' half a death to each of two service functions, and a whole death to its one
#' cause group. Both sets of columns therefore sum to the total, and the two
#' `pct_deaths` columns are directly comparable.
#'
#' @param geo An `sf` (or data frame) carrying the `<prefix>_*_obs` and
#'   `<prefix>_*_exp` columns, normally `smr_geo_full`.
#' @param prefix Column prefix: `"M"` for service functions, `"G"` for OECD
#'   cause groups, `"C"` for individual causes.
#' @param area Optional area key(s) to restrict to. `NULL` (default) uses the
#'   whole study area.
#' @param n_years Years the counts span. Default `3`.
#' @param min_deaths Strata with fewer than this many deaths are still returned
#'   but flagged `sparse`, since a share computed on four deaths is not a
#'   finding. Default `20`.
#'
#' @return A tibble, one row per stratum, ordered by deaths: `stratum`, `stem`,
#'   `deaths`, `deaths_yr`, `pct_deaths`, `expected`, `excess`, `excess_yr`,
#'   `pct_excess`, `smr`, `sparse`. `pct_excess` is the share of the total
#'   positive-and-negative excess and is signed.
#'
#' @examples
#' \dontrun{
#' stratum_shares(smr_geo_full, "G")                    # cause groups, study-wide
#' stratum_shares(smr_geo_full, "M", area = "015146_55")  # one NIL
#' }
#' @seealso [mechanism_label()], [stratum_label()]
#' @importFrom rlang .data
#' @export
stratum_shares <- function(geo, prefix = c("M", "G", "C"),
                           area       = NULL,
                           n_years    = 3,
                           min_deaths = 20) {

  prefix <- match.arg(prefix)
  tab    <- .plain(geo)

  if (!is.null(area)) {
    if (is.null(tab[["area"]])) {
      stop("`geo` has no `area` column to subset on.", call. = FALSE)
    }
    keep <- as.character(tab[["area"]]) %in% as.character(area)
    if (!any(keep)) {
      stop("None of the requested areas are in `geo`: ",
           paste(area, collapse = ", "), call. = FALSE)
    }
    tab <- tab[keep, , drop = FALSE]
  }

  obs_cols <- grep(sprintf("^%s_.*_obs$", prefix), names(tab), value = TRUE)
  if (!length(obs_cols)) {
    stop("No columns match '^", prefix, "_.*_obs$'. preprocess_smr() emits ",
         "these for each level of its `class_vars`; check the prefix.",
         call. = FALSE)
  }
  stems    <- sub("_obs$", "", obs_cols)
  exp_cols <- paste0(stems, "_exp")
  missing  <- setdiff(exp_cols, names(tab))
  if (length(missing)) {
    stop("Observed columns without a matching expected column: ",
         paste(sub("_exp$", "", missing), collapse = ", "), ".", call. = FALSE)
  }

  obs <- unname(vapply(obs_cols, function(v) sum(tab[[v]], na.rm = TRUE),
                       numeric(1)))
  exd <- unname(vapply(exp_cols, function(v) sum(tab[[v]], na.rm = TRUE),
                       numeric(1)))

  # The additivity the shares depend on, checked rather than assumed. A stratum
  # set that does not sum to the all-cause total means the classification is
  # not a partition of the deaths, and every percentage below is then a share
  # of the wrong denominator.
  if ("total_obs" %in% names(tab)) {
    tot <- sum(tab[["total_obs"]], na.rm = TRUE)
    if (abs(sum(obs) - tot) > 1e-6 * max(1, tot)) {
      warning("The '", prefix, "_' strata sum to ", format(sum(obs)),
              " deaths but total_obs is ", format(tot),
              ". The shares below are shares of the stratum sum, not of the ",
              "study total; the classification is not partitioning the deaths.",
              call. = FALSE)
    }
  }

  excess     <- obs - exd
  tot_obs    <- sum(obs)
  tot_excess <- sum(excess)

  out <- tibble::tibble(
    stratum    = stratum_label(stems),
    stem       = stems,
    deaths     = as.numeric(obs),
    deaths_yr  = as.numeric(obs) / n_years,
    pct_deaths = 100 * as.numeric(obs) / tot_obs,
    expected   = as.numeric(exd),
    excess     = as.numeric(excess),
    excess_yr  = as.numeric(excess) / n_years,
    pct_excess = if (abs(tot_excess) < 1e-9) NA_real_ else
      100 * as.numeric(excess) / tot_excess,
    smr        = ifelse(exd > 0, as.numeric(obs) / as.numeric(exd), NA_real_),
    sparse     = as.numeric(obs) < min_deaths
  )
  out[order(-out[["deaths"]]), , drop = FALSE]
}


#' Deaths and excess by stratum, as a chart
#'
#' Two panels side by side because the comparison is the finding: the share of
#' deaths a stratum carries, against the share of the excess it carries.
#'
#' @param shares Output of [stratum_shares()].
#' @param subtitle,caption Passed to `labs()`.
#' @return A ggplot.
#' @examples
#' \dontrun{
#' plot_stratum_shares(stratum_shares(smr_geo_full, "G"))
#' }
#' @export
plot_stratum_shares <- function(shares, subtitle = NULL, caption = NULL) {

  require_cols(shares, c("stratum", "pct_deaths", "pct_excess"), "shares")

  long <- shares |>
    dplyr::select(dplyr::all_of(c("stratum", "pct_deaths", "pct_excess",
                                  "sparse"))) |>
    tidyr::pivot_longer(cols = dplyr::all_of(c("pct_deaths", "pct_excess")),
                        names_to = "quantity", values_to = "pct") |>
    dplyr::mutate(quantity = factor(
      .data[["quantity"]],
      levels = c("pct_deaths", "pct_excess"),
      labels = c("% of avoidable deaths", "% of excess over expectation")
    ))

  ord <- shares[["stratum"]][order(shares[["pct_deaths"]])]
  long[["stratum"]] <- factor(long[["stratum"]], levels = ord)

  ggplot2::ggplot(long, ggplot2::aes(x = .data[["pct"]],
                                     y = .data[["stratum"]])) +
    ggplot2::geom_vline(xintercept = 0, linetype = "22", colour = "grey45") +
    ggplot2::geom_col(ggplot2::aes(alpha = .data[["sparse"]]),
                      fill = "#4682b4", width = 0.7) +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.35),
                                labels = c(`FALSE` = "", `TRUE` = "sparse"),
                                name = NULL, guide = "none") +
    ggplot2::facet_wrap(~ .data[["quantity"]], scales = "free_x") +
    ggplot2::labs(x = "%", y = NULL, subtitle = subtitle, caption = caption) +
    theme_atlas()
}
