# Model dispatcher -------------------------------------------------------------
#
# The thesis names seven models M1-M7. Before this file, the pipeline had
# targets called model_base, model_di, model_ivsm and model_stroke, none of
# which mapped onto that numbering, and three of the seven did not exist. A
# reader of Table 3 had no way to find the code that produced a row.
#
# fit_model() is the single entry point. Every model in the study is one row of
# the `model_specs` tibble in _targets.R: an id, a label, an engine and a
# right-hand side. Changing a specification means editing one string in one
# place, and the target name, the table row and the figure panel all follow.


#' Fit one model from a specification
#'
#' Dispatches to the appropriate \pkg{geostan} sampler and returns its fit. All
#' three engines use the same Poisson likelihood with \code{log(expected)} as an
#' offset, which is what makes the leave-one-out comparison across M1-M7
#' legitimate: the pointwise log-likelihood is reconstructed the same way for
#' every model in \code{\link{compare_bym2}}.
#'
#' \strong{Engines.}
#' \describe{
#'   \item{\code{"glm"}}{\code{geostan::stan_glm()} - no area random effect.
#'     M1 (intercept only) and M2 (covariates, no random effect). These exist to
#'     show what the covariates achieve before any smoothing, and to give the
#'     LOO comparison a floor.}
#'   \item{\code{"bym2"}}{\code{geostan::stan_icar(type = "bym2")} with a
#'     properly scaled ICAR component. M3-M5 and M7.}
#'   \item{\code{"esf"}}{\code{geostan::stan_esf()} - eigenvector spatial
#'     filtering under a regularised horseshoe prior. M6. This is the
#'     alternative specification the methods section argues for: the partition
#'     between a smooth covariate and a smooth random effect is weakly
#'     identified under BYM2, and ESF selects among the map's own eigenvectors
#'     instead, shrinking unsupported patterns to zero.}
#' }
#'
#' The formula is assembled rather than passed in, so the offset can never be
#' forgotten on one model and present on another - a difference that would make
#' the elpd values incomparable while looking perfectly reasonable in the code.
#'
#' @param geo An \code{sf}, one row per area, in the same order as \code{C}.
#' @param rhs Right-hand side of the formula as a string, excluding the offset.
#'   \code{"1"} for an intercept-only model.
#' @param engine One of \code{"bym2"}, \code{"glm"}, \code{"esf"}.
#' @param C Binary adjacency matrix. Required for \code{"bym2"} and
#'   \code{"esf"}.
#' @param scale_factor BYM2 scaling vector from
#'   \code{\link{compute_scale_factor}}. Computed from \code{C} when
#'   \code{NULL}.
#' @param obs_col,exp_col Observed-count and expected-count columns.
#' @param chains,iter,control Sampler settings. Defaults match those reported in
#'   the methods.
#' @param ... Passed to the underlying \pkg{geostan} function.
#'
#' @return The fitted \pkg{geostan} object, with the specification recorded in
#'   the \code{"atlaspm_spec"} attribute so downstream summaries can label
#'   themselves without being told twice.
#'
#' @examples
#' \dontrun{
#' m5 <- fit_model(smr_geo_full, "di_score_z + no2_z + gp_density_z",
#'                 engine = "bym2", C = C_matrix, scale_factor = scale_factor)
#' }
#' @seealso [fit_bym2()], [compare_bym2()], [residual_excess()]
#' @importFrom stats as.formula poisson
#' @export
fit_model <- function(geo,
                      rhs          = "1",
                      engine       = c("bym2", "glm", "esf"),
                      C            = NULL,
                      scale_factor = NULL,
                      obs_col      = "total_obs",
                      exp_col      = "total_exp",
                      chains       = 4,
                      iter         = 4000,
                      control      = list(adapt_delta = 0.97,
                                          max_treedepth = 12),
                      ...) {

  engine <- match.arg(engine)
  require_cols(geo, c(obs_col, exp_col), "geo")

  E <- geo[[exp_col]]
  if (any(is.na(E) | E <= 0)) {
    stop(sum(is.na(E) | E <= 0), " area(s) have expected counts <= 0 or NA, so ",
         "log(", exp_col, ") is undefined. A death-free area still has a ",
         "positive expected count; a zero here means something was zero-filled ",
         "upstream.", call. = FALSE)
  }

  # Fractional observed counts arise from the split causes. Round half UP, per
  # the methods; base round() is round-half-to-even and would disagree on
  # exactly the values that occur most often here.
  y_raw <- geo[[obs_col]]
  y_int <- round_half_up(y_raw)
  if (!isTRUE(all.equal(y_raw, y_int))) {
    message(sprintf(
      "  %s: rounded to integer (%d areas changed, total absolute shift %.1f deaths).",
      obs_col, sum(y_raw != y_int), sum(abs(y_raw - y_int))
    ))
  }
  geo[[obs_col]] <- y_int

  # Assert that the covariates named in `rhs` are actually present, BEFORE
  # handing the formula to geostan. Without this, a specification that failed
  # to substitute - or a covariate column that was renamed upstream - produces
  # a model that fits perfectly well and is silently intercept-only. That
  # failure is invisible until the coefficient table comes back empty, several
  # targets later, and looks like a plotting bug.
  covs <- spec_covariates(rhs)
  if (length(covs)) {
    require_cols(geo, covs, "geo")
    bad <- covs[vapply(covs, function(v) {
      x <- sf::st_drop_geometry(geo)[[v]]
      anyNA(x) || stats::sd(x, na.rm = TRUE) == 0
    }, logical(1))]
    if (length(bad)) {
      stop("Covariate(s) ", paste(bad, collapse = ", "),
           " are constant or contain NA across the ", nrow(geo),
           " modelled areas. A constant covariate is collinear with the ",
           "intercept and the model would not be identified.", call. = FALSE)
    }
  }

  f <- stats::as.formula(
    sprintf("%s ~ %s + offset(log(%s))", obs_col, rhs, exp_col)
  )

  if (engine %in% c("bym2", "esf") && is.null(C)) {
    stop("Engine '", engine, "' needs an adjacency matrix `C`.", call. = FALSE)
  }

  fit <- switch(
    engine,
    glm = geostan::stan_glm(
      formula = f, data = geo, family = stats::poisson(),
      chains = chains, iter = iter, control = control, ...
    ),
    bym2 = {
      if (is.null(scale_factor)) scale_factor <- compute_scale_factor(C)
      geostan::stan_icar(
        formula = f, data = geo, family = stats::poisson(),
        C = C, type = "bym2", scale_factor = scale_factor,
        chains = chains, iter = iter, control = control, ...
      )
    },
    esf = geostan::stan_esf(
      formula = f, data = geo, family = stats::poisson(),
      C = C, chains = chains, iter = iter, control = control, ...
    )
  )

  # Post-condition: the fit must contain one coefficient per covariate. If it
  # does not, the covariates were dropped somewhere between here and Stan, and
  # it is far better to know now than to discover it in a figure.
  if (length(covs)) {
    avail <- rownames(rstan::summary(fit$stanfit)$summary)
    if (!any(grepl("^(beta|gamma)[\\[_]", avail)) &&
        !all(covs %in% avail)) {
      stop("Fitted '", rhs, "' but the resulting model has no coefficient ",
           "parameters, i.e. it is intercept-only. The covariates did not ",
           "reach geostan.", call. = FALSE)
    }
  }

  attr(fit, "atlaspm_spec") <- list(
    rhs        = rhs,
    engine     = engine,
    obs_col    = obs_col,
    exp_col    = exp_col,
    covariates = spec_covariates(rhs),
    formula    = deparse(f)
  )
  fit
}


#' Extract covariate names from a formula right-hand side
#'
#' @param rhs Right-hand side string, e.g. \code{"di_score_z + no2_z"}.
#' @return Character vector of variable names; \code{character(0)} for
#'   \code{"1"}.
#' @examples
#' spec_covariates("di_score_z + no2_z")
#' spec_covariates("1")
#' @export
spec_covariates <- function(rhs) {
  if (is.null(rhs) || !nzchar(trimws(rhs)) || trimws(rhs) == "1") {
    return(character(0))
  }
  all.vars(stats::as.formula(paste("~", rhs)))
}


#' Recover the specification attached to a fit
#'
#' @param fit A model from [fit_model()].
#' @return The specification list.
#' @export
model_spec <- function(fit) {
  sp <- attr(fit, "atlaspm_spec")
  if (is.null(sp)) {
    stop("This fit carries no `atlaspm_spec`. It was not produced by ",
         "fit_model(), so its covariates cannot be recovered automatically; ",
         "pass them explicitly.", call. = FALSE)
  }
  sp
}
