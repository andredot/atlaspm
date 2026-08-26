# Tracer condition and control outcomes ----------------------------------------
#
# Objective 3. Before this file, `model_stroke` regressed `total_obs` - ALL
# avoidable mortality - on travel time to a stroke hub. That is not the tracer
# design the methods describe. A tracer works because the condition has a
# time-critical intervention; pooling every avoidable death together destroys
# exactly the property being exploited, and produces the confounded comparison
# the negative controls were introduced to guard against.
#
# The four models:
#
#   T1  tracer                    cerebrovascular ~ travel time to hub
#   T2  negative control OUTCOME  all cancer      ~ travel time to hub
#   T3  negative control EXPOSURE cerebrovascular ~ <undecided>
#   T4  positive control          lifestyle/NCD   ~ deprivation
#
# T2 works because cancer mortality reflects exposures accumulated over decades
# and is insensitive to the timeliness of emergency care, while sharing with
# cerebrovascular mortality the whole confounding structure of area type. If
# travel time predicts cancer as strongly as it predicts stroke, the exposure is
# measuring rurality rather than access to time-critical care.
#
# T3 is stubbed pending a decision on which exposure to use.


#' Attach the tracer and control outcome columns
#'
#' `preprocess_smr()` already emits an observed/expected pair for every cause,
#' group and mechanism; this selects the four pairs the control analysis needs
#' and gives them short, stable names.
#'
#' Column names are \emph{derived} through [smr_col()] rather than hardcoded,
#' because `janitor::make_clean_names()` transliterates in ways that are not
#' guessable - "Lifestyle and NCDs" becomes `lifestyle_and_nc_ds`, not
#' `lifestyle_and_ncds`. Hardcoding would work until the day someone edits a
#' label in the lookup.
#'
#' A note on the counts. Cerebrovascular disease is one of the seven causes
#' split across the preventable and treatable lists, so each death contributes
#' two rows of `weight = 0.5`. Both rows share the same `cause` label, so
#' `preprocess_smr()`'s cause-level aggregation sums them back to 1.0 and
#' `cvd_obs` is a whole death count, not a half one. The same holds for
#' `cancer_obs`, where colorectal and breast cancer are split across mechanisms
#' but share a group.
#'
#' @param geo The modelling `sf` carrying the wide SMR columns.
#' @param tracer_cause Cause label for the tracer. Default
#'   `"Cerebrovascular diseases"`.
#' @param nco_group Group label for the negative control outcome. Default
#'   `"Cancer"` (all cancer, not lung cancer alone).
#' @param positive_mechanism Mechanism label for the positive control. Default
#'   `"Lifestyle and NCDs"`.
#'
#' @return `geo` with `cvd_obs`/`cvd_exp`, `cancer_obs`/`cancer_exp` and
#'   `poscontrol_obs`/`poscontrol_exp` added.
#' @seealso [fit_controls()]
#' @export
attach_tracer_outcomes <- function(geo,
                                   tracer_cause       = "Cerebrovascular diseases",
                                   nco_group          = "Cancer",
                                   positive_mechanism = "Lifestyle and NCDs") {

  pairs <- list(
    cvd        = c(smr_col("C", tracer_cause, "_obs"),
                   smr_col("C", tracer_cause, "_exp")),
    cancer     = c(smr_col("G", nco_group, "_obs"),
                   smr_col("G", nco_group, "_exp")),
    poscontrol = c(smr_col("M", positive_mechanism, "_obs"),
                   smr_col("M", positive_mechanism, "_exp"))
  )

  needed <- unlist(pairs, use.names = FALSE)
  require_cols(geo, needed, "geo")

  for (nm in names(pairs)) {
    geo[[paste0(nm, "_obs")]] <- round_half_up(geo[[pairs[[nm]][1]]])
    geo[[paste0(nm, "_exp")]] <- geo[[pairs[[nm]][2]]]

    E <- geo[[paste0(nm, "_exp")]]
    if (any(is.na(E) | E <= 0)) {
      stop("Outcome '", nm, "' has ", sum(is.na(E) | E <= 0),
           " area(s) with expected <= 0 or NA, so it cannot carry a log ",
           "offset. This stratum is too sparse to model at this resolution.",
           call. = FALSE)
    }
  }

  attr(geo, "tracer_labels") <- c(
    cvd        = tracer_cause,
    cancer     = nco_group,
    poscontrol = positive_mechanism
  )
  geo
}


#' Fit the tracer and control models
#'
#' Fits T1, T2 and T4 (and T3 when `nce_exposure` is supplied) on a shared
#' adjacency graph and scale factor, so that differences between the
#' coefficients cannot be attributed to differences in the smoothing.
#'
#' Every model uses the same BYM2 engine and the same exposure column, so the
#' comparison across them is a comparison of outcomes, which is the only thing
#' the negative-control argument licenses.
#'
#' @param geo Output of [attach_tracer_outcomes()], with the exposure columns
#'   attached by [add_stroke_access()] and [add_covariate()].
#' @param C Adjacency matrix.
#' @param scale_factor Shared BYM2 scale factor.
#' @param exposure Travel-time exposure column. Default `"t_hub_mean_z"`.
#' @param deprivation Deprivation column for the positive control. Default
#'   `"di_score_z"`.
#' @param nce_exposure Negative control exposure column. `NULL` skips T3; see
#' @param adjust Covariates to adjust for, applied identically to the tracer and
#'   the negative control outcome. The DAG for Objective 3 has one backdoor
#'   path, `A <- hub siting <- U -> Y`, where `U` is the latent urban-rural
#'   position; these are its measured children on the outcome side, so
#'   adjusting for them is proxy control for `U`.
#' @param engine Estimator, passed to [fit_model()]. Defaults to `"glm"`. A
#'   spatial random effect is a nonparametric estimate of that same `U`, and
#'   because travel time is very nearly a deterministic function of position,
#'   conditioning on a smooth spatial surface conditions on the exposure. The
#'   BYM2 coefficient is then the effect of the exposure's non-spatial residual,
#'   which is a different estimand - fit it as a sensitivity bound, not here.
#' @param tracer_obs,tracer_exp Observed and expected counts for the tracer
#'   and, when fitted, the negative control exposure. Default to the
#'   cerebrovascular 0-74 outcome; the pipeline points them at the all-age
#'   cerebral-infarction outcome (`i63_obs` / `i63_exp`), which is the
#'   reference stroke model reported in the Results.
#'   the note below.
#' @param ... Passed to [fit_model()].
#'
#' @return A named list of fits: `tracer`, `nc_outcome`, `positive`, and
#'   `nc_exposure` when fitted.
#'
#' @section The negative control exposure:
#' T3 is deliberately not given a default. A negative control exposure has to be
#' chosen on substantive grounds - it must share the confounding structure of
#' travel-time-to-hub while being incapable of affecting stroke case fatality -
#' and picking one silently would give the analysis an unearned appearance of
#' completeness. Until it is chosen, `nce_exposure` stays `NULL` and the
#' control panel reports the slot as not yet filled.
#'
#' @export
fit_controls <- function(geo, C, scale_factor,
                         exposure     = "t_hub_mean_z",
                         deprivation  = "di_score_z",
                         adjust       = character(0),
                         engine       = "glm",
                         nce_exposure = NULL,
                         tracer_obs   = "cvd_obs",
                         tracer_exp   = "cvd_exp",
                         ...) {

  require_cols(geo, c(exposure, deprivation, adjust,
                      tracer_obs, tracer_exp,
                      "cancer_obs", "cancer_exp",
                      "poscontrol_obs", "poscontrol_exp"), "geo")

  # The tracer and the negative control outcome MUST share a right-hand side.
  # The inference is the contrast between them, and a contrast between two
  # differently-adjusted models measures the difference in adjustment as much
  # as the difference in outcome.
  rhs_exposed <- paste(c(exposure, setdiff(adjust, exposure)), collapse = " + ")

  # The positive control's exposure IS deprivation, so deprivation cannot also
  # sit in its adjustment set - it would be conditioning on the exposure.
  rhs_positive <- paste(c(deprivation, setdiff(adjust, deprivation)),
                        collapse = " + ")

  fits <- list(
    tracer = fit_model(geo, rhs = rhs_exposed, engine = engine, C = C,
                       scale_factor = scale_factor,
                       obs_col = tracer_obs, exp_col = tracer_exp, ...),
    nc_outcome = fit_model(geo, rhs = rhs_exposed, engine = engine, C = C,
                           scale_factor = scale_factor,
                           obs_col = "cancer_obs", exp_col = "cancer_exp", ...),
    positive = fit_model(geo, rhs = rhs_positive, engine = engine, C = C,
                         scale_factor = scale_factor,
                         obs_col = "poscontrol_obs",
                         exp_col = "poscontrol_exp", ...)
  )

  if (!is.null(nce_exposure)) {
    require_cols(geo, nce_exposure, "geo")
    fits$nc_exposure <- fit_model(
      geo, rhs = paste(c(nce_exposure, setdiff(adjust, nce_exposure)),
                       collapse = " + "),
      engine = engine, C = C, scale_factor = scale_factor,
      obs_col = tracer_obs, exp_col = tracer_exp, ...)
  } else {
    message("No negative control exposure supplied; T3 not fitted. ",
            "Set `nce_exposure` once the variable has been chosen.")
  }

  fits
}


#' Assemble the control panel for reporting
#'
#' Puts the four coefficients side by side with an explicit statement of what
#' each one is for and what its result would mean. The interpretation column is
#' written \emph{before} the numbers are seen, which is the point: it is a
#' pre-specification, not a post-hoc reading.
#'
#' @param fits Output of [fit_controls()].
#' @param labels Optional display labels for the exposures.
#' @param probs Credible-interval bounds.
#' @param tracer_outcome Display name of the tracer outcome, which must match
#'   the `tracer_obs` / `tracer_exp` pair given to [fit_controls()]. Hardcoding
#'   it here once meant the table kept saying "cerebrovascular mortality" after
#'   the reference outcome changed.
#' @param exposure Term name of the exposure under test. This is the *only* row
#'   in the tracer and negative-control-outcome arms that the pre-specification
#'   speaks to; the remaining rows in those arms are adjustment covariates and
#'   are deliberately left unjudged.
#' @param deprivation Term name of the positive control's exposure. Same logic:
#'   the "should be present" expectation was written for deprivation against
#'   lifestyle-preventable mortality, not for every covariate in that fit.
#' @param nce_exposure Term name of the negative control exposure, or `NULL`
#'   while T3 is unfitted.
#'
#' @return A tibble: `role`, `outcome`, `term`, `label`, `estimate`, `ci_low`,
#'   `ci_high`, `crosses_null`, `scored`, `verdict`, `expectation`. `scored` is
#'   `TRUE` only for the designated exposure in each arm; `verdict` is
#'   meaningful only on those rows.
#' @export
collect_controls <- function(fits, labels = NULL, probs = c(0.025, 0.975),
                             tracer_outcome = "I63 cerebral infarction, all ages",
                             exposure     = "t_hub_mean_z",
                             deprivation  = "di_score_z",
                             nce_exposure = NULL) {

  roles <- c(
    tracer      = "Tracer",
    nc_outcome  = "Negative control outcome",
    nc_exposure = "Negative control exposure",
    positive    = "Positive control"
  )
  outcomes <- c(
    tracer      = tracer_outcome,
    nc_outcome  = "All-cancer mortality",
    nc_exposure = tracer_outcome,
    positive    = "Lifestyle/NCD-preventable mortality"
  )
  expectation <- c(
    tracer      = "An association is consistent with access to time-critical care mattering.",
    nc_outcome  = "Should be null. An association of similar size indicates the exposure is measuring area character, not access.",
    nc_exposure = "Should be null. An association indicates residual confounding by area type.",
    positive    = "Should be present. Its absence would mean the design lacks the power or the measurement to detect anything, making the null results uninterpretable."
  )

  coefs <- collect_coefficients(fits, term_labels = labels, exponentiate = TRUE,
                                probs = probs)

  # No early return on an empty table. Returning `coefs` unchanged would hand
  # back a data frame with no `role` column at all, and the failure would then
  # appear as "Layer 1 is missing `role`" from ggplot - a message about the
  # figure, not about the models that produced nothing.
  coefs[["role"]]        <- unname(roles[coefs[["model"]]])
  coefs[["outcome"]]     <- unname(outcomes[coefs[["model"]]])
  coefs[["expectation"]] <- unname(expectation[coefs[["model"]]])

  # A control is only informative if its DIRECTION is checked, not merely
  # whether its interval clears the null. The positive control is
  # deprivation -> lifestyle-preventable mortality, which the literature says
  # is positive; an interval that excludes 1 on the PROTECTIVE side is not the
  # control passing, it is the control failing in the most informative way
  # available - it says the exposure is not measuring what it is supposed to.
  expected_dir <- c(positive = 1, tracer = NA, nc_outcome = 0, nc_exposure = 0)
  dir_obs      <- sign(log(coefs[["estimate"]]))
  want         <- unname(expected_dir[coefs[["model"]]])

  # Each arm makes a claim about ONE coefficient - the exposure it was built to
  # interrogate. The others in the same fit are the adjustment set, and judging
  # them against the arm's expectation is a category error: deprivation is
  # supposed to predict all-cancer mortality, so scoring it against the
  # negative control's "should be null" produces a FAILS that says nothing
  # about the control. Rows outside the designated term are marked not
  # applicable rather than silently dropped, so the reader can see that the
  # adjustment set was reported but deliberately not scored.
  expected_term <- c(
    tracer      = exposure,
    nc_outcome  = exposure,
    nc_exposure = if (is.null(nce_exposure)) NA_character_ else nce_exposure,
    positive    = deprivation
  )
  is_target <- !is.na(coefs[["term"]]) &
    coefs[["term"]] == unname(expected_term[coefs[["model"]]])
  is_target[is.na(is_target)] <- FALSE

  coefs[["scored"]] <- is_target
  want[!is_target]  <- NA

  coefs[["verdict"]] <- dplyr::case_when(
    !is_target                                    ~ "not applicable (adjustment covariate)",
    is.na(want)                                   ~ "no directional prediction",
    want == 0 &  coefs[["crosses_null"]]          ~ "as expected (null)",
    want == 0 & !coefs[["crosses_null"]]          ~ "FAILS: association where none expected",
    want == 1 &  coefs[["crosses_null"]]          ~ "FAILS: no association where one expected",
    want == 1 & dir_obs > 0                       ~ "as expected (positive)",
    want == 1 & dir_obs < 0                       ~ "FAILS: association is in the OPPOSITE direction",
    TRUE                                          ~ NA_character_
  )

  coefs <- coefs[order(match(coefs[["model"]], names(roles))), , drop = FALSE]

  failed <- coefs[["model"]][grepl("^FAILS", coefs[["verdict"]])]
  if (length(failed)) {
    warning("Control(s) not behaving as pre-specified: ",
            paste(failed, collapse = ", "),
            ". Read the verdict column before interpreting any null result ",
            "in this study.", call. = FALSE)
  }

  coefs[, c("role", "outcome", "term", "label", "estimate", "ci_low",
            "ci_high", "crosses_null", "scored", "verdict", "expectation")]
}


#' Travel-time descriptives for the stroke sub-model
#'
#' Range, median and population-weighted mean of travel time to the nearest
#' stroke centre and to the nearest thrombectomy hub, plus the share of the
#' population beyond the regional 45-minute centralisation threshold.
#'
#' The distinction the draft blurs: a *centre* is any accredited stroke unit,
#' including spokes that give thrombolysis but not thrombectomy; a *hub* is a
#' level-II centre that performs mechanical thrombectomy. They are different
#' networks with different catchments and the results should name which is
#' which.
#'
#' @param geo The modelling `sf` with the stroke access columns.
#'
#' @return A tibble, one row per travel-time measure.
#' @export
stroke_time_summary <- function(geo) {

  tab  <- sf::st_drop_geometry(geo)
  vars <- c(
    "Nearest stroke centre (any level), minutes"        = "t_centre_mean",
    "Nearest thrombectomy hub (level II), minutes"      = "t_hub_mean",
    "Nearest thrombectomy hub, 90th percentile, minutes" = "t_hub_p90",
    "Population share beyond 45 minutes from a hub"     = "pop_share_over_45min_hub"
  )
  vars <- vars[vars %in% names(tab)]

  if (!length(vars)) {
    stop("No stroke accessibility columns found. Was add_stroke_access() run?",
         call. = FALSE)
  }
  covariate_summary(geo, vars)
}
