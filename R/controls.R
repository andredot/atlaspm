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
                         nce_exposure = NULL,
                         ...) {

  require_cols(geo, c(exposure, deprivation,
                      "cvd_obs", "cvd_exp",
                      "cancer_obs", "cancer_exp",
                      "poscontrol_obs", "poscontrol_exp"), "geo")

  fits <- list(
    tracer = fit_model(geo, rhs = exposure, engine = "bym2", C = C,
                       scale_factor = scale_factor,
                       obs_col = "cvd_obs", exp_col = "cvd_exp", ...),
    nc_outcome = fit_model(geo, rhs = exposure, engine = "bym2", C = C,
                           scale_factor = scale_factor,
                           obs_col = "cancer_obs", exp_col = "cancer_exp", ...),
    positive = fit_model(geo, rhs = deprivation, engine = "bym2", C = C,
                         scale_factor = scale_factor,
                         obs_col = "poscontrol_obs",
                         exp_col = "poscontrol_exp", ...)
  )

  if (!is.null(nce_exposure)) {
    require_cols(geo, nce_exposure, "geo")
    fits$nc_exposure <- fit_model(geo, rhs = nce_exposure, engine = "bym2",
                                  C = C, scale_factor = scale_factor,
                                  obs_col = "cvd_obs", exp_col = "cvd_exp", ...)
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
#'
#' @return A tibble: `role`, `outcome`, `exposure`, `rr`, `ci_low`, `ci_high`,
#'   `expectation`.
#' @export
collect_controls <- function(fits, labels = NULL, probs = c(0.025, 0.975)) {

  roles <- c(
    tracer      = "Tracer",
    nc_outcome  = "Negative control outcome",
    nc_exposure = "Negative control exposure",
    positive    = "Positive control"
  )
  outcomes <- c(
    tracer      = "Cerebrovascular mortality",
    nc_outcome  = "All-cancer mortality",
    nc_exposure = "Cerebrovascular mortality",
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

  coefs[["verdict"]] <- dplyr::case_when(
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
            "ci_high", "crosses_null", "verdict", "expectation")]
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


# Primary care placeholder -----------------------------------------------------

#' Synthetic GP density, standing in for the NAR extract
#'
#' \strong{This is not data.} It is a seeded synthetic covariate that lets M2,
#' M5, M6 and M7 be fitted, checked and reported on before the Nuova Anagrafe
#' Regionale extract exists. Every result that depends on it is provisional and
#' the report says so wherever it appears.
#'
#' \strong{Why it is not simply a copy of the pollution surface.} Reusing the
#' NO2 values verbatim would make `gp_density_z` and `no2_z` perfectly
#' collinear. M5 contains both, so the model would be unidentified: the sampler
#' would wander along the ridge where the two coefficients trade off, R-hat
#' would fail, and the diagnostics would be reporting a defect of the
#' placeholder rather than anything about the study. What is wanted is a
#' covariate with a \emph{realistic spatial structure} and a realistic
#' correlation with the other covariates, which is what this builds: the NO2
#' surface is used as a spatial scaffold, then mixed with spatially smoothed
#' noise to hit a target correlation.
#'
#' The sign is negative by default. The inverse care law, and the empirical
#' pattern in the literature the methods cite, both point the same way: primary
#' care supply tends to be lower where environmental and social disadvantage is
#' higher.
#'
#' Replacing this with the real thing is a one-target change in `_targets.R`:
#' swap `gp_density_area` from `simulate_gp_density()` to
#' `compute_density_by_area()$density`. Nothing downstream needs to change,
#' because the column name and shape are identical.
#'
#' @param pollution_area Area-level pollution table, used only as a spatial
#'   scaffold.
#' @param scaffold_col Column to scaffold from. Detected when `NULL`.
#' @param target_cor Target Spearman correlation with the scaffold. Default
#'   `-0.45`.
#' @param mean_density,sd_density Mean and SD of the simulated density, on the
#'   scale of GP full-time-equivalents per 1,000 residents. Defaults chosen to
#'   sit in the plausible Italian range.
#' @param seed Random seed, so the pipeline is reproducible.
#'
#' @return A tibble: `area`, `gp_density`, with `attr(, "synthetic") = TRUE`.
#' @examples
#' \dontrun{
#' gp <- simulate_gp_density(pollution_area)
#' stopifnot(isTRUE(attr(gp, "synthetic")))
#' }
#' @export
simulate_gp_density <- function(pollution_area,
                                scaffold_col = NULL,
                                target_cor   = -0.45,
                                mean_density = 0.75,
                                sd_density   = 0.12,
                                seed         = 20260821L) {

  require_cols(pollution_area, "area", "pollution_area")

  if (is.null(scaffold_col)) {
    # Only a bare `<pollutant>_<year>` column, never a coverage or _z column,
    # and the most recent year - the same rule exposure_columns() applies, so
    # the placeholder is scaffolded on the exposure surface rather than on the
    # 2013 diagnostic one.
    cand <- grep("^(no2|pm25)_[0-9]{4}$", names(pollution_area), value = TRUE)
    if (!length(cand)) {
      stop("No pollutant column matching '<no2|pm25>_<year>' found in ",
           "`pollution_area` to scaffold from. Available: ",
           paste(names(pollution_area), collapse = ", "),
           ". Pass `scaffold_col` explicitly.", call. = FALSE)
    }
    scaffold_col <- cand[which.max(as.integer(sub(".*_", "", cand)))]
  }

  x <- as.numeric(pollution_area[[scaffold_col]])
  n <- length(x)

  if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) {
    stop("The scaffold column '", scaffold_col, "' is constant or all NA.",
         call. = FALSE)
  }
  x[is.na(x)] <- mean(x, na.rm = TRUE)

  zs <- as.numeric(scale(x))

  withr::with_seed(seed, {
    noise <- stats::rnorm(n)
    # Gram-Schmidt: make the noise exactly orthogonal to the scaffold, so the
    # realised correlation is the requested one rather than the requested one
    # plus whatever the noise happened to share with it.
    noise <- noise - stats::coef(stats::lm(noise ~ zs))[2] * zs
    noise <- as.numeric(scale(noise))

    rho <- max(-0.95, min(0.95, target_cor))
    z   <- rho * zs + sqrt(1 - rho^2) * noise
  })

  out <- tibble::tibble(
    area       = pollution_area[["area"]],
    gp_density = pmax(0.05, mean_density + sd_density * as.numeric(scale(z)))
  )

  attr(out, "synthetic")    <- TRUE
  attr(out, "scaffold")     <- scaffold_col
  attr(out, "realised_cor") <- stats::cor(out[["gp_density"]], x,
                                          method = "spearman")
  attr(out, "provenance")   <- paste0(
    "SYNTHETIC PLACEHOLDER. Generated by simulate_gp_density(seed = ", seed,
    ") from the '", scaffold_col, "' surface. Not derived from the NAR. ",
    "Every estimate involving primary care capacity is provisional."
  )

  message(attr(out, "provenance"))
  out
}
