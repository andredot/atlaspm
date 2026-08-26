# Fit -> tidy table collectors -------------------------------------------------
#
# One principle: a stanfit is large and expensive, a summary table is neither.
# The report reads only the tables, so rendering it never has to hold seven
# posterior sample arrays in memory at once. Every function here takes fits and
# returns a data frame.


#' Convergence and sampler diagnostics for a set of fits
#'
#' The table behind Results "Model comparison": *"All models converged:
#' rank-normalised R-hat below 1.01, bulk and tail effective sample sizes above
#' \[N\] for every parameter, \[N\] divergent transitions."*
#'
#' Runs [check_bym2_fit()] over a named list and stacks the results, so the
#' claim can be checked against every model at once rather than eyeballed one
#' at a time.
#'
#' @param fits Named list of fitted models.
#' @param labels Optional display labels, defaulting to `names(fits)`.
#'
#' @return A tibble, one row per model: `model`, `n_draws`, `rhat_max`,
#'   `ess_min`, `n_divergent`, `pct_divergent`, `treedepth_hits`, `ebfmi_min`,
#'   `rho`, `waic`.
#' @examples
#' \dontrun{
#' collect_diagnostics(list(M3 = fit_M3, M5 = fit_M5))
#' }
#' @export
collect_diagnostics <- function(fits, labels = names(fits)) {
  if (is.null(labels)) labels <- paste0("model_", seq_along(fits))

  rows <- Map(function(fit, lab) {
    d <- check_bym2_fit(fit, print = FALSE)
    tibble::tibble(
      model          = lab,
      n_draws        = d$n_draws,
      rhat_max       = d$rhat_max,
      ess_min        = d$neff_min,
      n_rhat_above   = d$n_rhat_above,
      n_ess_below    = d$n_ess_below,
      n_divergent    = d$divergences,
      pct_divergent  = d$pct_divergent,
      treedepth_hits = d$treedepth_hits,
      ebfmi_min      = min(d$ebfmi),
      rho            = d$rho_mean,
      waic           = d$waic
    )
  }, fits, labels)

  do.call(rbind, rows)
}


#' Coefficient table across models, on the relative-risk scale
#'
#' Covariates are standardised before fitting, so each coefficient is the
#' multiplicative change in relative risk per one standard deviation of the
#' covariate. Exponentiating is therefore the only presentation that means
#' anything to a reader, and it is the default.
#'
#' Term names are recovered from each fit's specification (see [fit_model()]),
#' so `beta[1]` becomes `di_score_z` automatically and the table cannot silently
#' mislabel a coefficient when a formula changes.
#'
#' @param fits Named list of fitted models.
#' @param term_labels Optional named character vector mapping variable names to
#'   display labels, e.g. `c(di_score_z = "Deprivation (per SD)")`.
#' @param exponentiate Return relative risks rather than log-scale
#'   coefficients. Default `TRUE`.
#' @param probs Credible-interval bounds.
#'
#' @return A tibble: `model`, `term`, `label`, `estimate`, `ci_low`, `ci_high`,
#'   `rhat`, `ess`, and `crosses_null` (whether the interval contains 1, or 0 on
#'   the log scale).
#' @examples
#' \dontrun{
#' collect_coefficients(list(M4 = fit_M4, M5 = fit_M5),
#'                      term_labels = c(no2_z = "NO2 (per SD)"))
#' }
#' @export
collect_coefficients <- function(fits,
                                 term_labels  = NULL,
                                 exponentiate = TRUE,
                                 probs        = c(0.025, 0.975)) {

  one <- function(fit, lab) {
    sp   <- tryCatch(model_spec(fit), error = function(e) NULL)
    covs <- if (is.null(sp)) character(0) else sp$covariates

    # An intercept-only model (M1, M3) legitimately has nothing to report.
    if (!length(covs)) return(NULL)

    betas <- .beta_pars(fit, covs)
    smry  <- rstan::summary(fit$stanfit, pars = betas, probs = probs)$summary
    term  <- covs

    lo <- sprintf("%g%%", 100 * probs[1])
    hi <- sprintf("%g%%", 100 * probs[2])

    tf <- if (exponentiate) exp else identity

    tibble::tibble(
      model    = lab,
      term     = term,
      label    = if (is.null(term_labels)) term else
        ifelse(term %in% names(term_labels), term_labels[term], term),
      estimate = tf(smry[betas, "mean"]),
      ci_low   = tf(smry[betas, lo]),
      ci_high  = tf(smry[betas, hi]),
      rhat     = smry[betas, "Rhat"],
      ess      = smry[betas, "n_eff"]
    )
  }

  out <- do.call(rbind, Map(one, fits, names(fits)))

  if (is.null(out) || !nrow(out)) {
    stop("No coefficients were recovered from any of: ",
         paste(names(fits), collapse = ", "),
         ". Either every model is intercept-only, or the fits carry no ",
         "`atlaspm_spec` (were they built by fit_model()?). Returning an ",
         "empty table here would only move the failure downstream into the ",
         "figures.", call. = FALSE)
  }

  null_value <- if (exponentiate) 1 else 0
  out[["crosses_null"]] <- out[["ci_low"]] < null_value &
                           out[["ci_high"]] > null_value
  out
}


#' Pareto-k diagnostics from a loo object, whatever the accessor is called
#'
#' `loo::pareto_k_values()` has moved and changed behaviour across loo
#' releases. Wrapping it in a `tryCatch()` that falls back to `numeric(0)` -
#' which is what this package did - converts an accessor failure into a count
#' of ZERO problematic observations. That is not a missing diagnostic, it is a
#' diagnostic that actively says the opposite of the truth, and it contradicted
#' loo's own "Some Pareto k diagnostic values are too high" warning without
#' anything erroring.
#'
#' Try the documented accessor, then the two places the values are stored, then
#' fail.
#'
#' @param lo A `loo` object.
#' @return Numeric vector of Pareto-k values, one per observation.
#' @noRd
.pareto_k <- function(lo) {
  k <- tryCatch(loo::pareto_k_values(lo), error = function(e) NULL)

  if (is.null(k) || !length(k)) k <- lo$diagnostics$pareto_k
  if (is.null(k) || !length(k)) {
    k <- tryCatch(lo$psis_object$diagnostics$pareto_k, error = function(e) NULL)
  }
  if (is.null(k) || !length(k)) {
    stop("Cannot read Pareto-k values from this loo object. Available ",
         "elements: ", paste(names(lo), collapse = ", "),
         ". Reporting zero problematic observations would be worse than ",
         "reporting none at all.", call. = FALSE)
  }
  as.numeric(k)
}


#' Leave-one-out comparison table, with the Pareto-k audit
#'
#' Wraps a `bym2_comparison` into the table Results asks for: models ordered by
#' expected log predictive density, elpd differences from the best with their
#' standard errors, and a verdict on whether each difference is larger than its
#' own uncertainty.
#'
#' Also surfaces the Pareto-k counts, which Results quotes and which the print
#' method buries. A high k means the leave-one-out approximation is unreliable
#' for that area, and in this study the offenders are predictably the least
#' populous units - worth reporting rather than hiding.
#'
#' @param comparison A `bym2_comparison` from [compare_bym2()].
#' @param k_threshold Pareto-k value above which an observation is flagged.
#'   Default `0.7`.
#'
#' @return A tibble: `model`, `elpd_loo`, `se_elpd`, `p_loo`, `looic`,
#'   `elpd_diff`, `se_diff`, `ratio`, `verdict`, `n_pareto_high`.
#' @export
collect_loo <- function(comparison, k_threshold = 0.7) {

  cmp <- comparison$loo_compare

  # loo has changed the shape of loo_compare() between versions. Older releases
  # return a matrix whose ROWNAMES are the model labels; newer ones return a
  # data frame with an explicit `model` column, integer row names, and extra
  # diagnostic columns. Reading rownames() unconditionally silently yields
  # "1".."7", which then fails to index comparison$loo and - because
  # tibble() DROPS NULL arguments rather than erroring - produces a table that
  # is simply missing elpd_loo, se_elpd, p_loo and looic. No error anywhere.
  cmp_df <- as.data.frame(cmp, stringsAsFactors = FALSE)

  models <- if ("model" %in% names(cmp_df)) {
    as.character(cmp_df[["model"]])
  } else {
    rownames(cmp_df)
  }

  if (is.null(models) || !length(models)) {
    stop("Cannot recover model labels from `loo_compare`. It has columns: ",
         paste(names(cmp_df), collapse = ", "), ".", call. = FALSE)
  }

  known <- names(comparison$loo)
  if (!all(models %in% known)) {
    stop("loo_compare labels (", paste(models, collapse = ", "),
         ") do not match the fitted models (", paste(known, collapse = ", "),
         "). The comparison and the fits have come apart.", call. = FALSE)
  }

  # Pull a column by name for a given model, whichever shape cmp is in.
  cell <- function(m, col) {
    if (!col %in% names(cmp_df)) return(NA_real_)
    v <- cmp_df[[col]][match(m, models)]
    if (length(v) != 1L) NA_real_ else as.numeric(v)
  }

  # Never let a missing estimate vanish: tibble() drops NULL, so coerce to NA.
  est <- function(lo, row) {
    if (is.null(lo) || is.null(lo$estimates) ||
        !row %in% rownames(lo$estimates)) return(c(NA_real_, NA_real_))
    c(lo$estimates[row, "Estimate"], lo$estimates[row, "SE"])
  }

  verdict_for <- function(r) {
    if (is.na(r))   "se = 0"
    else if (r < 2) "indistinguishable"
    else if (r < 4) "weak"
    else            "clear"
  }

  rows <- lapply(seq_along(models), function(i) {
    m  <- models[i]
    lo <- comparison$loo[[m]]
    k  <- .pareto_k(lo)

    e_loo <- est(lo, "elpd_loo")
    p_loo <- est(lo, "p_loo")
    looic <- est(lo, "looic")

    ed <- cell(m, "elpd_diff")
    se <- cell(m, "se_diff")
    r  <- if (!is.na(se) && se > 0) abs(ed) / se else NA_real_

    tibble::tibble(
      model         = m,
      elpd_loo      = e_loo[1],
      se_elpd       = e_loo[2],
      p_loo         = p_loo[1],
      looic         = looic[1],
      elpd_diff     = ed,
      se_diff       = se,
      ratio         = r,
      verdict       = if (i == 1L) "reference" else verdict_for(r),
      n_pareto_high = sum(k > k_threshold)
    )
  })

  out <- dplyr::bind_rows(rows)

  if (all(is.na(out[["elpd_loo"]]))) {
    stop("Every elpd_loo is missing. comparison$loo is named (",
         paste(known, collapse = ", "),
         ") but nothing matched; check that compare_bym2() labelled its ",
         "loo list and its comparison table the same way.", call. = FALSE)
  }
  out
}


#' Which areas the leave-one-out approximation struggles on
#'
#' Results says the high-Pareto-k areas "were the \[N\] least populous units,
#' where the leave-one-out approximation is least reliable". This checks that
#' claim rather than assuming it: it returns the flagged areas with their
#' expected counts and population, so the sentence can be written from evidence.
#'
#' @param comparison A `bym2_comparison`.
#' @param geo The `sf` the models were fitted to.
#' @param model Which model in the comparison to audit. Defaults to the first.
#' @param k_threshold Flag threshold. Default `0.7`.
#'
#' @return A tibble of flagged areas, ordered by expected count.
#' @export
pareto_k_summary <- function(comparison, geo, model = NULL, k_threshold = 0.7) {

  if (is.null(model)) model <- names(comparison$loo)[1]
  k <- .pareto_k(comparison$loo[[model]])

  if (length(k) != nrow(geo)) {
    stop("The loo object has ", length(k), " observations but `geo` has ",
         nrow(geo), " rows.", call. = FALSE)
  }

  tab <- sf::st_drop_geometry(geo)
  out <- tibble::tibble(
    model    = model,
    area     = if ("area" %in% names(tab)) tab[["area"]] else seq_len(nrow(tab)),
    pareto_k = k,
    expected = if ("total_exp" %in% names(tab)) tab[["total_exp"]] else NA_real_,
    population = if ("population" %in% names(tab)) tab[["population"]] else NA_real_
  )

  # The reporting layer needs the reference point, not just the flagged rows:
  # "median expected 109" is only interpretable next to the study-wide median.
  # Carrying it as an attribute keeps the table shape unchanged.
  med_all <- stats::median(out[["expected"]], na.rm = TRUE)

  out <- out[out[["pareto_k"]] > k_threshold, , drop = FALSE]
  out <- out[order(out[["expected"]]), , drop = FALSE]
  attr(out, "median_expected_all") <- med_all
  out
}


#' Per-mechanism model summary
#'
#' Results asks for "coefficients and mixing parameters" for each mechanism
#' stratum. The per-mechanism models are intercept-only BYM2 fits, so what there
#' is to report is the mixing parameter, the residual scale, the range of
#' smoothed relative risks, and how many areas show excess.
#'
#' @param fits Named list from [fit_bym2_mechanisms()].
#' @param geo The augmented `sf` from [augment_bym2_mechanisms()].
#' @param prob_cutoff Exceedance probability cutoff. Default `0.80`.
#' @param out_suffix Suffix used by [augment_bym2_mechanisms()].
#'
#' @return A tibble, one row per mechanism: deaths, rho, spatial scale, RR range
#'   and the count of areas above the cutoff.
#' @export
collect_mechanisms <- function(fits, geo, prob_cutoff = 0.80,
                               out_suffix = "_bym2") {

  tab <- sf::st_drop_geometry(geo)

  rows <- lapply(names(fits), function(stem) {
    s   <- rstan::summary(fits[[stem]]$stanfit)$summary
    rr  <- tab[[paste0(stem, out_suffix)]]
    exc <- tab[[paste0(stem, out_suffix, "_exc")]]
    obs <- tab[[paste0(stem, "_obs")]]

    grab <- function(par, col = "mean") {
      if (par %in% rownames(s)) s[par, col] else NA_real_
    }

    tibble::tibble(
      mechanism   = mechanism_label(stem),
      stem        = stem,
      deaths      = if (is.null(obs)) NA_real_ else sum(obs),
      rho         = grab("rho"),
      rho_low     = grab("rho", "2.5%"),
      rho_high    = grab("rho", "97.5%"),
      spatial_sd  = grab("spatial_scale"),
      rr_min      = if (is.null(rr)) NA_real_ else min(rr, na.rm = TRUE),
      rr_max      = if (is.null(rr)) NA_real_ else max(rr, na.rm = TRUE),
      n_exceed    = if (is.null(exc)) NA_integer_ else sum(exc > prob_cutoff)
    )
  })

  out <- do.call(rbind, rows)
  out[order(-out[["deaths"]]), , drop = FALSE]
}


#' Turn a wide-table column stem back into a readable mechanism name
#'
#' `M_lifestyle_and_nc_ds` is what `janitor::make_clean_names()` does to
#' "Lifestyle and NCDs". Nobody wants that in a thesis table.
#'
#' @param stem Column stem, e.g. `"M_tertiary_prevention"`.
#' @return A display label.
#' @examples
#' mechanism_label("M_lifestyle_and_nc_ds")
#' @export
mechanism_label <- function(stem) {
  known <- c(
    "Immunisation and Prophylaxis", "Lifestyle and NCDs",
    "Environment and Safety", "Gender and Reproductive Health",
    "Mental Health Services", "Screening", "Tertiary prevention"
  )
  key <- stats::setNames(known, smr_col("M", known))

  out <- unname(key[stem])
  # Fall back to a de-snaked version rather than NA, so an unexpected stem is
  # legible in the table instead of vanishing.
  ifelse(is.na(out),
         gsub("_", " ", sub("^[CGM]_", "", stem)),
         out)
}


#' Concordance of the per-mechanism risk surfaces
#'
#' Results asks which clusters persist across strata and which are
#' stratum-specific. A Spearman correlation matrix of the smoothed relative-risk
#' surfaces answers the first half directly: strata that pick out the same
#' places correlate, strata that pick out different places do not.
#'
#' Spearman rather than Pearson, because what matters is whether the strata
#' agree on the *ordering* of areas, not on the magnitude of the risk - the
#' magnitudes are on different scales because the strata have very different
#' expected counts.
#'
#' @param geo Augmented `sf` from [augment_bym2_mechanisms()].
#' @param pattern Regex selecting the smoothed-RR columns.
#'
#' @return A list: `matrix` (labelled Spearman matrix) and `long` (a tidy tibble
#'   of unique pairs, ordered by correlation).
#' @export
rank_concordance <- function(geo, pattern = "^M_.*_bym2$") {

  tab  <- sf::st_drop_geometry(geo)
  cols <- grep(pattern, names(tab), value = TRUE)
  cols <- cols[!grepl("_exc$", cols)]

  if (length(cols) < 2L) {
    stop("Fewer than two columns matched '", pattern,
         "'; there is nothing to compare.", call. = FALSE)
  }

  m <- stats::cor(as.matrix(tab[, cols, drop = FALSE]), method = "spearman")
  labs <- mechanism_label(sub("_bym2$", "", cols))
  dimnames(m) <- list(labs, labs)

  idx  <- which(upper.tri(m), arr.ind = TRUE)
  long <- tibble::tibble(
    mechanism_a = labs[idx[, "row"]],
    mechanism_b = labs[idx[, "col"]],
    rho         = m[idx]
  )

  list(matrix = m, long = long[order(-long[["rho"]]), , drop = FALSE])
}


#' Package versions actually used, for the software paragraph
#'
#' The methods list packages by hand, which drifts. This reads the versions from
#' the session so the paragraph can be regenerated.
#'
#' @param packages Packages to report.
#' @return A tibble: `package`, `version`.
#' @export
session_table <- function(packages = c(
  "dplyr", "tidyr", "purrr", "janitor", "readr", "readxl",
  "sf", "spdep", "terra", "exactextractr",
  "dodgr", "osmextract",
  "geostan", "rstan", "loo",
  "ggplot2", "gtsummary", "gt",
  "targets", "tarchetypes", "crew", "quarto",
  "arrow", "duckdb", "DBI", "tidygeocoder",
  "devtools", "roxygen2", "renv"
)) {
  v <- vapply(packages, function(p) {
    tryCatch(as.character(utils::packageVersion(p)),
             error = function(e) NA_character_)
  }, character(1))

  out <- tibble::tibble(package = packages, version = unname(v))
  out[!is.na(out[["version"]]), , drop = FALSE]
}


#' Stack per-model results, failing loudly on missing ones
#'
#' `dplyr::bind_rows()` on a list of `NULL`s returns a tibble with zero rows
#' **and zero columns**, which is the worst possible failure mode here: the
#' pipeline reports success, the target is stored, and the error only surfaces
#' much later inside the report as "Column `statistic` doesn't exist", pointing
#' at the report rather than at the model that actually failed.
#'
#' This makes the absence explicit and names the models responsible.
#'
#' @param x Named list of one-row data frames, one per model.
#' @param id Name of the identifier column. Default `"model"`.
#' @param require_all Error if any element is missing, rather than warning and
#'   stacking the rest. Default `FALSE`.
#'
#' @return A tibble with an `id` column, one row per non-empty element.
#' @examples
#' stack_by_model(list(M1 = data.frame(x = 1), M2 = data.frame(x = 2)))
#' @export
stack_by_model <- function(x, id = "model", require_all = FALSE) {

  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    stop("`x` must be a fully named list; the names become the `", id,
         "` column.", call. = FALSE)
  }

  empty <- names(x)[vapply(x, function(e) is.null(e) || NROW(e) == 0L,
                           logical(1))]

  if (length(empty) == length(x)) {
    stop("Every element of `x` is NULL or empty (", paste(empty, collapse = ", "),
         "). Those targets did not produce a result - check them with ",
         "targets::tar_meta(fields = error, complete_only = TRUE) - rather ",
         "than debugging the report, which is only where the symptom appears.",
         call. = FALSE)
  }

  if (length(empty)) {
    msg <- paste0(length(empty), " of ", length(x), " model(s) produced no ",
                  "result and are absent from the table: ",
                  paste(empty, collapse = ", "), ".")
    if (require_all) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }

  dplyr::bind_rows(x[setdiff(names(x), empty)], .id = id)
}
