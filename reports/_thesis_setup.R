# Formatting helpers for thesis_results.qmd ------------------------------------
#
# These live in the report rather than in R/ because they are presentation, not
# analysis: they decide how many decimal places a number gets, not what the
# number is. Nothing here computes anything.
#
# The one rule: every number that reaches the prose goes through one of these,
# so the thesis never has 3.4% in one sentence and 3.42% in the next.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Integer with thousands separators.
fmt_n <- function(x, digits = 0) {
  formatC(round(as.numeric(x), digits), format = "f", digits = digits,
          big.mark = ",")
}

#' Fixed-decimal number.
fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "\u2014",
         formatC(as.numeric(x), format = "f", digits = digits, big.mark = ","))
}

#' Percentage, one decimal place, with the sign.
fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "\u2014", paste0(fmt_num(x, digits), "%"))
}

#' p-value, never reported as exactly zero.
fmt_p <- function(p) {
  ifelse(is.na(p), "p = \u2014",
         ifelse(p < 0.001, "p < 0.001", paste0("p = ", fmt_num(p, 3))))
}

#' Estimate with credible interval.
fmt_ci <- function(est, low, high, digits = 2) {
  sprintf("%s (95%% CrI %s to %s)",
          fmt_num(est, digits), fmt_num(low, digits), fmt_num(high, digits))
}

# --- accessors into covariate_summary() ---------------------------------------

.cov_row <- function(tbl, variable) {
  i <- match(variable, tbl[["variable"]])
  if (is.na(i)) NULL else tbl[i, , drop = FALSE]
}

rng <- function(tbl, variable, digits = 2) {
  r <- .cov_row(tbl, variable)
  if (is.null(r)) return("\u2014")
  sprintf("%s to %s", fmt_num(r[["min"]], digits), fmt_num(r[["max"]], digits))
}

mean_of <- function(tbl, variable, digits = 2) {
  r <- .cov_row(tbl, variable); if (is.null(r)) "\u2014" else fmt_num(r[["mean"]], digits)
}
sd_of <- function(tbl, variable, digits = 2) {
  r <- .cov_row(tbl, variable); if (is.null(r)) "\u2014" else fmt_num(r[["sd"]], digits)
}
median_of <- function(tbl, variable, digits = 2) {
  r <- .cov_row(tbl, variable); if (is.null(r)) "\u2014" else fmt_num(r[["median"]], digits)
}
min_of <- function(tbl, variable, digits = 2) {
  r <- .cov_row(tbl, variable); if (is.null(r)) "\u2014" else fmt_num(r[["min"]], digits)
}
max_of <- function(tbl, variable, digits = 2) {
  r <- .cov_row(tbl, variable); if (is.null(r)) "\u2014" else fmt_num(r[["max"]], digits)
}

# Reciprocal of a covariate_summary() statistic, on a per-`scale` basis: the
# primary-care indicator is a supply density (GP-equivalents per 1,000), but the
# prose also quotes the caseload (residents per GP-equivalent), and the two are
# reciprocals. Inverting `min` gives the LARGEST caseload, so the arguments are
# deliberately crossed at the call site rather than here.
recip_of <- function(tbl, variable, stat, scale = 1000, digits = 0) {
  r <- .cov_row(tbl, variable)
  if (is.null(r) || !is.finite(r[[stat]]) || r[[stat]] == 0) return("\u2014")
  fmt_n(scale / r[[stat]], digits)
}

# Tolerant of the exposure-year suffix: cor_between(cc, "no2", "pm25") finds
# no2_2023 and pm25_2023 without the year having to be written into the prose.
cor_between <- function(cc, a, b) {
  m <- cc$matrix
  resolve <- function(nm) {
    if (nm %in% rownames(m)) return(nm)
    hit <- grep(paste0("^", nm, "(_[0-9]{4})?$"), rownames(m), value = TRUE)
    if (length(hit)) hit[which.max(nchar(hit))] else NA_character_
  }
  a <- resolve(a); b <- resolve(b)
  if (is.na(a) || is.na(b)) return(NA_real_)
  m[a, b]
}

# --- prose generators ---------------------------------------------------------
#
# These write the sentences whose SHAPE depends on the result: "was" versus
# "was not", "three of seven" versus "all seven". Writing them by hand is how a
# results section ends up asserting something the numbers contradict.

top_causes_sentence <- function(top) {
  paste(
    sprintf("%s (%s, %s)", top[["cause"]], fmt_n(top[["n"]]),
            fmt_pct(top[["pct"]])),
    collapse = "; "
  )
}

convergence_sentence <- function(diag) {
  worst_rhat <- max(diag[["rhat_max"]], na.rm = TRUE)
  min_ess    <- min(diag[["ess_min"]], na.rm = TRUE)
  tot_div    <- sum(diag[["n_divergent"]], na.rm = TRUE)
  bad        <- diag[["model"]][diag[["rhat_max"]] > 1.01]

  if (!length(bad) && tot_div == 0) {
    sprintf(paste0("All models converged: the largest rank-normalised R-hat ",
                   "across every monitored quantity in every model was %s, the ",
                   "smallest effective sample size was %s, and there were no ",
                   "divergent transitions."),
            fmt_num(worst_rhat, 4), fmt_n(min_ess))
  } else {
    sprintf(paste0("The largest rank-normalised R-hat was %s and the smallest ",
                   "effective sample size %s. %s%s"),
            fmt_num(worst_rhat, 4), fmt_n(min_ess),
            if (length(bad))
              sprintf("R-hat exceeded 1.01 in %s. ",
                      paste(bad, collapse = ", ")) else "",
            if (tot_div > 0)
              sprintf("There were %s divergent transitions in total, so these ",
                      fmt_n(tot_div)) else "")
  }
}

loo_sentence <- function(tbl) {
  best <- tbl[["model"]][1]
  others <- tbl[-1, , drop = FALSE]
  clear <- others[["model"]][others[["verdict"]] == "clear"]
  indist <- others[["model"]][others[["verdict"]] == "indistinguishable"]

  paste0(
    sprintf("%s had the highest expected log predictive density. ", best),
    if (length(indist))
      sprintf("It is indistinguishable from %s (elpd difference within two standard errors). ",
              paste(indist, collapse = ", ")) else "",
    if (length(clear))
      sprintf("%s %s clearly worse (elpd difference more than four standard errors). ",
              paste(clear, collapse = ", "),
              if (length(clear) == 1) "is" else "are") else ""
  )
}

# Writes the whole sentence, because its shape changes with the result: with no
# flagged areas there is nothing to describe, and "exceeded 0.7 for 0 areas;
# none of them, ..." is not a sentence anyone would write by hand.
pareto_sentence <- function(pa, model = "M5") {
  if (is.null(pa) || !nrow(pa)) {
    return(sprintf(paste0("No area had a Pareto-k above 0.7 under %s, so the ",
                          "leave-one-out approximation is reliable across the ",
                          "whole study area and no observation exerts undue ",
                          "influence on the comparison."), model))
  }
  sprintf(paste0("Pareto-k exceeded 0.7 for %d area(s) under %s. Their median ",
                 "expected count was %s, so these are the sparsest units, ",
                 "where the leave-one-out approximation is least reliable and ",
                 "each observation's influence on the fit is correspondingly ",
                 "greater."),
          nrow(pa), model,
          fmt_num(stats::median(pa[["expected"]], na.rm = TRUE), 1))
}

rr_range <- function(aug, value = "bym2_rr") {
  x <- sf::st_drop_geometry(aug)[[value]]
  sprintf("%s to %s", fmt_num(min(x, na.rm = TRUE), 2),
          fmt_num(max(x, na.rm = TRUE), 2))
}

n_exceeding <- function(aug, cutoff = 0.80, value = "bym2_exceed") {
  sum(sf::st_drop_geometry(aug)[[value]] > cutoff, na.rm = TRUE)
}

surface_cor <- function(a, b, value = "bym2_rr") {
  stats::cor(sf::st_drop_geometry(a)[[value]],
             sf::st_drop_geometry(b)[[value]], method = "spearman")
}

esf_verdict <- function(a, b, value = "bym2_rr") {
  rho <- surface_cor(a, b, value)
  if (rho > 0.95) {
    paste0("The two specifications are close to interchangeable at this ",
           "resolution, so the conclusions do not depend on how the spatial ",
           "term is parameterised.")
  } else if (rho > 0.85) {
    paste0("The two specifications largely agree, with some reordering among ",
           "areas whose posteriors are wide; the substantive conclusions are ",
           "unchanged.")
  } else {
    paste0("The two specifications disagree substantially. Given the open ",
           "debate on partitioning smooth covariates against smooth random ",
           "effects, this divergence is itself a finding and should be ",
           "reported as a limit on how firmly the covariate coefficients can ",
           "be interpreted.")
  }
}

mechanism_sentence <- function(mt) {
  modelled <- mt[!is.na(mt[["rho"]]), , drop = FALSE]
  with_excess <- modelled[modelled[["n_exceed"]] > 0, , drop = FALSE]

  paste0(
    sprintf("Of the %d strata, %d show at least one area above the 0.80 probability cutoff",
            nrow(mt), nrow(with_excess)),
    if (nrow(with_excess))
      sprintf(": %s. ",
              paste(sprintf("%s (%d areas)", with_excess[["mechanism"]],
                            with_excess[["n_exceed"]]), collapse = "; "))
    else ". ",
    sprintf("The mixing parameter ranged from %s to %s across strata, so the ",
            fmt_num(min(modelled[["rho"]], na.rm = TRUE), 2),
            fmt_num(max(modelled[["rho"]], na.rm = TRUE), 2)),
    "share of residual variation that is spatially structured differs ",
    "appreciably by service function."
  )
}

stroke_time_sentence <- function(st) {
  hub <- .cov_row(st, "Nearest thrombectomy hub (level II), minutes")
  ctr <- .cov_row(st, "Nearest stroke centre (any level), minutes")
  if (is.null(hub)) return("Travel times were not available.")

  sprintf(paste0("Population-weighted travel time to the nearest thrombectomy ",
                 "hub ranged from %s to %s minutes, with a median of %s; ",
                 "travel time to the nearest stroke centre of any level ranged ",
                 "from %s to %s minutes."),
          fmt_num(hub[["min"]], 1), fmt_num(hub[["max"]], 1),
          fmt_num(hub[["median"]], 1),
          if (is.null(ctr)) "\u2014" else fmt_num(ctr[["min"]], 1),
          if (is.null(ctr)) "\u2014" else fmt_num(ctr[["max"]], 1))
}

control_sentence <- function(ct) {
  get <- function(role) {
    r <- ct[ct[["role"]] == role, , drop = FALSE]
    if (!nrow(r)) NULL else r[1, ]
  }
  tr  <- get("Tracer")
  nco <- get("Negative control outcome")
  pos <- get("Positive control")

  parts <- character(0)
  if (!is.null(tr)) {
    parts <- c(parts, sprintf(
      "The coefficient for travel time to the nearest thrombectomy hub on %s was %s.",
      tolower(if (is.null(tr[["outcome"]])) "the tracer outcome" else tr[["outcome"]]),
      fmt_ci(tr[["estimate"]], tr[["ci_low"]], tr[["ci_high"]], 3)))
  }
  if (!is.null(nco)) {
    parts <- c(parts, sprintf(
      "The corresponding coefficient for the negative control outcome, all-cancer mortality, was %s.",
      fmt_ci(nco[["estimate"]], nco[["ci_low"]], nco[["ci_high"]], 3)))
  }
  if (!is.null(tr) && !is.null(nco)) {
    same_side <- sign(log(tr[["estimate"]])) == sign(log(nco[["estimate"]]))
    similar   <- abs(log(tr[["estimate"]]) - log(nco[["estimate"]])) < 0.05
    parts <- c(parts, if (similar && same_side) {
      paste0("The two are of similar magnitude and direction, which is what ",
             "would be seen if travel time were measuring the character of an ",
             "area rather than access to time-critical care. On this evidence ",
             "the tracer association should not be read as reflecting the ",
             "proposed mechanism.")
    } else if (isTRUE(nco[["crosses_null"]]) && !isTRUE(tr[["crosses_null"]])) {
      paste0("The tracer association is present where the negative control ",
             "outcome shows none, which is the pattern the design was built to ",
             "detect. This is consistent with, but does not establish, an ",
             "effect operating through access to time-critical care.")
    } else {
      paste0("Neither association is clearly distinguishable from the null, so ",
             "the comparison does not discriminate between the mechanism and ",
             "the confounding explanation.")
    })
  }
  if (!is.null(pos)) {
    # Three outcomes, not two. Checking only whether the interval clears the
    # null treats a strongly PROTECTIVE deprivation coefficient as the control
    # passing, when it is the most informative failure the design can produce.
    verdict <- if (isTRUE(pos[["crosses_null"]])) {
      paste0("which is not distinguishable from the null. Because this ",
             "association is well established in the literature, its absence ",
             "here indicates that the exposure measurement or the design's ",
             "power is inadequate, and the null results above are therefore ",
             "uninterpretable.")
    } else if (pos[["estimate"]] < 1) {
      paste0("i.e. deprivation is associated with LOWER lifestyle-preventable ",
             "mortality. This is the opposite of the pre-specified direction ",
             "and of the established literature, so the positive control has ",
             "not merely failed to fire - it has fired backwards. Until this ",
             "is explained, the deprivation exposure cannot be assumed to ",
             "measure deprivation, and no coefficient on it should be ",
             "interpreted substantively.")
    } else {
      paste0("in the expected direction, which supports the design having ",
             "enough power and measurement quality for the null results to ",
             "be interpretable.")
    }
    parts <- c(parts, sprintf(
      paste0("The positive control, deprivation against lifestyle-preventable ",
             "mortality, was %s, %s"),
      fmt_ci(pos[["estimate"]], pos[["ci_low"]], pos[["ci_high"]], 3),
      verdict))
  }
  paste(parts, collapse = " ")
}

software_paragraph <- function(st) {
  v <- stats::setNames(st[["version"]], st[["package"]])
  g <- function(p) if (p %in% names(v)) sprintf("%s %s", p, v[[p]]) else p
  paste0(
    "\n**Suggested wording for Methods, Software.** Analyses were conducted in ",
    R.version.string, " using a `targets` pipeline (", g("targets"),
    "), with data manipulation in ", g("dplyr"), " and ", g("tidyr"),
    "; spatial data handled by ", g("sf"), " and ", g("spdep"),
    "; environmental raster extraction by ", g("terra"), " and ",
    g("exactextractr"), "; road-network routing by ", g("dodgr"), " and ",
    g("osmextract"), "; Bayesian estimation by ", g("geostan"), " and ",
    g("rstan"), "; model comparison by ", g("loo"),
    "; and tables and figures by ", g("gtsummary"), " and ", g("ggplot2"),
    ". The analysis is built as an R package documented with ", g("roxygen2"),
    " and its dependencies pinned with ", g("renv"), ".\n"
  )
}

# --- the numbers card ---------------------------------------------------------

numbers_card <- function(desc, area, excess, varred, moran_crude,
                         moran_resid, flow) {
  m5 <- moran_resid[moran_resid[["model"]] == "M5", , drop = FALSE]

  tibble::tribble(
    ~Quantity,                                   ~Value,
    "Avoidable deaths, 2022-2024",               fmt_n(desc$n_total),
    "Avoidable deaths per year",                 fmt_n(desc$n_per_year),
    "Median age at death (IQR)",                 sprintf("%s (%s-%s)",
                                                         fmt_num(desc$median_age, 0),
                                                         fmt_num(desc$iqr_low, 0),
                                                         fmt_num(desc$iqr_high, 0)),
    "At least partially preventable",            sprintf("%s (%s)", fmt_n(desc$n_preventable), fmt_pct(desc$pct_preventable)),
    "At least partially treatable",              sprintf("%s (%s)", fmt_n(desc$n_treatable), fmt_pct(desc$pct_treatable)),
    "Both",                                      sprintf("%s (%s)", fmt_n(desc$n_both), fmt_pct(desc$pct_both)),
    "Men",                                       sprintf("%s (%s)", fmt_n(desc$n_men), fmt_pct(desc$pct_men)),
    "Women",                                     sprintf("%s (%s)", fmt_n(desc$n_women), fmt_pct(desc$pct_women)),
    "Areal units",                               fmt_n(area$n_areas),
    "  of which comuni",                         fmt_n(area$n_comuni),
    "  of which Milan NILs",                     fmt_n(area$n_nil),
    "Mean annual population",                    fmt_n(area$mean_annual_population),
    "Records in extract",                        fmt_n(flow$n[1]),
    "Records lost to age restriction",           fmt_n(flow$lost[2]),
    "Records not on the avoidable list",         fmt_n(flow$lost[3]),
    "Moran's I, crude SMR",                      sprintf("%s (%s)", fmt_num(moran_crude$statistic, 3), fmt_p(moran_crude$p_value)),
    "Moran's I, M5 residuals",                   if (nrow(m5)) sprintf("%s (%s)", fmt_num(m5$statistic, 3), fmt_p(m5$p_value)) else "\u2014",
    "Residual excess, deaths/year",              fmt_ci(excess$total$excess_per_year, excess$total$ci_low, excess$total$ci_high, 0),
    "Residual excess, % of total",               fmt_pct(excess$total$pct_of_total),
    "Areas flagged (P > 0.80, RR > 1.10)",       fmt_n(excess$n_flagged),
    "Between-area variance reduction M3 to M5",  fmt_pct(varred$reduction$pct_reduction)
  ) |>
    knitr::kable(caption = "Every scalar the Results prose uses.")
}

# One theme for every figure rendered in this document.
ggplot2::theme_set(theme_atlas())

#' Range of a numeric column in an sf or data frame, formatted.
rng_col <- function(x, col, digits = 2) {
  v <- sf::st_drop_geometry(x)[[col]]
  sprintf("%s to %s", fmt_num(min(v, na.rm = TRUE), digits),
          fmt_num(max(v, na.rm = TRUE), digits))
}

#' One value from a per-model table, or an em dash if that model is absent.
#'
#' `tbl$col[tbl$model == "M5"]` returns `numeric(0)` when M5 is missing, which
#' interpolates into the prose as nothing at all - a silently truncated
#' sentence. This returns a visible placeholder instead.
by_model <- function(tbl, model, col, digits = 3, fmt = fmt_num) {
  if (is.null(tbl) || !nrow(tbl) || !"model" %in% names(tbl)) return("\u2014")
  v <- tbl[[col]][tbl[["model"]] == model]
  if (!length(v) || is.na(v[1])) return("\u2014")

  # Not every formatter takes `digits`: fmt_p() decides its own precision,
  # because a p-value is reported as "p < 0.001" rather than to n places.
  # Pass digits only where it is accepted.
  if ("digits" %in% names(formals(fmt))) fmt(v[1], digits) else fmt(v[1])
}

#' "had" / "had not", or an explicit gap when the model is missing.
moran_verdict <- function(tbl, model, alpha = 0.05) {
  if (is.null(tbl) || !nrow(tbl) || !"model" %in% names(tbl)) {
    return("[NOT COMPUTED]")
  }
  p <- tbl[["p_value"]][tbl[["model"]] == model]
  if (!length(p) || is.na(p[1])) "[NOT COMPUTED]" else
    if (p[1] > alpha) "had" else "had not"
}


# --- stroke sub-model prose ---------------------------------------------------

#' One sentence on how much of the cerebrovascular burden the reference outcome
#' actually covers. The dilution arithmetic is deliberately explicit: an access
#' effect can only act on the thrombectomy-eligible subset, and stating that
#' fraction is what stops a null being read as "access does not matter".
stroke_subtype_sentence <- function(sub,
                                    lvo = 0.20, in_window = 0.50,
                                    mortality_reduction = 0.30) {
  if (is.null(sub) || !nrow(sub)) return("Subtype counts were not available.")

  tot  <- sum(sub[["deaths"]])
  i63  <- sum(sub[["deaths"]][grepl("^I63", sub[["subtype"]])])
  haem <- sum(sub[["deaths"]][grepl("^I60", sub[["subtype"]])])
  if (!tot) return("Subtype counts were not available.")

  reachable <- (i63 / tot) * lvo * in_window

  sprintf(paste0(
    "Of %s cerebrovascular deaths at all ages, %s were cerebral infarction ",
    "(%s) and %s were haemorrhagic (%s), which is a different care channel. ",
    "If the exposure acted through thrombectomy alone it could touch roughly ",
    "%s of the reference outcome (%s large-vessel occlusion, %s presenting in ",
    "window), so a %s reduction in mortality within that subset would move ",
    "total I63 mortality by about %s \u2014 below what a design of this size ",
    "can resolve."),
    fmt_n(tot), fmt_n(i63), fmt_pct(100 * i63 / tot),
    fmt_n(haem), fmt_pct(100 * haem / tot),
    fmt_pct(100 * reachable / (i63 / tot)),
    fmt_pct(100 * lvo), fmt_pct(100 * in_window),
    fmt_pct(100 * mortality_reduction),
    fmt_pct(100 * reachable * mortality_reduction))
}

#' One sentence pairing what the smoothing did with how many areas survive it.
i63_smoothing_sentence <- function(sm, exc) {
  val <- function(q) {
    i <- grep(q, sm[["Quantity"]], fixed = TRUE)
    if (!length(i)) "\u2014" else sm[["Value"]][i[1]]
  }
  n80 <- if ("p80" %in% names(exc)) exc[["p80"]] else NA_integer_
  n95 <- if ("p95" %in% names(exc)) exc[["p95"]] else NA_integer_

  sprintf(paste0(
    "With a median of %s expected deaths per area, the raw ratio ranged %s and ",
    "the BYM2-smoothed relative risk %s; %s of the residual variation was ",
    "spatially structured. %s of %s areal units carry at least an 80%% ",
    "posterior probability that all-age I63 mortality runs more than 20%% ",
    "above expectation, and %s reach 95%%."),
    val("Median expected deaths per area"),
    val("Range of raw SMR"),
    val("Range of smoothed relative risk"),
    val("Mixing parameter rho"),
    fmt_n(n80), fmt_n(exc[["n_areas"]]), fmt_n(n95))
}
