# Descriptive layer ------------------------------------------------------------
#
# Results "Descriptive epidemiology" is written as prose with roughly thirty
# numbers interpolated into it. Every one of them is produced here, as a named
# list, so the report can write `r desc$n_total` and the text can never disagree
# with the pipeline.
#
# The rule throughout: descriptives about DECEDENTS come from `deaths` (one row
# per person, from build_deaths()), never from `mort_count` (one row per death x
# avoidability arm). Counting rows of mort_count would double-count every one of
# the nine split causes.


#' Analysis flow counts
#'
#' The numbers behind a STROBE/RECORD flow diagram, and behind the bias
#' discussion: how many records the register supplied, how many survived each
#' restriction, and how many were lost at each step.
#'
#' The external-causes exclusion is worth quantifying rather than asserting.
#' Discussion says deaths from chapters V-Y are absent from ReNCaM and that this
#' understates avoidable mortality; the size of the gap is the difference
#' between records that failed to match the lookup and records that matched.
#'
#' @param mort_raw Output of [import_mortality()].
#' @param mort_count Output of [preprocess_mortality()].
#' @param mort_count_area `mort_count` restricted to modelled areas.
#' @param area_shp The modelling geography.
#' @param age_threshold Upper age bound used, for the label. Default `75`.
#'
#' @return A tibble: `step`, `n`, `lost`, with one row per stage.
#' @examples
#' \dontrun{
#' flow_counts(mort_raw, mort_count, mort_count_area, area_shp)
#' }
#' @export
flow_counts <- function(mort_raw, mort_count, mort_count_area, area_shp,
                        age_threshold = 75) {

  n_raw       <- nrow(mort_raw)
  n_premature <- sum(mort_raw[["eta"]] < age_threshold, na.rm = TRUE)
  n_avoidable <- dplyr::n_distinct(mort_count[["death_id"]])
  n_in_area   <- dplyr::n_distinct(mort_count_area[["death_id"]])

  out <- tibble::tibble(
    step = c(
      "Deaths in the ReNCaM extract",
      sprintf("Aged under %d at death", age_threshold),
      "Cause of death on the avoidable list",
      "Resident in a modelled area"
    ),
    n = c(n_raw, n_premature, n_avoidable, n_in_area)
  )
  out[["lost"]] <- c(NA_integer_, -diff(out[["n"]]))

  attr(out, "n_areas")     <- nrow(area_shp)
  attr(out, "n_rows_count") <- nrow(mort_count)
  out
}


#' Every scalar the Results prose interpolates
#'
#' One call, one named list. Adding a number to the thesis text means adding it
#' here, which keeps the text and the pipeline in lockstep.
#'
#' Preventable and treatable percentages are deliberately allowed to sum to more
#' than 100: a death from one of the seven causes split across the two lists is
#' *at least partially* both, which is how the draft phrases it. The
#' `avoidability` breakdown, which does sum to 100, is returned alongside.
#'
#' @param deaths One row per decedent, from [build_deaths()].
#' @param arms The long arm table, one row per (death, type, mechanism). Taken
#'   from the `"arms"` attribute of `deaths` by default. Required, because the
#'   weighted mechanism totals cannot be recovered from the collapsed mechanism
#'   columns: a breast cancer death reads "Screening + Tertiary prevention" at
#'   weight 1 there, when it is half a death in each of two real categories.
#' @param n_top Number of leading causes to return. Default `5`.
#' @param n_years Study duration in years, for annualised figures. Default `3`.
#'
#' @return A named list of scalars and small tibbles. See the examples for the
#'   naming convention.
#' @examples
#' \dontrun{
#' d <- desc_stats(deaths_area)
#' d$n_total; d$median_age; d$top_causes
#' }
#' @importFrom rlang .data
#' @export
desc_stats <- function(deaths, arms = attr(deaths, "arms"),
                       n_top = 5, n_years = 3) {

  require_cols(deaths, c("death_id", "eta", "preventable", "treatable",
                         "avoidability", "group", "cause"), "deaths")

  n <- nrow(deaths)
  pc <- function(x) 100 * x / n

  age_q <- stats::quantile(deaths[["eta"]], c(0.25, 0.5, 0.75), na.rm = TRUE,
                           names = FALSE)

  count_pct <- function(var) {
    tab <- deaths |>
      dplyr::count(.lvl = as.character(.data[[var]]), name = "n") |>
      dplyr::mutate(pct = pc(.data[["n"]])) |>
      dplyr::arrange(dplyr::desc(.data[["n"]]))
    names(tab)[1] <- var
    tab
  }

  by_group <- count_pct("group")
  by_cause <- count_pct("cause")

  # Named lookup so the prose can ask for a group by its label without knowing
  # its rank, and get 0 rather than an error if it is absent from this extract.
  grp <- function(label) {
    i <- match(label, by_group[["group"]])
    if (is.na(i)) list(n = 0L, pct = 0) else
      list(n = by_group[["n"]][i], pct = by_group[["pct"]][i])
  }

  sex_tab <- if ("sex" %in% names(deaths)) {
    dplyr::count(deaths, .lvl = as.character(.data[["sex"]]), name = "n")
  } else {
    dplyr::count(deaths, .lvl = as.character(.data[["sesso"]]), name = "n")
  }
  sex_n <- stats::setNames(sex_tab[["n"]], sex_tab[[".lvl"]])
  n_men   <- unname(sum(sex_n[names(sex_n) %in% c("Male", "1")]))
  n_women <- unname(sum(sex_n[names(sex_n) %in% c("Female", "2")]))

  # Mechanism totals are WEIGHTED and must come from the long arm table, never
  # from the collapsed mechanism columns of `deaths`. A female breast cancer
  # death carries the single label "Screening + Tertiary prevention" at weight
  # 1; counting that would give one death to a category that does not exist and
  # zero to the two that do. The arm table keeps the two halves apart.
  if (is.null(arms)) {
    stop("No arm table available. Pass `arms = build_death_arms(mort_count)`, ",
         "or use a `deaths` object straight from build_deaths(), which ",
         "attaches one. Weighted mechanism totals cannot be recovered from ",
         "the collapsed columns.", call. = FALSE)
  }

  mech <- arms |>
    dplyr::filter(!is.na(.data[["mechanism"]])) |>
    dplyr::count(mechanism = as.character(.data[["mechanism"]]),
                 wt = .data[["weight"]], name = "deaths") |>
    dplyr::mutate(pct = 100 * .data[["deaths"]] / n) |>
    dplyr::arrange(dplyr::desc(.data[["deaths"]]))

  # The whole point of weighting: this must hold.
  if (abs(sum(mech[["deaths"]]) - n) > 1e-6) {
    stop("Weighted mechanism totals sum to ", sum(mech[["deaths"]]),
         " but there are ", n, " deaths.", call. = FALSE)
  }

  mch <- function(label) {
    i <- match(label, mech[["mechanism"]])
    if (is.na(i)) list(n = 0, pct = 0) else
      list(n = mech[["deaths"]][i], pct = mech[["pct"]][i])
  }

  avoid <- deaths |>
    dplyr::count(.data[["avoidability"]], name = "n") |>
    dplyr::mutate(pct = pc(.data[["n"]]))

  list(
    n_total        = n,
    n_per_year     = n / n_years,
    n_years        = n_years,

    median_age     = age_q[2],
    iqr_low        = age_q[1],
    iqr_high       = age_q[3],

    n_preventable  = sum(deaths[["preventable"]] == "Yes"),
    pct_preventable = pc(sum(deaths[["preventable"]] == "Yes")),
    n_treatable    = sum(deaths[["treatable"]] == "Yes"),
    pct_treatable  = pc(sum(deaths[["treatable"]] == "Yes")),
    n_both         = sum(deaths[["avoidability"]] == "Preventable and treatable"),
    pct_both       = pc(sum(deaths[["avoidability"]] == "Preventable and treatable")),
    avoidability   = avoid,

    n_men          = n_men,
    pct_men        = pc(n_men),
    n_women        = n_women,
    pct_women      = pc(n_women),

    by_group       = by_group,
    by_cause       = by_cause,
    top_causes     = utils::head(by_cause, n_top),
    pct_top_causes = sum(utils::head(by_cause[["pct"]], n_top)),

    cancer         = grp("Cancer"),
    circulatory    = grp("Diseases of the circulatory system"),
    infectious     = grp("Infectious diseases"),
    respiratory    = grp("Diseases of the respiratory system"),
    alcohol_drug   = grp("Alcohol-related and drug-related deaths"),

    by_mechanism   = mech,
    lifestyle      = mch("Lifestyle and NCDs"),
    tertiary       = mch("Tertiary prevention"),
    screening      = mch("Screening"),
    environment    = mch("Environment and Safety"),
    mental_health  = mch("Mental Health Services"),
    immunisation   = mch("Immunisation and Prophylaxis"),
    reproductive   = mch("Gender and Reproductive Health")
  )
}


#' Table 1
#'
#' Decedent characteristics, overall and by avoidability category. One row per
#' person, so the denominators are people rather than death-by-arm records.
#'
#' `avoidability` rather than `preventable` as the stratifier: it has three
#' mutually exclusive levels that sum to the total, whereas "preventable
#' yes/no" and "treatable yes/no" overlap for the seven split causes and would
#' give a table whose columns do not add up.
#'
#' @param deaths One row per decedent, from [build_deaths()].
#' @param by Stratifying column. Default `"avoidability"`; `NULL` for an
#'   unstratified table.
#' @param n_causes Number of individual causes to show before lumping the rest
#'   into "Other". Default `15`.
#'
#' @return A `gtsummary` table.
#' @examples
#' \dontrun{
#' tbl_one(deaths_area)
#' }
#' @export
tbl_one <- function(deaths, by = "avoidability", n_causes = 15) {

  vars <- c("eta", "sex", "group", "cause", "mechanism_any")
  vars <- vars[vars %in% names(deaths)]

  d <- deaths |>
    dplyr::select(dplyr::all_of(unique(c(vars, by)))) |>
    dplyr::mutate(
      cause = forcats::fct_lump_n(factor(.data[["cause"]]), n_causes,
                                  other_level = "Other avoidable cause"),
      group = forcats::fct_infreq(factor(.data[["group"]]))
    )

  labs <- list(
    eta           = "Age at death (years)",
    sex           = "Sex",
    group         = "Cause group",
    cause         = "Cause of death",
    mechanism_any = "Service function that could have avoided the death"
  )
  labs <- labs[names(labs) %in% names(d)]

  tbl <- gtsummary::tbl_summary(
    d,
    by      = dplyr::all_of(by),
    label   = labs,
    missing = "no",
    type    = list(eta ~ "continuous"),
    statistic = list(
      gtsummary::all_continuous()  ~ "{median} ({p25}\u2013{p75})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(gtsummary::all_categorical() ~ c(0, 1))
  )

  if (!is.null(by)) tbl <- gtsummary::add_overall(tbl, last = TRUE)

  tbl |>
    gtsummary::bold_labels() |>
    gtsummary::modify_footnote(
      gtsummary::all_stat_cols() ~
        paste0("Median (IQR) or n (%). One row per decedent. Deaths from a ",
               "cause split between the preventable and treatable lists appear ",
               "in the 'Preventable and treatable' column, so the columns are ",
               "mutually exclusive and sum to the total.")
    )
}


#' Distribution of the model covariates
#'
#' The numbers behind Results "Covariates distribution and spatial structure":
#' range, mean, SD and median for every covariate and for the expected counts.
#'
#' @param geo The modelling `sf`.
#' @param vars Named character vector: names are display labels, values are
#'   column names.
#'
#' @return A tibble: `variable`, `n`, `n_missing`, `mean`, `sd`, `min`, `p25`,
#'   `median`, `p75`, `max`.
#' @export
covariate_summary <- function(geo, vars = NULL) {

  tab <- sf::st_drop_geometry(geo)

  if (is.null(vars)) {
    candidates <- c(
      "Expected deaths"            = "total_exp",
      "Observed deaths"            = "total_obs",
      "Standardised mortality ratio" = "total_smr",
      "Population (mean annual)"   = "population",
      "Deprivation Index"          = "di_score",
      "IVSM"                       = "ivsm",
      "GP-equivalents per 1,000"   = "gp_density",
      "Travel time to hub (min)"   = "t_hub_mean",
      "Travel time to centre (min)" = "t_centre_mean",
      "Population > 45 min from hub" = "pop_share_over_45min_hub"
    )
    vars <- candidates[candidates %in% names(tab)]

    # Pollution columns carry the exposure year (no2_2023, pm25_2023), so they
    # are found by pattern rather than by name. The display label deliberately
    # omits the year: the report refers to "NO2 (ug/m3)", and hardcoding a year
    # into the label would break that the moment the exposure year changed.
    vars <- c(vars, exposure_columns(tab))
  }

  rows <- lapply(seq_along(vars), function(i) {
    x <- as.numeric(tab[[vars[i]]])
    q <- stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    tibble::tibble(
      variable  = names(vars)[i],
      column    = unname(vars[i]),
      n         = sum(!is.na(x)),
      n_missing = sum(is.na(x)),
      mean      = mean(x, na.rm = TRUE),
      sd        = stats::sd(x, na.rm = TRUE),
      min       = min(x, na.rm = TRUE),
      p25       = q[1],
      median    = q[2],
      p75       = q[3],
      max       = max(x, na.rm = TRUE)
    )
  })

  do.call(rbind, rows)
}


#' Appendix B: the avoidable-cause coding table
#'
#' The full lookup as it was actually applied, with the local adaptations
#' highlighted. Built from `avoidable_lookup_v3.csv` itself rather than
#' maintained by hand, so it cannot fall out of step with the classification the
#' analysis used - which is the whole reason the appendix exists.
#'
#' ICD-10 keys are collapsed into ranges per cause, because 968 individual rows
#' is not an appendix anyone reads.
#'
#' @param lookup Path to the lookup CSV, or the data frame itself.
#' @param adapted_only Show only the rows that depart from the published
#'   OECD/Eurostat classification. Default `FALSE`.
#'
#' @return A tibble ready for `knitr::kable()`.
#' @export
tbl_lookup <- function(lookup, adapted_only = FALSE) {

  if (is.character(lookup)) {
    lookup <- readr::read_csv(lookup, show_col_types = FALSE)
  }
  require_cols(lookup, c("key", "group", "cause", "type", "mechanism",
                         "weight"), "lookup")

  flag_col <- if ("flag" %in% names(lookup)) "flag" else NULL

  out <- lookup |>
    dplyr::group_by(.data[["group"]], .data[["cause"]], .data[["type"]],
                    .data[["mechanism"]], .data[["weight"]]) |>
    dplyr::summarise(
      icd10   = collapse_codes(.data[["key"]]),
      n_codes = dplyr::n(),
      note    = if (is.null(flag_col)) NA_character_ else {
        u <- unique(stats::na.omit(.data[["flag"]]))
        u <- u[nzchar(u)]
        if (length(u)) paste(u, collapse = "; ") else NA_character_
      },
      .groups = "drop"
    ) |>
    dplyr::arrange(.data[["group"]], .data[["cause"]], .data[["type"]])

  if (adapted_only) {
    out <- dplyr::filter(out, !is.na(.data[["note"]]))
  }
  out
}


#' Collapse a set of ICD-10 keys into readable ranges
#'
#' `c("A00","A01","A02","A05")` becomes `"A00-A02, A05"`.
#'
#' @param keys Character vector of ICD-10 prefixes.
#' @return A length-1 string.
#' @examples
#' collapse_codes(c("A00", "A01", "A02", "A05"))
#' @export
collapse_codes <- function(keys) {
  k <- sort(unique(as.character(keys)))
  if (!length(k)) return(NA_character_)

  letter <- substr(k, 1, 1)
  num    <- suppressWarnings(as.numeric(substring(k, 2)))

  # Anything non-standard (4-character carve-outs, odd keys) is listed rather
  # than forced into a range, where it would be silently absorbed.
  if (anyNA(num)) return(paste(k, collapse = ", "))

  grp   <- cumsum(c(TRUE, diff(num) != 1 | letter[-1] != letter[-length(letter)]))
  parts <- vapply(split(k, grp), function(g) {
    if (length(g) == 1L) g else paste0(g[1], "\u2013", g[length(g)])
  }, character(1))

  paste(unname(parts), collapse = ", ")
}


#' Drop areas with no resident population
#'
#' Milan's NIL 8 (Parco Sempione) is a park. It appears in the NIL shapefile but
#' not in the population file, so it has no denominator, no expected count, and
#' no defined `log(expected)` offset. Left in, it propagates an `NA` into the
#' adjacency graph and every model fails on a row that was never meant to be
#' analysed.
#'
#' Removing it here, once, is what makes the analysis-unit count in the methods
#' true: 193 comuni, minus Milan, plus its 87 populated NILs, is 279.
#'
#' @param area_shp The full modelling geography.
#' @param pop_area_table Population table keyed on `area`.
#' @param min_pop Minimum population for an area to be retained. Default `1`.
#'
#' @return `area_shp` with unpopulated areas removed, and the dropped keys
#'   recorded in the `"dropped_areas"` attribute.
#' @export
drop_unpopulated_areas <- function(area_shp, pop_area_table, min_pop = 1) {

  require_cols(area_shp, "area", "area_shp")
  require_cols(pop_area_table, c("area", "numero"), "pop_area_table")

  pops <- stats::aggregate(
    list(pop = as.numeric(pop_area_table[["numero"]])),
    by = list(area = pop_area_table[["area"]]),
    FUN = sum, na.rm = TRUE
  )

  keep_keys <- pops[["area"]][pops[["pop"]] >= min_pop]
  dropped   <- setdiff(area_shp[["area"]], keep_keys)

  if (length(dropped)) {
    message("Dropping ", length(dropped), " unpopulated area(s): ",
            paste(dropped, collapse = ", "))
  }

  # st_as_sf() then dplyr::filter(), NOT `[`. targets' rds storage can hand back
  # an sf whose class or sf_column attribute no longer lines up with the
  # columns, and `[` propagates that damage instead of repairing it: the result
  # keeps the sf class but loses the geometry registration, so every downstream
  # st_transform() fails with either "sf_column does not point to a geometry
  # column" or "missing crs" - a long way from here, in stroke_aoi and
  # section_xwalk. st_as_sf() re-registers the geometry; filter() preserves it.
  out <- sf::st_as_sf(area_shp)
  out <- dplyr::filter(out, .data[["area"]] %in% keep_keys)

  if (is.na(sf::st_crs(out))) {
    stop("`area_shp` has no CRS after filtering. Every spatial target ",
         "downstream will fail; fix the CRS at build_area_shp().", call. = FALSE)
  }

  attr(out, "dropped_areas") <- dropped
  out
}


#' Study-area size, for the methods paragraph
#'
#' The counts the methods section states and which should never be typed by
#' hand: how many areal units, how many are comuni, how many are NILs, and the
#' resident population they cover.
#'
#' @param area_shp The modelling geography.
#' @param pop_area_table Population table keyed on `area`.
#' @param pop_year Year(s) to sum. A vector gives person-years, so the mean
#'   annual population is reported alongside.
#' @param milan_code Municipality code whose areas are NILs.
#'
#' @return A one-row tibble.
#' @export
summarise_study_area <- function(area_shp, pop_area_table,
                                 pop_year = 2022:2024,
                                 milan_code = "015146") {

  is_nil <- grepl(paste0("^", milan_code, "_"), area_shp[["area"]])

  pt <- pop_area_table[
    as.character(pop_area_table[["anno"]]) %in% as.character(pop_year), ,
    drop = FALSE
  ]
  person_years <- sum(as.numeric(pt[["numero"]]), na.rm = TRUE)
  n_years      <- length(unique(pop_year))

  areas_km2 <- tryCatch(
    as.numeric(sum(sf::st_area(sf::st_as_sf(area_shp)))) / 1e6,
    error = function(e) NA_real_
  )

  tibble::tibble(
    n_areas       = nrow(area_shp),
    n_comuni      = sum(!is_nil),
    n_nil         = sum(is_nil),
    n_comuni_total = sum(!is_nil) + 1L,   # + Milan, which the NILs replace
    person_years  = person_years,
    mean_annual_population = person_years / n_years,
    area_km2      = areas_km2,
    years         = paste(range(pop_year), collapse = "\u2013")
  )
}


#' Correlation matrix among the model covariates
#'
#' Results notes that NO2 and PM2.5 were "strongly correlated and mutually
#' attenuating when fitted jointly", and the pollutant selection rule turns on
#' collinearity with deprivation. Both claims need a number.
#'
#' Spearman by default: what matters for a covariate competing with a spatial
#' random effect is whether it orders the areas the same way, not whether the
#' relationship is linear.
#'
#' @param geo The modelling `sf`.
#' @param vars Columns to correlate. Sensible defaults are detected.
#' @param method Correlation method. Default `"spearman"`.
#'
#' @return A list: `matrix` and a tidy `long` tibble of unique pairs ordered by
#'   absolute correlation.
#' @export
covariate_correlation <- function(geo, vars = NULL, method = "spearman") {

  tab <- sf::st_drop_geometry(geo)

  if (is.null(vars)) {
    candidates <- c("di_score", "ivsm", "gp_density",
                    "t_hub_mean", "t_centre_mean", "total_smr")
    vars <- c(candidates[candidates %in% names(tab)],
              unname(exposure_columns(tab)))
  }
  if (length(vars) < 2L) {
    stop("Fewer than two covariates available to correlate.", call. = FALSE)
  }

  m <- stats::cor(as.matrix(tab[, vars, drop = FALSE]),
                  method = method, use = "pairwise.complete.obs")

  idx  <- which(upper.tri(m), arr.ind = TRUE)
  long <- tibble::tibble(
    var_a = rownames(m)[idx[, "row"]],
    var_b = colnames(m)[idx[, "col"]],
    rho   = m[idx]
  )

  list(matrix = m,
       long   = long[order(-abs(long[["rho"]])), , drop = FALSE],
       method = method)
}


#' Find the pollutant exposure columns, whatever year they carry
#'
#' `add_pollution()` names its output after the exposure year - `no2_2023`,
#' `pm25_2023` - and also attaches the 2013 columns, which are diagnostic and
#' must not appear in a covariate table as though they were exposures. This
#' picks the most recent year present for each pollutant and gives it a
#' year-free display label, so the report's prose can refer to "NO2 (ug/m3)"
#' without the exposure year being wired into the text.
#'
#' @param tab A data frame or dropped-geometry `sf`.
#' @param pollutants Named vector: names are display stems, values are column
#'   prefixes.
#'
#' @return A named character vector of column names, suitable for the `vars`
#'   argument of [covariate_summary()]. Empty if no pollutant column is present.
#' @examples
#' exposure_columns(data.frame(no2_2013 = 1, no2_2023 = 2, pm25_2023 = 3))
#' @export
exposure_columns <- function(tab,
                             pollutants = c("NO2" = "no2", "PM2.5" = "pm25")) {
  out <- character(0)

  for (i in seq_along(pollutants)) {
    hits <- grep(paste0("^", pollutants[i], "_[0-9]{4}$"), names(tab),
                 value = TRUE)
    if (!length(hits)) next
    # most recent year wins: the earlier one is the rank-stability diagnostic
    yr   <- as.integer(sub(".*_", "", hits))
    pick <- hits[which.max(yr)]
    out[paste0(names(pollutants)[i], " (ug/m3)")] <- pick
  }
  out
}
