# Deprivation index -----------------------------------------------------------
#
# Split out of preprocess_.R for two reasons. It had grown to the point where it
# was the only thing in that file with its own methodological argument, and the
# change of resolution (municipality -> modelling area) needs testing on its own
# rather than as a side effect of a mortality-preprocessing test.
#
# WHAT CHANGED AND WHY
#
# build_deprivation_proxy() aggregated census sections to PROCOM, i.e. to the
# municipality. expand_cov_to_area() then copied Milan's single municipal value
# onto all 87 of its NILs. Roughly a third of the modelling units therefore
# carried no within-city deprivation contrast at all, while the outcome varied
# freely across them. A covariate that is flat where the outcome is not cannot
# explain that variation, so the "deprivation association is minimal" finding
# was partly a property of the measurement rather than of the territory.
#
# build_deprivation() aggregates sections to the MODELLING AREA instead - NIL
# inside Milan, comune everywhere else - which is the resolution the rest of the
# pipeline works at.
#
# WHAT DID NOT CHANGE
#
# The reference distribution. z-scores and quintile cut points are still taken
# from the national distribution across all ~7,900 Italian municipalities, so
# a score of +1 means the same thing it meant before. Only the resolution of the
# study-area values improves. Mixing NIL-sized and comune-sized units into one
# national reference distribution would have changed the meaning of the scale,
# so it is deliberately not done.


#' Auto-detect the census-section identifier column
#'
#' ISTAT has used several names for the section code across releases
#' (`SEZ2011`, `SEZ21_ID`, `SEZ`, `IDSEZ`, `SEZIONE`), and the 2021 shapefile
#' and the 2023 permanent-census workbooks do not necessarily agree. Rather
#' than hardcode one and fail obscurely on the next release, look for any of
#' them and say what was found.
#'
#' @param x A data frame or `sf`.
#' @param candidates Character vector of acceptable column names, tried in
#'   order.
#' @param what Description of `x` for the error message.
#'
#' @return The name of the first matching column.
#' @export
detect_section_key <- function(x,
                               candidates = c("SEZ2021", "SEZ21_ID", "SEZ2011",
                                              "IDSEZ", "SEZ", "SEZIONE",
                                              "COD_SEZ", "sez_code"),
                               what = "census table") {
  hit <- candidates[candidates %in% names(x)]
  if (!length(hit)) {
    stop("No census-section identifier found in `", what, "`.\n",
         "  Tried: ", paste(candidates, collapse = ", "), "\n",
         "  Available: ", paste(names(x), collapse = ", "), "\n",
         "  Pass the correct name via the `sez_key` argument.",
         call. = FALSE)
  }
  hit[1]
}


#' Crosswalk census sections to modelling areas
#'
#' Assigns every census section polygon to the modelling area that contains it,
#' by spatial join on the section's point-on-surface. The same device
#' [build_section_points()] uses for the stroke origins, and for the same
#' reason: reconstructing the composite `area` key from ISTAT codes would break
#' the moment the key definition changes, whereas a spatial join cannot drift
#' out of step with `area_shp`.
#'
#' Point-on-surface rather than centroid: a centroid can fall outside a concave
#' section, which for a section wrapping a city block would place it in the
#' neighbouring NIL.
#'
#' @param sez_shp `sf` of ISTAT census section polygons.
#' @param area_shp `sf` of modelling areas carrying an `area` column.
#' @param sez_key Name of the section identifier column in `sez_shp`. Detected
#'   automatically when `NULL`.
#'
#' @return A tibble with `sez_code` and `area`, one row per section that falls
#'   inside a modelling area.
#' @seealso [build_deprivation()], [build_section_points()]
#' @export
build_section_area_xwalk <- function(sez_shp, area_shp, sez_key = NULL) {

  if (is.null(sez_key)) {
    sez_key <- detect_section_key(sez_shp, what = "sez_shp")
  }
  require_cols(area_shp, "area", "area_shp")

  # Re-register both layers as sf before touching them. targets' rds storage
  # can return an sf whose geometry registration is stale, and st_transform()
  # then fails with "sf_column does not point to a geometry column".
  sez_shp  <- sf::st_as_sf(sez_shp)
  area_shp <- sf::st_as_sf(area_shp)

  crs_m <- 32632L

  sez <- sez_shp |>
    sf::st_transform(crs_m) |>
    sf::st_make_valid()

  pts <- suppressWarnings(sf::st_point_on_surface(sez))

  pts <- sf::st_join(
    pts,
    area_shp |>
      sf::st_transform(crs_m) |>
      sf::st_make_valid() |>
      dplyr::select("area"),
    join = sf::st_within
  )

  out <- tibble::tibble(
    sez_code = as.character(sf::st_drop_geometry(pts)[[sez_key]]),
    area     = pts[["area"]]
  ) |>
    dplyr::filter(!is.na(.data[["area"]]))

  n_area <- dplyr::n_distinct(out[["area"]])
  message(nrow(out), " census sections mapped onto ", n_area,
          " modelling areas (of ", nrow(area_shp), " in area_shp).")

  if (n_area < nrow(area_shp)) {
    missing <- setdiff(area_shp[["area"]], out[["area"]])
    warning(length(missing), " modelling area(s) contain no census section ",
            "and will have a missing deprivation score: ",
            paste(utils::head(missing, 10), collapse = ", "),
            if (length(missing) > 10) ", ..." else "", call. = FALSE)
  }

  out
}


#' Form the four deprivation indicators from aggregated census counts
#'
#' Internal. Kept separate from the aggregation so the same arithmetic is
#' provably applied to the national reference units and to the study areas -
#' if the two ever diverged, the z-scores would be meaningless.
#'
#' All four are oriented so that a higher value means more disadvantage, which
#' is what licenses summing the z-scores without weights.
#'
#' @param x A data frame of summed census counts.
#' @return `x` with `edu_low`, `nonemp`, `foreign`, `crowd` and `pop` added.
#' @noRd
.deprivation_indicators <- function(
    x,
    vintage = c("2023", "2011"),
    bounds = list(
      edu_low = c(0, 1),      # proportions by construction
      nonemp  = c(0, 1),
      foreign = c(0, 1),
      crowd   = c(1, 8)       # occupants per occupied dwelling
    )) {

  vintage <- match.arg(vintage)

  # ONE function computes both census vintages, deliberately. The 2011-vs-2023
  # comparison is a test of whether the TERRITORY changed; if the two years
  # went through two different implementations, any difference would confound
  # that with the two implementations disagreeing, and the comparison would be
  # worthless. Only the variable codes differ below - the arithmetic, the
  # denominator guards, the clamping and the standardisation are shared.
  #
  # Variable mapping, permanent census 2023 -> general census 2011:
  #
  #   edu_low  P86+P87+P88 / P83     ->  P49+P50+P51+P52 / P46
  #            (lower secondary or less; the 2023 numerator gives a median of
  #            0.44, which matches "at most licenza media" rather than "at most
  #            elementary", so the 2011 numerator includes P49)
  #   nonemp   (P17:P26 - P101)      ->  (P17:P26 - P61)
  #            (the 15-64 age bands carry the SAME codes P17-P26 in both
  #            releases. P61 is employed aged 15+, so a small number of working
  #            over-65s sit in the numerator; this understates non-employment
  #            by roughly a percentage point and is documented rather than
  #            silently corrected)
  #   foreign  ST1 / P1              ->  ST1 / P1        (identical)
  #   crowd    P1 / A2               ->  P1 / A2         (identical)

  band_cols <- paste0("P", 17:26)   # 15-19 ... 60-64, same codes both vintages

  spec <- if (vintage == "2023") {
    list(edu_num = c("P86", "P87", "P88"), edu_den = "P83", emp = "P101")
  } else {
    list(edu_num = c("P49", "P50", "P51", "P52"), edu_den = "P46", emp = "P61")
  }

  require_cols(x, c("P1", spec$edu_den, spec$edu_num, spec$emp, "ST1", "A2",
                    band_cols), paste("census counts,", vintage))

  pop_1564 <- rowSums(dplyr::select(x, dplyr::all_of(band_cols)), na.rm = TRUE)

  # Guard every denominator. A NIL that is almost entirely park or industrial
  # land can have A2 == 0 or P83 == 0; x/0 would give Inf, which then poisons
  # the mean and SD of the whole national distribution.
  safe <- function(num, den) dplyr::if_else(den > 0, num / den, NA_real_)

  edu_num <- rowSums(dplyr::select(x, dplyr::all_of(spec$edu_num)),
                     na.rm = TRUE)

  out <- dplyr::mutate(
    x,
    pop     = .data[["P1"]],
    edu_low = safe(edu_num, .data[[spec$edu_den]]),
    nonemp  = safe(pop_1564 - .data[[spec$emp]], pop_1564),
    foreign = safe(.data[["ST1"]], .data[["P1"]]),
    crowd   = safe(.data[["P1"]], .data[["A2"]])
  )

  # A POSITIVE denominator is not the same as a PLAUSIBLE one, and the
  # difference is not cosmetic. An area with 951 residents and one recorded
  # occupied dwelling - a barracks, a care home, a student residence, or a
  # section where dwellings simply were not captured - yields crowd = 951
  # against a national median of 2.3. Summed into the score, that single area
  # reached di_score = 4126 while every other area sat between -6 and 19.
  # Because add_covariate() then standardises across areas, the whole real
  # deprivation range compressed into 0.1 SD and `di_score_z` became, in
  # effect, a dummy variable for that one NIL. Nothing errored; the only
  # symptom was a deprivation coefficient with the wrong sign.
  #
  # The bounds are substantive, not distributional - the same argument as
  # clamp_supply() for GP list sizes. Three of the four indicators are
  # proportions and cannot leave [0, 1] except through a broken denominator.
  # Occupants per dwelling above `crowd_max` is not a crowded area, it is an
  # area whose dwelling count is wrong.
  for (v in names(bounds)) {
    lo <- bounds[[v]][1]; hi <- bounds[[v]][2]
    bad <- !is.na(out[[v]]) & (out[[v]] < lo | out[[v]] > hi)
    if (any(bad)) {
      warning(sum(bad), " value(s) of `", v, "` fell outside the plausible ",
              "range [", lo, ", ", hi, "] (observed up to ",
              signif(max(out[[v]][bad]), 4),
              ") and were clamped to it. This indicates a broken denominator, ",
              "not an extreme area; check the census codes for those units.",
              call. = FALSE)
    }
    out[[paste0(v, "_clamped")]] <- bad
    out[[v]] <- pmin(pmax(out[[v]], lo), hi)
  }

  out
}


#' Build the deprivation index at modelling-area resolution
#'
#' The permanent-census reformulation of the Italian Deprivation Index
#' (Rosano et al. 2020), computed for every modelling area - Milan's NILs
#' individually, every other comune as a whole - and standardised against the
#' national distribution across municipalities.
#'
#' \strong{This is a proxy, not the validated 2011 five-indicator index.} Two of
#' the original indicators, non-home-ownership and single-parent family, are not
#' released in the permanent census and are replaced by the share of foreign
#' residents; overcrowding is approximated by occupants per occupied dwelling
#' and unemployment by non-employment among the 15-64 population. The four
#' indicators used are therefore low education, non-employment, foreign
#' residents, and crowding. The methods section must say so.
#'
#' \strong{Two-stage standardisation.} The reference distribution is formed at
#' municipal level across every Italian municipality in \code{census}: the four
#' indicators are computed from municipal counts, and their national means and
#' standard deviations are recorded. Those constants are then applied to the
#' study-area indicators. Quintile cut points come from the same national
#' municipal distribution, population-weighted.
#'
#' Doing it this way keeps the scale interpretable. Pooling 87 NILs in with
#' 7,900 municipalities to form one distribution would have changed what a
#' z-score of 1 means, and would have made the index incomparable with every
#' published application of it.
#'
#' @param census Section-level table from [import_census_2023()], carrying a
#'   section identifier and `PROCOM`.
#' @param xwalk Section-to-area crosswalk from [build_section_area_xwalk()].
#' @param sez_key Section identifier column in `census`. Detected when `NULL`.
#' @param vintage Census release the variable codes belong to: `"2023"` (the
#'   permanent census) or `"2011"` (the general census). The indicators, the
#'   guards, the clamping and the two-stage standardisation are identical
#'   across vintages; only the source variable codes differ. Running both years
#'   through one function is what makes a temporal comparison a statement about
#'   the territory rather than about two implementations.
#'
#' @return A tibble, one row per modelling area: `area`, `population`, the four
#'   raw indicators, `di_score` (sum of nationally standardised indicators) and
#'   `di_quintile` (1 = least deprived, on national population-weighted cut
#'   points). The national means and SDs used are attached as the
#'   `"national_reference"` attribute.
#'
#' @examples
#' \dontrun{
#' xwalk <- build_section_area_xwalk(sez_shp, area_shp)
#' di    <- build_deprivation(census_2023, xwalk)
#' attr(di, "national_reference")
#' }
#' @seealso [build_section_area_xwalk()], [wtd_quantile_group()]
#' @importFrom dplyr across all_of group_by left_join mutate select summarise
#' @importFrom rlang .data
#' @export
build_deprivation <- function(census, xwalk, sez_key = NULL,
                              vintage = c("2023", "2011")) {

  vintage <- match.arg(vintage)

  if (is.null(sez_key)) {
    sez_key <- detect_section_key(census, what = "census")
  }
  require_cols(census, c(sez_key, "PROCOM"), "census")
  require_cols(xwalk, c("sez_code", "area"), "xwalk")

  count_cols <- setdiff(names(census), c(sez_key, "PROCOM"))

  # ---- 1. national reference: aggregate sections to municipality -----------
  nat <- census |>
    dplyr::group_by(comune = .data[["PROCOM"]]) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(count_cols),
                                   ~ sum(.x, na.rm = TRUE)),
                     .groups = "drop") |>
    .deprivation_indicators(vintage = vintage)

  ind_names <- c("edu_low", "nonemp", "foreign", "crowd")
  ref <- lapply(ind_names, function(v) {
    c(mean = mean(nat[[v]], na.rm = TRUE), sd = stats::sd(nat[[v]], na.rm = TRUE))
  })
  names(ref) <- ind_names

  if (any(vapply(ref, function(r) !is.finite(r[["sd"]]) || r[["sd"]] == 0,
                 logical(1)))) {
    stop("A national indicator has zero or non-finite standard deviation; ",
         "the census table is probably not national in scope.", call. = FALSE)
  }

  zref <- function(v, nm) (v - ref[[nm]][["mean"]]) / ref[[nm]][["sd"]]

  # national quintile cut points, population-weighted, on the reference units
  nat <- dplyr::mutate(
    nat,
    di_score = zref(.data[["edu_low"]], "edu_low") +
               zref(.data[["nonemp"]],  "nonemp")  +
               zref(.data[["foreign"]], "foreign") +
               zref(.data[["crowd"]],   "crowd")
  )
  ok  <- !is.na(nat[["di_score"]]) & !is.na(nat[["pop"]])
  cuts <- stats::quantile(
    rep(nat[["di_score"]][ok],
        times = pmax(1L, round(nat[["pop"]][ok] / 100))),
    probs = seq(0.2, 0.8, by = 0.2), na.rm = TRUE, names = FALSE
  )

  # ---- 2. study areas: aggregate the SAME sections to `area` ---------------
  area_counts <- census |>
    dplyr::mutate(sez_code = as.character(.data[[sez_key]])) |>
    dplyr::inner_join(xwalk, by = "sez_code") |>
    dplyr::group_by(.data[["area"]]) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(count_cols),
                                   ~ sum(.x, na.rm = TRUE)),
                     .groups = "drop") |>
    .deprivation_indicators(vintage = vintage)

  if (!nrow(area_counts)) {
    stop("No census section matched the crosswalk. The section identifier in ",
         "`census` (", sez_key, ") and in `sez_shp` are probably different ",
         "vintages; pass `sez_key` explicitly on both sides.", call. = FALSE)
  }

  out <- area_counts |>
    dplyr::mutate(
      di_score = zref(.data[["edu_low"]], "edu_low") +
                 zref(.data[["nonemp"]],  "nonemp")  +
                 zref(.data[["foreign"]], "foreign") +
                 zref(.data[["crowd"]],   "crowd"),
      di_quintile = as.integer(cut(.data[["di_score"]],
                                   breaks = c(-Inf, cuts, Inf),
                                   labels = FALSE))
    ) |>
    dplyr::select("area", population = "pop",
                  dplyr::all_of(ind_names), "di_score", "di_quintile")

  check_deprivation_plausibility(out)

  n_na <- sum(is.na(out[["di_score"]]))
  if (n_na) {
    warning(n_na, " area(s) have no deprivation score, because at least one ",
            "indicator had a zero denominator (typically an unpopulated NIL). ",
            "They will be dropped by the model.", call. = FALSE)
  }

  message("Deprivation built for ", nrow(out), " modelling areas; national ",
          "reference from ", nrow(nat), " municipalities.")

  attr(out, "national_reference") <- ref
  attr(out, "national_quintile_cuts") <- cuts
  attr(out, "vintage") <- vintage
  out
}


#' Compare the area-level and municipality-level deprivation indices
#'
#' Quantifies what the change of resolution bought. Reports the number of
#' distinct values inside Milan under each construction, and the standard
#' deviation of the score across NILs - which is exactly zero under the
#' municipal version, since every NIL inherited the same municipal figure.
#'
#' Intended for the methods appendix: it turns "we improved the resolution"
#' into a number.
#'
#' @param di_area Output of [build_deprivation()].
#' @param di_comune Municipality-level index expanded to areas, e.g. via
#'   [expand_cov_to_area()]. Optional.
#' @param milan_code Municipality code whose areas are the NILs.
#'
#' @return A one-row tibble of diagnostics.
#' @export
check_deprivation_resolution <- function(di_area, di_comune = NULL,
                                         milan_code = "015146") {
  in_milan <- grepl(paste0("^", milan_code, "_"), di_area[["area"]])

  tibble::tibble(
    n_areas            = nrow(di_area),
    n_nil              = sum(in_milan),
    nil_distinct_area  = dplyr::n_distinct(round(di_area[["di_score"]][in_milan], 6)),
    nil_sd_area        = stats::sd(di_area[["di_score"]][in_milan], na.rm = TRUE),
    nil_range_area     = diff(range(di_area[["di_score"]][in_milan], na.rm = TRUE)),
    nil_distinct_comune = if (is.null(di_comune)) NA_integer_ else {
      m <- grepl(paste0("^", milan_code, "_"), di_comune[["area"]])
      dplyr::n_distinct(round(di_comune[["di_score"]][m], 6))
    },
    overall_sd         = stats::sd(di_area[["di_score"]], na.rm = TRUE)
  )
}


#' Sanity-check the deprivation indicators against known national magnitudes
#'
#' Every indicator here is a proportion or a ratio whose plausible range is
#' known from published national statistics. Checking the magnitudes catches
#' the failure mode that no amount of code review will: a census variable code
#' that points at the wrong column. The arithmetic is then perfectly correct
#' and the result is silently inverted.
#'
#' This is not hypothetical. The four-indicator proxy in this package produced
#' a deprivation index that correlated at 0.08 with ISTAT's validated IVSM over
#' the same 279 areas, and gave a strongly protective association with
#' avoidable mortality, because \code{edu_low} - which drives roughly 70\% of
#' the summed score - was suspected of counting the wrong education
#' categories. Nothing errored. The only symptom was a coefficient with the
#' wrong sign, three stages downstream.
#'
#' @param di Output of [build_deprivation()].
#' @param bounds Named list of `c(low, high)` plausible ranges.
#' @param error Whether to stop rather than warn. Default `FALSE`, because a
#'   genuinely unusual territory should not block a run - but the warning must
#'   be read.
#'
#' @return `di`, invisibly.
#' @examples
#' \dontrun{
#' check_deprivation_plausibility(build_deprivation(census, xwalk))
#' }
#' @seealso [build_deprivation()]
#' @export
check_deprivation_plausibility <- function(
    di,
    bounds = list(
      # share aged 9+ without upper-secondary education. Around half the
      # Italian adult population; a median near 0.2 suggests the numerator is
      # counting graduates instead.
      edu_low = c(0.25, 0.75),
      # non-employment among the 15-64 population. Italian employment rate is
      # roughly 60-65%, so non-employment is roughly 0.35-0.40.
      nonemp  = c(0.20, 0.60),
      # foreign residents as a share of the population
      foreign = c(0.01, 0.40),
      # occupants per occupied dwelling
      crowd   = c(1.5, 3.5)
    ),
    error = FALSE) {

  problems <- character(0)

  for (v in names(bounds)) {
    if (!v %in% names(di)) next
    med <- stats::median(di[[v]], na.rm = TRUE)
    if (is.na(med)) next
    lo <- bounds[[v]][1]; hi <- bounds[[v]][2]
    if (med < lo || med > hi) {
      problems <- c(problems, sprintf(
        "  %s: median %.3f, outside the plausible range %.2f-%.2f",
        v, med, lo, hi))
    }
  }

  # Checking medians alone is not enough, and this is not a hypothetical
  # refinement: on the ATS Milano data every median was plausible while a
  # single `crowd` value of 951 was quietly turning the standardised score
  # into a dummy variable for one NIL. Screen the score itself for the shape
  # that does the damage - one unit far outside the rest.
  if ("di_score" %in% names(di)) {
    v   <- stats::na.omit(di[["di_score"]])
    med <- stats::median(v)
    mad <- stats::mad(v)
    if (mad > 0) {
      far <- abs(v - med) / mad
      if (max(far) > 20) {
        i <- which.max(abs(di[["di_score"]] - med))
        problems <- c(problems, sprintf(
          paste0("  di_score: %s is %.0f MADs from the median (%.2f vs %.2f). ",
                 "After standardisation this one area will absorb almost the ",
                 "entire scale and the covariate becomes an indicator for it."),
          if ("area" %in% names(di)) di[["area"]][i] else paste("row", i),
          max(far), di[["di_score"]][i], med))
      }
    }
  }

  if (length(problems)) {
    msg <- paste0(
      "Deprivation indicator(s) have implausible magnitudes:\n",
      paste(problems, collapse = "\n"), "\n",
      "  The arithmetic may be correct while the census variable codes point ",
      "at the wrong columns. Verify the codes in .deprivation_indicators() ",
      "against the ISTAT tracciato before using this index, and cross-check ",
      "the result against IVSM: two measures of area disadvantage over the ",
      "same areas should correlate strongly and positively."
    )
    if (error) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }

  invisible(di)
}
