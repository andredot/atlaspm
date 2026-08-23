# 2011 census and the deprivation stability check -----------------------------
#
# Appendix C. Answers the assumption the Discussion rests on: that an area's
# relative deprivation is stable enough over a decade for a 2023 measurement to
# stand in for the accumulated exposure that actually causes premature death.
#
# WHAT THIS IS AND IS NOT
#
# This is NOT the validated Italian Deprivation Index for 2011. That index
# (Rosano et al. 2020, Epidemiol Prev 44(2-3):162-170) was built from
# INDIVIDUAL census microdata, which is what allowed its authors to restrict
# low education to ages 15-60 and single-parent status to households with
# minor children. Neither restriction is reconstructable from the public
# section-level release: it has no age-by-education crosstab and no family
# typology at all (PF1-PF9 are household sizes only).
#
# What this file computes is the SAME four-indicator proxy used for 2023,
# applied to 2011 data. That is the right instrument for the question being
# asked. A stability test needs the method held constant so that any change is
# attributable to the territory; comparing a five-indicator microdata index
# against a four-indicator aggregate proxy would confound real change with the
# difference between two constructions, and a low correlation would not say
# which.
#
# The cost is that this tells you nothing about how good the proxy is. That is
# a separate question, needing a separate comparison, and it is not answered
# here.


#' Import the 2011 general census at section level
#'
#' Reads the ISTAT regional section files and keeps the variables needed to
#' reproduce the four-indicator proxy. Same shape as
#' [import_census_2023()], so the same downstream code consumes both.
#'
#' @param path Directory of regional CSVs, or a vector of file paths.
#' @param pattern Regex selecting the files when `path` is a directory.
#' @param delim Field delimiter. ISTAT ships these semicolon-separated.
#' @param cols Variables to retain. The defaults are the 2011 codes
#'   corresponding to the 2023 proxy; see [.deprivation_indicators()] for the
#'   mapping.
#'
#' @return A tibble, one row per census section: `SEZ2011`, `PROCOM` and the
#'   selected counts as numeric.
#'
#' @examples
#' \dontrun{
#' census_2011 <- import_census_2011("data-raw/geodata/census_2011")
#' }
#' @seealso [build_deprivation()], [check_deprivation_stability()]
#' @export
import_census_2011 <- function(path,
                               pattern = "\\.csv$",
                               delim   = ";",
                               cols    = c("PROCOM", "SEZ2011",
                                           "P1", "P46",
                                           "P49", "P50", "P51", "P52",
                                           "P60", "P61",
                                           paste0("P", 17:26),
                                           "ST1", "A2")) {

  files <- if (length(path) == 1L && dir.exists(path)) {
    list.files(path, pattern = pattern, full.names = TRUE)
  } else {
    path
  }
  if (!length(files)) {
    stop("No 2011 census files found at: ", paste(path, collapse = ", "),
         call. = FALSE)
  }

  pad <- function(x) sprintf("%06d", as.integer(x))

  read_one <- function(f) {
    d <- readr::read_delim(f, delim = delim, show_col_types = FALSE,
                           col_types = readr::cols(.default = readr::col_character()))
    missing <- setdiff(cols, names(d))
    if (length(missing)) {
      stop("`", basename(f), "` is missing: ", paste(missing, collapse = ", "),
           ".\n  Available: ", paste(utils::head(names(d), 40), collapse = ", "),
           call. = FALSE)
    }
    d <- dplyr::select(d, dplyr::all_of(cols))
    count_cols <- setdiff(cols, c("PROCOM", "SEZ2011"))

    dplyr::mutate(
      d,
      PROCOM  = pad(.data[["PROCOM"]]),
      SEZ2011 = as.character(.data[["SEZ2011"]]),
      dplyr::across(dplyr::all_of(count_cols),
                    ~ suppressWarnings(as.numeric(.x)))
    )
  }

  out <- dplyr::bind_rows(lapply(files, read_one))

  message("2011 census: ", nrow(out), " sections from ", length(files),
          " file(s), ", dplyr::n_distinct(out[["PROCOM"]]), " municipalities.")
  out
}


#' Is area deprivation stable enough to treat 2023 as a decade-long exposure?
#'
#' The Discussion argues that deprivation acts over decades and that a
#' contemporary measurement can stand in for accumulated exposure. That is only
#' true if areas hold their relative position. This quantifies it.
#'
#' Reports the rank correlation of the overall score, the same for each
#' component indicator, and the agreement between quintile assignments. The
#' per-indicator breakdown matters more than the headline: a stable index built
#' from unstable components is a different situation from a stable index built
#' from stable ones, and only the second supports the argument being made.
#'
#' Quintile agreement is reported as the share of areas staying in the same
#' quintile and the share moving by two or more. Areas moving two quintiles have
#' materially changed their position in the national distribution, and a
#' meaningful count of them undercuts the exposure interpretation regardless of
#' what the overall correlation says.
#'
#' @param di_old,di_new Outputs of [build_deprivation()] for the two vintages.
#' @param by_layer Report comuni and Milan NILs separately as well as pooled.
#'   Default `TRUE`; the two layers can behave very differently, and pooling
#'   hides it.
#' @param milan_code Municipality whose areas are NILs.
#'
#' @return A list with `overall`, `by_indicator` and `quintiles`.
#'
#' @examples
#' \dontrun{
#' check_deprivation_stability(deprivation_2011, deprivation_area)
#' }
#' @seealso [build_deprivation()], [check_aq_ranks()]
#' @export
check_deprivation_stability <- function(di_old, di_new,
                                        by_layer   = TRUE,
                                        milan_code = "015146") {

  ind <- c("edu_low", "nonemp", "foreign", "crowd")

  j <- dplyr::inner_join(
    dplyr::select(di_old, "area", "di_score", "di_quintile",
                  dplyr::all_of(ind)),
    dplyr::select(di_new, "area", "di_score", "di_quintile",
                  dplyr::all_of(ind)),
    by = "area", suffix = c("_old", "_new")
  )

  if (!nrow(j)) {
    stop("No areas in common between the two vintages. The 2011 and 2021 ",
         "section geographies differ, so the crosswalks must both resolve to ",
         "the SAME `area` keys - check that build_section_area_xwalk() was ",
         "given the matching shapefile for each year.", call. = FALSE)
  }

  j[["layer"]] <- ifelse(grepl(paste0("^", milan_code, "_"), j[["area"]]),
                         "Milan NILs", "Comuni")

  rho <- function(d, a, b) {
    if (nrow(d) < 3L) return(NA_real_)
    stats::cor(d[[a]], d[[b]], method = "spearman",
               use = "complete.obs")
  }

  groups <- list(All = j)
  if (by_layer) {
    for (l in unique(j[["layer"]])) groups[[l]] <- j[j[["layer"]] == l, ]
  }

  overall <- do.call(rbind, lapply(names(groups), function(g) {
    d <- groups[[g]]
    tibble::tibble(
      layer   = g,
      n_areas = nrow(d),
      rho     = rho(d, "di_score_old", "di_score_new"),
      pearson = if (nrow(d) < 3L) NA_real_ else
        stats::cor(d[["di_score_old"]], d[["di_score_new"]],
                   use = "complete.obs")
    )
  }))

  by_indicator <- do.call(rbind, lapply(ind, function(v) {
    tibble::tibble(
      indicator = v,
      rho       = rho(j, paste0(v, "_old"), paste0(v, "_new")),
      mean_2011 = mean(j[[paste0(v, "_old")]], na.rm = TRUE),
      mean_2023 = mean(j[[paste0(v, "_new")]], na.rm = TRUE)
    )
  }))

  shift <- abs(j[["di_quintile_new"]] - j[["di_quintile_old"]])
  quintiles <- do.call(rbind, lapply(names(groups), function(g) {
    d  <- groups[[g]]
    sh <- abs(d[["di_quintile_new"]] - d[["di_quintile_old"]])
    sh <- sh[!is.na(sh)]
    tibble::tibble(
      layer          = g,
      n_areas        = length(sh),
      pct_same       = 100 * mean(sh == 0),
      pct_moved_1    = 100 * mean(sh == 1),
      pct_moved_2plus = 100 * mean(sh >= 2)
    )
  }))

  list(overall = overall, by_indicator = by_indicator, quintiles = quintiles,
       joined = j)
}
