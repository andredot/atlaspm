# nar_density_db.R
# Person-time general practitioner density from a NAR-like event history —
# database-friendly version, aligned with the SAS export contract.
#
# EXPORT CONTRACT (what the SAS extraction guarantees — see nar_export.sas):
#   For each person residing in ATS Milan, the event table contains
#     (a) AT LEAST ONE ANCHOR ROW: the person's most recent NAR event strictly
#         before `study_start`, so their status AT study start is known by
#         carrying that event forward (the NAR records changes, so the state
#         at any time is the last event before it);
#     (b) all status-change events between `study_start` and `study_end`
#         (zero or more per person);
#     (c) DEATH as an explicit event: death ends the "assisted by a GP"
#         status, and here it also ends the person's OBSERVATION — after
#         death a person contributes neither assigned nor unassigned
#         person-time. A resident-month counts in the denominator only while
#         the person is alive (on the first day of the month).
#
# STATED APPROXIMATIONS (all parameterised or documented, none silent):
#   - BOUNDARY CLOSURE: no GP practising outside ATS Milan is assumed to
#     serve residents inside it, and list sizes are computed from ATS Milan
#     residents only. For GPs near the boundary whose list includes external
#     patients this understates L (overstates supply) slightly; with choice-
#     based registration concentrated locally this is a small approximation,
#     but it is an approximation and is reported as such.
#   - PLAUSIBLE LIST SIZES: see clamp_supply() — artifact rule (< 2 patients)
#     and clamping of effective L to [l_min, l_max] = [500, 2000] by default,
#     so each assigned patient contributes between 1/2000 and 1/500 of a GP.
#     Every clamp is counted and returned in `clamp_report`.
#   - MONTHLY ATTRIBUTION: a month belongs to the state (and vital status)
#     holding on its first day. Sub-monthly transitions, including mid-month
#     deaths, are resolved to that day; the error is bounded by half a month
#     per transition.
#
# STATUS: the pipeline does NOT run this estimator. The area-level result is
# computed upstream against the register and read in by import_gp_density();
# what follows is the reference implementation of that computation, kept so the
# definition of the indicator is auditable in one place and so the estimator can
# be re-run here if the event-level export is ever loaded directly.
#
# ARCHITECTURE (each function = one candidate {targets} node):
#   build_spells()            -> assignment intervals from events
#   build_observation()       -> per-resident observation end (death/censor)
#   clamp_supply()            -> the visible list-size assumptions
#   compute_density_by_area() -> month-by-month density; DB-friendly, with
#                                progress reporting
# No residents x months grid is ever materialised: each month touches at most
# one row per living resident and collects only per-area aggregates.

# NOTE: this file previously attached dplyr, lubridate, purrr and tibble with
# library() at the top. That is a package-level side effect: it modifies the
# search path of whoever loads atlaspm, and R CMD check rejects it. Every call
# below is now namespaced, and the packages are declared in DESCRIPTION.

#' @importFrom dplyr filter mutate select group_by summarise arrange left_join
#' @importFrom dplyr n distinct bind_rows if_else coalesce
#' @importFrom lubridate floor_date year
#' @importFrom tibble tibble
#' @name nar-imports
#' @keywords internal
NULL

# ---------------------------------------------------------------------------

#' Reconstruct assignment spells from the NAR event history
#'
#' The NAR records changes, so a person's state at any time is their most
#' recent event before it: the anchor row exported before `study_start`
#' therefore fixes the initial status by simply carrying forward. Any event
#' — a new assignment, a revocation, or a death — terminates the open spell;
#' death terminates it definitively. Spells never extend past `study_end`.
#'
#' @param events Tibble: `resident_id`, `event_date`, `event_type`
#'   (`"assignment"`/`"revocation"`/`"death"`), `gp_id`.
#' @param study_end Date. Administrative censoring for still-open spells.
#' @return Tibble: `resident_id`, `gp_id`, `spell_start`, `spell_end`
#'   (inclusive).
#' @export
build_spells <- function(events, study_end) {
  events |>
    dplyr::arrange(resident_id, event_date) |>
    dplyr::group_by(resident_id) |>
    dplyr::mutate(spell_end = dplyr::coalesce(dplyr::lead(event_date) - lubridate::days(1), study_end)) |>
    dplyr::ungroup() |>
    dplyr::filter(event_type == "assignment") |>
    dplyr::transmute(resident_id, gp_id, spell_start = event_date,
              spell_end = pmin(spell_end, study_end)) |>
    dplyr::filter(spell_end >= spell_start)
}

# ---------------------------------------------------------------------------

#' Derive each resident's observation window from the event history
#'
#' Death ends observation: from the month after the month of death (deaths on
#' the 1st excluded from that month too, per the first-day attribution rule)
#' a resident contributes no person-time of any kind — they exit the
#' denominator cleanly rather than accumulating spurious "unassigned" time.
#' Residents without a death event are censored at `study_end`.
#'
#' Emigration out of ATS Milan, if exported as its own event type, slots in
#' here identically (treat it as an observation end); it is not simulated.
#'
#' @param residents Tibble: `resident_id`, `area_id`.
#' @param events Event tibble (only `death` rows are used).
#' @param study_end Date. Censoring date for survivors.
#' @return `residents` with `obs_end` (Date): last day of observation.
#' @export
build_observation <- function(residents, events, study_end) {
  deaths <- events |>
    dplyr::filter(event_type == "death") |>
    dplyr::group_by(resident_id) |>
    dplyr::summarise(death_date = min(event_date), .groups = "drop")

  residents |>
    dplyr::left_join(deaths, by = "resident_id") |>
    dplyr::mutate(obs_end = pmin(dplyr::coalesce(death_date, study_end), study_end)) |>
    dplyr::select(resident_id, area_id, obs_end)
}

# ---------------------------------------------------------------------------

#' Apply the plausible-list-size assumptions to one month of list sizes
#'
#' This function IS the assumptions, isolated so they are visible, testable,
#' and reported.
#'
#' @section Assumptions:
#' \describe{
#'   \item{artifact_min_patients (default 2)}{A GP observed with fewer
#'     patients than this in a month is an administrative artefact or a
#'     non-practising registration; the assignment is disregarded and their
#'     patients count as unassigned that month.}
#'   \item{l_min (default 500)}{Lists below 500 (trainees, part-time,
#'     ramp-up) are raised to 500: part of that GP's time is assumed not
#'     devoted to patients, so no patient contributes more than 1/500 of a
#'     GP-month.}
#'   \item{l_max (default 2000)}{Lists above 2000 are capped: marginal
#'     capacity is assumed exhausted, so no assigned patient contributes
#'     less than 1/2000 of a GP-month.}
#'   \item{Boundary closure}{L is computed from ATS Milan residents only and
#'     external GPs are absent from the export; both directions of
#'     cross-boundary registration are assumed negligible.}
#' }
#'
#' @param list_sizes Tibble: `gp_id`, `list_size` (one month).
#' @param l_min,l_max Numeric clamping bounds for effective list size.
#' @param artifact_min_patients Integer artifact threshold (strictly below).
#' @return Input plus `is_artifact`, `l_clamped`, `clamp_status`
#'   ("artifact"/"raised_to_min"/"capped_at_max"/"kept"), `supply`
#'   (= 1/l_clamped, NA for artifacts).
#' @export
clamp_supply <- function(list_sizes,
                         l_min = 500,
                         l_max = 2000,
                         artifact_min_patients = 2L) {
  stopifnot(l_min <= l_max, artifact_min_patients >= 0)
  list_sizes |>
    dplyr::mutate(
      is_artifact = list_size < artifact_min_patients,
      l_clamped   = pmin(pmax(list_size, l_min), l_max),
      clamp_status = dplyr::case_when(
        is_artifact       ~ "artifact",
        list_size < l_min ~ "raised_to_min",
        list_size > l_max ~ "capped_at_max",
        TRUE              ~ "kept"
      ),
      supply = dplyr::if_else(is_artifact, NA_real_, 1 / l_clamped)
    )
}

# ---------------------------------------------------------------------------

#' Person-time GP density by area, computed month by month
#'
#' For each month m (attributed by its first day):
#'   1. denominator: residents under observation at m (`obs_end >= m`) —
#'      the dead have exited cleanly;
#'   2. numerator inputs: spells active at m (<= 1 row per living resident);
#'   3. observed monthly list sizes -> [clamp_supply()] (visible assumptions);
#'   4. per-area aggregates: assigned/artifact person-months, summed supply,
#'      list-size components; unassigned = alive - assigned.
#'
#' `spells` and `residents` may be tibbles or dbplyr lazy tables (DuckDB,
#' Postgres, SQLite): every verb inside the loop translates to SQL, so only
#' small aggregates cross into R each month. NOTE on database date columns:
#' store dates as native DATE (DuckDB, Postgres) or ISO-8601 text (SQLite,
#' which has no DATE type); integer epoch-days will compare wrongly against
#' the Date literals used here. Verified to give identical results on
#' tibbles and on a SQL backend.
#'
#' Progress: a cli progress bar (message() fallback) advances per month; a
#' closing summary reports the person-time affected by each assumption.
#'
#' @param spells Tibble/lazy table: `resident_id`, `gp_id`, `spell_start`,
#'   `spell_end`.
#' @param residents Tibble/lazy table: `resident_id`, `area_id`, `obs_end`
#'   (from [build_observation()]).
#' @param study_start,study_end Dates bounding the monthly grid.
#' @param l_min,l_max,artifact_min_patients Passed to [clamp_supply()].
#' @param gps Optional tibble `gp_id`, `birth_year` for the 65+ measure.
#' @param progress Logical.
#'
#' @return list(density, by_month, clamp_report, assumptions):
#'   \describe{
#'     \item{density}{Per `area_id`: `person_months` (alive person-time),
#'       `gp_density`, `gp_density_per_100k`, `mean_list_size_assigned`
#'       (observed, unclamped), `prop_unassigned` (share of ALIVE person-time
#'       unassigned; artifact months included), `prop_artifact`, `gps_fte`,
#'       optional `prop_gp_65plus`.}
#'     \item{by_month}{Same per `area_id` x `month`, incl. `n_alive` — use
#'       it to inspect retirement-churn windows and mortality-driven
#'       denominator decline.}
#'     \item{clamp_report}{Per month x `clamp_status`: GPs and person-months
#'       affected — the audit trail of the assumptions.}
#'     \item{assumptions}{One-row echo of all thresholds and the implied
#'       supply bounds, plus the boundary-closure flag.}
#'   }
#' @export
compute_density_by_area <- function(spells,
                                    residents,
                                    study_start,
                                    study_end,
                                    l_min = 500,
                                    l_max = 2000,
                                    artifact_min_patients = 2L,
                                    gps = NULL,
                                    progress = TRUE) {
  months_grid <- seq(lubridate::floor_date(study_start, "month"),
                     lubridate::floor_date(study_end, "month"), by = "month")

  has_cli <- requireNamespace("cli", quietly = TRUE)
  if (progress && has_cli) {
    cli::cli_alert_info(paste0(
      "Assumptions: L clamped to [{l_min}, {l_max}] ",
      "(supply in [1/{l_max}, 1/{l_min}]); GPs with < {artifact_min_patients} ",
      "patients = artifacts; boundary closure (no cross-ATS registration); ",
      "death exits the denominator."))
    bar <- cli::cli_progress_bar(
      "Computing monthly density",
      total = length(months_grid),
      format = "{cli::pb_bar} {cli::pb_percent} | month {cli::pb_current}/{cli::pb_total} ({format(months_grid[cli::pb_current], '%Y-%m')}) | ETA {cli::pb_eta}",
      .envir = environment())
  }

  monthly <- vector("list", length(months_grid))
  clamps  <- vector("list", length(months_grid))

  for (i in seq_along(months_grid)) {
    m <- months_grid[[i]]

    # 1. Denominator: residents alive (under observation) on day 1 of m.
    alive_m <- residents |>
      dplyr::filter(obs_end >= m) |>
      dplyr::count(area_id, name = "n_alive") |>
      dplyr::collect()

    # 2. Assignment state on day 1 of m. Spells cannot outlive obs_end by
    # construction (death is an event), so no extra filter is needed here.
    active <- spells |>
      dplyr::filter(spell_start <= m, spell_end >= m)

    # 3. Observed monthly list sizes -> visible assumptions.
    lists_m <- active |>
      dplyr::count(gp_id, name = "list_size") |>
      dplyr::collect() |>
      clamp_supply(l_min = l_min, l_max = l_max,
                   artifact_min_patients = artifact_min_patients)

    if (!is.null(gps)) {
      lists_m <- lists_m |>
        dplyr::left_join(gps, by = "gp_id") |>
        dplyr::mutate(gp_65plus = (lubridate::year(m) - birth_year) >= 65)
    } else {
      lists_m <- dplyr::mutate(lists_m, gp_65plus = NA)
    }

    clamps[[i]] <- lists_m |>
      dplyr::group_by(clamp_status) |>
      dplyr::summarise(n_gps = dplyr::n(), person_months = sum(list_size),
                .groups = "drop") |>
      dplyr::mutate(month = m)

    # 4. Per-area aggregation of assigned person-time.
    assigned_m <- active |>
      dplyr::inner_join(residents, by = "resident_id") |>
      dplyr::count(area_id, gp_id, name = "n_patients") |>
      dplyr::collect() |>
      dplyr::left_join(lists_m, by = "gp_id") |>
      dplyr::group_by(area_id) |>
      dplyr::summarise(
        assigned_pm = sum(n_patients * !is_artifact),
        artifact_pm = sum(n_patients * is_artifact),
        supply_sum  = sum(n_patients * dplyr::coalesce(supply, 0)),
        list_pm_sum = sum(n_patients * list_size * !is_artifact),
        gp65_pm     = sum(n_patients * dplyr::coalesce(as.numeric(gp_65plus), 0) *
                          !is_artifact),
        .groups = "drop")

    monthly[[i]] <- alive_m |>
      dplyr::left_join(assigned_m, by = "area_id") |>
      dplyr::mutate(dplyr::across(-c(area_id, n_alive), ~ dplyr::coalesce(.x, 0)),
             month = m,
             unassigned_pm = n_alive - assigned_pm)

    if (progress && has_cli) cli::cli_progress_update(id = bar)
    else if (progress) message(sprintf("[%d/%d] %s done", i,
                                       length(months_grid), format(m, "%Y-%m")))
  }
  if (progress && has_cli) cli::cli_progress_done(id = bar)

  by_month <- purrr::list_rbind(monthly) |>
    dplyr::mutate(
      gp_density              = supply_sum / n_alive,
      prop_unassigned         = unassigned_pm / n_alive,
      mean_list_size_assigned = dplyr::if_else(assigned_pm > 0,
                                        list_pm_sum / assigned_pm, NA_real_))

  density <- by_month |>
    dplyr::group_by(area_id) |>
    dplyr::summarise(
      person_months           = sum(n_alive),
      gp_density              = sum(supply_sum) / sum(n_alive),
      gp_density_per_100k     = gp_density * 1e5,
      mean_list_size_assigned = sum(list_pm_sum) / sum(assigned_pm),
      prop_unassigned         = sum(unassigned_pm) / sum(n_alive),
      prop_artifact           = sum(artifact_pm) / sum(n_alive),
      gps_fte                 = sum(supply_sum) / 12,
      prop_gp_65plus          = if (!is.null(gps))
                                  sum(gp65_pm) / sum(assigned_pm)
                                else NA_real_,
      .groups = "drop")

  clamp_report <- purrr::list_rbind(clamps) |>
    dplyr::select(month, clamp_status, n_gps, person_months)

  if (progress && has_cli) {
    tot <- clamp_report |>
      dplyr::group_by(clamp_status) |>
      dplyr::summarise(pm = sum(person_months), .groups = "drop")
    denom <- sum(tot$pm)
    for (s in c("raised_to_min", "capped_at_max", "artifact")) {
      pm <- tot$pm[tot$clamp_status == s]
      if (length(pm))
        cli::cli_alert_warning(
          "Assumption '{s}': {pm} assigned person-months affected ({round(100*pm/denom, 1)}%).")
    }
  }

  list(
    density  = density,
    by_month = by_month |>
      dplyr::select(area_id, month, n_alive, gp_density, prop_unassigned,
             mean_list_size_assigned, assigned_pm, artifact_pm),
    clamp_report = clamp_report,
    assumptions  = tibble::tibble(
      l_min = l_min, l_max = l_max,
      artifact_min_patients = artifact_min_patients,
      supply_lower = 1 / l_max, supply_upper = 1 / l_min,
      boundary_closure = TRUE,
      death_exits_denominator = TRUE)
  )
}

# ---------------------------------------------------------------------------
# {targets} sketch, DuckDB-backed, on the SAS export:
#
# # _targets.R
# library(targets)
# tar_source("R/nar_density_db.R")
# study_start <- as.Date("2022-01-01"); study_end <- as.Date("2024-12-31")
# list(
#   tar_target(export_file, "data/nar_export.parquet", format = "file"),
#   tar_target(events,    arrow::read_parquet(export_file)),   # SAS output
#   tar_target(res_file,  "data/residents.parquet", format = "file"),
#   tar_target(residents0, arrow::read_parquet(res_file)),
#   tar_target(spells,    build_spells(events, study_end)),
#   tar_target(residents, build_observation(residents0, events, study_end)),
#   tar_target(db_path, {
#     con <- DBI::dbConnect(duckdb::duckdb(), "nar.duckdb")
#     on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
#     DBI::dbWriteTable(con, "spells",    spells,    overwrite = TRUE)
#     DBI::dbWriteTable(con, "residents", residents, overwrite = TRUE)
#     "nar.duckdb"
#   }, format = "file"),
#   tar_target(density, {
#     con <- DBI::dbConnect(duckdb::duckdb(), db_path, read_only = TRUE)
#     on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
#     compute_density_by_area(tbl(con, "spells"),
#                             tbl(con, "residents"),
#                             study_start, study_end,
#                             l_min = 500, l_max = 2000,
#                             artifact_min_patients = 2L, gps = gps_table)
#   })
# )
# Sensitivity branches: tarchetypes::tar_map() over (l_min, l_max), with
# (1, Inf) recovering the unclamped estimator.
