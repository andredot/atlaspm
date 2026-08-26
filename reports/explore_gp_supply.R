# ---------------------------------------------------------------------------
# Exploratory look at the primary-care supply indicator (2023)
#
# The indicator is a pipeline target now, not a loose CSV: import_gp_density()
# reads and validates it, check_gp_density() audits it against the modelling
# geography, and add_geo() attaches the polygons. This script only *looks* at
# those targets - it deliberately re-derives nothing, so the numbers here and
# the numbers in the thesis cannot drift apart.
#
#   gp_density_area  area, gp_density (burden-adjusted GP-equivalents per
#                    1,000 patients), gp_caseload, gp_supply, gp_population,
#                    gp_density_unweighted, gp_supply_unweighted,
#                    gp_weight_ratio
#   gp_density_geo   the same, joined to area_shp (keep = "geometry", so every
#                    modelled area is present and unmatched ones draw grey)
#   gp_density_audit the coverage audit
#
# TWO INDICATORS, NOT ONE
#
# The 2023 extract weights each patient by the national age-sex care-burden
# schedule, so an area whose patients are older shows lower effective supply
# for the same headcount of doctors. `gp_density` is that weighted measure and
# is the covariate the models use. `gp_density_unweighted` is the previous
# definition - a plain count of GP-equivalents per patient - carried forward so
# the effect of the re-weighting can be measured rather than asserted.
#
# The two are not a rescaling of each other. Sections 3, 5 and 6 quantify the
# difference and map where it falls, which is the evidence needed before
# writing that the weighting mattered.
#
# Output: nine figures in output/figures/ + the audit and comparison, printed.
#
# Run with the package loaded (devtools::load_all()) or after library(atlaspm),
# and after tar_make() has built the primary-care targets.
# ---------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sf)
targets::tar_source()

FIG_DIR <- "output/figures"

targets::tar_load(c(gp_density_area, gp_density_geo, gp_density_audit))

has_unweighted <- "gp_density_unweighted" %in% names(gp_density_area)
if (!has_unweighted) {
  message("`gp_density_unweighted` is absent: the indicator file predates the ",
          "care-burden weighting. The comparison sections will be skipped.")
}

# --- 1. Coverage -------------------------------------------------------------

cat("\n-- coverage against the modelling geography ------------------------\n")
print(gp_density_audit)
cat("\nModelled areas with no indicator: ",
    paste(c(attr(gp_density_audit, "area_no_indicator"), "(none)"),
          collapse = ", "), "\n",
    "Indicator rows with no modelled area: ",
    paste(c(attr(gp_density_audit, "indicator_no_area"), "(none)"),
          collapse = ", "), "\n", sep = "")

# --- 2. Distribution, by geography type --------------------------------------
# The study area mixes three very different unit types, and pooling them hides
# that: a Lodi comune of 1,500 people and a Milan NIL of 60,000 are one row
# each. Splitting them is the first thing worth looking at.

supply <- gp_density_area |>
  mutate(
    stratum = case_when(
      grepl("_", area)            ~ "Milan (NIL)",
      substr(area, 1, 3) == "015" ~ "Milan province (comuni)",
      substr(area, 1, 3) == "098" ~ "Lodi province (comuni)",
      TRUE                        ~ "other"
    )
  )

STRATUM_PAL <- c("Milan (NIL)"             = "#b5651d",
                 "Milan province (comuni)" = "#5a8a7d",
                 "Lodi province (comuni)"  = "#dca678",
                 "other"                   = "grey50")

wmean_density <- weighted.mean(supply$gp_density, supply$gp_population)
wmean_unw <- if (has_unweighted) {
  weighted.mean(supply$gp_density_unweighted, supply$gp_population)
} else NA_real_

cat("\n-- burden-adjusted indicator, by stratum --------------------------\n")
supply |>
  group_by(stratum) |>
  summarise(
    n            = n(),
    pop          = sum(gp_population),
    density_med  = median(gp_density),
    caseload_p10 = quantile(gp_caseload, 0.10),
    caseload_med = median(gp_caseload),
    caseload_p90 = quantile(gp_caseload, 0.90),
    .groups = "drop"
  ) |>
  print()

cat(sprintf(
  "\nPopulation-weighted mean: %.3f burden-adjusted GP-equivalents per 1,000 (%.0f patients per GP-equivalent)\n",
  wmean_density, 1000 / wmean_density))

# --- 3. What did the care-burden weighting do? -------------------------------
# The question is not "are the two indicators different" - a different
# numerator guarantees that - but "do they identify the SAME areas as
# underserved". If they rank areas alike, the weighting is a rescaling and
# nothing downstream changes. If they do not, the choice of indicator is
# substantive and has to be defended in the methods.

if (has_unweighted) {

  cmp <- supply |>
    mutate(
      diff      = gp_density - gp_density_unweighted,
      pct_diff  = 100 * diff / gp_density_unweighted,
      rank_w    = rank(gp_density),
      rank_u    = rank(gp_density_unweighted),
      rank_move = rank_u - rank_w   # positive = looks WORSE once weighted
    )

  rho_s <- cor(cmp$gp_density, cmp$gp_density_unweighted, method = "spearman")
  rho_p <- cor(cmp$gp_density, cmp$gp_density_unweighted)

  cat("\n-- burden-adjusted vs unweighted ----------------------------------\n")
  cat(sprintf("Spearman rank correlation : %.4f\n", rho_s))
  cat(sprintf("Pearson correlation       : %.4f\n", rho_p))
  cat(sprintf("Weight ratio (adjusted supply / unweighted supply): %.3f to %.3f, median %.3f\n",
              min(cmp$gp_weight_ratio, na.rm = TRUE),
              max(cmp$gp_weight_ratio, na.rm = TRUE),
              median(cmp$gp_weight_ratio, na.rm = TRUE)))
  cat(sprintf("Population-weighted mean  : %.3f adjusted vs %.3f unweighted\n",
              wmean_density, wmean_unw))
  cat(sprintf("Areas moving >25 rank positions: %d of %d\n",
              sum(abs(cmp$rank_move) > 25), nrow(cmp)))

  cat("\n-- looks WORSE once care burden is counted ------------------------\n")
  cmp |>
    arrange(desc(rank_move)) |>
    select(area, stratum, gp_density_unweighted, gp_density,
           gp_weight_ratio, pct_diff, rank_move) |>
    head(10) |>
    print()

  cat("\n-- looks BETTER once care burden is counted -----------------------\n")
  cmp |>
    arrange(rank_move) |>
    select(area, stratum, gp_density_unweighted, gp_density,
           gp_weight_ratio, pct_diff, rank_move) |>
    head(10) |>
    print()

  cat("\n-- weight ratio by stratum ----------------------------------------\n")
  cmp |>
    group_by(stratum) |>
    summarise(n = n(),
              ratio_med = median(gp_weight_ratio, na.rm = TRUE),
              ratio_p10 = quantile(gp_weight_ratio, 0.10, na.rm = TRUE),
              ratio_p90 = quantile(gp_weight_ratio, 0.90, na.rm = TRUE),
              .groups = "drop") |>
    print()
}

# --- 4. Histograms -----------------------------------------------------------
# One bar per AREA, not per resident, so every unit carries equal weight
# regardless of size. The dashed line is the population-weighted mean, which is
# the honest central value and sits somewhere else entirely.

p_hist <- ggplot(supply, aes(x = gp_density)) +
  geom_histogram(binwidth = 0.015, boundary = 0,
                 fill = "#5a8a7d", colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = wmean_density,
             linetype = "22", colour = "#8c2d04", linewidth = 0.6) +
  labs(
    title    = "Primary-care supply across the modelling areas, 2023",
    subtitle = sprintf("Burden-adjusted GP-equivalents per 1,000 patients; dashed line = population-weighted mean (%.3f)",
                       wmean_density),
    x = "Burden-adjusted GP-equivalents per 1,000 patients",
    y = "Number of areas",
    caption = "One bar per area (comune, or NIL within Milan), unweighted by population."
  ) +
  theme_atlas()

p_hist_strata <- p_hist +
  facet_wrap(~ stratum, ncol = 1, scales = "free_y") +
  labs(subtitle = "By geography type; note the different y-scales")

# Both indicators overlaid: same bins, same axis, so a shift in location and a
# change in spread are separable by eye.

if (has_unweighted) {
  long <- bind_rows(
    supply |> transmute(area, stratum, value = gp_density,
                        indicator = "Burden-adjusted (current)"),
    supply |> transmute(area, stratum, value = gp_density_unweighted,
                        indicator = "Unweighted count (previous)")
  )

  p_hist_both <- ggplot(long, aes(x = value, fill = indicator)) +
    geom_histogram(binwidth = 0.015, boundary = 0, position = "identity",
                   alpha = 0.55, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("Burden-adjusted (current)"   = "#b5651d",
                                 "Unweighted count (previous)" = "#5a8a7d"),
                      name = NULL) +
    labs(
      title    = "What the care-burden weighting does to the distribution",
      subtitle = sprintf("Population-weighted means: %.3f adjusted, %.3f unweighted",
                         wmean_density, wmean_unw),
      x = "GP-equivalents per 1,000 patients", y = "Number of areas",
      caption = "Same bins and axis for both. A shift in location alone would be a rescaling; a change in spread or in ordering is not."
    ) +
    theme_atlas()
}

# --- 5. Do the two agree about which areas are underserved? ------------------
# The diagonal is the answer. Points on it mean the weighting changed the level
# but not the ordering; scatter around it means areas trade places, and the
# choice of indicator decides who looks underserved.

if (has_unweighted) {
  p_scatter <- ggplot(cmp, aes(x = gp_density_unweighted, y = gp_density)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey45") +
    geom_point(aes(colour = stratum, size = gp_population), alpha = 0.6) +
    scale_size_area(name = "Patients", max_size = 5) +
    scale_colour_manual(values = STRATUM_PAL, name = NULL) +
    coord_equal() +
    labs(
      title    = "Burden-adjusted against unweighted supply, by area",
      subtitle = sprintf("Spearman rank correlation %.3f. Points off the diagonal change rank once burden is counted.",
                         rho_s),
      x = "Unweighted GP-equivalents per 1,000 patients",
      y = "Burden-adjusted GP-equivalents per 1,000 patients"
    ) +
    theme_atlas()

  # Rank movement against area size. A weighting effect concentrated in small
  # areas is a different problem from one that moves the large ones.
  p_rankmove <- ggplot(cmp, aes(x = gp_population, y = rank_move)) +
    geom_hline(yintercept = 0, linetype = "22", colour = "grey45") +
    geom_point(aes(colour = stratum), alpha = 0.6, size = 1.9) +
    scale_x_log10(labels = scales::comma) +
    scale_colour_manual(values = STRATUM_PAL, name = NULL) +
    labs(
      title    = "Which areas the weighting reclassifies",
      subtitle = "Positive = the area looks worse served once care burden is counted",
      x = "Patients (log scale)",
      y = "Change in rank (unweighted \u2212 adjusted)"
    ) +
    theme_atlas()
}

# --- 6. Maps -----------------------------------------------------------------
# 6a. Continuous caseload. This is the same figure the pipeline builds as
#     fig_gp_density_map; it is rebuilt here so the script stands alone.

p_map_caseload <- plot_travel_time(
  gp_density_geo,
  value        = "gp_caseload",
  legend_title = "Patients per\nGP-equivalent",
  caption      = "Reciprocal of the burden-adjusted indicator; darker = each GP-equivalent carries more weighted caseload."
) +
  labs(title    = "Primary-care caseload, 2023",
       subtitle = "Burden-adjusted patients per GP-equivalent, by area")

# 6b. Binned and anchored on the ATS average, reusing the SMR map. Oriented as
#     a STRAIN ratio (local caseload / average caseload) rather than as supply,
#     because the palette encodes "high = worse than expected" and inverting
#     the measure is what keeps the colours honest.

strain_geo <- gp_density_geo |>
  mutate(gp_strain = wmean_density / gp_density)

p_map_strain <- plot_smr_map(
  strain_geo,
  value    = "gp_strain",
  title    = "Primary-care strain relative to the ATS average, 2023",
  subtitle = "Caseload ratio: 1 = the population-weighted average caseload",
  caption  = "Above 1 = each GP-equivalent carries more weighted caseload than the ATS average."
)
p_map_strain$labels$fill <- "Caseload ratio\n(local / average)"

# 6c. WHERE the weighting bites. This maps gp_weight_ratio rather than the
#     difference in densities, because the ratio isolates the age-sex
#     composition effect on its own, independent of how many doctors an area
#     happens to have. Above 1 means the area's patients are heavier than the
#     national schedule's average, so the same headcount represents less
#     effective supply. Reused from the SMR map because the question has the
#     same shape - above or below a meaningful 1 - so the palette keeps its
#     meaning. Breaks are tighter than the SMR default because the ratio spans
#     a narrower range.

if (has_unweighted) {
  ratio_geo <- gp_density_geo |>
    left_join(select(gp_density_area, area, gp_weight_ratio), by = "area")

  p_map_ratio <- plot_smr_map(
    ratio_geo,
    value    = "gp_weight_ratio.x",
    breaks   = c(-Inf, 0.9, 0.95, 1.05, 1.1, Inf),
    title    = "Where the care-burden weighting changes the picture, 2023",
    subtitle = "Ratio of burden-adjusted to unweighted supply",
    caption  = "Above 1 = this area's patients are heavier than the national age-sex average, so the same number of GPs represents less effective supply."
  )
  p_map_ratio$labels$fill <- "Adjusted /\nunweighted"

  # 6d. Rank movement: the decision-relevant version. Not how much the number
  #     moved, but whether the area changed place in the queue.
  move_geo <- gp_density_geo |>
    left_join(select(cmp, area, rank_move), by = "area")

  p_map_move <- ggplot(sf::st_as_sf(move_geo)) +
    geom_sf(aes(fill = rank_move), colour = NA) +
    scale_fill_gradient2(
      low = "#5a8a7d", mid = "#f2e8d5", high = "#b5651d", midpoint = 0,
      na.value = "grey88", name = "Rank change"
    ) +
    labs(
      title    = "Areas reclassified by the care-burden weighting",
      subtitle = "Terracotta = looks worse served once burden is counted; sage = looks better",
      caption  = "Change in rank position between the unweighted and burden-adjusted indicators."
    ) +
    theme_atlas(map = TRUE)
}

# --- 7. Save -----------------------------------------------------------------
save_figure(p_hist,         "gp_supply_hist",        dir = FIG_DIR, height = 110)
save_figure(p_hist_strata,  "gp_supply_hist_strata", dir = FIG_DIR, height = 180)
save_figure(p_map_caseload, "gp_caseload_map",       dir = FIG_DIR)
save_figure(p_map_strain,   "gp_strain_map",         dir = FIG_DIR)

if (has_unweighted) {
  save_figure(p_hist_both, "gp_supply_hist_compare", dir = FIG_DIR, height = 110)
  save_figure(p_scatter,   "gp_supply_scatter",      dir = FIG_DIR, height = 170)
  save_figure(p_rankmove,  "gp_supply_rank_move",    dir = FIG_DIR, height = 120)
  save_figure(p_map_ratio, "gp_weight_ratio_map",    dir = FIG_DIR)
  save_figure(p_map_move,  "gp_rank_move_map",       dir = FIG_DIR)
}

# --- 8. Note on the covariate ------------------------------------------------
# smr_geo_full carries `gp_density` (and `gp_density_z`) via add_covariate(),
# and that is the BURDEN-ADJUSTED measure. `gp_density_unweighted` travels
# alongside it but is not a second covariate and is not modelled; it exists so
# the weighting can be audited. `gp_strain`, `gp_weight_ratio` and `rank_move`
# above are presentational transforms for these maps only.
#
# If section 3 shows the two indicators ranking areas differently, that belongs
# in the methods as a stated choice with its consequence quantified, not as an
# unremarked switch between extracts.
#
# Note also that this file is the density component alone: the proposal's "GP
# Strain Index" also folds in the share of residents without a registered GP
# and the share of GPs over 62, which the upstream extract can supply as
# prop_unassigned / prop_gp_65plus (see compute_density_by_area()).
