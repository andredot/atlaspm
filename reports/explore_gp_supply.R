# ---------------------------------------------------------------------------
# Exploratory look at the primary-care supply indicator (2023)
#
# The indicator is a pipeline target now, not a loose CSV: import_gp_density()
# reads and validates it, check_gp_density() audits it against the modelling
# geography, and add_geo() attaches the polygons. This script only *looks* at
# those targets - it deliberately re-derives nothing, so the numbers here and
# the numbers in the thesis cannot drift apart.
#
#   gp_density_area  area, gp_density (GP-equivalents per 1,000 residents),
#                    gp_caseload (residents per GP-equivalent), gp_supply,
#                    gp_population
#   gp_density_geo   the same, joined to area_shp (keep = "geometry", so every
#                    modelled area is present and unmatched ones draw grey)
#   gp_density_audit the coverage audit
#
# Output: four figures in output/figures/ + the audit, printed.
#
# Run with the package loaded (devtools::load_all()) or after library(atlaspm),
# and after tar_make() has built the primary-care targets.
# ---------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sf)
library(atlaspm)

FIG_DIR <- "output/figures"

targets::tar_load(c(gp_density_area, gp_density_geo, gp_density_audit))

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

wmean_density <- weighted.mean(supply$gp_density, supply$gp_population)

cat("\n-- indicator, by stratum ------------------------------------------\n")
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
  "\nPopulation-weighted mean: %.2f GP-equivalents per 1,000 (%.0f residents per GP-equivalent)\n",
  wmean_density, 1000 / wmean_density))

# --- 3. Histogram ------------------------------------------------------------
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
    subtitle = sprintf("GP-equivalents per 1,000 residents; dashed line = population-weighted mean (%.2f)",
                       wmean_density),
    x = "GP-equivalents per 1,000 residents",
    y = "Number of areas",
    caption = "One bar per area (comune, or NIL within Milan), unweighted by population."
  ) +
  theme_atlas()

p_hist_strata <- p_hist +
  facet_wrap(~ stratum, ncol = 1, scales = "free_y") +
  labs(subtitle = "By geography type; note the different y-scales")

# --- 4. Maps -----------------------------------------------------------------
# 4a. Continuous caseload. This is the same figure the pipeline builds as
#     fig_gp_density_map; it is rebuilt here so the script stands alone.

p_map_caseload <- plot_travel_time(
  gp_density_geo,
  value        = "gp_caseload",
  legend_title = "Residents per\nGP-equivalent",
  caption      = "Reciprocal of the supply indicator; darker = each GP-equivalent covers more residents."
) +
  labs(title    = "Primary-care caseload, 2023",
       subtitle = "Residents per GP-equivalent, by area")

# 4b. Binned and anchored on the ATS average, reusing the SMR map. Oriented as
#     a STRAIN ratio (local caseload / average caseload) rather than as supply,
#     because the palette encodes "high = worse than expected" and inverting
#     the measure is what keeps the colours honest.

strain_geo <- gp_density_geo |>
  mutate(gp_strain = wmean_density / gp_density)

p_map_strain <- plot_smr_map(
  strain_geo,
  value    = "gp_strain",
  breaks   = c(-Inf, 0.9, 0.95, 1.05, 1.1, Inf),
  title    = "Primary-care strain relative to the ATS average, 2023",
  subtitle = "Caseload ratio: 1 = the population-weighted average caseload",
  caption  = "Above 1 = each GP-equivalent covers more residents than the ATS average."
)
p_map_strain$labels$fill <- "Caseload ratio\n(local / average)"

# --- 5. Save -----------------------------------------------------------------
save_figure(p_hist,         "gp_supply_hist",        dir = FIG_DIR, height = 110)
save_figure(p_hist_strata,  "gp_supply_hist_strata", dir = FIG_DIR, height = 180)
save_figure(p_map_caseload, "gp_caseload_map",       dir = FIG_DIR)
save_figure(p_map_strain,   "gp_strain_map",         dir = FIG_DIR)

# --- 6. Note on the covariate ------------------------------------------------
# smr_geo_full carries `gp_density` (and `gp_density_z`) via add_covariate();
# `gp_strain` above is a presentational transform for this map only and is not
# a second covariate. Note also that this file is the density component alone:
# the proposal's "GP Strain Index" also folds in the share of residents without
# a registered GP and the share of GPs over 62, which the upstream extract can
# supply as prop_unassigned / prop_gp_65plus (see compute_density_by_area()).
