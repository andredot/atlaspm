# ---------------------------------------------------------------------------
# Exploratory look at the primary-care supply indicator (2023)
#
# Input : indicatore_assistenza_2023.csv, one row per modelling `area`
#           area         mixed comune / Milan-NIL key ("015011", "015146_79")
#           numeratore   summed GP supply  = sum over residents of 1/list_size
#           denominatore person-time / residents in the denominator
#           indicatore   numeratore / denominatore  = GP-equivalents per resident
#         i.e. the `gp_density` column of compute_density_by_area().
#
# Output: three figures in output/figures/ + a printed join audit.
#
# Run with the package loaded (devtools::load_all()) or after library(atlaspm).
# ---------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sf)
library(atlaspm)

INDICATOR_CSV <- "indicatore_assistenza_2023.csv"   # adjust / wrap in get_input_data_path()
FIG_DIR       <- "output/figures"

# --- 1. Read -----------------------------------------------------------------
# `area` MUST be character: read as numeric it loses the leading zero and the
# NIL keys ("015146_79") collapse to NA.

supply_raw <- readr::read_csv(
  get_input_data_path(INDICATOR_CSV),
  col_types = readr::cols(area = readr::col_character(),
                          .default = readr::col_double())
)

# --- 2. Derive the three readable scalings -----------------------------------
# `indicatore` is ~7e-4, which no axis renders usefully. Three equivalent views:
#   gp_per_100k     supply per 100,000 residents      (higher = better)
#   patients_per_gp 1 / indicatore, the caseload       (higher = worse)
#   gp_strain       caseload relative to the ATS-wide, population-weighted mean
#                   (1 = as expected, >1 = thinner cover than average)
#
# gp_strain is the one that plugs into the package's SMR machinery: it is
# anchored on 1.0 and oriented so that "high = bad", which is the semantics the
# earthy diverging palette encodes.

wmean_supply <- weighted.mean(supply_raw$indicatore, supply_raw$denominatore)

supply <- supply_raw |>
  mutate(
    gp_per_100k     = indicatore * 1e5,
    patients_per_gp = 1 / indicatore,
    gp_strain       = wmean_supply / indicatore,
    stratum = case_when(
      grepl("_", area)          ~ "Milan (NIL)",
      substr(area, 1, 3) == "015" ~ "Milan province (comuni)",
      substr(area, 1, 3) == "098" ~ "Lodi province (comuni)",
      TRUE                      ~ "other"
    )
  )

cat("\n-- indicator, by stratum ------------------------------------------\n")
supply |>
  group_by(stratum) |>
  summarise(
    n            = n(),
    pop          = sum(denominatore),
    gp_per_100k  = median(gp_per_100k),
    caseload_p10 = quantile(patients_per_gp, 0.10),
    caseload_med = median(patients_per_gp),
    caseload_p90 = quantile(patients_per_gp, 0.90),
    .groups = "drop"
  ) |>
  print()

cat(sprintf("\nPopulation-weighted mean: %.1f GP-equivalents per 100,000 (%.0f residents per GP)\n",
            wmean_supply * 1e5, 1 / wmean_supply))

# --- 3. Geometry + join audit ------------------------------------------------
# area_shp is the NIL-aware layer (Milan split into NILs, everything else at
# comune level). Both keys are already zero-padded strings and one side holds
# non-numeric values, so pad_keys = FALSE is mandatory here: the default would
# push "015146_79" through as.integer() and collapse every NIL to a single NA.

area_shp <- targets::tar_read(area_shp)   # or: sf::st_read(...) |> build_area_shp()

shp_keys  <- area_shp$area
data_keys <- supply$area

cat("\n-- join audit -----------------------------------------------------\n")
cat(sprintf("areas in shapefile : %d\nrows in indicator  : %d\nmatched            : %d\n",
            length(shp_keys), length(data_keys),
            length(intersect(shp_keys, data_keys))))

orphan_rows <- setdiff(data_keys, shp_keys)   # indicator with nowhere to draw
empty_areas <- setdiff(shp_keys, data_keys)   # polygons with no indicator

if (length(orphan_rows)) {
  cat("\nIndicator rows with no polygon (dropped from the map):\n")
  supply |> filter(area %in% orphan_rows) |>
    select(area, numeratore, denominatore, indicatore) |> print()
}
if (length(empty_areas)) {
  cat("\nPolygons with no indicator (drawn grey):\n"); print(empty_areas)
}

supply_geo <- add_geo(
  supply, area_shp,
  data_key = "area", shp_key = "area",
  keep     = "geometry",   # keep every polygon, so gaps show as grey, not as holes
  pad_keys = FALSE
)

# --- 4. Histogram ------------------------------------------------------------
# One bar per AREA, not per resident, so the Milan NILs (~5k-65k people each)
# and the small Lodi comuni (~1.5k) carry equal weight. The dashed line is the
# population-weighted mean, which is the honest central value.

p_hist <- ggplot(supply, aes(x = gp_per_100k)) +
  geom_histogram(binwidth = 1.5, boundary = 0,
                 fill = "#5a8a7d", colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = wmean_supply * 1e5,
             linetype = "22", colour = "#8c2d04", linewidth = 0.6) +
  labs(
    title    = "Primary-care supply across the modelling areas, 2023",
    subtitle = sprintf("GP-equivalents per 100,000 residents; dashed line = population-weighted mean (%.1f)",
                       wmean_supply * 1e5),
    x = "GP-equivalents per 100,000 residents",
    y = "Number of areas",
    caption = "One bar per area (comune, or NIL within Milan), unweighted by population."
  ) +
  theme_atlas()

p_hist_strata <- p_hist +
  facet_wrap(~ stratum, ncol = 1, scales = "free_y") +
  labs(subtitle = "By geography type; note the different y-scales")

# --- 5. Maps -----------------------------------------------------------------
# 5a. Continuous caseload, reusing the package's covariate map. cividis
#     reversed puts the dark end on the high values, i.e. on the thinnest cover.

p_map_caseload <- plot_travel_time(
  supply_geo,
  value        = "patients_per_gp",
  legend_title = "Residents per\nGP-equivalent",
  caption      = "Reciprocal of the supply indicator; darker = each GP-equivalent covers more residents."
) +
  labs(title    = "Primary-care caseload, 2023",
       subtitle = "Residents per GP-equivalent, by area")

# 5b. Binned, anchored on the ATS average, reusing the SMR map. Same palette
#     and same reading as the mortality maps: terracotta = worse than expected.

p_map_strain <- plot_smr_map(
  supply_geo,
  value    = "gp_strain",
  breaks   = c(-Inf, 0.9, 0.95, 1.05, 1.1, Inf),
  title    = "Primary-care strain relative to the ATS average, 2023",
  subtitle = "Caseload ratio: 1 = the population-weighted average caseload",
  caption  = "Above 1 = each GP-equivalent covers more residents than the ATS average."
)
p_map_strain$labels$fill <- "Caseload ratio\n(local / average)"

# --- 6. Save -----------------------------------------------------------------
save_figure(p_hist,         "gp_supply_hist",        dir = FIG_DIR, height = 110)
save_figure(p_hist_strata,  "gp_supply_hist_strata", dir = FIG_DIR, height = 180)
save_figure(p_map_caseload, "gp_caseload_map",       dir = FIG_DIR)
save_figure(p_map_strain,   "gp_strain_map",         dir = FIG_DIR)

# --- 7. Next step ------------------------------------------------------------
# To take this into the model, gp_strain (or a z-scored gp_per_100k) is the
# beta_2 "GP Strain Index" term:
#   add_covariate(smr_geo, sf::st_drop_geometry(supply_geo), var = "gp_strain", by = "area")
# Note this file only carries the density component; the proposal's index also
# folds in the share of unassigned residents and the share of GPs over 62,
# which compute_density_by_area() returns as prop_unassigned / prop_gp_65plus.
