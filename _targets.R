library(targets)
library(tarchetypes)
library(crew)  # parallel computing

controller <- crew::crew_controller_local(
  name = "atlaspm_controller",
  workers = 1
)

# Set target-specific options such as packages.
tar_option_set(
  # error handling
  error = "abridge", # "continue" (do other), "null" (NULL if error)
  workspace_on_error = TRUE,
  format = "rds",
  # parallel computing
  storage = "worker",
  retrieval = "worker",
  controller = controller
)

# Define custom functions and other global objects.
# This is where you write source(\"R/functions.R\")
# if you keep your functions in external scripts.
tar_source()


# End this file with a list of target objects.
list(
  # LOOK UP TABLES
  tar_target(lookup_causes, get_input_data_path("avoidable_lookup_v2.csv")),
  tar_target(pop_table, get_input_data_path("pop_finale.csv") |>
               readr::read_delim(delim = ";", show_col_types = FALSE)),
  tar_target(pop_shp,
             get_input_data_path("geodata/Com01012025_g/Com01012025_g_WGS84.shp") |>
               sf::st_read(quiet = TRUE) |> sf::st_make_valid()),

  # IMPORT
  tar_target(mort_path, get_input_data_path("mort_2023.csv")),
  tar_target(mort_raw, import_mortality(mort_path)),
  tar_target(ivsm_path, get_input_data_path("Indicatori_Regione_Lombardia.csv")),
  tar_target(ivsm_raw, import_ivsm(ivsm_path)),
  tar_target(census_2023, get_input_data_path("census_2023") |>
               import_census_2023()),
  tar_target(deprivation, build_deprivation_proxy(census_2023, mort_raw)),

  # PREPROCESSING
  tar_target(mort_count, preprocess_mortality(mort_raw,
                                              lookup_causes,
                                              code_col = "causa",
                                              age_col = "eta")),
  tar_target(mort_crude, preprocess_cmr(mort_count,
                                        pop_table)),
  tar_target(mort_smr, preprocess_smr(mort_count,
                                      pop_table)),
  tar_target(smr_geo, add_geo(mort_smr, pop_shp, data_key = "comune")),
  tar_target(C_matrix, build_adjacency(smr_geo)),

  tar_target(smr_geo_ivsm, add_covariate(smr_geo, ivsm_raw, var = "ivsm")),
  tar_target(smr_geo_di, add_covariate(smr_geo, deprivation, var = "di_score")),

  # BYM MODEL ------------------------------------------------------------
  tar_target(scale_factor, compute_scale_factor(C_matrix)),

  ## BASE
  tar_target(
    model_base,
    fit_bym2(
      smr_geo, C_matrix,
      formula      = total_obs ~ offset(log(total_exp)),
      scale_factor = scale_factor,                       # <- pass it in
      cores        = 4,
      refresh      = 0
    )
  ),
  tar_target(diag_base, check_bym2_fit(model_base, print = FALSE)),  # stored metrics
  tar_target(smr_geo_bym2, augment_bym2(smr_geo, model_base, threshold = 1.10)),
  tar_target(
    map_smr_bym2,
    plot_smr_map(
      smr_geo_bym2,
      value    = "bym2_rr",
      title    = "BYM2-smoothed preventable mortality, by comune",
      subtitle = "ICAR-smoothed relative risk (base model, no covariates)"
    )
  ),
  tar_target(
    map_exceedance,
    plot_exceedance_map(smr_geo_bym2)        # label auto-derived from stored threshold
  ),

  ## IVSM
  tar_target(
    model_ivsm,
    fit_bym2(
      smr_geo_ivsm, C_matrix,
      formula      = total_obs ~ ivsm_z + offset(log(total_exp)),
      scale_factor = scale_factor,           # <- same factor, same graph
      cores        = 4,
      refresh      = 0
    )
  ),
  tar_target(diag_ivsm, check_bym2_fit(model_ivsm, print = FALSE)),
  tar_target(smr_geo_ivsm_bym2, augment_bym2(smr_geo_ivsm, model_ivsm, threshold = 1.10)),
  tar_target(
    map_smr_ivsm,
    plot_smr_map(
      smr_geo_ivsm_bym2,
      value    = "bym2_rr",
      title    = "BYM2-smoothed preventable mortality, adjusted for deprivation (IVSM)",
      subtitle = "ICAR-smoothed relative risk, IVSM covariate model"
    )
  ),
  tar_target(
    map_exceedance_ivsm,
    plot_exceedance_map(smr_geo_ivsm_bym2)
  ),

  ## Deprivation Index
  tar_target(
    model_di,
    fit_bym2(
      smr_geo_di, C_matrix,
      formula      = total_obs ~ di_score_z + offset(log(total_exp)),
      scale_factor = scale_factor,           # <- same factor, same graph
      cores        = 4,
      refresh      = 0
    )
  ),
  tar_target(diag_di, check_bym2_fit(model_di, print = FALSE)),
  tar_target(smr_geo_di_bym2, augment_bym2(smr_geo_di, model_di, threshold = 1.10)),
  tar_target(
    map_smr_di,
    plot_smr_map(
      smr_geo_di_bym2,
      value    = "bym2_rr",
      title    = "BYM2-smoothed preventable mortality, adjusted for deprivation (DI)",
      subtitle = "ICAR-smoothed relative risk, DI covariate model"
    )
  ),
  tar_target(
    map_exceedance_di,
    plot_exceedance_map(smr_geo_di_bym2)
  ),

  tar_target(
    bym2_comparison,
    compare_bym2(
      fits = list(base = model_base, ivsm = model_ivsm, di = model_di),
      data = list(base = smr_geo,    ivsm = smr_geo_ivsm, di = smr_geo_di),
      param_labels = c("beta[1]" = "index_z")
    )
  ),

  # SCATTER
  tar_target(scatter_cmr_isr_overall,   plot_cmr_isr(mort_crude, mort_smr)),
  tar_target(scatter_cmr_isr_mechanism, plot_cmr_isr_facets(mort_crude, mort_smr)),
  tar_target(scatter_smr_ivsm, plot_scatter_smr_index(
    mort_smr, ivsm_raw,
    index_col = "ivsm",
    ref_line  = 100,
    xlab      = "IVSM (social & material vulnerability index)",
    title     = "Standardised mortality vs social/material vulnerability, by comune",
    subtitle  = "Each point a comune; x = IVSM (national average = 100), y = indirectly standardised rate")),

  tar_target(scatter_smr_di, plot_scatter_smr_index(
    mort_smr, deprivation,
    index_col = "di_score",
    ref_line  = 0,
    xlab      = "Italian Deprivation Index (sum of national z-scores)",
    title     = "Standardised mortality vs Deprivation Index, by comune",
    subtitle  = "Each point a comune; x = Deprivation Index, y = indirectly standardised rate")),

  # MAPS
  tar_target(map_smr_overall, mort_smr |> add_geo(pop_shp, data_key = "comune") |>
               plot_smr_map()),
  tar_target(map_smr_mechanism, mort_smr |> add_geo(pop_shp, data_key = "comune") |>
               plot_smr_facets()),

  # REPORT
  tar_quarto(explore_mort_count, path = "reports\\mortality_explore.qmd"),


  tar_target(JustDontCareLastComma, NULL)
)
