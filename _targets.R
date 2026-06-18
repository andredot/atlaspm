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

  # BYM MODEL
  ## BASE
  tar_target(
    model_base,
    fit_bym2(
      smr_geo, C_matrix,
      formula = total_obs ~ offset(log(total_exp)),   # base: spatial effect only
      chains  = 4,
      cores = 4,
      iter    = 4000,
      control = list(max_treedepth = 12, adapt_delta = 0.95),
      refresh = 0
    )
  ),

  tar_target(smr_geo_bym2, augment_bym2(smr_geo, model_base, threshold = 1.10)),

  tar_target(
    map_smr_bym2,
    plot_smr_map(
      smr_geo_bym2,                       # <- not model_base
      value    = "bym2_rr",
      title    = "BYM2-smoothed preventable mortality, by comune",
      subtitle = "ICAR-smoothed relative risk (base model, no covariates)"
    )
  ),

  tar_target(
    map_exceedance,
    plot_exceedance_map(smr_geo_bym2,     # <- not model_base
                        threshold_label = "+20% (RR > 1.20)")
  ),

  ## IVSM
  tar_target(
    model_ivsm,
    fit_bym2(
      smr_geo_ivsm, C_matrix,
      formula = total_obs ~ ivsm_z + offset(log(total_exp)),   # covariate + spatial effect
      chains  = 4,
      iter    = 4000,
      control = list(max_treedepth = 12, adapt_delta = 0.95),
      refresh = 0
    )
  ),

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
    plot_exceedance_map(smr_geo_ivsm_bym2, threshold_label = "+20% (RR > 1.20)")
  ),

  # SCATTER
  tar_target(scatter_cmr_isr_overall,   plot_cmr_isr(mort_crude, mort_smr)),
  tar_target(scatter_cmr_isr_mechanism, plot_cmr_isr_facets(mort_crude, mort_smr)),
  tar_target(scatter_smr_ivsm, plot_scatter_smr_ivsm(mort_smr, ivsm_raw)),

    # MAPS
  tar_target(map_smr_overall, mort_smr |> add_geo(pop_shp, data_key = "comune") |>
               plot_smr_map()),
  tar_target(map_smr_mechanism, mort_smr |> add_geo(pop_shp, data_key = "comune") |>
               plot_smr_facets()),

  # REPORT
  tar_quarto(explore_mort_count, path = "reports\\mortality_explore.qmd"),


  tar_target(JustDontCareLastComma, NULL)
)
