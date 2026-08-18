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
  tar_target(lookup_causes, get_input_data_path("avoidable_lookup_v3.csv")),
  tar_target(pop_table, get_input_data_path("pop_finale.csv") |>
               readr::read_delim(delim = ";", show_col_types = FALSE)),
  tar_target(pop_nil, get_input_data_path("pop_nil.csv") |>
               readr::read_delim(delim = ",", show_col_types = FALSE)),
  tar_target(pop_shp,
             get_input_data_path("geodata/Com01012025_g/Com01012025_g_WGS84.shp") |>
               sf::st_read(quiet = TRUE) |> sf::st_make_valid()),
  tar_target(nil_shp,
             get_input_data_path("geodata/ds964_nil_wm/NIL_WM.shp") |>
               sf::st_read(quiet = TRUE) |> sf::st_make_valid()),
  tar_target(area_shp, build_area_shp(nil_shp, pop_shp)),

  # IMPORT
  tar_target(mort_path, get_input_data_path("mort.csv")),
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
  tar_target(mort_count_area,
             dplyr::filter(mort_count, area_residenza %in% area_shp$area)),
  # one row per decedent - Table 1 and the reported N come from here,
  # never from nrow(mort_count)
  tar_target(deaths, build_deaths(mort_count)),
  tar_target(deaths_area, build_deaths(mort_count_area)),
  tar_target(layer_sizes,
             dplyr::count(mort_count_area, mechanism, wt = weight,
                          name = "deaths")),
  tar_target(pop_area_table,
             build_pop_area_table(pop_table, pop_nil, mort_count)),
  tar_target(mort_crude, preprocess_cmr(mort_count_area,
                                        pop_area_table,
                                        group_var = "area",
                                        mort_col  = "area_residenza",
                                        pop_col   = "area",
                                        pad_area  = FALSE)),
  tar_target(mort_smr, preprocess_smr(mort_count_area,
                                      pop_area_table,
                                      group_var = "area",
                                      mort_col  = "area_residenza",
                                      pop_col   = "area",
                                      pad_area  = FALSE)),
  tar_target(smr_geo, add_geo(mort_smr, area_shp,
                              data_key = "area",
                              shp_key  = "area",
                              pad_keys = FALSE)),
  tar_target(C_matrix, build_adjacency(smr_geo)),

  # POLLUTION (S8) -------------------------------------------------------
  tar_target(aq_dir, get_input_data_path("eea_aq"), format = "file"),
  tar_target(aq_manifest, discover_aq_files(aq_dir)),
  tar_target(pollution_area,
             build_pollution_area(smr_geo, manifest = aq_manifest)),
  tar_target(aq_ranks, check_aq_ranks(pollution_area, years = c(2013L, 2023L))),
  tar_target(
    pollutant_selection,
    pollutant_selection_table(
      dplyr::left_join(pollution_area, deprivation_area, by = "area"),
      year   = 2023L,
      di_col = "di_score",
      z_low  = NULL,
      ranks  = aq_ranks
    )
  ),

  # STROKE NETWORK ACCESS ----------------------------------------
  tar_target(stroke_centres_path,
             get_input_data_path("stroke_centres_dgr7473.csv"),
             format = "file"),
  tar_target(stroke_centres_raw, import_stroke_centres(stroke_centres_path)),

  tar_target(
    stroke_centres,
    geocode_stroke_centres(
      stroke_centres_raw,
      registry = NULL,
      # Populate from the failures reported by check_stroke_centres()
      overrides = list(
        "Ospedale di Circolo di Varese" = c(45.80989, 8.8391),
        "Ospedale di Circolo Desio" = c(45.62642, 9.19632),
        "Policlinico San Marco Zingonia" = c(45.60409, 9.5911),
        "Istituto Clinico S. Anna" = c(45.55361, 10.18027),
        "Fondazione IRCCS Policlinico San Matteo" = c(45.19622, 9.14884),
        "Ospedale G. Salvini" = c(45.58284, 9.09504)
      )
    )
  ),

  # tar_target(stroke_centres_ok,
  #            check_stroke_centres(stroke_centres, pop_shp,
  #                                 name_col = "COMUNE")),

  # Origins: ISTAT census sections. Column names differ between the 2011 and
  # 2021 releases -- build_section_points() fails loudly with the available
  # names rather than guessing.
  tar_target(sez_shp,
             get_input_data_path("geodata/R03_21/R03_21_WGS84.shp") |>
               sf::st_read(quiet = TRUE) |> sf::st_make_valid()),
  tar_target(section_points, build_section_points(sez_shp, area_shp,
                                                  pop_col = "POP21")),
  tar_target(urban_mask, build_urban_mask(sez_shp, tipo_loc_col = "TIPO_LOC")),

  # Routing area: the modelled areas plus a buffer wide enough that routes
  # leaving the study area are not truncated. 40 km is generous for ATS
  # Milano; raise it if the study area ever extends into alpine comuni.
  tar_target(stroke_aoi,
             smr_geo |> sf::st_transform(32632L) |> sf::st_union() |>
               sf::st_buffer(40000)),

  # The slow target: 10-40 min and 8-16 GB. Cached by targets thereafter.
  # speed_model = "areu" reproduces the 33/60/90 km/h assumptions behind the
  # DGR's own centralisation maps, so the output is commensurable with the
  # regional 45-minute threshold. "osm" models a private car instead.
  tar_target(stroke_network,
             build_stroke_network(stroke_aoi, urban_mask,
                                  speed_model = "areu",
                                  osm_dir = get_input_data_path("osm"))),

  tar_target(stroke_times,
             build_stroke_times(stroke_network, stroke_centres,
                                section_points)),
  tar_target(stroke_area, build_stroke_area(stroke_times)),
  tar_target(diag_stroke, check_stroke_access(stroke_times, stroke_area,
                                              print = FALSE)),

  tar_target(smr_geo_stroke, add_stroke_access(smr_geo, stroke_area)),

  # Model. pop_share_over_45min_hub is the alternative worth fitting alongside
  # t_hub_mean_z: bounded on [0,1], and anchored to the DGR's own decision
  # rule rather than to an arbitrary scale.
  tar_target(
    model_stroke,
    fit_bym2(
      smr_geo_stroke, C_matrix,
      formula      = total_obs ~ t_hub_mean_z + offset(log(total_exp)),
      scale_factor = scale_factor,
      cores        = 4,
      refresh      = 0
    )
  ),
  tar_target(diag_stroke_fit, check_bym2_fit(model_stroke, print = FALSE)),
  tar_target(smr_geo_stroke_bym2,
             augment_bym2(smr_geo_stroke, model_stroke, threshold = 1.10)),
  tar_target(
    map_smr_stroke,
    plot_smr_map(
      smr_geo_stroke_bym2,
      value    = "bym2_rr",
      title    = "BYM2-smoothed preventable mortality, adjusted for stroke network access",
      subtitle = "ICAR-smoothed relative risk, travel time to nearest level-II hub"
    )
  ),

  # covariates attached to the modelling sf by the same route as IVSM/DI
  tar_target(smr_geo_poll, add_pollution(smr_geo, pollution_area)),

  tar_target(ivsm_area,        expand_cov_to_area(ivsm_raw,   area_shp$area, by = "comune")),
  tar_target(deprivation_area, expand_cov_to_area(deprivation, area_shp$area, by = "comune")),
  tar_target(smr_geo_ivsm, add_covariate(smr_geo, ivsm_area,        var = "ivsm",     by = "area")),
  tar_target(smr_geo_di,   add_covariate(smr_geo, deprivation_area, var = "di_score", by = "area")),

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

  ## PER-MECHANISM BYM2 (seven models -> one faceted smoothed map)
  tar_target(
    models_mechanism,
    fit_bym2_mechanisms(smr_geo, C_matrix, scale_factor)   # reuses shared scale_factor
  ),
  tar_target(
    smr_geo_mech_bym2,
    augment_bym2_mechanisms(smr_geo, models_mechanism, threshold = 1.10)
  ),
  tar_target(
    map_smr_facets_bym2,
    plot_smr_facets(
      smr_geo_mech_bym2,
      cols         = dplyr::matches("^M_.*_bym2$"),
      breaks       = c(-Inf, 0.90, 0.95, 1.05, 1.10, Inf),
      strip_suffix = "_bym2$",
      title    = "BYM2-smoothed preventable mortality by mechanism, by comune",
      subtitle = "ICAR-smoothed relative risk (model SMR); 1 = matches the age-sex expectation",
      caption  = "Per-mechanism BYM2 on shared adjacency; bins shared across panels."
      # title    = "Mortalità prevenibile per comune, suddivisa per funzione",
      # subtitle = "Rischio Relativo; 1 = corrisponde a quanto atteso per sesso ed età",
      # caption  = "BYM2 su matrice di adiacenza; raggruppamento unico per territorio"
    )
  ),
  tar_target(
    map_exceedance_facets_bym2,
    plot_exceedance_facets(smr_geo_mech_bym2)   # threshold read from the stored attribute
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
    mort_smr, ivsm_area,
    index_col = "ivsm",
    ref_line  = 100,
    xlab      = "IVSM (social & material vulnerability index)",
    title     = "Standardised mortality vs social/material vulnerability, by comune",
    subtitle  = "Each point a comune; x = IVSM (national average = 100), y = indirectly standardised rate")),

  tar_target(scatter_smr_di, plot_scatter_smr_index(
    mort_smr, deprivation_area,
    index_col = "di_score",
    ref_line  = 0,
    xlab      = "Italian Deprivation Index (sum of national z-scores)",
    title     = "Standardised mortality vs Deprivation Index, by comune",
    subtitle  = "Each point a comune; x = Deprivation Index, y = indirectly standardised rate")),

  # MAPS
  # area_shp with pad_keys = FALSE, not pop_shp with the default padding:
  # pop_shp is keyed on PRO_COM_T and the default pad coerces through
  # as.integer(), so every NIL key ("015146_79") would become NA and all ~80
  # Milan NILs would collapse into one row.
  tar_target(map_smr_overall,
             mort_smr |>
               add_geo(area_shp, data_key = "area", shp_key = "area",
                       pad_keys = FALSE) |>
               plot_smr_map()),
  tar_target(map_smr_mechanism,
             mort_smr |>
               add_geo(area_shp, data_key = "area", shp_key = "area",
                       pad_keys = FALSE) |>
               plot_smr_facets()),

  # REPORT
  tar_quarto(explore_mort_count, path = "reports\\mortality_explore.qmd"),


  tar_target(JustDontCareLastComma, NULL)
)
