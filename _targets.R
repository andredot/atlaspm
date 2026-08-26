library(targets)
library(tarchetypes)
library(crew)  # parallel computing

controller <- crew::crew_controller_local(
  name = "atlaspm_controller",
  workers = 2
)

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

tar_source()

# ---------------------------------------------------------------------------
# STUDY CONSTANTS
# ---------------------------------------------------------------------------
STUDY_YEARS <- 2022:2024
N_YEARS     <- length(STUDY_YEARS)

# Exceedance: the relative-risk threshold, and the posterior probability above
# which an area is called. 0.80 follows Richardson et al. (2004).
RR_THRESHOLD <- 1.10
RR_THRESHOLD_STROKE <- 1.20
PROB_CUTOFF  <- 0.80

AQ_YEAR <- 2023L
NO2_Z   <- sprintf("no2_%d_z",  AQ_YEAR)
PM25_Z  <- sprintf("pm25_%d_z", AQ_YEAR)

RHS_FULL <- paste("di_score_z", NO2_Z,  "gp_density_z", sep = " + ")
RHS_PM25 <- paste("di_score_z", PM25_Z, "gp_density_z", sep = " + ")

# Objective 3. The DAG has one backdoor path from travel time to mortality:
# A <- hub siting <- U -> Y, where U is the latent urban-rural position. U is
# unmeasured, so these three - its measured children on the outcome side - are
# the proxy adjustment set. Reperfusion timing and case fatality are mediators
# and are deliberately absent.
STROKE_ADJUST  <- c("di_score_z", NO2_Z, "gp_density_z")
STROKE_EXPOSURE <- "t_hub_mean_z"
STROKE_RHS <- paste(c(STROKE_EXPOSURE, STROKE_ADJUST), collapse = " + ")

# ---------------------------------------------------------------------------
# MODEL SPECIFICATIONS
# ---------------------------------------------------------------------------
model_specs <- tibble::tribble(
  ~id,  ~label,                                     ~engine, ~rhs,
  "M1", "Intercept only",                           "glm",   "1",
  "M2", "Covariates, no random effect",             "glm",   RHS_FULL,
  "M3", "BYM2, no covariates",                      "bym2",  "1",
  "M4", "BYM2 + deprivation",                       "bym2",  "di_score_z",
  "M5", "BYM2 + deprivation, NO2, primary care",    "bym2",  RHS_FULL,
  "M6", "ESF + deprivation, NO2, primary care",     "esf",   RHS_FULL,
  "M7", "BYM2 + deprivation, PM2.5, primary care",  "bym2",  RHS_PM25)

# Sensitivity fits that sit outside the M1-M7 sequence.
sens_specs <- tibble::tribble(
  ~id,      ~label,                                          ~engine, ~rhs,
  "SIVSM",  "BYM2 + IVSM (alternative deprivation measure)",  "bym2",  "ivsm_z"
)

TERM_LABELS <- c(
  di_score_z      = "Deprivation Index (per SD)",
  ivsm_z          = "IVSM (per SD)",
  gp_density_z    = "Primary care density (per SD)",
  t_hub_mean_z    = "Travel time to thrombectomy hub (per SD)",
  t_centre_mean_z = "Travel time to nearest stroke centre (per SD)"
)
TERM_LABELS[NO2_Z]  <- "NO2 (per SD)"
TERM_LABELS[PM25_Z] <- "PM2.5 (per SD)"

# ---------------------------------------------------------------------------

list(

  # === LOOK-UP TABLES ======================================================
  tar_target(lookup_path, get_input_data_path("avoidable_lookup_v3.csv"),
             format = "file"),
  tar_target(lookup_causes, readr::read_csv(lookup_path,
                                            show_col_types = FALSE)),
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

  tar_target(area_shp_all, build_area_shp(nil_shp, pop_shp)),
  tar_target(area_shp, drop_unpopulated_areas(area_shp_all, pop_area_table)),
  tar_target(study_area_summary, summarise_study_area(area_shp, pop_area_table,
                                                      pop_year = STUDY_YEARS)),

  # === IMPORT ==============================================================
  tar_target(mort_path, get_input_data_path("mort.csv"), format = "file"),
  tar_target(mort_raw, import_mortality(mort_path)),
  tar_target(ivsm_path, get_input_data_path("Indicatori_Regione_Lombardia.csv")),
  tar_target(ivsm_raw, import_ivsm(ivsm_path)),
  tar_target(census_2023, get_input_data_path("census_2023") |>
               import_census_2023()),

  # === PREPROCESSING =======================================================
  tar_target(mort_count, preprocess_mortality(mort_raw,
                                              lookup_causes,
                                              code_col = "causa",
                                              age_col  = "eta")),
  tar_target(mort_count_area,
             dplyr::filter(mort_count, area_residenza %in% area_shp$area)),

  # One row per decedent. Table 1 and every reported N come from here, never
  # from nrow(mort_count), which holds one row per death x avoidability arm.
  tar_target(deaths, build_deaths(mort_count)),
  tar_target(deaths_area, build_deaths(mort_count_area)),
  tar_target(death_arms, build_death_arms(mort_count_area)),
  tar_target(layer_sizes,
             dplyr::count(mort_count_area, mechanism, wt = weight,
                          name = "deaths")),

  tar_target(pop_area_table,
             build_pop_area_table(pop_table, pop_nil, mort_count)),

  tar_target(mort_crude, preprocess_cmr(mort_count_area, pop_area_table,
                                        group_var = "area",
                                        mort_col  = "area_residenza",
                                        pop_col   = "area",
                                        pad_area  = FALSE,
                                        pop_year  = STUDY_YEARS)),
  tar_target(mort_smr, preprocess_smr(mort_count_area, pop_area_table,
                                      group_var = "area",
                                      mort_col  = "area_residenza",
                                      pop_col   = "area",
                                      pad_area  = FALSE,
                                      pop_year  = STUDY_YEARS)),
  tar_target(smr_geo, add_geo(mort_smr, area_shp,
                              data_key = "area", shp_key = "area",
                              pad_keys = FALSE)),
  tar_target(crude_geo, add_geo(mort_crude, area_shp,
                                data_key = "area", shp_key = "area",
                                pad_keys = FALSE)),
  tar_target(C_matrix, build_adjacency(smr_geo)),
  tar_target(scale_factor, compute_scale_factor(C_matrix)),

  tar_target(flow, flow_counts(mort_raw, mort_count, mort_count_area,
                               area_shp)), # STROBE/RECORD flow

  # === DEPRIVATION (at modelling-area resolution) ==========================
  tar_target(sez_shp,
             get_input_data_path("geodata/R03_21/R03_21_WGS84.shp") |>
               sf::st_read(quiet = TRUE) |> sf::st_make_valid()),
  tar_target(section_xwalk, build_section_area_xwalk(sez_shp, area_shp)),
  tar_target(deprivation_area, build_deprivation(census_2023, section_xwalk)),
  tar_target(deprivation_resolution,
             check_deprivation_resolution(deprivation_area)),

  tar_target(ivsm_area, expand_cov_to_area(ivsm_raw, area_shp$area,
                                           by = "comune")),

  # === POLLUTION ===========================================================
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

  # === PRIMARY CARE ========================================================
  tar_target(gp_density_path,
             get_input_data_path("indicatore_assistenza_2023.csv"),
             format = "file"),
  tar_target(gp_density_area, import_gp_density(gp_density_path)),
  tar_target(gp_density_audit, check_gp_density(gp_density_area, area_shp)),
  tar_target(gp_density_geo,
             add_geo(gp_density_area, area_shp,
                     data_key = "area", shp_key = "area",
                     keep = "geometry", pad_keys = FALSE)),

  # === STROKE NETWORK ACCESS ===============================================
  tar_target(stroke_centres_path,
             get_input_data_path("stroke_centres_dgr7473.csv"),
             format = "file"),
  tar_target(stroke_centres_raw, import_stroke_centres(stroke_centres_path)),
  tar_target(
    stroke_centres,
    geocode_stroke_centres(
      stroke_centres_raw,
      registry = NULL,
      overrides = list(
        "Ospedale di Circolo di Varese"           = c(45.80989, 8.8391),
        "Ospedale di Circolo Desio"               = c(45.62642, 9.19632),
        "Policlinico San Marco Zingonia"          = c(45.60409, 9.5911),
        "Istituto Clinico S. Anna"                = c(45.55361, 10.18027),
        "Fondazione IRCCS Policlinico San Matteo" = c(45.19622, 9.14884),
        "Ospedale G. Salvini"                     = c(45.58284, 9.09504),
        "Ospedale di Chiari"                      = c(45.53827, 9.93333),
        "Istituto Clinico Humanitas"              = c(45.37363,	9.16496)
      )
    )
  ),
  # --- ROUTING: computed once, then read back -------------------------------
  #
  # The network build is 10-40 minutes and 8-16 GB and needs an OSM extract on
  # disk, and nothing downstream of `stroke_times` touches the network object.
  # The routing has been run and its output saved, so the pipeline now resumes
  # from the table.
  #
  # The functions below are unchanged and still tested; only the pipeline's
  # entry point has moved. TO RE-RUN THE ROUTING - after a change to the
  # centre list, the speed model, or the section origins - comment out the
  # `stroke_times` target below and uncomment this block, then save the result
  # back to the shared folder:
  #
  #   saveRDS(targets::tar_read(stroke_times),
  #           file.path(Sys.getenv("PRJ_SHARED_PATH"),
  #                     Sys.getenv("INPUT_DATA_FOLDER"), "stroke_times.rds"))
  #
  # tar_target(section_points, build_section_points(sez_shp, area_shp,
  #                                                 pop_col = "POP21")),
  # tar_target(urban_mask, build_urban_mask(sez_shp, tipo_loc_col = "TIPO_LOC")),
  # tar_target(stroke_aoi,
  #            smr_geo |> sf::st_transform(32632L) |> sf::st_union() |>
  #              sf::st_buffer(40000)),
  # tar_target(stroke_network,
  #            build_stroke_network(stroke_aoi, urban_mask,
  #                                 speed_model = "areu",
  #                                 osm_dir = get_input_data_path("osm"))),
  # tar_target(stroke_times,
  #            build_stroke_times(stroke_network, stroke_centres,
  #                               section_points)),
  tar_target(stroke_times_path,
             get_input_data_path("stroke_times.csv"),
             format = "file"),
  tar_target(stroke_times,
             import_stroke_times(stroke_times_path, areas = area_shp$area)),

  tar_target(stroke_area, build_stroke_area(stroke_times)),
  tar_target(diag_stroke, check_stroke_access(stroke_times, stroke_area,
                                              print = FALSE)),

  # === THE MODELLING FRAME =================================================
  tar_target(
    smr_geo_full,
    smr_geo |>
      add_covariate(deprivation_area, var = "di_score",   by = "area") |>
      add_covariate(ivsm_area,        var = "ivsm",       by = "area") |>
      add_covariate(gp_density_area,  var = "gp_density", by = "area") |>
      add_pollution(pollution_area) |>
      add_stroke_access(stroke_area)
  ),
  tar_target(covariate_table, covariate_summary(smr_geo_full)),
  tar_target(covariate_correlations, covariate_correlation(smr_geo_full)),

  # === PRE-MODEL SPATIAL STRUCTURE =========================================
  tar_target(moran_crude, moran_test_raw(smr_geo$total_smr, C_matrix)),

  # === MODELS M1-M7 ========================================================
  tar_map(
    values = model_specs,
    names  = "id",
    tar_target(fit,   fit_model(smr_geo_full, rhs = rhs, engine = engine,
                                C = C_matrix, scale_factor = scale_factor,
                                refresh = 0)),
    tar_target(diag,  check_bym2_fit(fit, print = FALSE)),
    tar_target(aug,   augment_bym2(smr_geo_full, fit,
                                   pop_col   = "person_years",
                                   threshold = RR_THRESHOLD)),
    tar_target(moran, moran_test(fit, smr_geo_full, C_matrix))
  ),
  tar_map(
    values = sens_specs,
    names  = "id",
    tar_target(fit,  fit_model(smr_geo_full, rhs = rhs, engine = engine,
                               C = C_matrix, scale_factor = scale_factor,
                               refresh = 0)),
    tar_target(diag, check_bym2_fit(fit, print = FALSE)),
    tar_target(aug,  augment_bym2(smr_geo_full, fit,
                                  pop_col   = "person_years",
                                  threshold = RR_THRESHOLD))
  ),

  tar_target(model_fits, list(M1 = fit_M1, M2 = fit_M2, M3 = fit_M3,
                              M4 = fit_M4, M5 = fit_M5, M6 = fit_M6,
                              M7 = fit_M7)),
  tar_target(model_diagnostics, collect_diagnostics(model_fits)),
  tar_target(model_coefficients,
             collect_coefficients(model_fits[c("M2", "M4", "M5", "M6", "M7")],
                                  term_labels = TERM_LABELS)),
  tar_target(sensitivity_coefficients,
             collect_coefficients(list(M4 = fit_M4, SIVSM = fit_SIVSM),
                                  term_labels = TERM_LABELS)),

  tar_target(loo_comparison,
             compare_bym2(fits = model_fits,
                          data = rep(list(smr_geo_full),
                                     length(model_fits)))),
  tar_target(tbl_loo, collect_loo(loo_comparison)),
  tar_target(pareto_audit, pareto_k_summary(loo_comparison, smr_geo_full,
                                            model = "M5")),
  tar_target(moran_residual,
             stack_by_model(list(M1 = moran_M1, M2 = moran_M2, M3 = moran_M3,
                                 M4 = moran_M4, M5 = moran_M5, M6 = moran_M6,
                                 M7 = moran_M7))),

  # === THE ESTIMAND ========================================================
  tar_target(excess, residual_excess(fit_M5, smr_geo_full,
                                     n_years     = N_YEARS,
                                     threshold   = RR_THRESHOLD,
                                     prob_cutoff = PROB_CUTOFF)),
  tar_target(variance_reduction,
             variance_decomposition(fit_M3, fit_M5, smr_geo_full,
                                    labels = c("M3", "M5"))),

  # Testable implication of the DAG: given deprivation and the spatial term,
  # primary care capacity should be uncorrelated with the residuals of M4.
  tar_target(cond_independence,
             test_conditional_independence(fit_M4, smr_geo_full,
                                           var = "gp_density_z")),

  # === PER-MECHANISM MODELS (all seven strata) =============================
  tar_target(models_mechanism,
             fit_bym2_mechanisms(smr_geo_full, C_matrix, scale_factor)),
  tar_target(smr_geo_mech_bym2,
             augment_bym2_mechanisms(smr_geo_full, models_mechanism,
                                     threshold = RR_THRESHOLD)),
  tar_target(mechanism_table,
             collect_mechanisms(models_mechanism, smr_geo_mech_bym2,
                                prob_cutoff = PROB_CUTOFF)),
  tar_target(mechanism_diagnostics, collect_diagnostics(models_mechanism)),
  tar_target(mechanism_concordance, rank_concordance(smr_geo_mech_bym2)),

  # === TRACER AND CONTROLS =================================================
  tar_target(smr_geo_tracer, attach_tracer_outcomes(smr_geo_full)),
  tar_target(control_fits,
             fit_controls(smr_geo_allage, C_matrix, scale_factor,
                          exposure     = STROKE_EXPOSURE,
                          deprivation  = "di_score_z",
                          adjust       = STROKE_ADJUST,
                          engine       = "glm",
                          nce_exposure = NULL,
                          tracer_obs   = "i63_obs",
                          tracer_exp   = "i63_exp",
                          refresh      = 0)),
  tar_target(control_table,
             collect_controls(control_fits, labels = TERM_LABELS,
                              tracer_outcome = "I63 cerebral infarction, all ages",
                              exposure       = STROKE_EXPOSURE,
                              deprivation    = "di_score_z",
                              nce_exposure   = NULL)),
  tar_target(control_diagnostics, collect_diagnostics(control_fits)),
  tar_target(stroke_times_table, stroke_time_summary(smr_geo_full)),

  # === FIGURES =============================================================
  tar_target(fig_cmr_map,
             plot_cmr_map(crude_geo,
                          caption = "Crude rate, all avoidable causes, 2022-2024, per 100,000 person-years.")),
  tar_target(fig_smr_map, plot_smr_map(smr_geo, title = NULL)),
  tar_target(fig_rr_map,
             plot_smr_map(aug_M3, value = "bym2_rr", title = NULL,
                          subtitle = "BYM2-smoothed relative risk (M3, no covariates)")),
  tar_target(fig_exceedance, plot_exceedance_map(aug_M3)),
  tar_target(fig_excess_map, plot_excess_map(excess, smr_geo_full)),
  tar_target(fig_excess_flagged,
             plot_excess_map(excess, smr_geo_full, flagged_only = TRUE)),
  tar_target(fig_forest, plot_forest(model_coefficients, facet = "model")),
  tar_target(fig_controls, plot_forest(control_table, facet = "role")),
  tar_target(fig_rr_compare,
             plot_rr_compare(aug_M5, aug_M6,
                             labels = c("BYM2 (M5)", "ESF (M6)"))),
  # label_fun = mechanism_label turns "M_lifestyle_and_nc_ds" into
  # "Lifestyle and NCDs"; ncol = 3 gives seven panels room to be maps rather
  # than thumbnails.
  tar_target(fig_mech_facets,
             plot_smr_facets(smr_geo_mech_bym2,
                             cols   = dplyr::matches("^M_.*_bym2$"),
                             strip_suffix = "_bym2$",
                             label_fun = mechanism_label,
                             ncol = 3,
                             title = NULL)),
  tar_target(fig_mech_exceedance,
             plot_exceedance_facets(smr_geo_mech_bym2,
                                    label_fun = mechanism_label,
                                    ncol = 3)),
  tar_target(fig_concordance, plot_concordance(mechanism_concordance)),
  tar_target(fig_gp_density_map,
             plot_travel_time(
               gp_density_geo,
               value        = "gp_caseload",
               legend_title = "Residents per\nGP-equivalent",
               caption      = paste("Reciprocal of the primary-care density indicator;",
                                    "darker = each GP-equivalent covers more residents."))),
  tar_target(scatter_smr_gp,
             plot_scatter_smr_index(
               mort_smr, gp_density_area,
               index_col = "gp_density",
               xlab      = "GP-equivalents per 1,000 residents",
               title     = NULL,
               subtitle  = "Each point one areal unit")),
  tar_target(fig_pollution_pair, plot_pollution_pair(smr_geo_full)),
  tar_target(fig_pollution_shared,
             plot_pollution_pair(smr_geo_full, shared_scale = TRUE)),

  tar_target(scatter_cmr_isr_overall, plot_cmr_isr(mort_crude, mort_smr)),
  tar_target(scatter_smr_di,
             plot_scatter_smr_index(
               mort_smr, deprivation_area,
               index_col = "di_score",
               ref_line  = 0,
               xlab      = "Italian Deprivation Index (sum of national z-scores)",
               title     = NULL,
               subtitle  = "Each point one areal unit")),

  # === DESCRIPTIVES AND TABLES =============================================
  tar_target(desc_numbers, desc_stats(deaths_area, arms = death_arms,
                                      n_years = N_YEARS)),
  tar_target(table_one, tbl_one(deaths_area)),
  tar_target(appendix_b, tbl_lookup(lookup_causes)),
  tar_target(appendix_b_adapted, tbl_lookup(lookup_causes,
                                            adapted_only = TRUE)),
  tar_target(session_tbl, session_table()),

  # === REPORT ==============================================================
  tar_quarto(thesis_results, path = file.path("reports", "thesis_results.qmd")),

  # ==========================================================================
  # APPENDIX C: DEPRIVATION STABILITY, 2011 vs 2023
  # ==========================================================================
  tar_target(census_2011,
             get_input_data_path("census_2011") |> import_census_2011()),
  tar_target(sez_shp_2011,
             get_input_data_path("geodata/R03_11/R03_11_WGS84.shp") |>
               sf::st_read(quiet = TRUE) |> sf::st_make_valid()),
  tar_target(section_xwalk_2011,
             build_section_area_xwalk(sez_shp_2011, area_shp,
                                      sez_key = "SEZ2011")),
  tar_target(deprivation_2011,
             build_deprivation(census_2011, section_xwalk_2011,
                               sez_key = "SEZ2011", vintage = "2011")),
  tar_target(deprivation_stability,
             check_deprivation_stability(deprivation_2011, deprivation_area)),

  # ==========================================================================
  # STROKE SUB-MODEL (Objective 3)
  # ==========================================================================
  # The 0-74 cerebrovascular comparison keeps its own BYM2 fit. It used to
  # borrow control_fits$tracer, which now carries the all-age I63 reference
  # outcome - the two tables would silently have been comparing different
  # outcomes row for row.
  tar_target(fit_cvd_hub_bym2,
             fit_model(smr_geo_tracer, rhs = "t_hub_mean_z", engine = "bym2",
                       C = C_matrix, scale_factor = scale_factor,
                       obs_col = "cvd_obs", exp_col = "cvd_exp", refresh = 0)),
  tar_target(
    tracer_exposures,
    list(
      `Hub (SU II): thrombectomy and neurosurgery` = fit_cvd_hub_bym2,
      `Any accredited stroke centre (SU I or II)`  = fit_model(
        smr_geo_tracer, rhs = "t_centre_mean_z", engine = "bym2",
        C = C_matrix, scale_factor = scale_factor,
        obs_col = "cvd_obs", exp_col = "cvd_exp", refresh = 0)
    )
  ),
  tar_target(tracer_exposure_table,
             collect_coefficients(tracer_exposures, term_labels = TERM_LABELS)),

  # Specification sensitivity: the exposure is a smooth spatial surface and so
  # is the BYM2 random effect. When two terms describe the same variation the
  # coefficient is pulled toward zero whether or not an effect exists.
  tar_target(
    tracer_engines,
    list(
      `No spatial term (GLM)` = fit_model(
        smr_geo_tracer, rhs = "t_hub_mean_z", engine = "glm",
        obs_col = "cvd_obs", exp_col = "cvd_exp", refresh = 0),
      `BYM2` = fit_cvd_hub_bym2,
      `ESF`  = fit_model(
        smr_geo_tracer, rhs = "t_hub_mean_z", engine = "esf", C = C_matrix,
        obs_col = "cvd_obs", exp_col = "cvd_exp", refresh = 0)
    )
  ),
  tar_target(tracer_engine_table,
             collect_coefficients(tracer_engines, term_labels = TERM_LABELS)),

  tar_target(exposure_contrast_table,
             exposure_contrast(smr_geo_full, C = C_matrix)),
  tar_target(tracer_per_10min,
             mde_per_unit(control_fits$tracer, smr_geo_allage,
                          "t_hub_mean", per = 10)),
  tar_target(i63_all_ages,
             build_icd_outcome(mort_raw, pop_area_table, area_shp$area,
                               prefixes = "I63", label = "i63",
                               pop_year = STUDY_YEARS)),
  tar_target(cvd_all_ages,
             build_icd_outcome(mort_raw, pop_area_table, area_shp$area,
                               prefixes = paste0("I6", 0:9), label = "cvdall",
                               pop_year = STUDY_YEARS)),
  tar_target(haem_under75,
             build_icd_outcome(mort_raw, pop_area_table, area_shp$area,
                               prefixes = c("I60", "I61", "I62"),
                               age_max = 74, label = "haem",
                               pop_year = STUDY_YEARS)),

  # Feasibility BEFORE fitting. A BYM2 on a surface that is zero across most
  # units estimates the prior, not the data.
  tar_target(outcome_feasibility_table,
             dplyr::bind_rows(
               outcome_feasibility(i63_all_ages, "i63"),
               outcome_feasibility(cvd_all_ages, "cvdall"),
               outcome_feasibility(haem_under75, "haem")
             )),

  tar_target(
    smr_geo_allage,
    sf::st_as_sf(smr_geo_tracer) |>
      dplyr::left_join(i63_all_ages, by = "area") |>
      dplyr::left_join(cvd_all_ages, by = "area") |>
      dplyr::left_join(haem_under75, by = "area")
  ),

  tar_target(
    allage_fits,
    list(
      `I63 all ages ~ hub`    = fit_model(
        smr_geo_allage, rhs = "t_hub_mean_z", engine = "bym2",
        C = C_matrix, scale_factor = scale_factor,
        obs_col = "i63_obs", exp_col = "i63_exp", refresh = 0),
      `I63 all ages ~ centre` = fit_model(
        smr_geo_allage, rhs = "t_centre_mean_z", engine = "bym2",
        C = C_matrix, scale_factor = scale_factor,
        obs_col = "i63_obs", exp_col = "i63_exp", refresh = 0),
      `All cerebrovascular, all ages ~ hub` = fit_model(
        smr_geo_allage, rhs = "t_hub_mean_z", engine = "bym2",
        C = C_matrix, scale_factor = scale_factor,
        obs_col = "cvdall_obs", exp_col = "cvdall_exp", refresh = 0)
    )
  ),
  tar_target(allage_table,
             collect_coefficients(allage_fits, term_labels = TERM_LABELS)),
  tar_target(allage_diagnostics, collect_diagnostics(allage_fits)),

  tar_target(fit_cvd_bym2,
             fit_model(smr_geo_tracer, rhs = "1", engine = "bym2",
                       C = C_matrix, scale_factor = scale_factor,
                       obs_col = "cvd_obs", exp_col = "cvd_exp", refresh = 0)),
  tar_target(aug_cvd,
             augment_bym2(smr_geo_tracer, fit_cvd_bym2,
                          exp_col = "cvd_exp", pop_col = "person_years",
                          threshold = RR_THRESHOLD_STROKE)),
  tar_target(diag_cvd_bym2, check_bym2_fit(fit_cvd_bym2, print = FALSE)),

  tar_target(fit_i63_bym2,
             fit_model(smr_geo_allage, rhs = "1", engine = "bym2",
                       C = C_matrix, scale_factor = scale_factor,
                       obs_col = "i63_obs", exp_col = "i63_exp", refresh = 0)),
  tar_target(aug_i63,
             augment_bym2(smr_geo_allage, fit_i63_bym2,
                          exp_col = "i63_exp", pop_col = "person_years",
                          threshold = RR_THRESHOLD_STROKE)),
  tar_target(diag_i63_bym2, check_bym2_fit(fit_i63_bym2, print = FALSE)),

  tar_target(fit_haem_bym2,
             fit_model(smr_geo_allage, rhs = "1", engine = "bym2",
                       C = C_matrix, scale_factor = scale_factor,
                       obs_col = "haem_obs", exp_col = "haem_exp",
                       refresh = 0)),
  tar_target(aug_haem,
             augment_bym2(smr_geo_allage, fit_haem_bym2,
                          exp_col = "haem_exp", pop_col = "person_years",
                          threshold = RR_THRESHOLD_STROKE)),
  # The specification ladder, in increasing order of how much spatial variation
  # is removed from the exposure. Rung two is the reference: it closes the one
  # backdoor path the DAG identifies, using measured proxies for U, and leaves
  # the exposure's spatial contrast intact. Rungs three and four remove that
  # contrast along with the confounding, and are reported as a bound.
  tar_target(
    i63_engines,
    list(
      `GLM, exposure only` = fit_model(
        smr_geo_allage, rhs = STROKE_EXPOSURE, engine = "glm",
        obs_col = "i63_obs", exp_col = "i63_exp", refresh = 0),
      # Identical to the control panel's tracer arm, so it is reused rather
      # than sampled a second time.
      `GLM + adjustment set (reference)` = control_fits$tracer,
      `BYM2 + adjustment set` = fit_model(
        smr_geo_allage, rhs = STROKE_RHS, engine = "bym2", C = C_matrix,
        scale_factor = scale_factor,
        obs_col = "i63_obs", exp_col = "i63_exp", refresh = 0),
      `ESF + adjustment set` = fit_model(
        smr_geo_allage, rhs = STROKE_RHS, engine = "esf", C = C_matrix,
        obs_col = "i63_obs", exp_col = "i63_exp", refresh = 0)
    )
  ),
  tar_target(i63_engine_table,
             collect_coefficients(i63_engines, term_labels = TERM_LABELS)),

  # --- reporting summaries for the reference outcome ------------------------
  tar_target(i63_smoothing,
             smoothing_summary(smr_geo_allage, aug_i63, diag_i63_bym2,
                               obs_col = "i63_obs", exp_col = "i63_exp")),
  tar_target(i63_exceedance, exceedance_count(aug_i63)),
  tar_target(stroke_subtypes, stroke_subtype_table(mort_raw)),
  tar_target(stroke_subtypes_under75,
             stroke_subtype_table(mort_raw, age_max = 74)),

  # --- figures for the stroke sub-model -------------------------------------
  tar_target(fig_hub_time,
             plot_travel_time(
               smr_geo_full, "t_hub_mean", centres = stroke_centres,
               caption = paste("Routed on the OSM network with AREU emergency",
                               "speeds, from population-weighted",
                               "census-section origins."))),
  tar_target(fig_i63_pair,
             map_pair(
               overlay_hubs(plot_smr_map(smr_geo_allage, value = "i63_smr",
                                         title = NULL, subtitle = NULL,
                                         caption = NULL), stroke_centres),
               overlay_hubs(plot_smr_map(aug_i63, value = "bym2_rr",
                                         title = NULL, subtitle = NULL,
                                         caption = NULL), stroke_centres),
               titles = c("Raw SMR (observed / expected)",
                          "BYM2-smoothed relative risk"))),
  tar_target(fig_i63_exceedance,
             overlay_hubs(plot_exceedance_map(aug_i63, title = NULL,
                                              subtitle = NULL),
                          stroke_centres)),
  tar_target(fig_i63_raw,
             plot_tracer_raw(smr_geo_allage, "t_hub_mean",
                             obs_col = "i63_obs", exp_col = "i63_exp")),
  # Haemorrhagic stroke is NOT a clean negative control: level-II centres do
  # neurosurgery as well as thrombectomy, so the exposure has a plausible path
  # to this outcome too. It stays in the stroke sub-report as a descriptive
  # companion and is out of the thesis Results.
  tar_target(fig_haem_pair,
             map_pair(
               overlay_hubs(plot_smr_map(aug_haem, value = "bym2_rr",
                                         title = NULL, subtitle = NULL,
                                         caption = NULL), stroke_centres),
               overlay_hubs(plot_exceedance_map(aug_haem, title = NULL,
                                                subtitle = NULL),
                            stroke_centres),
               titles = c("Smoothed relative risk",
                          "Posterior P(RR > 1.20)"))),
  # --- the stroke report ----------------------------------------------------
  tar_quarto(stroke_results, path = file.path("reports", "stroke_results.qmd")),

  tar_target(JustDontCareLastComma, NULL)
)
