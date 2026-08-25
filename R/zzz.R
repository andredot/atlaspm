# Non-standard evaluation declarations ----------------------------------------
#
# R CMD check reads the package with a static analyser that has no concept of a
# data mask. Every bare column name inside dplyr::mutate(), tidyselect, or a
# formula therefore looks like a reference to an undefined global variable, and
# the check emits one NOTE per occurrence - around a hundred of them here.
#
# WHAT THIS FILE IS FOR, AND WHAT IT IS NOT FOR
#
# Declaring these names silences a false positive: the analyser is wrong, the
# code is right. That is legitimate.
#
# It is NOT a way to silence a real one. Two of the entries the check reported
# were genuine bugs, not NSE at all - `transmute`, `collect` and `year` had
# been left unqualified when the library() calls were removed from
# R/nar_density.R, so they would have failed at run time with "could not find
# function". Adding those to a globalVariables() list would have hidden a
# crash. They were fixed properly, by namespacing them.
#
# The rule: only DATA COLUMNS belong in this list. If the check reports
# "no visible global function definition", that is never an NSE artefact and
# must be fixed in the code.
#
# The better long-term fix is `.data[["col"]]`, which needs no declaration and
# fails loudly when a column is missing. Newer code in this package already
# uses it. This list covers the older functions, which are not worth churning
# just to satisfy a static analyser.

utils::globalVariables(c(

  # -- internal temporaries created inside a pipeline -------------------------
  # Prefixed with a dot to keep them clear of real columns.
  ".area", ".age", ".sex", ".pop", ".lvl", ".arm", ".w", ".row", ".sig",

  # -- mortality register (R/import_.R, R/preprocess_.R) ----------------------
  "causa_finale", "comune_residenza", "asst", "nil", "distretto",
  "eta", "anno", "sesso", "death_id", "mechanism", "weight", "type",
  "area_residenza", "cause", "group", "flag",

  # -- population and census (R/preprocess_.R, R/import_.R) -------------------
  "Codice comune", "Comune", "Eta", "numero", "area",
  "nil_label", "nil_num", "genere", "eta_label",
  "PROCOM", "PRO_COM_T", "ID_NIL",

  # -- standardisation intermediates (R/preprocess_.R) ------------------------
  "std_deaths", "std_pop", "std_rate", "observed", "expected",
  "smr", "isr", "total_smr", "total_isr", "deaths",
  "person_years", "population",

  # -- primary-care indicator (R/import_.R) -----------------------------------
  "numeratore", "denominatore", "indicatore",

  # -- NAR density reference implementation (R/nar_density.R) -----------------
  "resident_id", "gp_id", "event_date", "event_type", "death_date",
  "birth_year", "area_id", "obs_end", "spell_start", "spell_end",
  "list_size", "is_artifact", "l_clamped", "clamp_status",
  "n_patients", "n_alive", "n_gps", "supply", "supply_sum",
  "assigned_pm", "unassigned_pm", "artifact_pm", "gp65_pm",
  "list_pm_sum", "person_months", "prop_unassigned",
  "mean_list_size_assigned", "gp_65plus", "gp_density", "month",

  # -- targets referenced inside tar_map() specifications ---------------------
  # These are supplied by tarchetypes at pipeline build time, not by this
  # package, but roxygen sees the _targets.R symbols when the file is sourced.
  "rhs", "engine", "id"
))
