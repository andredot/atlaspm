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
  # fast data formats
  format = "qs",
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

  # Import your file from custom (shared) location, and preprocess them
  # tar_files_input(
  #   db_raw_path,
  #   # get_input_data_path("db_raw.csv"), # single
  #   get_input_data_path() |>
  #     list.files(pattern = "\\.csv$", full.names = TRUE) # multiple
  # ),

  # LOOK UP TABLES
  tar_target(lookup_causes, get_input_data_path("avoidable_lookup.csv")),

  # IMPORT

  tar_target(mort_path, get_input_data_path("mort_2023.csv")),
  tar_target(mort_raw, import_mortality(mort_path)),

  # PREPROCESSING

  tar_target(mort_count, preprocess_mortality(mort_raw,
                                              lookup_causes,
                                              code_col = "causa",
                                              age_col = "eta")),


  # # Call your custom functions as needed.
  # tar_target(relevantResult, relevant_computation(db)),
  #
  # # compile yor report
  tar_quarto(explore_mort_count, path = "reports\\mortality_explore.qmd"),
  #
  #
  # # Decide what to share with other, and do it in a standard RDS format
  # tar_target(
  #   objectToShare,
  #   list(
  #     relevant_result = relevantResult
  #   )
  # ),
  # tar_target(
  #   shareOutput,
  #   share_objects(objectToShare),
  #   format = "file",
  #   pattern = map(objectToShare)
  # ),


  tar_target(JustDontCareLastComma, NULL)
)
