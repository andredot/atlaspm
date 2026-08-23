
# Project packages (TO BE UPDATED EVERY NEW PACKAGE USED) ----------
meta_pkgs <- c("tidyverse")  # e.g., tidyverse, tidymodels, ...
renv::install(meta_pkgs)

prj_pkgs <- c(
  # data manipulation
  "dplyr", "tidyr", "purrr", "tibble", "readr", "readxl",
  "stringr", "forcats", "janitor", "lubridate", "fs", "rlang", "withr",
  # spatial
  "sf", "spdep", "terra", "exactextractr", "Matrix",
  # routing (stroke sub-model)
  "dodgr", "osmextract", "tidygeocoder", "stringdist",
  # NAR primary-care pipeline
  # "arrow", "duckdb", "DBI",
  # estimation
  "geostan", "rstan", "loo",
  # output
  "ggplot2", "gtsummary", "knitr",
  # pipeline
  "targets", "tarchetypes", "crew", "quarto", "here", "cli", "methods"
)
renv::install(prj_pkgs)
purrr::walk(prj_pkgs, usethis::use_package)

dev_pkgs <- c(
  "checkmate", "covr", "devtools", "distill", "fs", "here", "htmltools",
  "knitr", "lintr", "lubridate", "purrr", "rstudioapi",
  "spelling", "stringr", "targets", "tarchetypes", "testthat",
  "usethis", "withr"
)
renv::install(dev_pkgs)
purrr::walk(dev_pkgs, usethis::use_package, type = "Suggests")

usethis::use_tidy_description()
devtools::document()
renv::status()
# renv::snapshot()

# Functions definitions -------------------------------------------

## if you need more structure respect to include your functions inside
## `R/functions.R`, you can create other couple of test/function-script
## by running the following lines of code as needed.

# usethis::use_test("<my_fun>")
# usethis::use_r(<"my_fun">)
