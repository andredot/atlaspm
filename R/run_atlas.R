#' Launch the exploratory atlas
#'
#' Opens the Shiny front end over the pipeline's outputs. The app reads the
#' bundle written by the `app_data` target and never touches a stanfit, so it
#' starts in about a second and holds a few megabytes rather than a few hundred.
#'
#' The bundle is a pipeline product, not a side script: `targets::tar_make()`
#' builds it and rebuilds it whenever a fit, a covariate or the geography
#' changes. If the target is out of date this function says so rather than
#' opening an app that renders the previous run's results perfectly.
#'
#' @param data_path Path to the bundle. When `NULL` (default) it is resolved
#'   from the `app_data` target, falling back to the copy sitting next to the
#'   app.
#' @param check Whether to ask `targets` if the bundle is out of date before
#'   launching. Costs a graph walk; set `FALSE` to skip it.
#' @param ... Passed to [shiny::runApp()], e.g. `port`, `launch.browser`.
#'
#' @return Invisibly, whatever [shiny::runApp()] returns; called for the app.
#'
#' @examples
#' \dontrun{
#' targets::tar_make()
#' run_atlas()
#' run_atlas(launch.browser = FALSE, port = 4321)
#' }
#' @seealso [build_app_data()]
#' @export
run_atlas <- function(data_path = NULL, check = TRUE, ...) {

  for (pkg in c("shiny", "leaflet", "htmltools", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is needed to run the atlas app. ",
           "install.packages(\"", pkg, "\")", call. = FALSE)
    }
  }

  app_dir <- system.file("shiny", "atlas", package = "atlaspm")
  if (!nzchar(app_dir)) {
    # devtools::load_all() does not install, so fall back to the source tree.
    # Deliberately not R/: _targets.R calls tar_source(), which evaluates every
    # .R file under R/ as part of loading the pipeline, and an app has no
    # business being loaded by the pipeline.
    candidates <- file.path(c("inst/shiny/atlas", "shiny/atlas",
                              "inst/shiny", "shiny"), "app.R")
    hit <- candidates[file.exists(candidates)]
    app_dir <- if (length(hit)) dirname(hit[1]) else ""
  }
  if (!nzchar(app_dir) || !file.exists(file.path(app_dir, "app.R"))) {
    stop("Cannot find app.R. Expected it at inst/shiny/atlas/app.R.\n",
         "  Note it must NOT live under R/: tar_source() would evaluate it ",
         "every time the pipeline loads.", call. = FALSE)
  }

  # Resolve the bundle from the pipeline first. tar_read() on a format = "file"
  # target returns the path it wrote, which is the authoritative answer; the
  # path beside the app is only a fallback for a store that is not there.
  if (is.null(data_path) && requireNamespace("targets", quietly = TRUE)) {
    data_path <- tryCatch(
      as.character(targets::tar_read_raw("app_data"))[1],
      error = function(e) NULL
    )
  }
  if (is.null(data_path) || !nzchar(data_path)) {
    data_path <- file.path(app_dir, "app_data.rds")
  }
  if (!file.exists(data_path)) {
    stop("No app data at '", data_path, "'.\n",
         "  The bundle is built by the `app_data` target:\n",
         "    targets::tar_make(app_data)",
         call. = FALSE)
  }

  # A bundle older than the fits renders perfectly and shows the previous run.
  # targets already knows whether that is so, so ask it rather than comparing
  # file mtimes and guessing.
  if (isTRUE(check) && requireNamespace("targets", quietly = TRUE)) {
    stale <- tryCatch(
      "app_data" %in% targets::tar_outdated(names = "app_data",
                                            reporter = "silent"),
      error = function(e) FALSE
    )
    if (isTRUE(stale)) {
      warning("The `app_data` target is out of date: the app will show the ",
              "previous run's results. Run targets::tar_make(app_data).",
              call. = FALSE)
    }
  }

  withr::with_envvar(
    c(ATLASPM_APP_DATA = normalizePath(data_path)),
    shiny::runApp(app_dir, ...)
  )
}
