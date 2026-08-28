# Numbered tables and figures for Draft_Atlas_v3 -------------------------------
#
# Like _thesis_setup.R, this is presentation rather than analysis: it decides
# which pipeline object becomes "Figure 4" in the thesis and at what size it is
# written, not what any number is. Nothing here computes anything.
#
# The one rule: the mapping from thesis number to pipeline target lives in the
# manifests below and nowhere else, so renumbering a figure is a one-line edit
# and the DOCX, the PNG filenames and the caption list can never disagree.

# --- internal formatting ------------------------------------------------------
#
# Deliberately a local copy of fmt_n() from _thesis_setup.R rather than a call
# to it. This file has to work from either R/ (sourced by tar_source, before
# any report helper exists) or reports/ (sourced after them), and a function
# whose enclosing environment lacks fmt_n fails at call time with nothing to
# suggest why. Same output, no load-order dependency.

.fmt_int <- function(x) {
  formatC(round(as.numeric(x), 0), format = "f", digits = 0, big.mark = ",")
}

# --- manifests ----------------------------------------------------------------
#
# `target` is the name of the tar_load()ed object. Multi-panel figures carry
# several targets separated by "|" and one `panel` label per target.
#
# Sizes are millimetres, matching save_figure(). 190 is full text width.

fig_manifest <- tibble::tribble(
  ~number, ~stem,                  ~target,                            ~panel,                                  ~width, ~height, ~caption,

  1L, "fig01_flow",           "flow",                                  NA_character_,                            160, 210,
  paste("Derivation of the analysis dataset, main analysis and stroke",
        "sub-analysis. The stroke arm branches from the full extract rather",
        "than from the age-restricted population, because the tracer outcome",
        "is defined for all ages."),

  2L, "fig02_pollution_pair", "fig_pollution_shared",                  NA_character_,                            190, 110,
  paste("Nitrogen dioxide and PM2.5 across the study area, on a shared",
        "scale. The shared legend is the honest version of the selection",
        "argument: PM2.5 varies far less across the territory than NO2 does",
        "and therefore carries less contrast for a model to use."),

  3L, "fig03_excess",         "fig_excess_map|fig_excess_flagged",     "All areal units|Areas above the 0.80 cutoff", 190, 150,
  paste("Residual excess avoidable deaths per year under M5. Panel A shows",
        "all 279 areal units; panel B restricts to areas exceeding the 0.80",
        "posterior probability cutoff."),

  4L, "fig04_rr_m3",          "fig_rr_map|fig_exceedance",             "Smoothed relative risk|Exceedance probability (RR > 1.10)", 190, 150,
  paste("Smoothed relative risk and exceedance probability, M3 (no",
        "covariates). Panel A is the BYM2-smoothed relative risk surface;",
        "panel B is the posterior probability that an area's relative risk",
        "exceeds 1.10."),

  5L, "fig05_forest",         "fig_forest",                            NA_character_,                            190, 200,
  "Covariate associations across specifications, relative risk per standard deviation.",

  6L, "fig06_mech_rr",        "fig_mech_facets",                       NA_character_,                            190, 200,
  "BYM2-smoothed relative risk by service function.",

  7L, "fig07_mech_exceed",    "fig_mech_exceedance",                   NA_character_,                            190, 200,
  "Posterior probability that relative risk exceeds 1.10, by service function.",

  8L, "fig08_hub_time",       "fig_hub_time",                          NA_character_,                            160, 150,
  "Population-weighted travel time to the nearest level-II thrombectomy hub.",

  9L, "fig09_i63_exceedance", "fig_i63_exceedance",                    NA_character_,                            160, 150,
  "Posterior probability that all-age I63 relative risk exceeds 1.20.",

  10L, "fig10_i63_raw",       "fig_i63_raw",                           NA_character_,                            160, 130,
  paste("All-age I63 standardised mortality ratio against travel time to the",
        "nearest thrombectomy hub. Point size is the expected count, so the",
        "large points are the areas the model learns from."),

  # --- suggested additions, numbered after the ten the thesis already has ---

  11L, "fig11_rr_compare",    "fig_rr_compare",                        NA_character_,                            160, 140,
  paste("BYM2 (M5) against ESF (M6) smoothed relative risk. Supports the",
        "Spearman rank correlation quoted in Sensitivity to model",
        "specifications, which currently has no figure."),

  12L, "fig12_concordance",   "fig_concordance",                       NA_character_,                            160, 140,
  paste("Rank concordance among the seven mechanism risk surfaces. Supports",
        "the correlation between the lifestyle/NCD and tertiary prevention",
        "surfaces that the Discussion turns on."),

  13L, "fig13_controls",      "fig_controls",                          NA_character_,                            190, 130,
  paste("Tracer and control coefficients, relative risk per standard",
        "deviation. Makes the control panel legible at a glance.")
)

# Figures the pipeline builds that no thesis slot currently uses. Listed so the
# report can print them as an inventory rather than leaving them invisible.
fig_unused <- c(
  "fig_cmr_map", "fig_smr_map", "scatter_smr_di", "scatter_smr_gp",
  "fig_gp_density_map", "fig_pollution_pair", "fig_i63_pair", "fig_haem_pair",
  "scatter_cmr_isr_overall"
)

# --- figure export ------------------------------------------------------------

#' Write every manifest figure to `dir` as a 300 dpi PNG.
#'
#' Uses [save_figure()] so the thesis PNGs are produced by the same code path
#' as every other figure in the project. Multi-panel entries are written as
#' `<stem>a`, `<stem>b`, ... because Word handles two placed images better than
#' one wide composite.
#'
#' @param env Environment holding the tar_load()ed figure objects.
#' @param dir Output directory.
#' @param manifest Defaults to [fig_manifest].
#' @return A tibble of number, panel, target and path, invisibly.
export_thesis_figures <- function(env = parent.frame(),
                                  dir = here::here("output", "thesis"),
                                  manifest = fig_manifest) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- c("a", "b", "c", "d")

  out <- lapply(seq_len(nrow(manifest)), function(i) {
    row     <- manifest[i, ]
    targets <- strsplit(row[["target"]], "|", fixed = TRUE)[[1]]
    panels  <- if (is.na(row[["panel"]])) NA_character_ else
      strsplit(row[["panel"]], "|", fixed = TRUE)[[1]]

    lapply(seq_along(targets), function(j) {
      nm <- targets[[j]]

      # Figure 1 is drawn here rather than in the pipeline, because `flow` is a
      # table target and no plot target consumes it.
      p <- if (nm == "flow") {
        if (exists("flow", envir = env)) {
          plot_flow_diagram(get("flow", envir = env),
                            stroke = if (exists("stroke_subtypes", envir = env))
                              get("stroke_subtypes", envir = env) else NULL)
        } else {
          NULL
        }
      } else if (exists(nm, envir = env)) {
        get(nm, envir = env)
      } else {
        NULL
      }

      if (is.null(p)) {
        cli::cli_warn("Target {.val {nm}} not loaded; figure {row$number} skipped.")
        return(NULL)
      }

      stem <- if (length(targets) > 1L)
        paste0(row[["stem"]], suffix[[j]]) else row[["stem"]]

      path <- save_figure(p, stem, dir = dir,
                          width = row[["width"]], height = row[["height"]])

      tibble::tibble(number = row[["number"]],
                     panel  = if (length(panels) >= j) panels[[j]] else NA_character_,
                     target = nm, path = path)
    })
  })

  res <- dplyr::bind_rows(unlist(out, recursive = FALSE))
  cli::cli_alert_success("Wrote {nrow(res)} figure file{?s} to {.path {dir}}.")
  invisible(res)
}

# --- Figure 1: analysis flow --------------------------------------------------

#' STROBE-style derivation diagram, drawn from the `flow` target.
#'
#' The counts are read from `flow` rather than typed, so the diagram cannot
#' drift from the pipeline. The stroke arm is drawn only when `stroke` is
#' supplied, and its I63 count is taken from the subtype table so that the
#' exclusions shown actually sum to the total.
#'
#' @param flow The `flow` target, from [flow_counts()]: columns `step`, `n`,
#'   `lost`.
#' @param stroke Optionally the `stroke_subtypes` target: columns `subtype`,
#'   `deaths`.
#' @return A ggplot.
plot_flow_diagram <- function(flow, stroke = NULL) {

  # flow_counts() returns step/n/lost. Validating here turns a wrong column
  # name into a message naming the columns that do exist, rather than a
  # zero-length sprintf() and a tibble recycling error several frames down.
  require_cols(flow, c("step", "n", "lost"), "flow")
  if (!is.null(stroke) && nrow(stroke)) {
    require_cols(stroke, c("subtype", "deaths"), "stroke")
  }

  # `arm` separates the two chains. An arrow is drawn from a box to the one
  # below it only within an arm: the stroke arm branches from the extract, not
  # from the main analysis, so a connector between them would assert a
  # derivation that did not happen.
  main <- tibble::tibble(
    y     = seq(from = 0, by = -1, length.out = nrow(flow)),
    arm   = "main",
    label = sprintf("%s\nn = %s", flow[["step"]], .fmt_int(flow[["n"]]))
  )
  drops <- which(!is.na(flow[["lost"]]))
  excl <- tibble::tibble(
    y     = main[["y"]][drops] + 0.5,
    label = sprintf("Excluded at this stage\nn = %s", .fmt_int(flow[["lost"]][drops]))
  )

  if (!is.null(stroke) && nrow(stroke)) {
    tot  <- sum(stroke[["deaths"]], na.rm = TRUE)
    i63  <- sum(stroke[["deaths"]][grepl("^I63", stroke[["subtype"]])], na.rm = TRUE)
    base <- min(main[["y"]]) - 1.8
    main <- dplyr::bind_rows(main, tibble::tibble(
      y = c(base, base - 1),
      arm = "stroke",
      label = c(sprintf("Cerebrovascular deaths, all ages (I60-I69)\nn = %s",
                        .fmt_int(tot)),
                sprintf("TRACER ANALYSIS\nI63 cerebral infarction, all ages\nn = %s",
                        .fmt_int(i63)))))
    excl <- dplyr::bind_rows(excl, tibble::tibble(
      y = base - 0.5,
      label = sprintf("Excluded: not cerebral infarction\nn = %s",
                      .fmt_int(tot - i63))))
  }

  # Within-arm connectors only.
  nxt  <- c(main[["arm"]][-1], NA_character_)
  link <- main[!is.na(nxt) & nxt == main[["arm"]], , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = link,
      ggplot2::aes(x = 0, xend = 0, y = .data$y - 0.18, yend = .data$y - 0.82),
      linewidth = 0.3,
      arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm"), type = "closed")) +
    ggplot2::geom_segment(
      data = excl,
      ggplot2::aes(x = 0, xend = 0.55, y = .data$y, yend = .data$y),
      linewidth = 0.3,
      arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm"), type = "closed")) +
    ggplot2::geom_label(data = main,
                        ggplot2::aes(x = 0, y = .data$y, label = .data$label),
                        size = 2.9, label.padding = ggplot2::unit(2, "mm"),
                        label.r = ggplot2::unit(1, "mm"), lineheight = 1.1) +
    ggplot2::geom_label(data = excl,
                        ggplot2::aes(x = 0.60, y = .data$y, label = .data$label),
                        size = 2.6, hjust = 0, colour = "grey30",
                        label.padding = ggplot2::unit(1.6, "mm"),
                        label.r = ggplot2::unit(1, "mm"), lineheight = 1.1) +
    ggplot2::scale_x_continuous(limits = c(-0.75, 1.65)) +
    ggplot2::theme_void()
}

# --- table helper -------------------------------------------------------------

#' One numbered thesis table.
#'
#' Wraps [knitr::kable()] so every table in the assets document is titled the
#' same way, and so the source target is recorded in the caption. `x` may be a
#' data frame or a `gtsummary` object.
#'
#' @param x Data frame or gtsummary table.
#' @param number Thesis table number, e.g. `"3b"`.
#' @param title Caption text, without the "Table N." prefix.
#' @param source Name of the pipeline target the table comes from.
#' @param col_names Passed to [knitr::kable()].
#' @return A kable.
thesis_table <- function(x, number, title, source, col_names = NULL) {
  cap <- sprintf("Table %s. %s (Source: %s.)", number, title, source)
  if (inherits(x, "gtsummary")) {
    return(gtsummary::as_kable(x, caption = cap))
  }
  if (is.null(col_names)) {
    knitr::kable(x, caption = cap)
  } else {
    knitr::kable(x, col.names = col_names, caption = cap)
  }
}
