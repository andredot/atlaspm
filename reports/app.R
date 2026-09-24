# Atlas of Preventable Mortality - exploratory Shiny app -----------------------
#
# Reads one pre-built bundle (dev/build_app_data.R) and nothing else. No Stan
# object is loaded here and nothing is refitted: everything on screen was
# computed by the pipeline or by the build script.
#
# TWO CONVENTIONS THAT RUN THROUGH THE WHOLE FILE
#
# 1. Base maps are rendered once, empty, and repainted with leafletProxy. If a
#    selector rebuilt the whole map the view would snap back to the full extent
#    on every change, which makes the search box useless: you would find an
#    area, change the cause, and lose it.
#
# 2. Where a quantity does not exist for the current selection - a smoothed
#    surface for a cause group, an ICAR/IID split for the ESF model - the map
#    is painted grey and says why, in place, rather than silently showing the
#    previous selection or a blank panel.

library(shiny)
library(leaflet)
library(sf)
library(dplyr)
library(ggplot2)
library(htmltools)

`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) y else x
}

# ---- data ----------------------------------------------------------------------
#
# Loaded on first use, never when this file is sourced.
#
# `_targets.R` calls tar_source(), which recursively evaluates every .R file
# under R/, and covr, lintr and a stray source() all walk the project the same
# way. A script that reads a file and stop()s at the top level turns every one
# of those into a pipeline failure that has nothing to do with the pipeline.
# Nothing below touches the disk until Shiny serves a request.

.atlas_cache <- new.env(parent = emptyenv())

#' Walk up from `start` looking for the project root.
#'
#' The app can be launched from its own directory (shiny::runApp sets the
#' working directory to the app folder), from the project root, or from
#' anywhere via run_atlas(). The bundle lives at a path relative to the root,
#' so the root has to be found rather than assumed.
.project_root <- function(start = getwd(), max_up = 6L) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(max_up)) {
    if (any(file.exists(file.path(d, c("_targets.R", "DESCRIPTION"))))) return(d)
    up <- dirname(d)
    if (identical(up, d)) break
    d <- up
  }
  NULL
}

#' Ask the pipeline where it put the bundle.
#'
#' `app_data` is a format = "file" target, so tar_read() returns the path it
#' wrote. That is the authoritative answer and it survives someone changing the
#' `path` argument in _targets.R, which guessing at fixed locations does not.
.bundle_from_targets <- function(root) {
  if (is.null(root) || !requireNamespace("targets", quietly = TRUE)) return(NULL)
  if (!dir.exists(file.path(root, "_targets"))) return(NULL)
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  p <- tryCatch(as.character(targets::tar_read_raw("app_data"))[1],
                error = function(e) NULL)
  if (is.null(p) || is.na(p) || !nzchar(p)) return(NULL)
  if (!file.exists(p)) p <- file.path(root, p)
  if (file.exists(p)) normalizePath(p, winslash = "/") else NULL
}

#' Every place the bundle could reasonably be, in order of authority.
.bundle_candidates <- function(root = .project_root()) {
  rel <- c("app_data.rds",
           file.path("inst", "shiny", "atlas", "app_data.rds"),
           file.path("reports", "app_data.rds"),
           file.path("output", "app_data.rds"))
  paths <- file.path(getwd(), rel)
  if (!is.null(root)) paths <- c(paths, file.path(root, rel))
  unique(normalizePath(paths, winslash = "/", mustWork = FALSE))
}

atlas_data <- function(path = NULL) {

  if (!is.null(.atlas_cache$d)) return(.atlas_cache$d)

  root  <- .project_root()
  tried <- character(0)

  if (is.null(path)) {
    env <- Sys.getenv("ATLASPM_APP_DATA", unset = "")
    if (nzchar(env)) {
      tried <- c(tried, paste0(env, "   (ATLASPM_APP_DATA)"))
      if (file.exists(env)) path <- env
    }
  }
  if (is.null(path)) {
    p <- .bundle_from_targets(root)
    if (!is.null(p)) path <- p else
      tried <- c(tried, "the `app_data` target   (not built, or no _targets/ store here)")
  }
  if (is.null(path)) {
    cand <- .bundle_candidates(root)
    hit  <- cand[file.exists(cand)]
    tried <- c(tried, cand)
    if (length(hit)) path <- hit[1]
  }

  if (is.null(path)) {
    root_msg <- if (is.null(root)) {
      paste0("not found (no _targets.R or DESCRIPTION at or above ",
             getwd(), ")")
    } else root
    stop("No app data bundle found.\n",
         "  Project root: ", root_msg, "\n",
         "  Looked in:\n    ", paste(tried, collapse = "\n    "), "\n",
         "  Build it with:  targets::tar_make(app_data)\n",
         "  Then launch with run_atlas(), which sets the path for you, ",
         "or set ATLASPM_APP_DATA.",
         call. = FALSE)
  }

  A     <- readRDS(path)
  GEO   <- A$geo
  TAB   <- sf::st_drop_geometry(GEO)
  AREAS <- as.character(TAB$area)
  LABS  <- as.character(TAB$label)

  # Per-area bounding boxes, computed once: the search box needs them on every
  # keystroke, and st_bbox on a polygon is not free at 279 areas x 20 maps.
  BBOX        <- lapply(sf::st_geometry(GEO), sf::st_bbox)
  names(BBOX) <- AREAS

  .atlas_cache$d <- list(
    A = A, GEO = GEO, TAB = TAB, AREAS = AREAS, LABS = LABS,
    BBOX = BBOX, META = A$meta,
    AREA_CHOICES = stats::setNames(AREAS, LABS)[order(LABS)]
  )
  .atlas_cache$d
}

# ---- palettes ------------------------------------------------------------------
#
# The SMR breaks and colours are the ones plot_smr_map() uses, so a map here and
# a figure in the thesis mean the same thing by the same colour.

SMR_BREAKS <- c(0.8, 0.9, 1.1, 1.2)
SMR_COLS   <- c("#5a8a7d", "#a3c4b5", "#f2e8d5", "#dca678", "#b5651d")

EXC_BREAKS <- c(-1e-9, 0.5, 0.8, 0.95, 1 + 1e-9)
EXC_COLS   <- c("#EAEAEA", "#FDB863", "#D94801", "#8C2D04")
EXC_LABS   <- c("Unlikely (< 50%)", "Possible (50-80%)",
                "Likely (80-95%)", "Almost certain (>= 95%)")

DIV_COLS <- c("#2c7bb6", "#abd9e9", "#f7f7f7", "#fdae61", "#d7191c")
SEQ_COLS <- c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45", "#00441b")

#' Finite bins around a set of interior cut points.
#'
#' leaflet::colorBin will not take Inf, so the outer edges are pushed just past
#' the data. Done here rather than inline so that an empty or all-NA layer -
#' which happens whenever a selection has no surface - cannot produce
#' non-monotonic bins and a cryptic error from colorBin.
finite_bins <- function(x, cuts) {
  rng <- suppressWarnings(range(x, na.rm = TRUE))
  if (!all(is.finite(rng))) rng <- c(min(cuts) - 1, max(cuts) + 1)
  c(min(rng[1], min(cuts)) - 1e-9, cuts, max(rng[2], max(cuts)) + 1e-9)
}

pal_smr <- function(x) {
  colorBin(SMR_COLS, domain = x, bins = finite_bins(x, SMR_BREAKS),
           na.color = "#d9d9d9")
}
pal_exceed <- function(x) {
  colorBin(EXC_COLS, domain = c(0, 1), bins = EXC_BREAKS,
           na.color = "#d9d9d9")
}
pal_seq <- function(x) {
  colorNumeric(SEQ_COLS, domain = x, na.color = "#d9d9d9")
}
#' Diverging palette, forced symmetric about `centre`.
#'
#' A diverging scale whose midpoint drifts with the data is a lie: the pale
#' middle band stops meaning "no excess" and starts meaning "average excess".
pal_div <- function(x, centre = 0) {
  m <- suppressWarnings(max(abs(x - centre), na.rm = TRUE))
  if (!is.finite(m) || m == 0) m <- 1
  colorNumeric(DIV_COLS, domain = c(centre - m, centre + m),
               na.color = "#d9d9d9")
}

fmt <- function(x, d = 2) ifelse(is.na(x), "\u2014", formatC(x, format = "f", digits = d))

# ---- map plumbing --------------------------------------------------------------

#' An empty base map: tiles, extent, and optionally the stroke centres.
#' The CARTO key, if one is configured.
#'
#' Since late August 2026 CARTO serves watermarked tiles ("API KEY REQUIRED")
#' to any request without a key. The key is free - request it at
#' https://carto.com/basemaps/apikey/ - and is read from the environment so it
#' never lands in the repository. Put it in your USER-level ~/.Renviron:
#'
#'   CARTO_API_KEY=your_key_here
#'
#' and restart R. Not the project .Renviron: that file is committed.
.carto_key <- function() {
  k <- Sys.getenv("CARTO_API_KEY", unset = "")
  if (nzchar(k)) k else NULL
}

CARTO_ATTR <- paste0(
  '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> ',
  'contributors &copy; <a href="https://carto.com/attributions">CARTO</a>'
)

#' Add basemaps, with a keyless fallback.
#'
#' With a key: CARTO Positron, which is the quiet grey a choropleth wants.
#' Without one: Esri World Gray Canvas, which needs no key and is visually the
#' closest substitute, rather than serving watermarked CARTO tiles. OSM is
#' always offered as a second option because it is the provider least likely
#' to change terms underneath the app.
#'
#' The watermarked tile comes back as an ordinary HTTP 200 image, so leaflet
#' cannot detect it and fall back by itself; the choice has to be made here.
add_basemaps <- function(m) {
  key <- .carto_key()

  if (!is.null(key)) {
    carto <- function(style) {
      paste0("https://{s}.basemaps.cartocdn.com/", style,
             "/{z}/{x}/{y}{r}.png?key=", utils::URLencode(key, reserved = TRUE))
    }
    opts <- tileOptions(subdomains = "abcd", maxZoom = 20)
    m <- m |>
      addTiles(urlTemplate = carto("light_nolabels"), attribution = CARTO_ATTR,
               group = "Light", options = opts) |>
      addTiles(urlTemplate = carto("light_all"), attribution = CARTO_ATTR,
               group = "Light + labels", options = opts)
    groups <- c("Light", "Light + labels")
  } else {
    m <- addProviderTiles(m, "Esri.WorldGrayCanvas", group = "Grey")
    groups <- "Grey"
  }

  m <- addTiles(m, group = "OpenStreetMap")
  addLayersControl(m, baseGroups = c(groups, "OpenStreetMap"),
                   options = layersControlOptions(collapsed = TRUE))
}

base_map <- function(centres = NULL) {
  GEO <- atlas_data()$GEO
  bb  <- st_bbox(GEO)
  m <- leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
    add_basemaps() |>
    fitBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])

  if (!is.null(centres)) {
    hub   <- centres[centres$is_hub %in% TRUE, , drop = FALSE]
    spoke <- centres[!(centres$is_hub %in% TRUE), , drop = FALSE]
    if (nrow(spoke)) {
      m <- addCircleMarkers(
        m, lng = spoke$lon, lat = spoke$lat, radius = 5,
        color = "#333333", weight = 1.5, fillColor = "#ffffff",
        fillOpacity = 0.95, group = "Stroke centres",
        label = lapply(paste0("<b>", spoke$facility, "</b><br/>Stroke unit (SU I)"), HTML)
      )
    }
    if (nrow(hub)) {
      m <- addCircleMarkers(
        m, lng = hub$lon, lat = hub$lat, radius = 8,
        color = "#111111", weight = 2, fillColor = "#8C2D04",
        fillOpacity = 0.95, group = "Stroke centres",
        label = lapply(paste0("<b>", hub$facility,
                              "</b><br/>Thrombectomy hub (SU II)"), HTML)
      )
    }
  }
  m
}

#' Repaint a map's polygon layer.
#'
#' `values` is a numeric (or factor) vector in the row order of GEO. `pal` is a
#' leaflet palette function; when `values` is a factor a categorical palette is
#' built here. `note` replaces the whole layer with a grey surface and an
#' explanation, which is what an unavailable selection gets.
paint <- function(id, values, pal = NULL, legend_title = "", digits = 2,
                  note = NULL, extra = NULL) {

  D <- atlas_data()
  GEO <- D$GEO; AREAS <- D$AREAS; LABS <- D$LABS

  p <- leafletProxy(id, data = GEO)
  p <- clearGroup(p, "poly")
  # NOT clearControls(): that would also remove the basemap switcher, which
  # this function does not own. Remove only the two controls it adds.
  p <- removeControl(p, "atlas_legend")
  p <- removeControl(p, "atlas_note")

  if (!is.null(note)) {
    p <- addPolygons(
      p, group = "poly", layerId = AREAS,
      fillColor = "#e8e8e8", fillOpacity = 0.7,
      color = "white", weight = 0.6,
      label = lapply(LABS, HTML)
    )
    p <- addControl(p, layerId = "atlas_note", HTML(paste0(
      "<div style='background:#fff;padding:8px 10px;max-width:260px;",
      "border-left:3px solid #b5651d;font-size:12px;line-height:1.35'>", note,
      "</div>")), position = "topright")
    return(invisible(p))
  }

  is_factor <- is.factor(values) || is.character(values)
  if (is_factor) {
    values <- as.factor(values)
    pal <- colorFactor("Set2", domain = levels(values), na.color = "#d9d9d9")
    lab_txt <- paste0("<b>", LABS, "</b><br/>", legend_title, ": ",
                      ifelse(is.na(values), "\u2014", as.character(values)))
  } else {
    if (is.null(pal)) pal <- pal_seq(values)
    lab_txt <- paste0("<b>", LABS, "</b><br/>", legend_title, ": ",
                      fmt(values, digits))
  }
  if (!is.null(extra)) lab_txt <- paste0(lab_txt, "<br/>", extra)

  p <- addPolygons(
    p, group = "poly", layerId = AREAS,
    fillColor   = pal(values),
    fillOpacity = 0.85,
    color = "white", weight = 0.6, smoothFactor = 0.3,
    label = lapply(lab_txt, HTML),
    highlightOptions = highlightOptions(weight = 3, color = "#222222",
                                        bringToFront = TRUE)
  )

  if (identical(legend_title, "Probability")) {
    p <- addLegend(p, position = "bottomright", layerId = "atlas_legend",
                   colors = EXC_COLS, labels = EXC_LABS,
                   title = legend_title, opacity = 0.85)
  } else {
    p <- addLegend(p, position = "bottomright", layerId = "atlas_legend",
                   pal = pal, values = values,
                   title = legend_title, opacity = 0.85,
                   na.label = "no data")
  }
  invisible(p)
}

#' Centre every map on one area and outline it.
focus <- function(ids, area) {
  D <- atlas_data()
  BBOX <- D$BBOX; GEO <- D$GEO; AREAS <- D$AREAS
  if (is.null(area) || !nzchar(area) || is.null(BBOX[[area]])) return(invisible())
  bb <- BBOX[[area]]
  g  <- GEO[AREAS == area, ]
  for (id in ids) {
    try({
      leafletProxy(id) |>
        clearGroup("found") |>
        addPolygons(data = g, group = "found", fill = FALSE,
                    color = "#111111", weight = 3, opacity = 1) |>
        flyToBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
    }, silent = TRUE)
  }
  invisible()
}

MAPS <- list(
  covariates = c("cov_dep", "cov_ivsm", "cov_gpw", "cov_gpu",
                 "cov_no2", "cov_pm25"),
  outcome    = c("out_cmr", "out_smr", "out_rr", "out_exc"),
  models     = c("mod_rr", "mod_resid", "mod_count", "mod_decomp"),
  stroke     = c("str_hub", "str_centre", "str_unadj_rr", "str_unadj_exc",
                 "str_adj_rr", "str_adj_exc")
)

mapbox <- function(id, title, height = 330) {
  div(class = "atlas-mapbox",
      div(class = "atlas-maptitle", title),
      leafletOutput(id, height = height))
}

# ---- UI --------------------------------------------------------------------------

CSS <- HTML("
  .atlas-maptitle { font-weight: 600; font-size: 13px; margin: 2px 0 4px 2px;
                    color: #333; }
  .atlas-mapbox   { margin-bottom: 14px; }
  .atlas-colhead  { font-weight: 700; font-size: 13px; text-transform: uppercase;
                    letter-spacing: .04em; color: #555; border-bottom: 2px solid #ddd;
                    padding-bottom: 4px; margin-bottom: 10px; }
  .atlas-note     { font-size: 12px; color: #666; line-height: 1.4;
                    background:#fafafa; border-left:3px solid #ccc;
                    padding:8px 10px; margin-bottom:12px; }
  .navbar { margin-bottom: 12px; }
")

build_ui <- function() {

  D <- atlas_data()
  GEO <- D$GEO; META <- D$META; AREA_CHOICES <- D$AREA_CHOICES

  navbarPage(
  title = "Atlas of Preventable Mortality",
  id = "tab",
  header = tagList(
    tags$head(tags$style(CSS)),
    fluidRow(
      column(5, selectizeInput(
        "find", NULL, choices = c("Find an area..." = "", AREA_CHOICES),
        width = "100%",
        options = list(placeholder = "Find an area (type a name)"))),
      column(7, div(style = "padding-top:7px;font-size:12px;color:#666;",
                    textOutput("selected_txt", inline = TRUE)))
    )
  ),

  # ---- 1. covariates ---------------------------------------------------------
  tabPanel(
    "Covariates",
    fluidRow(
      column(3,
        div(class = "atlas-note",
            "Each row is one measurement choice. The left column is the ",
            "variable carried into M5; the right is the alternative that was ",
            "not selected. On the z-score view the two maps in a row share ",
            "one scale, so they can be compared directly; on the raw view ",
            "they cannot, because the units differ."),
        radioButtons("cov_scale", "Scale",
                     c("Raw units" = "raw", "z-score (comparable in a row)" = "z"),
                     selected = "z")
      ),
      column(9,
        fluidRow(
          column(6, div(class = "atlas-colhead", "Selected (in M5)")),
          column(6, div(class = "atlas-colhead", "Not chosen"))
        ),
        fluidRow(
          column(6, mapbox("cov_dep",  "Italian Deprivation Index")),
          column(6, mapbox("cov_ivsm", "IVSM (regional vulnerability index)"))
        ),
        fluidRow(
          column(6, mapbox("cov_gpw", "Primary care density, weighted")),
          column(6, mapbox("cov_gpu", "Primary care density, unweighted"))
        ),
        fluidRow(
          column(6, mapbox("cov_no2",  "NO2")),
          column(6, mapbox("cov_pm25", "PM2.5"))
        )
      )
    )
  ),

  # ---- 2. outcome ------------------------------------------------------------
  tabPanel(
    "Avoidable mortality",
    fluidRow(
      column(3,
        radioButtons("out_class", "Classification",
                     c("All avoidable causes" = "total",
                       "OECD cause group"     = "group",
                       "Specific cause"       = "cause",
                       "Service function"     = "mechanism")),
        uiOutput("out_level_ui"),
        div(class = "atlas-note",
            "The crude rate and the SMR are available for every level. The ",
            "smoothed surface and the exceedance probability come from a ",
            "fitted model, and models were fitted for all avoidable causes ",
            "(M3) and for the seven service functions only. Cause groups and ",
            "individual causes have no smoothed surface and the lower two ",
            "panels say so rather than showing something else."),
        uiOutput("out_summary")
      ),
      column(9,
        fluidRow(
          column(6, mapbox("out_cmr", "Crude mortality rate (per 100,000 py)")),
          column(6, mapbox("out_smr", "SMR (observed / expected)"))
        ),
        fluidRow(
          column(6, mapbox("out_rr",  "Smoothed relative risk (BYM2)")),
          column(6, mapbox("out_exc", "Exceedance probability"))
        )
      )
    )
  ),

  # ---- 3. adjusted models -----------------------------------------------------
  tabPanel(
    "Adjusted models",
    fluidRow(
      column(3,
        radioButtons("mod_which", "Model",
                     c("M5 - BYM2 + deprivation, NO2, primary care" = "M5",
                       "M6 - ESF + deprivation, NO2, primary care"  = "M6")),
        div(class = "atlas-note",
            "Residual relative risk and residual excess are measured against ",
            "a prediction that already contains the covariates, so they are ",
            "what the covariates could not explain. The decomposition map is ",
            "measured against a different reference - an area at the study ",
            "average on every covariate - and splits the whole departure ",
            "from it, which is why the two are not the same quantity."),
        div(class = "atlas-note",
            "Click any area to decompose it in the panel below.")
      ),
      column(9,
        fluidRow(
          column(6, mapbox("mod_rr",     "Smoothed relative risk")),
          column(6, mapbox("mod_resid",  "Residual relative risk (covariate-adjusted)"))
        ),
        fluidRow(
          column(6, mapbox("mod_count",  "Residual excess deaths per year")),
          column(6, mapbox("mod_decomp", "Dominant contributor to the excess"))
        )
      )
    ),
    fluidRow(
      column(12, hr(),
        h4(textOutput("decomp_title")),
        plotOutput("decomp_plot", height = 280),
        div(class = "atlas-note", textOutput("decomp_note"))
      )
    )
  ),

  # ---- 4. stroke ---------------------------------------------------------------
  tabPanel(
    "Stroke sub-model",
    fluidRow(
      column(3,
        div(class = "atlas-note",
            "Outcome is I63 cerebral infarction, all ages. The unadjusted ",
            "surface is an intercept-only BYM2; the adjusted one adds ",
            "deprivation, NO2 and primary care density, which are the ",
            "measured proxies for the latent urban-rural position the DAG ",
            "identifies as the one backdoor path. Reperfusion timing and ",
            "case fatality are mediators and are deliberately absent."),
        div(class = "atlas-note",
            paste0("Exceedance here is P(RR > ",
                   formatC(META$rr_threshold_str, format = "f", digits = 2),
                   "), not the ",
                   formatC(META$rr_threshold, format = "f", digits = 2),
                   " used for avoidable mortality.")),
        checkboxInput("str_centres", "Show stroke centres", TRUE)
      ),
      column(9,
        fluidRow(
          column(6, mapbox("str_hub",    "Travel time to nearest thrombectomy hub (SU II), minutes")),
          column(6, mapbox("str_centre", "Travel time to nearest accredited centre (SU I or II), minutes"))
        ),
        fluidRow(
          column(6, mapbox("str_unadj_rr",  "Unadjusted: smoothed relative risk")),
          column(6, mapbox("str_unadj_exc", "Unadjusted: exceedance probability"))
        ),
        fluidRow(
          column(6, mapbox("str_adj_rr",  "Adjusted: smoothed relative risk")),
          column(6, mapbox("str_adj_exc", "Adjusted: exceedance probability"))
        )
      )
    )
  ),

  tabPanel(
    "About",
    fluidRow(column(8, offset = 1,
      h4("What this is"),
      p("An exploratory front end over the outputs of the atlaspm pipeline. ",
        "Everything shown was computed by targets or by dev/build_app_data.R; ",
        "nothing is fitted here."),
      h4("What it is not"),
      tags$ul(
        tags$li("Not a causal display. The maps show a fitted model, and the ",
                "split between a covariate and the spatial term is partly a ",
                "modelling choice - it moves between M5, M6 and M7."),
        tags$li("Not evenly reliable across the map. Colour is drawn with the ",
                "same confidence in an area with three expected deaths as in ",
                "one with three hundred. The exceedance maps are where the ",
                "uncertainty is; read them alongside every other panel."),
        tags$li("Not a substitute for the reported intervals. The ",
                "decomposition chart shows posterior means only.")
      ),
      p(tags$small(if (is.null(.carto_key()))
        HTML(paste0("Basemap: Esri grey canvas (no CARTO_API_KEY set). For ",
                    "the CARTO Positron style, request a free key at ",
                    "<a href='https://carto.com/basemaps/apikey/' target='_blank'>",
                    "carto.com/basemaps/apikey</a>, add CARTO_API_KEY to your ",
                    "user ~/.Renviron and restart R."))
        else "Basemap: CARTO Positron (CARTO_API_KEY set).")),
      p(tags$small(paste("Bundle built", format(META$built_at, "%Y-%m-%d %H:%M"),
                         "|", nrow(GEO), "areas |",
                         META$n_years, "years of deaths.")))
    ))
  )
)
}

# ---- server ------------------------------------------------------------------------

server <- function(input, output, session) {

  D <- atlas_data()
  A <- D$A; GEO <- D$GEO; TAB <- D$TAB
  AREAS <- D$AREAS; LABS <- D$LABS; META <- D$META

  selected <- reactiveVal(NULL)

  # --- search and click both set the selection --------------------------------
  observeEvent(input$find, {
    if (nzchar(input$find %||% "")) selected(input$find)
  }, ignoreInit = TRUE)

  lapply(unlist(MAPS, use.names = FALSE), function(id) {
    observeEvent(input[[paste0(id, "_shape_click")]], {
      ev <- input[[paste0(id, "_shape_click")]]
      if (!is.null(ev$id)) {
        selected(ev$id)
        updateSelectizeInput(session, "find", selected = ev$id)
      }
    }, ignoreInit = TRUE)
  })

  # Recentre only the maps on the tab in view. Flying maps the user cannot see
  # costs nothing visually and a great deal in redraws.
  observeEvent(list(selected(), input$tab), {
    a <- selected()
    if (is.null(a)) return()
    ids <- switch(input$tab %||% "",
                  "Covariates"        = MAPS$covariates,
                  "Avoidable mortality" = MAPS$outcome,
                  "Adjusted models"   = MAPS$models,
                  "Stroke sub-model"  = MAPS$stroke,
                  character(0))
    focus(ids, a)
  }, ignoreInit = TRUE)

  output$selected_txt <- renderText({
    a <- selected()
    if (is.null(a)) "No area selected. Type above, or click a map."
    else paste0("Selected: ", LABS[match(a, AREAS)], "  (", a, ")")
  })

  # --- base maps ---------------------------------------------------------------
  for (id in c(MAPS$covariates, MAPS$outcome, MAPS$models)) {
    local({
      i <- id
      output[[i]] <- renderLeaflet(base_map())
    })
  }
  for (id in MAPS$stroke) {
    local({
      i <- id
      output[[i]] <- renderLeaflet(base_map(centres = A$centres))
    })
  }

  observeEvent(input$str_centres, {
    for (id in MAPS$stroke) {
      p <- leafletProxy(id)
      if (isTRUE(input$str_centres)) showGroup(p, "Stroke centres")
      else hideGroup(p, "Stroke centres")
    }
  })

  # ============================ TAB 1: covariates ==============================

  observe({
    zed <- identical(input$cov_scale, "z")
    sfx <- if (zed) "_z" else ""

    rows <- list(
      list(ids = c("cov_dep", "cov_ivsm"),
           cols = c(paste0("di_score", sfx), paste0("ivsm", sfx)),
           titles = c("Deprivation index", "IVSM")),
      list(ids = c("cov_gpw", "cov_gpu"),
           cols = c(paste0("gp_density", sfx),
                    paste0("gp_density_unweighted", sfx)),
           titles = c("GP density (weighted)", "GP density (unweighted)")),
      list(ids = c("cov_no2", "cov_pm25"),
           cols = c(paste0(META$no2, sfx), paste0(META$pm25, sfx)),
           titles = c("NO2", "PM2.5"))
    )

    for (r in rows) {
      present <- r$cols %in% names(TAB)
      # On the z view the pair shares one scale, which only means anything if
      # both members exist.
      shared <- if (zed && all(present)) {
        pal_div(c(TAB[[r$cols[1]]], TAB[[r$cols[2]]]), centre = 0)
      } else NULL

      for (k in seq_along(r$ids)) {
        if (!present[k]) {
          paint(r$ids[k], NULL,
                note = paste0("<b>", r$titles[k], "</b> is not in the bundle. ",
                              "Column <code>", r$cols[k], "</code> was not ",
                              "found; check that build_app_data.R joined it."))
          next
        }
        v <- TAB[[r$cols[k]]]
        pal <- if (!is.null(shared)) shared else
          if (zed) pal_div(v, 0) else pal_seq(v)
        paint(r$ids[k], v, pal,
              legend_title = paste0(r$titles[k], if (zed) " (z)" else ""),
              digits = if (zed) 2 else 2)
      }
    }
  })

  # ============================ TAB 2: outcome =================================

  out_menu <- reactive({
    cls <- input$out_class %||% "total"
    if (cls == "total") NULL else A$menus[[cls]]
  })

  output$out_level_ui <- renderUI({
    m <- out_menu()
    if (is.null(m) || !nrow(m)) return(NULL)
    ch <- setNames(m$stem, sprintf("%s  (%s deaths)", m$label,
                                   formatC(m$deaths, format = "f", digits = 0,
                                           big.mark = ",")))
    selectizeInput("out_level", "Level", choices = ch, width = "100%")
  })

  out_stem <- reactive({
    cls <- input$out_class %||% "total"
    if (cls == "total") "total" else input$out_level
  })

  observe({
    stem <- out_stem()
    req(!is.null(stem), nzchar(stem))
    cls  <- input$out_class %||% "total"
    lab  <- if (stem == "total") "All avoidable causes" else
      A$menus[[cls]]$label[match(stem, A$menus[[cls]]$stem)]

    # --- crude rate
    crude_col <- if (stem == "total") "total" else stem
    cv <- if (crude_col %in% names(A$crude)) {
      A$crude[[crude_col]][match(AREAS, as.character(A$crude$area))]
    } else NULL
    if (is.null(cv)) {
      paint("out_cmr", NULL, note = paste0(
        "No crude-rate column for <b>", lab, "</b>."))
    } else {
      paint("out_cmr", cv, pal_seq(cv), legend_title = "per 100,000 py",
            digits = 1)
    }

    # --- SMR
    smr_col <- paste0(stem, "_smr")
    sv <- if (smr_col %in% names(TAB)) TAB[[smr_col]] else NULL
    obs <- if (paste0(stem, "_obs") %in% names(TAB)) TAB[[paste0(stem, "_obs")]] else NA
    exp <- if (paste0(stem, "_exp") %in% names(TAB)) TAB[[paste0(stem, "_exp")]] else NA
    if (is.null(sv)) {
      paint("out_smr", NULL, note = paste0("No SMR column for <b>", lab, "</b>."))
    } else {
      paint("out_smr", sv, pal_smr(sv), legend_title = "SMR",
            extra = paste0("observed ", fmt(obs, 1), " / expected ", fmt(exp, 1)))
    }

    # --- smoothed surface and exceedance, where a model exists
    surf <- NULL
    if (stem == "total") {
      surf <- list(rr = A$surfaces$M3_rr[match(AREAS, A$surfaces$area)],
                   ex = A$surfaces$M3_exceed[match(AREAS, A$surfaces$area)],
                   what = "M3 (BYM2, no covariates)")
    } else if (cls == "mechanism") {
      rc <- paste0(stem, "_bym2"); ec <- paste0(stem, "_bym2_exc")
      if (all(c(rc, ec) %in% names(A$mech_surface))) {
        i <- match(AREAS, A$mech_surface$area)
        surf <- list(rr = A$mech_surface[[rc]][i],
                     ex = A$mech_surface[[ec]][i],
                     what = paste0("per-mechanism BYM2 (", lab, ")"))
      }
    }

    if (is.null(surf)) {
      msg <- paste0(
        "No smoothed surface for <b>", lab, "</b>. Models were fitted for all ",
        "avoidable causes and for the seven service functions. Fitting one ",
        "BYM2 per cause group or per individual cause would mostly estimate ",
        "the prior: at this level the counts are a handful of deaths spread ",
        "over ", nrow(GEO), " areas.")
      paint("out_rr",  NULL, note = msg)
      paint("out_exc", NULL, note = msg)
    } else {
      paint("out_rr", surf$rr, pal_smr(surf$rr),
            legend_title = "Relative risk",
            extra = paste0("<i>", surf$what, "</i>"))
      paint("out_exc", surf$ex, pal_exceed(surf$ex),
            legend_title = "Probability",
            extra = paste0("P(RR > ", fmt(META$rr_threshold, 2), ")"))
    }
  })

  output$out_summary <- renderUI({
    stem <- out_stem()
    req(!is.null(stem), nzchar(stem))
    oc <- paste0(stem, "_obs"); ec <- paste0(stem, "_exp")
    if (!all(c(oc, ec) %in% names(TAB))) return(NULL)
    o <- sum(TAB[[oc]], na.rm = TRUE); e <- sum(TAB[[ec]], na.rm = TRUE)
    tot <- sum(TAB$total_obs, na.rm = TRUE)
    div(class = "atlas-note",
        HTML(sprintf(
          paste0("<b>%s deaths</b> over %d years (%s per year), %s%% of all ",
                 "avoidable deaths. Study-wide SMR %s."),
          formatC(o, format = "f", digits = 1, big.mark = ","),
          META$n_years,
          formatC(o / META$n_years, format = "f", digits = 1, big.mark = ","),
          formatC(100 * o / tot, format = "f", digits = 1),
          formatC(o / e, format = "f", digits = 2))))
  })

  # ============================ TAB 3: adjusted models =========================

  mod_excess <- reactive({
    m <- input$mod_which %||% "M5"
    e <- A$excess[A$excess$model == m, , drop = FALSE]
    e[match(AREAS, e$area), , drop = FALSE]
  })

  mod_attr <- reactive({
    m <- input$mod_which %||% "M5"
    A$attribution[A$attribution$model == m, , drop = FALSE]
  })

  observe({
    m  <- input$mod_which %||% "M5"
    ex <- mod_excess()

    rr_col <- paste0(m, "_rr"); ec_col <- paste0(m, "_exceed")
    rr <- if (rr_col %in% names(A$surfaces))
      A$surfaces[[rr_col]][match(AREAS, A$surfaces$area)] else NULL

    if (is.null(rr)) {
      paint("mod_rr", NULL, note = paste0("No smoothed surface stored for ", m, "."))
    } else {
      paint("mod_rr", rr, pal_smr(rr), legend_title = "Relative risk")
    }

    paint("mod_resid", ex$rr_resid, pal_smr(ex$rr_resid),
          legend_title = "Residual RR",
          extra = paste0("P(resid RR > ", fmt(META$rr_threshold, 2), ") = ",
                         fmt(ex$p_exceed, 2),
                         ifelse(ex$flagged %in% TRUE,
                                " <b>&mdash; flagged</b>", "")))

    paint("mod_count", ex$excess, pal_div(ex$excess, 0),
          legend_title = "Excess deaths / yr",
          extra = paste0("95% CrI ", fmt(ex$excess_low, 1), " to ",
                         fmt(ex$excess_high, 1)))

    at <- mod_attr()
    if (!nrow(at)) {
      paint("mod_decomp", NULL, note = paste0(
        "No decomposition stored for ", m, "."))
    } else {
      dom <- at[at$dominant %in% TRUE, c("area", "label"), drop = FALSE]
      dom <- dom[!duplicated(dom$area), , drop = FALSE]
      tot  <- tapply(at$deaths_yr, at$area, sum)
      totv <- as.numeric(tot[AREAS])
      totv[is.na(totv)] <- 0
      lab  <- dom$label[match(AREAS, dom$area)]
      # An area whose every component is protective has no dominant driver and
      # must not be given the least-protective one by default.
      lab[is.na(lab) | totv <= 0] <- "No material excess"
      paint("mod_decomp", factor(lab), legend_title = "Dominant driver",
            extra = paste0("total excess ", fmt(totv, 1), "/yr"))
    }
  })

  output$decomp_title <- renderText({
    a <- selected()
    if (is.null(a)) "Decomposition - click an area on a map above"
    else paste0("Decomposition of ", LABS[match(a, AREAS)],
                " under ", input$mod_which %||% "M5")
  })

  output$decomp_plot <- renderPlot({
    a <- selected()
    validate(need(!is.null(a),
                  "Click an area on any map above, or use the search box."))
    d <- mod_attr()
    d <- d[d$area == a, , drop = FALSE]
    validate(need(nrow(d) > 0, "No decomposition available for this area."))

    d$label <- factor(d$label, levels = d$label[order(d$deaths_yr)])
    ggplot(d, aes(x = deaths_yr, y = label,
                  fill = deaths_yr > 0)) +
      geom_vline(xintercept = 0, linetype = "22", colour = "grey45") +
      geom_col(width = 0.65) +
      scale_fill_manual(values = c(`TRUE` = "#b5651d", `FALSE` = "#5a8a7d"),
                        guide = "none") +
      labs(x = "Attributable avoidable deaths per year (posterior mean)",
           y = NULL) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank())
  })

  output$decomp_note <- renderText({
    m <- input$mod_which %||% "M5"
    base <- paste(
      "Shapley attribution of the area's excess over a reference area of the",
      "same age and sex structure at the study-area average on every",
      "covariate. A negative bar is a determinant holding the excess down.",
      "Posterior means only - the intervals are in the thesis tables."
    )
    if (identical(m, "M6")) {
      paste(base,
            "M6 is an eigenvector spatial filter, which has no ICAR/IID split",
            "to make, so its spatial variation appears as one combined term.")
    } else {
      paste(base,
            "The split between a covariate and the spatial term is partly a",
            "modelling choice: these two components describe overlapping",
            "variation and the division between them moves between",
            "specifications.")
    }
  })

  # ============================ TAB 4: stroke ===================================

  observe({
    s <- A$stroke[match(AREAS, A$stroke$area), , drop = FALSE]

    for (spec in list(
      list(id = "str_hub",    col = "t_hub_mean",
           title = "Minutes to hub"),
      list(id = "str_centre", col = "t_centre_mean",
           title = "Minutes to nearest centre")
    )) {
      v <- if (spec$col %in% names(TAB)) TAB[[spec$col]] else NULL
      if (is.null(v)) {
        paint(spec$id, NULL, note = paste0(
          "Column <code>", spec$col, "</code> is not in the bundle."))
      } else {
        paint(spec$id, v, pal_seq(v), legend_title = spec$title, digits = 1)
      }
    }

    paint("str_unadj_rr", s$unadj_rr, pal_smr(s$unadj_rr),
          legend_title = "Relative risk",
          extra = paste0("observed ", fmt(s$i63_obs, 0), " / expected ",
                         fmt(s$i63_exp, 1)))
    paint("str_unadj_exc", s$unadj_exceed, pal_exceed(s$unadj_exceed),
          legend_title = "Probability",
          extra = paste0("P(RR > ", fmt(META$rr_threshold_str, 2), ")"))
    paint("str_adj_rr", s$adj_rr, pal_smr(s$adj_rr),
          legend_title = "Relative risk",
          extra = "adjusted for hub travel time + deprivation, NO2, primary care")
    paint("str_adj_exc", s$adj_exceed, pal_exceed(s$adj_exceed),
          legend_title = "Probability",
          extra = paste0("P(RR > ", fmt(META$rr_threshold_str, 2), ")"))
  })
}

# `ui` as a function is evaluated per request, so sourcing this file builds an
# app object and reads nothing.
shinyApp(ui = function(req) build_ui(), server = server)
