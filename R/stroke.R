#' @keywords internal
#' @noRd
#' Reverse every directed edge in a dodgr graph.
#'
#' Must be applied to an UNCONTRACTED graph. Contraction attaches hash
#' attributes tying the graph to dodgr's internal cache of its vertex table;
#' swapping endpoints afterwards leaves those attributes describing a topology
#' that no longer exists, and the C++ router reads the stale cache.
reverse_edges <- function(graph) {
  out <- graph
  out$from_id  <- graph$to_id
  out$to_id    <- graph$from_id
  out$from_lon <- graph$to_lon
  out$to_lon   <- graph$from_lon
  out$from_lat <- graph$to_lat
  out$to_lat   <- graph$from_lat
  out
}

#' @keywords internal
#' @noRd
#' Population-weighted quantile.
#'
#' @param x Numeric vector of values.
#' @param w Numeric vector of non-negative weights.
#' @param p Probability in `[0, 1]`.
#' @return Scalar double, or `NA_real_` if no finite weighted observation.
weighted_quantile_pop <- function(x, w, p = 0.5) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  o <- order(x)
  x <- x[o]
  w <- w[o]
  x[which(cumsum(w) / sum(w) >= p)[1]]
}

#' @keywords internal
#' @noRd
weighted_mean_pop <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

#' @keywords internal
#' @noRd
weighted_share_above <- function(x, w, thr) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(w[ok][x[ok] > thr]) / sum(w[ok])
}

#' @keywords internal
#' @noRd
#' Normalise a facility name for fuzzy matching.
#'
#' Lowercases, transliterates to ASCII, and strips legal forms and generic
#' hospital vocabulary, which otherwise dominate the string distance and make
#' every "Presidio Ospedaliero" look like every other one.
normalise_facility_name <- function(x) {
  x |>
    stringr::str_to_lower() |>
    iconv(to = "ASCII//TRANSLIT") |>
    stringr::str_replace_all("[^a-z0-9 ]", " ") |>
    stringr::str_replace_all(
      paste0(
        "\\b(po|presidio|ospedaliero|ospedale|osp|azienda|asst|irccs|",
        "fondazione|spa|srl|casa|di|cura|istituto|clinico|clinica|",
        "privata|accreditata|unico)\\b"
      ),
      " "
    ) |>
    stringr::str_squish()
}


# ---- 1. Import ---------------------------------------------------------------

#' Import the Lombardy stroke network node list
#'
#' Reads the node list transcribed from DGR Regione Lombardia n. XI/7473 del
#' 30/11/2022, Allegato *"Rete Stroke di Regione Lombardia"*, Fig. 3, and keeps
#' only the nodes that can actually treat a stroke.
#'
#' The DGR defines three tiers: 16 Stroke Unit di II livello (hub, mechanical
#' thrombectomy and neurosurgery h24), 25 Stroke Unit di I livello (spoke, IV
#' thrombolysis h24), and 55 hospitals with an emergency department but no
#' stroke unit. The third tier is dropped: those sites cannot recanalise, and
#' including them would measure proximity to a building rather than to
#' treatment.
#'
#' @param path Path to the node list CSV. Expected columns: `facility`,
#'   `level`, `unita_funzionale`, `ente`, `comune`, `provincia`.
#' @param levels Character vector of `level` values to retain. Defaults to both
#'   stroke unit tiers.
#'
#' @return A tibble with one row per stroke centre and a `centre_id` key.
#' @export
import_stroke_centres <- function(path,
                                  levels = c("HUB_SU_II", "SPOKE_SU_I")) {
  centres <- readr::read_csv(path, show_col_types = FALSE)

  required <- c("facility", "level", "unita_funzionale", "comune", "provincia")
  missing <- setdiff(required, names(centres))
  if (length(missing)) {
    stop(
      "Stroke centre file is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  centres <- centres |>
    dplyr::filter(.data$level %in% levels) |>
    dplyr::mutate(
      centre_id = sprintf("SC%02d", dplyr::row_number()),
      is_hub    = .data$level == "HUB_SU_II"
    ) |>
    dplyr::relocate("centre_id")

  n_hub <- sum(centres$is_hub)
  message(
    "Imported ", nrow(centres), " stroke centres (",
    n_hub, " hub / ", nrow(centres) - n_hub, " spoke)."
  )

  centres
}


# ---- 2. Geocode --------------------------------------------------------------

#' Geocode stroke centres
#'
#' Resolves each centre to a coordinate pair, preferring the Regione Lombardia
#' registry of accredited facilities over Nominatim.
#'
#' Fuzzy name matching against the registry is constrained to candidates in the
#' comune the DGR assigns to the centre. Name similarity alone produces
#' confident nonsense here: the region contains several facilities called
#' "Ospedale Civile", and a 0.9 Jaro-Winkler score between two of them is not
#' evidence of anything.
#'
#' @param centres Tibble from [import_stroke_centres()].
#' @param registry Optional tibble from the regional accredited-facilities
#'   dataset, with columns `reg_name`, `reg_comune`, `reg_lat`, `reg_lon`. When
#'   `NULL`, every centre goes to Nominatim.
#' @param overrides Named list mapping `facility` strings to `c(lat, lon)`.
#'   Applied last and unconditionally. This is where the handful of cases
#'   automated geocoding cannot resolve are pinned down by hand.
#' @param sim_min Minimum Jaro-Winkler similarity to accept a registry match.
#'
#' @return An `sf` POINT object in EPSG:4326 with a `geocode_source` column.
#' @export
geocode_stroke_centres <- function(centres,
                                   registry = NULL,
                                   overrides = list(),
                                   sim_min = 0.85) {
  centres <- centres |>
    dplyr::mutate(
      lat = NA_real_,
      lon = NA_real_,
      geocode_source = NA_character_
    )

  # -- registry match ----------------------------------------------------------
  if (!is.null(registry)) {
    reg <- registry |>
      dplyr::filter(is.finite(.data$reg_lat), is.finite(.data$reg_lon)) |>
      dplyr::mutate(
        key        = normalise_facility_name(.data$reg_name),
        comune_key = normalise_facility_name(.data$reg_comune)
      )

    cen_key <- normalise_facility_name(centres$facility)
    cen_com <- normalise_facility_name(centres$comune)

    hits <- purrr::map_dfr(seq_len(nrow(centres)), function(i) {
      cand <- dplyr::filter(reg, .data$comune_key == cen_com[i])
      if (!nrow(cand)) {
        return(tibble::tibble(lat = NA_real_, lon = NA_real_, sim = NA_real_))
      }
      d <- stringdist::stringdist(cen_key[i], cand$key, method = "jw", p = 0.1)
      j <- which.min(d)
      tibble::tibble(lat = cand$reg_lat[j], lon = cand$reg_lon[j], sim = 1 - d[j])
    })

    keep <- !is.na(hits$sim) & hits$sim >= sim_min
    centres$lat[keep] <- hits$lat[keep]
    centres$lon[keep] <- hits$lon[keep]
    centres$geocode_source[keep] <- "regione_lombardia"

    message(sum(keep), " / ", nrow(centres), " centres matched to the registry.")
  }

  # -- Nominatim fallback ------------------------------------------------------
  todo <- which(is.na(centres$lat))
  if (length(todo)) {
    query <- paste0(
      centres$facility[todo], ", ",
      # Multi-site presidi carry a "/" in comune; geocode the first site and
      # handle the second through `overrides`.
      stringr::str_remove(centres$comune[todo], "\\s*/.*$"), ", ",
      centres$provincia[todo], ", Italy"
    )
    nom <- tidygeocoder::geo(
      address = query, method = "osm",
      lat = "lat_n", long = "lon_n", quiet = TRUE
    )
    centres$lat[todo] <- nom$lat_n
    centres$lon[todo] <- nom$lon_n
    centres$geocode_source[todo] <- "nominatim"
    message(sum(!is.na(nom$lat_n)), " / ", length(todo), " resolved by Nominatim.")
  }

  # -- manual overrides --------------------------------------------------------
  for (nm in names(overrides)) {
    i <- which(centres$facility == nm)
    if (!length(i)) {
      warning("Override facility not found: ", nm, call. = FALSE)
      next
    }
    centres$lat[i] <- overrides[[nm]][1]
    centres$lon[i] <- overrides[[nm]][2]
    centres$geocode_source[i] <- "manual"
  }

  unresolved <- centres |> dplyr::filter(is.na(.data$lat) | is.na(.data$lon))
  if (nrow(unresolved)) {
    stop(
      nrow(unresolved), " centres could not be geocoded:\n",
      paste0("  - ", unresolved$facility, collapse = "\n"),
      "\nSupply coordinates through the `overrides` argument.",
      call. = FALSE
    )
  }

  sf::st_as_sf(centres, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
}


#' Validate geocoded stroke centres against administrative boundaries
#'
#' A gate, not a report. One hub displaced by a few kilometres silently
#' corrupts every travel time in its catchment, and the corruption is invisible
#' downstream because the numbers stay plausible.
#'
#' @param centres `sf` object from [geocode_stroke_centres()].
#' @param comuni_shp `sf` polygons carrying a comune name column.
#' @param name_col Name of the comune name column in `comuni_shp`.
#' @param stop_on_fail Whether to error (default) or warn when a centre falls
#'   outside its declared comune.
#'
#' @return The input `centres`, invisibly, when all checks pass.
#' @export
check_stroke_centres <- function(centres,
                                 comuni_shp,
                                 name_col = "COMUNE",
                                 stop_on_fail = TRUE) {
  crs_m <- 32632L # UTM 32N, metric, correct for Lombardy

  pts <- sf::st_transform(centres, crs_m)
  poly <- comuni_shp |>
    sf::st_transform(crs_m) |>
    sf::st_make_valid() |>
    dplyr::select(dplyr::all_of(name_col))

  joined <- sf::st_join(pts, poly, join = sf::st_within) |>
    sf::st_drop_geometry() |>
    dplyr::mutate(
      # Split presidi legitimately fail a single-comune test.
      multi_site = stringr::str_detect(.data$comune, "/"),
      declared   = normalise_facility_name(.data$comune),
      found      = normalise_facility_name(.data[[name_col]]),
      pass       = .data$multi_site |
        (!is.na(.data$found) & .data$declared == .data$found)
    )

  fails <- dplyr::filter(joined, !.data$pass)
  if (nrow(fails)) {
    msg <- paste0(
      nrow(fails), " stroke centres are not inside their declared comune:\n",
      paste0(
        "  - ", fails$facility, ": declared ", fails$comune,
        ", geocoded into ", fails[[name_col]],
        " (source: ", fails$geocode_source, ")",
        collapse = "\n"
      ),
      "\nResolve each through the `overrides` argument of ",
      "geocode_stroke_centres()."
    )
    if (stop_on_fail) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  } else {
    message("All ", nrow(joined), " stroke centres validated against comune boundaries.")
  }

  invisible(centres)
}


# ---- 3. Origins --------------------------------------------------------------

#' Build population-weighted origin points from census sections
#'
#' Census sections are the finest geography ISTAT publishes population for.
#' Routing from them, rather than from an area centroid, is what makes the
#' population weighting meaningful: travel time is convex over a road network,
#' so the time *from* a weighted centroid systematically understates the
#' population-weighted *mean* of times. The gap is negligible in compact
#' comuni and material wherever settlement is strung out, which is precisely
#' where an accessibility covariate is supposed to carry signal.
#'
#' Sections are assigned to modelling areas by spatial join against `area_shp`
#' rather than by reconstructing the composite `area` key. This guarantees the
#' origin layer partitions exactly the same areas the model is fitted on, and
#' survives any future change to how that key is built.
#'
#' @param sez_shp `sf` polygons of ISTAT census sections.
#' @param area_shp `sf` polygons of modelling areas, carrying an `area` column.
#' @param pop_col Name of the resident population column in `sez_shp`.
#' @param drop_unpopulated Drop sections with zero residents. They contribute
#'   nothing to a population-weighted statistic and are typically a third of
#'   all sections, so this is the cheapest speed-up in the pipeline.
#'
#' @return A tibble with one row per populated section: `sez_id`, `area`,
#'   `pop`, `lon`, `lat`.
#' @export
build_section_points <- function(sez_shp,
                                 area_shp,
                                 pop_col = "P1",
                                 drop_unpopulated = TRUE) {
  if (!pop_col %in% names(sez_shp)) {
    stop(
      "Population column '", pop_col, "' not found in census sections.\n",
      "Available: ", paste(names(sez_shp), collapse = ", "), "\n",
      "Column names differ between the 2011 and 2021 ISTAT releases.",
      call. = FALSE
    )
  }
  if (!"area" %in% names(area_shp)) {
    stop("`area_shp` must carry an `area` column.", call. = FALSE)
  }

  crs_m <- 32632L

  sez <- sez_shp |>
    sf::st_transform(crs_m) |>
    sf::st_make_valid() |>
    dplyr::mutate(
      sez_id = as.character(dplyr::row_number()),
      pop    = as.numeric(.data[[pop_col]])
    ) |>
    dplyr::filter(!is.na(.data$pop))

  if (drop_unpopulated) {
    n0 <- nrow(sez)
    sez <- dplyr::filter(sez, .data$pop > 0)
    message("Dropped ", n0 - nrow(sez), " unpopulated sections; ",
            nrow(sez), " retained.")
  }

  # point_on_surface, not centroid: a centroid can fall outside a concave
  # section (common where sections wrap a block or follow a watercourse),
  # which then snaps to the wrong side of the road network.
  pts <- suppressWarnings(sf::st_point_on_surface(sez))

  pts <- sf::st_join(
    pts,
    area_shp |> sf::st_transform(crs_m) |> sf::st_make_valid() |>
      dplyr::select("area"),
    join = sf::st_within
  )

  orphan <- sum(is.na(pts$area))
  if (orphan) {
    message(orphan, " sections (", format(sum(pts$pop[is.na(pts$area)]),
                                          big.mark = ","),
            " residents) fall outside every modelling area and are dropped.")
  }
  pts <- dplyr::filter(pts, !is.na(.data$area))

  ll <- pts |> sf::st_transform(4326) |> sf::st_coordinates()

  tibble::tibble(
    sez_id = pts$sez_id,
    area   = pts$area,
    pop    = pts$pop,
    lon    = ll[, 1],
    lat    = ll[, 2]
  )
}


#' Build the urban mask for the AREU speed model
#'
#' AREU distinguishes urban from extra-urban roads. ISTAT's `TIPO_LOC` already
#' encodes exactly that distinction (`1` = *centro abitato*), so the mask is
#' derived from the census sections rather than guessed from OSM tags.
#'
#' @param sez_shp `sf` polygons of ISTAT census sections.
#' @param tipo_loc_col Name of the locality-type column.
#'
#' @return An `sfc` geometry in EPSG:32632.
#' @export
build_urban_mask <- function(sez_shp, tipo_loc_col = "TIPO_LOC") {
  if (!tipo_loc_col %in% names(sez_shp)) {
    stop("Locality-type column '", tipo_loc_col, "' not found.", call. = FALSE)
  }
  sez_shp |>
    sf::st_transform(32632L) |>
    dplyr::filter(as.integer(.data[[tipo_loc_col]]) == 1L) |>
    sf::st_geometry() |>
    sf::st_union() |>
    sf::st_make_valid()
}


# ---- 4. Network --------------------------------------------------------------

#' Build a routable street network with the AREU speed model
#'
#' Downloads a Geofabrik OSM extract and builds a contracted `dodgr` graph. No
#' routing server and no Docker: everything runs in-process.
#'
#' `speed_model = "areu"` imposes the assumptions Regione Lombardia used to
#' draw the DGR's own centralisation maps (Figs. 4-6 of the Allegato): 33 km/h
#' urban, 60 km/h extra-urban, 90 km/h motorway. Use it if the resulting
#' covariate is to be commensurable with the regional 45-minute mothership
#' threshold. `"osm"` keeps dodgr's stock motorcar profile, which models a
#' private car in free flow and is *not* comparable to the DGR figures.
#'
#' Neither model includes congestion, dispatch delay, or on-scene time. What
#' comes out is network accessibility, not a realised prehospital interval, and
#' should be described that way.
#'
#' @param aoi `sf` or `sfc` polygon bounding the area to route within. Must
#'   extend well beyond the study area or routes that legitimately leave it
#'   will be truncated.
#' @param urban_mask `sfc` from [build_urban_mask()]. Required when
#'   `speed_model = "areu"`.
#' @param speed_model One of `"areu"` or `"osm"`.
#' @param osm_dir Directory for the cached `.pbf`.
#' @param place Geofabrik place name passed to [osmextract::oe_get()].
#' @param speeds Named numeric vector of km/h for the AREU model.
#'
#' @return A contracted `dodgr` graph (a data frame).
#' @export
build_stroke_network <- function(aoi,
                                 urban_mask = NULL,
                                 speed_model = c("areu", "osm"),
                                 osm_dir = "data/osm",
                                 place = "Lombardia",
                                 speeds = c(urban = 33, extra_urban = 60,
                                            motorway = 90)) {
  speed_model <- match.arg(speed_model)
  if (speed_model == "areu" && is.null(urban_mask)) {
    stop("`urban_mask` is required when speed_model = 'areu'.", call. = FALSE)
  }
  dir.create(osm_dir, showWarnings = FALSE, recursive = TRUE)

  drivable <- c(
    "motorway", "motorway_link", "trunk", "trunk_link",
    "primary", "primary_link", "secondary", "secondary_link",
    "tertiary", "tertiary_link", "unclassified", "residential",
    "living_street", "service", "road"
  )

  roads <- osmextract::oe_get(
    place = place,
    provider = "geofabrik",
    layer = "lines",
    download_directory = osm_dir,
    extra_tags = c("maxspeed", "oneway", "access", "service"),
    # SELECT * so the geometry column always survives, and so the columns
    # osmextract adds via extra_tags are actually present. The WHERE clause
    # still pushes the row filter down into GDAL, which is where the saving is.
    query = paste0(
      "SELECT * FROM lines WHERE highway IN ('",
      paste(drivable, collapse = "','"), "')"
    ),
    boundary = sf::st_transform(sf::st_geometry(aoi), 4326),
    boundary_type = "spat",
    quiet = TRUE
  )

  if (!inherits(roads, "sf")) {
    stop(
      "OSM read returned an object of class ",
      paste(class(roads), collapse = "/"),
      " rather than sf. The geometry column was dropped by the GDAL query.\n",
      "Columns returned: ", paste(names(roads), collapse = ", "),
      call. = FALSE
    )
  }

  if ("access" %in% names(roads)) {
    roads <- roads[is.na(roads$access) |
                     !roads$access %in% c("private", "no", "customers"), ]
  }
  if ("service" %in% names(roads)) {
    roads <- roads[is.na(roads$service) |
                     !roads$service %in% c("parking_aisle", "driveway"), ]
  }
  message("Retained ", nrow(roads), " drivable ways.")
  graph <- dodgr::weight_streetnet(
    roads,
    wt_profile = "motorcar",
    type_col = "highway",
    id_col = "osm_id"
  )

  # Keep the largest connected component. Islands and mis-tagged fragments
  # otherwise yield infinite travel times that are artefacts, not findings.
  graph <- dodgr::dodgr_components(graph)
  sizes <- sort(table(graph$component), decreasing = TRUE)
  graph <- graph[graph$component == as.integer(names(sizes)[1]), ]
  message("Network: ", format(nrow(graph), big.mark = ","), " edges (",
          round(100 * sizes[1] / sum(sizes), 1), "% of the largest component).")

  if (speed_model == "areu") {
    mid <- sf::st_as_sf(
      data.frame(
        x = (graph$from_lon + graph$to_lon) / 2,
        y = (graph$from_lat + graph$to_lat) / 2
      ),
      coords = c("x", "y"), crs = 4326
    ) |> sf::st_transform(32632L)

    in_urban <- lengths(
      sf::st_intersects(mid, sf::st_sf(geometry = urban_mask))
    ) > 0
    rm(mid)
    invisible(gc())

    is_mw <- graph$highway %in% c("motorway", "motorway_link",
                                  "trunk", "trunk_link")
    kmh <- ifelse(is_mw, speeds[["motorway"]],
                  ifelse(in_urban, speeds[["urban"]], speeds[["extra_urban"]]))

    message("Speeds -- motorway: ", sum(is_mw),
            " | urban: ", sum(!is_mw & in_urban),
            " | extra-urban: ", sum(!is_mw & !in_urban))

    # dodgr distances are metres and times seconds; dodgr_times() routes on
    # `time_weighted` and reports `time`, so both must be set.
    graph$time <- graph$d / (kmh / 3.6)
    graph$time_weighted <- graph$time
  }

  dodgr::dodgr_contract_graph(graph)
}

# ---- 5. Travel times ---------------------------------------------------------

#' Compute travel time from every census section to the nearest stroke centre
#'
#' Routes one-to-many from each of the 41 stroke centres across the network,
#' which is 41 Dijkstra runs rather than one per census section.
#'
#' Because the graph is not reversed, what is computed is the time *from* each
#' centre *to* each section, and it is used as an estimate of the journey in
#' the opposite direction. The two differ only where one-way restrictions make
#' the return leg longer. Against a speed model with no congestion term and
#' 33/60/90 km/h class speeds, that discrepancy sits well below the resolution
#' of the measure, but it is an approximation and should be described as one.
#'
#' @param network Contracted `dodgr` graph from [build_stroke_network()].
#' @param centres `sf` from [geocode_stroke_centres()], carrying `is_hub`,
#'   `lon` and `lat`.
#' @param sections Tibble from [build_section_points()], carrying `lon`, `lat`,
#'   `pop` and `area`.
#' @param max_snap_m Warn when a centre snaps further than this to the network.
#'
#' @return A tibble: one row per section, with `t_centre_min`, `t_hub_min`, the
#'   identity of the nearest centre and hub, and the `area` key.
#' @export
build_stroke_times <- function(network, centres, sections, max_snap_m = 500) {
  cen <- sf::st_drop_geometry(centres)

  # dodgr expects a plain numeric matrix with columns x/y. Passing a tibble
  # slice such as `cen[, c("lon", "lat")]` fails inside the C++ layer with
  # "Index out of bounds: [index='x']", because dodgr's column lookup does not
  # resolve against a tibble the way it does against a base data frame.
  xy <- function(lon, lat) {
    m <- as.matrix(data.frame(x = as.numeric(lon), y = as.numeric(lat)))
    storage.mode(m) <- "double"
    m
  }

  verts <- dodgr::dodgr_vertices(network)
  i_cen <- dodgr::match_points_to_verts(verts, xy(cen$lon, cen$lat))
  i_sez <- dodgr::match_points_to_verts(verts, xy(sections$lon, sections$lat))

  # How far each point had to move to reach the network. A centre dragged
  # hundreds of metres is a geocoding error, not a routing subtlety.
  snap_m <- function(lon, lat, idx) {
    a <- sf::st_as_sf(data.frame(lon = lon, lat = lat),
                      coords = c("lon", "lat"), crs = 4326) |>
      sf::st_transform(32632L)
    b <- sf::st_as_sf(data.frame(lon = verts$x[idx], lat = verts$y[idx]),
                      coords = c("lon", "lat"), crs = 4326) |>
      sf::st_transform(32632L)
    as.numeric(sf::st_distance(a, b, by_element = TRUE))
  }

  d_cen <- snap_m(cen$lon, cen$lat, i_cen)
  far <- which(d_cen > max_snap_m)
  if (length(far)) {
    warning(
      "Stroke centres snapped >", max_snap_m, " m to the road network. ",
      "This is a geocoding error, not a routing subtlety:\n",
      paste0("  - ", cen$facility[far], " (", round(d_cen[far]), " m)",
             collapse = "\n"),
      call. = FALSE
    )
  }
  message("Centre snap distance: median ", round(stats::median(d_cen)),
          " m, max ", round(max(d_cen)), " m.")

  message("Routing ", nrow(cen), " centres to ",
          format(nrow(sections), big.mark = ","), " sections...")
  tmat <- dodgr::dodgr_times(
    network,
    from = verts$id[i_cen],
    to   = verts$id[i_sez]
  ) / 60 # dodgr returns seconds

  stopifnot(nrow(tmat) == nrow(cen), ncol(tmat) == nrow(sections))

  safe_min   <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_which <- function(x) if (all(is.na(x))) NA_integer_ else which.min(x)

  hub <- tmat[cen$is_hub, , drop = FALSE]

  out <- sections |>
    dplyr::mutate(
      t_centre_min   = apply(tmat, 2, safe_min),
      t_hub_min      = apply(hub, 2, safe_min),
      nearest_centre = cen$facility[apply(tmat, 2, safe_which)],
      nearest_hub    = cen$facility[cen$is_hub][apply(hub, 2, safe_which)]
    )

  n_bad <- sum(is.na(out$t_centre_min))
  if (n_bad) {
    warning(
      n_bad, " sections (",
      format(sum(out$pop[is.na(out$t_centre_min)]), big.mark = ","),
      " residents) reached no centre. Almost always a disconnected network ",
      "fragment rather than a finding -- inspect before treating as missing.",
      call. = FALSE
    )
  }

  out
}

# ---- 6. Aggregate to modelling areas -----------------------------------------

#' Aggregate section travel times to modelling areas
#'
#' Population weighting is applied to the *times*, not to the geometry: each
#' section carries its own travel time and its own resident count, and the
#' weighted statistics are formed at this step.
#'
#' Section populations need not reconcile with `pop_area_table`; they are used
#' only as relative weights within an area, so a different census vintage
#' changes nothing material.
#'
#' `pop_share_over_45min_hub` is the covariate worth reaching for first. It is
#' bounded on `[0, 1]`, which behaves better in a BYM2 linear predictor than a
#' raw mean in minutes, and it maps directly onto the DGR's own mothership
#' criterion rather than onto an arbitrary threshold.
#'
#' @param stroke_times Tibble from [build_stroke_times()].
#' @param thresholds Named numeric vector of minute thresholds. Defaults to the
#'   DGR's own: 45 minutes for mothership centralisation, 60 for helicopter.
#'
#' @return A tibble with one row per `area`.
#' @export
build_stroke_area <- function(stroke_times,
                              thresholds = c(mothership = 45, helicopter = 60)) {
  stroke_times |>
    dplyr::group_by(.data$area) |>
    dplyr::summarise(
      n_sections    = dplyr::n(),
      pop_sections  = sum(.data$pop, na.rm = TRUE),

      t_centre_mean = weighted_mean_pop(.data$t_centre_min, .data$pop),
      t_centre_p50  = weighted_quantile_pop(.data$t_centre_min, .data$pop, 0.50),
      t_centre_p90  = weighted_quantile_pop(.data$t_centre_min, .data$pop, 0.90),

      t_hub_mean    = weighted_mean_pop(.data$t_hub_min, .data$pop),
      t_hub_p50     = weighted_quantile_pop(.data$t_hub_min, .data$pop, 0.50),
      t_hub_p90     = weighted_quantile_pop(.data$t_hub_min, .data$pop, 0.90),

      pop_share_over_45min_hub = weighted_share_above(
        .data$t_hub_min, .data$pop, thresholds[["mothership"]]
      ),
      pop_share_over_60min_hub = weighted_share_above(
        .data$t_hub_min, .data$pop, thresholds[["helicopter"]]
      ),

      dominant_centre = {
        tt <- tapply(.data$pop, .data$nearest_centre, sum)
        if (!length(tt)) NA_character_ else names(which.max(tt))
      },
      dominant_hub = {
        tt <- tapply(.data$pop, .data$nearest_hub, sum)
        if (!length(tt)) NA_character_ else names(which.max(tt))
      },
      .groups = "drop"
    )
}
#' Attach stroke accessibility covariates to the modelling sf
#'
#' Standardisation happens here, not in [build_stroke_area()], and deliberately
#' so. `stroke_area` covers every Lombardy area the census sections touch,
#' while the model is fitted on a subset. Z-scoring upstream would centre the
#' covariate on the regional mean rather than on the modelled one, leaving it
#' off-centre in the linear predictor and giving the coefficient a unit that no
#' other regressor in the model shares.
#'
#' @param shp `sf` of modelling areas (typically `smr_geo`).
#' @param stroke_area Tibble from [build_stroke_area()].
#' @param vars Raw covariate columns to attach.
#' @param scale_vars Subset of `vars` to additionally attach as `_z` columns,
#'   standardised across the modelled areas only.
#'
#' @return `shp` with the requested columns joined on `area`.
#' @export
add_stroke_access <- function(shp,
                              stroke_area,
                              vars = c("t_centre_mean", "t_hub_mean",
                                       "t_hub_p90", "pop_share_over_45min_hub"),
                              scale_vars = c("t_centre_mean", "t_hub_mean")) {
  missing_vars <- setdiff(c(vars, scale_vars), names(stroke_area))
  if (length(missing_vars)) {
    stop("Columns absent from `stroke_area`: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  # Join on the bare table and reattach geometry. dplyr's sf methods are not
  # reliably registered in a crew worker, and when they are not, left_join()
  # dispatches to the tibble method and errors on the sfc column.
  geom <- sf::st_geometry(shp)
  tab  <- sf::st_drop_geometry(shp)

  tab <- dplyr::left_join(
    tab,
    dplyr::select(stroke_area, dplyr::all_of(c("area", union(vars, scale_vars)))),
    by = "area"
  )

  unmatched <- sum(is.na(tab[[vars[1]]]))
  if (unmatched) {
    warning(
      unmatched, " of ", nrow(tab), " modelling areas received no stroke ",
      "accessibility value. Check that `area_shp` and the census sections ",
      "cover the same territory.",
      call. = FALSE
    )
  }

  for (v in scale_vars) {
    tab[[paste0(v, "_z")]] <- as.numeric(scale(tab[[v]]))
  }

  sf::st_set_geometry(tab, geom)
}

# ---- 7. Diagnostics ----------------------------------------------------------

#' Diagnose the stroke accessibility covariate before modelling
#'
#' Two questions this answers, both of which should be settled before the
#' covariate enters a model.
#'
#' *How much would the centroid shortcut have cost?* The convexity gap between
#' the population-weighted mean of section times and the time at the area's
#' median section. Worth reporting in the methods as justification for routing
#' at section level rather than merely asserting that it matters.
#'
#' *Is there variance to model?* Milan holds four hubs and three spokes, so
#' every NIL sits minutes from a stroke centre. A covariate with a
#' sub-3-minute standard deviation across the NIL layer carries essentially no
#' signal, and a null coefficient there is a property of the geography rather
#' than evidence about access.
#'
#' @param stroke_times Tibble from [build_stroke_times()].
#' @param stroke_area Tibble from [build_stroke_area()].
#' @param print Whether to message the summary.
#'
#' @return A list with `centroid_bias`, `layer_variance`, and `headline`.
#' @export
check_stroke_access <- function(stroke_times, stroke_area, print = TRUE) {
  centroid_proxy <- stroke_times |>
    dplyr::group_by(.data$area) |>
    dplyr::summarise(
      t_hub_median_section = stats::median(.data$t_hub_min, na.rm = TRUE),
      .groups = "drop"
    )

  bias <- stroke_area |>
    dplyr::select("area", "pop_sections", "t_hub_mean") |>
    dplyr::left_join(centroid_proxy, by = "area") |>
    dplyr::mutate(
      centroid_bias_min = .data$t_hub_mean - .data$t_hub_median_section
    ) |>
    dplyr::arrange(dplyr::desc(abs(.data$centroid_bias_min)))

  # NIL keys carry an underscore; comune keys do not.
  layer_variance <- stroke_area |>
    dplyr::mutate(
      layer = ifelse(stringr::str_detect(.data$area, "_"), "nil", "comune")
    ) |>
    dplyr::group_by(.data$layer) |>
    dplyr::summarise(
      n       = dplyr::n(),
      mean    = mean(.data$t_hub_mean, na.rm = TRUE),
      sd      = stats::sd(.data$t_hub_mean, na.rm = TRUE),
      min     = min(.data$t_hub_mean, na.rm = TRUE),
      max     = max(.data$t_hub_mean, na.rm = TRUE),
      .groups = "drop"
    )

  headline <- list(
    pop_covered = sum(stroke_times$pop, na.rm = TRUE),
    mean_to_centre = weighted_mean_pop(stroke_times$t_centre_min,
                                       stroke_times$pop),
    mean_to_hub = weighted_mean_pop(stroke_times$t_hub_min, stroke_times$pop),
    share_over_45 = weighted_share_above(stroke_times$t_hub_min,
                                         stroke_times$pop, 45),
    max_centroid_bias = max(abs(bias$centroid_bias_min), na.rm = TRUE)
  )

  if (print) {
    message(
      "\nStroke accessibility diagnostics\n",
      "  population covered:            ",
      format(headline$pop_covered, big.mark = ","), "\n",
      "  pop-weighted mean to centre:   ",
      round(headline$mean_to_centre, 1), " min\n",
      "  pop-weighted mean to hub:      ",
      round(headline$mean_to_hub, 1), " min\n",
      "  residents >45 min from a hub:  ",
      round(100 * headline$share_over_45, 2), "%\n",
      "  largest centroid bias:         ",
      round(headline$max_centroid_bias, 1), " min"
    )
    nil <- dplyr::filter(layer_variance, .data$layer == "nil")
    if (nrow(nil) && !is.na(nil$sd) && nil$sd < 3) {
      message(
        "  NOTE: NIL-layer SD is ", round(nil$sd, 2), " min. Milan contains ",
        "four hubs; this covariate has almost no sub-municipal variance and ",
        "should not be expected to explain NIL-level mortality."
      )
    }
  }

  list(centroid_bias = bias, layer_variance = layer_variance,
       headline = headline)
}


#' Read a previously computed section-level travel-time table
#'
#' Restores the output of [build_stroke_times()] from a file, so the pipeline
#' can resume from the routing results without rebuilding the road network.
#'
#' \strong{Why this exists.} [build_stroke_network()] takes 10-40 minutes and
#' 8-16 GB, and needs an OSM extract on disk. Nothing downstream of
#' `stroke_times` depends on the network object itself, so once the routing has
#' been run and its output saved, every subsequent analysis can start from the
#' table. The routing functions stay in the package and stay tested; only the
#' pipeline's entry point moves.
#'
#' \strong{What it validates.} The columns [build_stroke_area()] consumes, and
#' the two failure modes that would otherwise pass silently: an `area` key that
#' does not match the modelling geography (aggregation would drop those
#' sections), and travel times that are not plausibly minutes. A file saved
#' from an older run with seconds, or with kilometres, would aggregate happily
#' and produce an accessibility surface wrong by a constant factor.
#'
#' @param file_path Path to the saved table. `.rds`, `.csv` and `.parquet` are
#'   recognised by extension.
#' @param areas Optional character vector of valid modelling-area keys, for
#'   checking. Normally `area_shp$area`.
#' @param max_minutes Upper bound for a plausible travel time, used only to
#'   catch a units mismatch. Default `600`.
#'
#' @return A tibble with the columns [build_stroke_area()] expects: `area`,
#'   `pop`, `t_centre_min`, `t_hub_min`, `nearest_centre`, `nearest_hub`, plus
#'   whatever else the file carries.
#'
#' @examples
#' \dontrun{
#' stroke_times <- import_stroke_times(
#'   get_input_data_path("stroke_times.rds"), areas = area_shp$area
#' )
#' }
#' @seealso [build_stroke_times()], [build_stroke_area()]
#' @export
import_stroke_times <- function(file_path, areas = NULL, max_minutes = 600) {

  ext <- tolower(tools::file_ext(file_path))
  out <- switch(
    ext,
    rds     = readRDS(file_path),
    csv     = readr::read_csv(file_path, show_col_types = FALSE,
                              col_types = readr::cols(
                                area           = readr::col_character(),
                                nearest_centre = readr::col_character(),
                                nearest_hub    = readr::col_character(),
                                .default       = readr::col_double())),
    parquet = {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("Reading .parquet needs the arrow package.", call. = FALSE)
      }
      tibble::as_tibble(arrow::read_parquet(file_path))
    },
    stop("Unrecognised extension '", ext, "'. Expected rds, csv or parquet.",
         call. = FALSE)
  )

  needed <- c("area", "pop", "t_centre_min", "t_hub_min",
              "nearest_centre", "nearest_hub")
  require_cols(out, needed, "saved stroke_times")

  out[["area"]] <- as.character(out[["area"]])

  # Units. A table saved in seconds would aggregate without complaint and give
  # an accessibility surface wrong by a factor of sixty.
  for (v in c("t_centre_min", "t_hub_min")) {
    x <- out[[v]][is.finite(out[[v]])]
    if (length(x) && max(x) > max_minutes) {
      stop("`", v, "` reaches ", round(max(x)), ", above the plausible bound ",
           "of ", max_minutes, " minutes. The saved file is probably in ",
           "seconds rather than minutes.", call. = FALSE)
    }
    if (length(x) && any(x < 0)) {
      stop("`", v, "` contains negative values.", call. = FALSE)
    }
  }

  n_unreached <- sum(!is.finite(out[["t_hub_min"]]))

  if (!is.null(areas)) {
    missing_area <- setdiff(unique(out[["area"]]), areas)
    empty_area   <- setdiff(areas, unique(out[["area"]]))
    if (length(missing_area)) {
      warning(length(missing_area), " area key(s) in the saved routing are not ",
              "in the modelling geography and will be dropped on aggregation: ",
              paste(utils::head(missing_area, 8), collapse = ", "),
              call. = FALSE)
    }
    if (length(empty_area)) {
      warning(length(empty_area), " modelling area(s) have no routed section ",
              "and will get a missing travel time: ",
              paste(utils::head(empty_area, 8), collapse = ", "),
              call. = FALSE)
    }
  }

  message(sprintf(
    "stroke_times restored: %s sections, %s areas, %s population, %s unreached.",
    format(nrow(out), big.mark = ","),
    format(dplyr::n_distinct(out[["area"]]), big.mark = ","),
    format(round(sum(out[["pop"]], na.rm = TRUE)), big.mark = ","),
    format(n_unreached, big.mark = ",")))

  tibble::as_tibble(out)
}
