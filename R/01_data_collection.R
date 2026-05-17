# =============================================================================
# PharmOptim-Rennes | Script 01 — Collecte & Préparation des données
# Auteur : Hafiz | Statistique Spatiale — Magistère 2
# =============================================================================
# Unité statistique : IRIS (IGN 2024 — Lambert-93)
# Sources :
#   - Contours IRIS   : IGN CONTOURS-IRIS 3.0 (2024)
#   - Population IRIS : INSEE Recensement 2022
#   - POI             : OpenStreetMap via osmdata
# =============================================================================

library(osmdata)
library(sf)
library(tidyverse)

dir.create("data", showWarnings = FALSE)

BBOX_RENNES <- getbb("Rennes, France")

cat("========================================\n")
cat("  PharmOptim-Rennes | Collecte données  \n")
cat("========================================\n\n")

# =============================================================================
# 1. CONTOUR COMMUNE DE RENNES
# =============================================================================

cat("[1/6] Contour de la commune de Rennes...\n")

contour_raw <- opq(bbox = BBOX_RENNES) %>%
  add_osm_feature(key = "boundary",    value = "administrative") %>%
  add_osm_feature(key = "admin_level", value = "8") %>%
  add_osm_feature(key = "name",        value = "Rennes") %>%
  osmdata_sf()

rennes_poly <- contour_raw$osm_multipolygons %>%
  filter(name == "Rennes") %>%
  st_transform(2154)

if (nrow(rennes_poly) == 0) {
  rennes_poly <- contour_raw$osm_polygons %>%
    filter(name == "Rennes") %>%
    st_transform(2154)
}

cat(sprintf("  OK — %d polygone(s) commune\n\n", nrow(rennes_poly)))
saveRDS(rennes_poly, "data/rennes_contour.rds")

# =============================================================================
# 2. CONTOURS IRIS DE RENNES — IGN 2024
# =============================================================================

cat("[2/6] Chargement des contours IRIS IGN 2024...\n")

shp_iris <- paste0(
  "data/raw/IRIS/",
  "CONTOURS-IRIS_3-0__SHP_LAMB93_FXX_2024-01-01/CONTOURS-IRIS/",
  "1_DONNEES_LIVRAISON_2024-12-00164/",
  "CONTOURS-IRIS_3-0_SHP_LAMB93_FXX-ED2024-01-01/",
  "CONTOURS-IRIS.shp"
)

if (!file.exists(shp_iris)) stop(paste("Fichier IRIS introuvable :\n", shp_iris))

iris_rennes <- st_read(shp_iris, quiet = TRUE) %>%
  filter(INSEE_COM == "35238") %>%
  st_transform(2154)

cat(sprintf("  OK — %d IRIS chargés pour Rennes\n\n", nrow(iris_rennes)))

# =============================================================================
# 3. POPULATION PAR IRIS — INSEE Recensement 2022
# =============================================================================

cat("[3/6] Chargement des données de population INSEE RP 2022...\n")

pop_path <- "data/raw/base-ic-evol-struct-pop-2022.CSV"
if (!file.exists(pop_path)) stop(paste("Fichier population introuvable :", pop_path))

# Détection automatique du séparateur
first_line <- readLines(pop_path, n = 1, encoding = "latin1")
delim_used <- ifelse(str_count(first_line, ";") > str_count(first_line, ","), ";", ",")
cat(sprintf("  Séparateur détecté : '%s'\n", delim_used))

pop_raw <- read_delim(
  pop_path,
  delim          = delim_used,
  locale         = locale(encoding = "latin1"),
  col_types      = cols(.default = "c"),
  show_col_types = FALSE
)

pop_rennes <- pop_raw %>%
  filter(str_starts(IRIS, "35238")) %>%
  transmute(
    CODE_IRIS  = IRIS,
    POP_TOTAL  = as.numeric(P22_POP),
    POP_0_14   = as.numeric(P22_POP0014),
    POP_15_29  = as.numeric(P22_POP1529),
    POP_65P    = as.numeric(P22_POP65P),
    PART_0_14  = POP_0_14  / POP_TOTAL,
    PART_15_29 = POP_15_29 / POP_TOTAL,
    PART_65P   = POP_65P   / POP_TOTAL,
    POP_ETR    = as.numeric(P22_POP_ETR),
    PART_ETR   = POP_ETR   / POP_TOTAL
  ) %>%
  mutate(across(starts_with("PART_"),
                ~ifelse(is.nan(.x) | is.infinite(.x), 0, .x)))

cat(sprintf("  OK — %d IRIS de Rennes\n", nrow(pop_rennes)))
cat(sprintf("  Population totale : %s habitants\n",
            format(sum(pop_rennes$POP_TOTAL, na.rm = TRUE), big.mark = " ")))

# Jointure contours IRIS + population INSEE
iris_rennes <- iris_rennes %>%
  left_join(pop_rennes, by = "CODE_IRIS") %>%
  mutate(
    SUPERFICIE_KM2      = as.numeric(st_area(geometry)) / 1e6,
    DENSITE_POP         = POP_TOTAL / SUPERFICIE_KM2,
    SCORE_VULNERABILITE = (replace_na(PART_65P,  0) * 0.6 +
                             replace_na(PART_0_14, 0) * 0.4) * 100
  )

cat(sprintf("  Jointure : %d / %d IRIS avec population\n\n",
            sum(!is.na(iris_rennes$POP_TOTAL)), nrow(iris_rennes)))

# =============================================================================
# 4. POINTS D'INTÉRÊT — OpenStreetMap
# =============================================================================

cat("[4/6] Téléchargement des points d'intérêt (OSM)...\n")

# -----------------------------------------------------------------------
# get_osm_points : pour les POI de type point (pharmacies, médecins, bus)
# Prend uniquement les points OSM + centroïdes polygones si pas de points
# -----------------------------------------------------------------------
get_osm_points <- function(key, value, type_label) {
  Sys.sleep(3)
  
  raw <- tryCatch({
    opq(bbox = BBOX_RENNES, timeout = 60) %>%
      add_osm_feature(key = key, value = value) %>%
      osmdata_sf()
  }, error = function(e) {
    message(sprintf("    Erreur Overpass (%s=%s) : %s", key, value, e$message))
    return(NULL)
  })
  
  empty_sf <- st_sf(osm_id   = character(),
                    name     = character(),
                    type     = character(),
                    geometry = st_sfc(crs = 2154))
  
  if (is.null(raw)) return(empty_sf)
  
  pts <- raw$osm_points
  if (is.null(pts) || nrow(pts) == 0) {
    if (!is.null(raw$osm_polygons) && nrow(raw$osm_polygons) > 0) {
      pts <- raw$osm_polygons %>% suppressWarnings(st_centroid())
    } else return(empty_sf)
  }
  
  pts %>%
    mutate(
      name = if ("name" %in% names(.)) name else NA_character_,
      type = type_label
    ) %>%
    select(osm_id, name, type, geometry) %>%
    st_transform(2154) %>%
    st_filter(rennes_poly) %>%
    distinct(osm_id, .keep_all = TRUE)
}

# -----------------------------------------------------------------------
# get_osm_polygons : pour les grands équipements (hôpitaux)
# Prend multipolygones + polygones nommés — ignore les points individuels
# -----------------------------------------------------------------------
get_osm_polygons <- function(key, value, type_label) {
  Sys.sleep(3)
  
  raw <- tryCatch({
    opq(bbox = BBOX_RENNES, timeout = 60) %>%
      add_osm_feature(key = key, value = value) %>%
      osmdata_sf()
  }, error = function(e) {
    message(sprintf("    Erreur Overpass (%s=%s) : %s", key, value, e$message))
    return(NULL)
  })
  
  empty_sf <- st_sf(osm_id   = character(),
                    name     = character(),
                    type     = character(),
                    geometry = st_sfc(crs = 2154))
  
  if (is.null(raw)) return(empty_sf)
  
  extraire_poly <- function(sf_obj) {
    if (is.null(sf_obj) || nrow(sf_obj) == 0) return(NULL)
    sf_obj %>%
      filter(!is.na(name)) %>%
      suppressWarnings(st_centroid()) %>%
      mutate(type = type_label) %>%
      select(osm_id, name, type, geometry)
  }
  
  result <- bind_rows(
    extraire_poly(raw$osm_multipolygons),
    extraire_poly(raw$osm_polygons)
  )
  
  if (is.null(result) || nrow(result) == 0) return(empty_sf)
  
  result %>%
    st_transform(2154) %>%
    st_filter(rennes_poly) %>%
    distinct(osm_id, .keep_all = TRUE)
}

# --- Pharmacies ---
cat("  Pharmacies...\n")
pharmacies <- get_osm_points("amenity", "pharmacy", "pharmacie")
cat(sprintf("    => %d pharmacies\n", nrow(pharmacies)))

# --- Hôpitaux (multipolygones + polygones nommés uniquement) ---
cat("  Hôpitaux & cliniques...\n")
hopitaux <- get_osm_polygons("amenity", "hospital", "hopital")
cat(sprintf("    => %d hôpitaux/cliniques\n", nrow(hopitaux)))

# --- Médecins ---
cat("  Médecins...\n")
medecins <- get_osm_points("amenity", "doctors", "medecin")
cat(sprintf("    => %d médecins\n", nrow(medecins)))

# --- Arrêts bus ---
cat("  Arrêts bus...\n")
bus <- get_osm_points("highway", "bus_stop", "bus")
cat(sprintf("    => %d arrêts bus\n", nrow(bus)))

# --- Stations métro (points nommés uniquement, dédoublonnés par nom) ---
cat("  Stations métro...\n")
raw_metro <- tryCatch({
  Sys.sleep(3)
  opq(bbox = BBOX_RENNES, timeout = 60) %>%
    add_osm_feature(key = "station", value = "subway") %>%
    osmdata_sf()
}, error = function(e) NULL)

if (!is.null(raw_metro)) {
  metro <- raw_metro$osm_points %>%
    filter(!is.na(name)) %>%
    mutate(type = "metro") %>%
    select(osm_id, name, type, geometry) %>%
    st_transform(2154) %>%
    st_filter(rennes_poly) %>%
    distinct(name, .keep_all = TRUE)   # dédoublonnage par nom de station
} else {
  metro <- st_sf(osm_id = character(), name = character(),
                 type = character(), geometry = st_sfc(crs = 2154))
}
cat(sprintf("    => %d stations métro\n", nrow(metro)))

transport <- bind_rows(bus, metro)
cat(sprintf("    => %d arrêts transport total\n", nrow(transport)))

# --- EHPAD ---
cat("  EHPAD...\n")
ehpad <- get_osm_polygons("amenity", "nursing_home", "ehpad")
cat(sprintf("    => %d EHPAD\n", nrow(ehpad)))

saveRDS(pharmacies, "data/pharmacies.rds")
saveRDS(hopitaux,   "data/hopitaux.rds")
saveRDS(medecins,   "data/medecins.rds")
saveRDS(transport,  "data/transport.rds")
saveRDS(ehpad,      "data/ehpad.rds")

# =============================================================================
# 5. JOINTURE SPATIALE — Comptages par IRIS
# =============================================================================

cat("\n[5/6] Agrégation des POI par IRIS...\n")

# Fonction robuste : initialise à 0, puis met à jour par matching
ajouter_comptage <- function(iris_sf, points_sf, col_name) {
  iris_sf[[col_name]] <- 0L
  
  if (is.null(points_sf) || nrow(points_sf) == 0) return(iris_sf)
  
  counts <- st_join(points_sf, iris_sf["CODE_IRIS"], join = st_within) %>%
    st_drop_geometry() %>%
    filter(!is.na(CODE_IRIS)) %>%
    count(CODE_IRIS, name = col_name)
  
  idx <- match(iris_sf$CODE_IRIS, counts$CODE_IRIS)
  iris_sf[[col_name]][!is.na(idx)] <- counts[[col_name]][idx[!is.na(idx)]]
  
  iris_sf
}

iris_rennes <- iris_rennes %>%
  ajouter_comptage(pharmacies, "nb_pharmacies") %>%
  ajouter_comptage(hopitaux,   "nb_hopitaux")   %>%
  ajouter_comptage(medecins,   "nb_medecins")   %>%
  ajouter_comptage(transport,  "nb_transport")  %>%
  ajouter_comptage(ehpad,      "nb_ehpad")

# Distance minimale à la pharmacie la plus proche (centroïde IRIS → pharmacie)
centroid_iris <- suppressWarnings(st_centroid(iris_rennes))

if (nrow(pharmacies) > 0) {
  dist_matrix <- st_distance(centroid_iris, pharmacies)
  iris_rennes$dist_min_pharmacie_m <- as.numeric(apply(dist_matrix, 1, min))
} else {
  iris_rennes$dist_min_pharmacie_m <- NA_real_
}

# Distance au centre-ville + flags (calculés APRÈS les comptages)
iris_rennes <- iris_rennes %>%
  mutate(
    dist_centre_m = as.numeric(
      st_distance(
        suppressWarnings(st_centroid(geometry)),
        suppressWarnings(st_centroid(st_union(rennes_poly)))
      )
    ),
    zone_sans_pharmacie = as.integer(nb_pharmacies == 0),
    desert_pharma       = as.integer(dist_min_pharmacie_m > 500)
  )

saveRDS(iris_rennes, "data/iris_rennes.rds")

# =============================================================================
# 6. VÉRIFICATION FINALE
# =============================================================================

cat("\n[6/6] Vérification...\n")

vars_check <- c("POP_TOTAL", "PART_65P", "PART_0_14", "DENSITE_POP",
                "nb_pharmacies", "dist_min_pharmacie_m")

for (v in vars_check) {
  vals <- iris_rennes[[v]]
  cat(sprintf("  %-28s | min: %7.1f | max: %8.1f | NA: %d\n",
              v, min(vals, na.rm=TRUE), max(vals, na.rm=TRUE), sum(is.na(vals))))
}

cat("\n========================================\n")
cat("        COLLECTE TERMINÉE\n")
cat("========================================\n")
cat(sprintf("  IRIS Rennes            : %d\n",  nrow(iris_rennes)))
cat(sprintf("  Population totale      : %s hab.\n",
            format(sum(iris_rennes$POP_TOTAL, na.rm=TRUE), big.mark=" ")))
cat(sprintf("  Pharmacies             : %d\n",  nrow(pharmacies)))
cat(sprintf("  Hôpitaux/Cliniques     : %d\n",  nrow(hopitaux)))
cat(sprintf("  Médecins               : %d\n",  nrow(medecins)))
cat(sprintf("  Arrêts transport       : %d\n",  nrow(transport)))
cat(sprintf("  EHPAD                  : %d\n",  nrow(ehpad)))
cat(sprintf("  IRIS sans pharmacie    : %d\n",
            sum(iris_rennes$zone_sans_pharmacie)))
cat(sprintf("  Déserts pharmaceutiques: %d\n",
            sum(iris_rennes$desert_pharma, na.rm=TRUE)))
cat("========================================\n")
