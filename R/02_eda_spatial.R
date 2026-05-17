# =============================================================================
# PharmOptim-Rennes | Script 02 — Analyse Exploratoire Spatiale (EDA)
# Auteur : Hafiz | Statistique Spatiale — Magistère 2
# =============================================================================
# Objectif : Explorer la distribution spatiale des variables clés
#   - Cartes choroplèthes (population, densité, vulnérabilité)
#   - Distribution des pharmacies existantes
#   - Identification visuelle des déserts pharmaceutiques
#   - Statistiques descriptives par IRIS
# =============================================================================

library(sf)
library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 0. CHARGEMENT DES DONNÉES
# =============================================================================

cat("Chargement des données...\n")

iris_rennes  <- readRDS("data/iris_rennes.rds")
pharmacies   <- readRDS("data/pharmacies.rds")
hopitaux     <- readRDS("data/hopitaux.rds")
medecins     <- readRDS("data/medecins.rds")
transport    <- readRDS("data/transport.rds")
rennes_poly  <- readRDS("data/rennes_contour.rds")

cat(sprintf("  %d IRIS chargés\n", nrow(iris_rennes)))
cat(sprintf("  %d pharmacies chargées\n\n", nrow(pharmacies)))

# =============================================================================
# 1. STATISTIQUES DESCRIPTIVES GLOBALES
# =============================================================================

cat("=== STATISTIQUES DESCRIPTIVES ===\n\n")

cat("--- Population par IRIS ---\n")
cat(sprintf("  Total        : %s hab.\n",
            format(sum(iris_rennes$POP_TOTAL, na.rm=TRUE), big.mark=" ")))
cat(sprintf("  Moyenne/IRIS : %.0f hab.\n",
            mean(iris_rennes$POP_TOTAL, na.rm=TRUE)))
cat(sprintf("  Min          : %.0f hab.\n",
            min(iris_rennes$POP_TOTAL, na.rm=TRUE)))
cat(sprintf("  Max          : %.0f hab.\n",
            max(iris_rennes$POP_TOTAL, na.rm=TRUE)))
cat(sprintf("  Médiane      : %.0f hab.\n",
            median(iris_rennes$POP_TOTAL, na.rm=TRUE)))

cat("\n--- Densité de population ---\n")
cat(sprintf("  Moyenne      : %.0f hab/km²\n",
            mean(iris_rennes$DENSITE_POP, na.rm=TRUE)))
cat(sprintf("  Max          : %.0f hab/km²\n",
            max(iris_rennes$DENSITE_POP, na.rm=TRUE)))

cat("\n--- Pharmacies ---\n")
cat(sprintf("  Total                  : %d\n", nrow(pharmacies)))
cat(sprintf("  IRIS avec pharmacie    : %d / %d\n",
            sum(iris_rennes$nb_pharmacies > 0), nrow(iris_rennes)))
cat(sprintf("  IRIS sans pharmacie    : %d / %d (%.1f%%)\n",
            sum(iris_rennes$zone_sans_pharmacie),
            nrow(iris_rennes),
            100 * mean(iris_rennes$zone_sans_pharmacie)))
cat(sprintf("  Déserts pharmaceutiques: %d / %d (%.1f%%)\n",
            sum(iris_rennes$desert_pharma, na.rm=TRUE),
            nrow(iris_rennes),
            100 * mean(iris_rennes$desert_pharma, na.rm=TRUE)))
cat(sprintf("  Distance moyenne/pharma: %.0f m\n",
            mean(iris_rennes$dist_min_pharmacie_m, na.rm=TRUE)))
cat(sprintf("  Distance max/pharma    : %.0f m\n",
            max(iris_rennes$dist_min_pharmacie_m, na.rm=TRUE)))

cat("\n--- Population en désert pharmaceutique ---\n")
pop_desert <- iris_rennes %>%
  filter(desert_pharma == 1) %>%
  pull(POP_TOTAL) %>%
  sum(na.rm=TRUE)
cat(sprintf("  Population concernée : %s hab. (%.1f%% de Rennes)\n",
            format(pop_desert, big.mark=" "),
            100 * pop_desert / sum(iris_rennes$POP_TOTAL, na.rm=TRUE)))

# =============================================================================
# 2. THÈME GRAPHIQUE COMMUN
# =============================================================================

theme_pharma <- function() {
  theme_void() +
    theme(
      plot.title    = element_text(size = 13, face = "bold",
                                   margin = margin(b=5), color = "#1a1a2e"),
      plot.subtitle = element_text(size = 9, color = "#555555",
                                   margin = margin(b=10)),
      plot.caption  = element_text(size = 7, color = "#888888",
                                   margin = margin(t=8)),
      legend.position  = "bottom",
      legend.title     = element_text(size = 8, face = "bold"),
      legend.text      = element_text(size = 7),
      legend.key.width = unit(1.5, "cm"),
      legend.key.height= unit(0.3, "cm"),
      plot.margin      = margin(10, 10, 10, 10),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

# =============================================================================
# 3. CARTE 1 — DENSITÉ DE POPULATION PAR IRIS
# =============================================================================

cat("\n[1/6] Carte densité de population...\n")

p1 <- ggplot(iris_rennes) +
  geom_sf(aes(fill = DENSITE_POP), color = "white", linewidth = 0.2) +
  geom_sf(data = pharmacies, color = "red", size = 1,
          alpha = 0.7, shape = 16) +
  scale_fill_viridis_c(
    option   = "plasma",
    name     = "Habitants/km²",
    labels   = scales::comma_format(big.mark = " "),
    na.value = "grey90"
  ) +
  labs(
    title    = "Densité de population par IRIS",
    subtitle = "Points rouges = pharmacies existantes",
    caption  = "Sources : INSEE RP 2022 | OSM 2024"
  ) +
  theme_pharma()

# =============================================================================
# 4. CARTE 2 — PART DES 65 ANS ET PLUS
# =============================================================================

cat("[2/6] Carte population senior...\n")

p2 <- ggplot(iris_rennes) +
  geom_sf(aes(fill = PART_65P * 100), color = "white", linewidth = 0.2) +
  geom_sf(data = pharmacies, color = "red", size = 1,
          alpha = 0.8, shape = 16) +
  scale_fill_viridis_c(
    option   = "mako",
    name     = "Part 65+ (%)",
    direction = -1,
    na.value = "grey90"
  ) +
  labs(
    title    = "Part des 65 ans et plus par IRIS",
    subtitle = "Points verts = pharmacies existantes",
    caption  = "Sources : INSEE RP 2022 | OSM 2024"
  ) +
  theme_pharma()

# =============================================================================
# 5. CARTE 3 — DISTANCE À LA PHARMACIE LA PLUS PROCHE
# =============================================================================

cat("[3/6] Carte distance aux pharmacies...\n")

p3 <- ggplot(iris_rennes) +
  geom_sf(aes(fill = dist_min_pharmacie_m), color = "white", linewidth = 0.2) +
  geom_sf(data = pharmacies, color = "red", size = 1.5,
          alpha = 0.9, shape = 3) +
  scale_fill_viridis_c(
    option   = "inferno",
    name     = "Distance (m)",
    direction = -1,
    na.value = "grey90"
  ) +
  labs(
    title    = "Distance à la pharmacie la plus proche",
    subtitle = "Croix blanches = pharmacies | Zones sombres = zones bien desservies",
    caption  = "Sources : IGN 2024 | OSM 2024"
  ) +
  theme_pharma()

# =============================================================================
# 6. CARTE 4 — DÉSERTS PHARMACEUTIQUES
# =============================================================================

cat("[4/6] Carte déserts pharmaceutiques...\n")

iris_rennes <- iris_rennes %>%
  mutate(
    statut_pharma = case_when(
      nb_pharmacies > 0              ~ "IRIS avec pharmacie",
      desert_pharma == 1             ~ "Désert pharmaceutique (>500m)",
      zone_sans_pharmacie == 1       ~ "Sans pharmacie (<500m d'une autre)",
      TRUE                           ~ "Autre"
    ),
    statut_pharma = factor(statut_pharma, levels = c(
      "IRIS avec pharmacie",
      "Sans pharmacie (<500m d'une autre)",
      "Désert pharmaceutique (>500m)"
    ))
  )

p4 <- ggplot(iris_rennes) +
  geom_sf(aes(fill = statut_pharma), color = "white", linewidth = 0.3) +
  geom_sf(data = pharmacies, color = "#1d3557", size = 1.2,
          alpha = 0.9, shape = 16) +
  scale_fill_manual(
    name   = "Statut",
    values = c(
      "IRIS avec pharmacie"                  = "#52b788",
      "Sans pharmacie (<500m d'une autre)"   = "#ffd166",
      "Désert pharmaceutique (>500m)"        = "#e63946"
    ),
    na.value = "grey90"
  ) +
  labs(
    title    = "Déserts pharmaceutiques à Rennes",
    subtitle = "Points bleus = pharmacies existantes",
    caption  = "Seuil désert : distance > 500m à la pharmacie la plus proche"
  ) +
  theme_pharma() +
  theme(legend.position = "bottom",
        legend.direction = "vertical")

# =============================================================================
# 7. CARTE 5 — SCORE DE VULNÉRABILITÉ
# =============================================================================

cat("[5/6] Carte score de vulnérabilité...\n")

p5 <- ggplot(iris_rennes) +
  geom_sf(aes(fill = SCORE_VULNERABILITE), color = "white", linewidth = 0.2) +
  geom_sf(data = iris_rennes %>% filter(desert_pharma == 1),
          fill = NA, color = "#e63946", linewidth = 0.8) +
  scale_fill_viridis_c(
    option   = "rocket",
    name     = "Score (0-100)",
    direction = -1,
    na.value = "grey90"
  ) +
  labs(
    title    = "Score de vulnérabilité par IRIS",
    subtitle = "Contours rouges = déserts pharmaceutiques",
    caption  = "Score = 60% part 65+ | 40% part 0-14 ans"
  ) +
  theme_pharma()

# =============================================================================
# 8. GRAPHIQUE 6 — DISTRIBUTION DES VARIABLES CLÉS
# =============================================================================

cat("[6/6] Graphiques de distribution...\n")

# Distribution de la distance aux pharmacies
p6a <- ggplot(iris_rennes, aes(x = dist_min_pharmacie_m)) +
  geom_histogram(aes(fill = after_stat(x) > 500),
                 bins = 25, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 500, color = "#e63946",
             linewidth = 1, linetype = "dashed") +
  scale_fill_manual(values = c("FALSE" = "#52b788", "TRUE" = "#e63946"),
                    guide = "none") +
  annotate("text", x = 520, y = Inf, label = "Seuil désert\n500m",
           hjust = 0, vjust = 1.5, size = 3, color = "#e63946") +
  scale_x_continuous(labels = scales::comma_format(suffix = "m")) +
  labs(
    title = "Distribution des distances à la pharmacie",
    x = "Distance (m)", y = "Nombre d'IRIS"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11))

# Scatter : densité population vs distance pharmacie
p6b <- ggplot(iris_rennes,
              aes(x = DENSITE_POP, y = dist_min_pharmacie_m,
                  color = statut_pharma, size = POP_TOTAL)) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = 500, color = "#e63946",
             linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(
    name   = "Statut",
    values = c(
      "IRIS avec pharmacie"                = "#52b788",
      "Sans pharmacie (<500m d'une autre)" = "#ffd166",
      "Désert pharmaceutique (>500m)"      = "#e63946"
    )
  ) +
  scale_size_continuous(name = "Population", range = c(1, 6),
                        labels = scales::comma_format()) +
  scale_x_continuous(labels = scales::comma_format(suffix = " hab/km²")) +
  scale_y_continuous(labels = scales::comma_format(suffix = "m")) +
  labs(
    title = "Densité de population vs distance à la pharmacie",
    x = "Densité de population (hab/km²)",
    y = "Distance min. à une pharmacie (m)"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        legend.direction = "vertical")

# =============================================================================
# 9. SAUVEGARDE DES FIGURES
# =============================================================================

cat("\nSauvegarde des figures...\n")

ggsave("outputs/figures/01_densite_population.png",  p1, width=8, height=7, dpi=200)
ggsave("outputs/figures/02_population_senior.png",   p2, width=8, height=7, dpi=200)
ggsave("outputs/figures/03_distance_pharmacie.png",  p3, width=8, height=7, dpi=200)
ggsave("outputs/figures/04_deserts_pharma.png",      p4, width=8, height=7, dpi=200)
ggsave("outputs/figures/05_score_vulnerabilite.png", p5, width=8, height=7, dpi=200)
ggsave("outputs/figures/06a_distribution_distance.png", p6a, width=7, height=5, dpi=200)
ggsave("outputs/figures/06b_scatter_densite.png",    p6b, width=8, height=6, dpi=200)

# Planche combinée (cartes 1-4)
combined <- (p1 | p2) / (p3 | p4)
ggsave("outputs/figures/00_planche_eda.png", combined,
       width=14, height=12, dpi=200)

cat("  OK — figures sauvegardées dans outputs/figures/\n")

# =============================================================================
# 10. TABLEAU RÉCAPITULATIF — TOP IRIS PRIORITAIRES
# =============================================================================

cat("\n=== TOP 10 IRIS PRIORITAIRES (déserts pharmaceutiques) ===\n")

top_iris <- iris_rennes %>%
  st_drop_geometry() %>%
  filter(desert_pharma == 1) %>%
  arrange(desc(POP_TOTAL)) %>%
  select(CODE_IRIS, NOM_IRIS, POP_TOTAL, PART_65P,
         dist_min_pharmacie_m, nb_medecins, nb_transport) %>%
  mutate(
    PART_65P             = round(PART_65P * 100, 1),
    dist_min_pharmacie_m = round(dist_min_pharmacie_m, 0)
  ) %>%
  head(10)

print(top_iris)

saveRDS(iris_rennes, "data/iris_rennes.rds")

cat("\n=== EDA TERMINÉE ===\n")
cat("Figures disponibles dans outputs/figures/\n")
cat("Données enrichies sauvegardées dans data/iris_rennes.rds\n")