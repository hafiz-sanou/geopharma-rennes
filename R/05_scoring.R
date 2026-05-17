# =============================================================================
# PharmOptim-Rennes | Script 05 — Scoring Composite & Recommandation
# Auteur : Hafiz | Statistique Spatiale — Magistère 2
# =============================================================================
# Objectif (Tâche 3) :
#   - Construire un score composite multicritère par IRIS
#   - Identifier les zones optimales d'implantation
#   - Tableau comparatif IRIS optimal vs reste de Rennes
#   - Justification statistique du choix final
#
# =============================================================================

library(sf)
library(tidyverse)
library(ggplot2)
library(viridis)
library(scales)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/stats",   recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 0. CHARGEMENT DES DONNÉES
# =============================================================================

cat("Chargement des données...\n")

# On lit depuis la version BASE (propre, sans colonnes scores)
iris_rennes <- readRDS("data/iris_rennes_base.rds")
pharmacies  <- readRDS("data/pharmacies.rds")
hopitaux    <- readRDS("data/hopitaux.rds")
rennes_poly <- readRDS("data/rennes_contour.rds")

cat(sprintf("  %d IRIS chargés (%d colonnes)\n", nrow(iris_rennes), ncol(iris_rennes)))

# =============================================================================
# 1. FILTRE POPULATION >= 500 habitants
# =============================================================================

iris_scoring <- iris_rennes %>%
  filter(POP_TOTAL >= 500)

cat(sprintf("  IRIS retenus pour scoring (pop >= 500) : %d / %d\n\n",
            nrow(iris_scoring), nrow(iris_rennes)))

# =============================================================================
# 2. NORMALISATION ET SCORE COMPOSITE
# =============================================================================

cat("=== CONSTRUCTION DU SCORE COMPOSITE ===\n\n")

normalize <- function(x, reverse = FALSE) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  n <- (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 100
  if (reverse) n <- 100 - n
  n
}

iris_scoring <- iris_scoring %>%
  mutate(
    # Scores individuels 0-100
    score_distance  = normalize(dist_min_pharmacie_m),
    score_densite   = normalize(DENSITE_POP),
    score_seniors   = normalize(PART_65P),
    score_medecins  = normalize(nb_medecins),
    score_transport = normalize(nb_transport),
    score_jeunes    = normalize(PART_0_14),
    
    # Score composite pondéré
    # 30% distance | 25% densité | 20% seniors
    # 10% médecins | 10% transport | 5% jeunes
    SCORE_COMPOSITE = (
      score_distance  * 0.30 +
        score_densite   * 0.25 +
        score_seniors   * 0.20 +
        score_medecins  * 0.10 +
        score_transport * 0.10 +
        score_jeunes    * 0.05
    )
  )

cat("Scores calculés :\n")
cat(sprintf("  Score distance  : min=%.1f | max=%.1f | moy=%.1f\n",
            min(iris_scoring$score_distance),
            max(iris_scoring$score_distance),
            mean(iris_scoring$score_distance)))
cat(sprintf("  Score densité   : min=%.1f | max=%.1f | moy=%.1f\n",
            min(iris_scoring$score_densite),
            max(iris_scoring$score_densite),
            mean(iris_scoring$score_densite)))
cat(sprintf("  Score composite : min=%.1f | max=%.1f | moy=%.1f\n\n",
            min(iris_scoring$SCORE_COMPOSITE),
            max(iris_scoring$SCORE_COMPOSITE),
            mean(iris_scoring$SCORE_COMPOSITE)))

# =============================================================================
# 3. CLASSEMENT ET TOP 10
# =============================================================================

cat("=== TOP ZONES CANDIDATES ===\n\n")

iris_ranked <- iris_scoring %>%
  st_drop_geometry() %>%
  arrange(desc(SCORE_COMPOSITE)) %>%
  mutate(rang = row_number())

top10 <- iris_ranked %>%
  select(rang, CODE_IRIS, NOM_IRIS,
         SCORE_COMPOSITE, score_distance, score_densite,
         score_seniors, POP_TOTAL, PART_65P,
         dist_min_pharmacie_m, nb_medecins, nb_transport,
         desert_pharma, zone_sans_pharmacie) %>%
  head(10) %>%
  mutate(
    SCORE_COMPOSITE      = round(SCORE_COMPOSITE, 1),
    score_distance       = round(score_distance, 1),
    score_densite        = round(score_densite, 1),
    score_seniors        = round(score_seniors, 1),
    PART_65P             = round(PART_65P * 100, 1),
    dist_min_pharmacie_m = round(dist_min_pharmacie_m, 0)
  )

cat("Top 10 IRIS prioritaires :\n")
print(top10 %>% select(rang, NOM_IRIS, SCORE_COMPOSITE,
                       dist_min_pharmacie_m, POP_TOTAL, PART_65P))

saveRDS(top10, "outputs/stats/top10_iris.rds")
write.csv(top10, "outputs/stats/top10_iris.csv", row.names = FALSE)

# =============================================================================
# 4. IRIS OPTIMAL
# =============================================================================

cat("\n=== ANALYSE DÉTAILLÉE — IRIS OPTIMAL ===\n\n")

iris_optimal <- iris_scoring %>%
  arrange(desc(SCORE_COMPOSITE)) %>%
  slice(1)

cat(sprintf("IRIS retenu : %s (CODE : %s)\n",
            iris_optimal$NOM_IRIS, iris_optimal$CODE_IRIS))
cat(sprintf("  Score composite      : %.1f / 100\n", iris_optimal$SCORE_COMPOSITE))
cat(sprintf("  Distance pharmacie   : %.0f m\n",     iris_optimal$dist_min_pharmacie_m))
cat(sprintf("  Population totale    : %.0f hab.\n",  iris_optimal$POP_TOTAL))
cat(sprintf("  Part 65 ans et +     : %.1f%%\n",     iris_optimal$PART_65P * 100))
cat(sprintf("  Part 0-14 ans        : %.1f%%\n",     iris_optimal$PART_0_14 * 100))
cat(sprintf("  Densité population   : %.0f hab/km²\n", iris_optimal$DENSITE_POP))
cat(sprintf("  Médecins à proximité : %d\n",         iris_optimal$nb_medecins))
cat(sprintf("  Arrêts transport     : %d\n",         iris_optimal$nb_transport))
cat(sprintf("  Désert pharmaceutique: %s\n",
            ifelse(iris_optimal$desert_pharma == 1, "OUI ✓", "NON")))

# =============================================================================
# 5. TABLEAU COMPARATIF
# =============================================================================

cat("\n=== TABLEAU COMPARATIF (Tâche 3 du prof) ===\n\n")

autres_iris <- iris_scoring %>%
  st_drop_geometry() %>%
  filter(CODE_IRIS != iris_optimal$CODE_IRIS)

vars_comparaison <- c(
  "POP_TOTAL", "DENSITE_POP", "PART_65P", "PART_0_14",
  "dist_min_pharmacie_m", "nb_medecins", "nb_transport",
  "SCORE_COMPOSITE"
)

labels_vars <- c(
  "Population totale",
  "Densité (hab/km²)",
  "Part 65+ (%)",
  "Part 0-14 (%)",
  "Distance pharmacie (m)",
  "Nb médecins IRIS",
  "Nb arrêts transport",
  "Score composite"
)

tableau_comp <- data.frame(
  Variable     = labels_vars,
  IRIS_optimal = NA_real_,
  Moyenne      = NA_real_,
  Mediane      = NA_real_,
  Min          = NA_real_,
  Max          = NA_real_,
  Ecart_type   = NA_real_
)

for (i in seq_along(vars_comparaison)) {
  v       <- vars_comparaison[i]
  val_opt <- as.numeric(iris_optimal[[v]])
  vals    <- as.numeric(autres_iris[[v]])
  vals    <- vals[!is.na(vals)]
  if (v %in% c("PART_65P", "PART_0_14")) {
    val_opt <- val_opt * 100
    vals    <- vals * 100
  }
  tableau_comp[i, "IRIS_optimal"] <- round(val_opt, 2)
  tableau_comp[i, "Moyenne"]      <- round(mean(vals),   2)
  tableau_comp[i, "Mediane"]      <- round(median(vals), 2)
  tableau_comp[i, "Min"]          <- round(min(vals),    2)
  tableau_comp[i, "Max"]          <- round(max(vals),    2)
  tableau_comp[i, "Ecart_type"]   <- round(sd(vals),     2)
}

cat("Tableau comparatif :\n\n")
print(tableau_comp, row.names = FALSE)

saveRDS(tableau_comp, "outputs/stats/tableau_comparatif.rds")
write.csv(tableau_comp, "outputs/stats/tableau_comparatif.csv", row.names = FALSE)

# =============================================================================
# 6. JOINTURE SCORES → iris_rennes (pour les cartes)
# =============================================================================

scores_a_joindre <- iris_scoring %>%
  st_drop_geometry() %>%
  select(CODE_IRIS, score_distance, score_densite, score_seniors,
         score_medecins, score_transport, score_jeunes, SCORE_COMPOSITE)

iris_rennes <- iris_rennes %>%
  left_join(scores_a_joindre, by = "CODE_IRIS")

cat(sprintf("\nJointure OK — SCORE_COMPOSITE présent : %s\n",
            "SCORE_COMPOSITE" %in% names(iris_rennes)))

# Codes top 3
top3_codes <- iris_ranked %>% head(3) %>% pull(CODE_IRIS)

# =============================================================================
# 7. CARTE SCORE COMPOSITE
# =============================================================================

cat("\n[1/3] Carte score composite...\n")

p_score <- ggplot(iris_rennes) +
  geom_sf(aes(fill = SCORE_COMPOSITE),
          color = "white", linewidth = 0.3) +
  geom_sf(
    data     = iris_rennes %>% filter(CODE_IRIS %in% top3_codes),
    fill     = NA,
    color    = "#ffd166",
    linewidth = 1.5
  ) +
  geom_sf(
    data         = pharmacies,
    color        = "white",
    size         = 1,
    shape        = 3,
    inherit.aes  = FALSE
  ) +
  geom_sf_label(
    data = iris_rennes %>%
      filter(CODE_IRIS %in% top3_codes) %>%
      left_join(iris_ranked %>% select(CODE_IRIS, rang), by = "CODE_IRIS"),
    aes(label = paste0("#", rang)),
    size          = 3,
    fontface      = "bold",
    color         = "#1a1a2e",
    fill          = "#ffd166",
    label.padding = unit(0.15, "cm")
  ) +
  scale_fill_viridis_c(
    option    = "rocket",
    direction = -1,
    name      = "Score (0-100)",
    na.value  = "grey85"
  ) +
  labs(
    title    = "Score composite d'opportunité pharmaceutique par IRIS",
    subtitle = "Contours jaunes = Top 3 zones candidates | Croix = pharmacies existantes",
    caption  = paste0(
      "Score = 30% distance + 25% densité + 20% seniors + ",
      "10% médecins + 10% transport + 5% jeunes\n",
      "IRIS pop < 500 hab exclus du scoring | ",
      "Sources : INSEE RP 2022 | IGN 2024 | OSM 2024"
    )
  ) +
  theme_void() +
  theme(
    plot.title       = element_text(face = "bold", size = 13, color = "#1a1a2e"),
    plot.subtitle    = element_text(size = 9,  color = "#555555"),
    plot.caption     = element_text(size = 7,  color = "#888888"),
    legend.position  = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.background  = element_rect(fill = "white", color = NA)
  )

ggsave("outputs/figures/17_score_composite.png", p_score,
       width = 9, height = 8, dpi = 200)
cat("  Sauvegardé : 17_score_composite.png\n")

# =============================================================================
# 8. GRAPHIQUE PROFIL IRIS OPTIMAL
# =============================================================================

cat("[2/3] Graphique profil IRIS optimal...\n")

df_radar <- data.frame(
  critere = c(
    "Distance\npharmacie", "Densité\npopulation",
    "Population\nseniors",  "Présence\nmédecins",
    "Accessibilité\ntransport", "Population\njeune"
  ),
  iris_optimal = c(
    iris_optimal$score_distance,
    iris_optimal$score_densite,
    iris_optimal$score_seniors,
    iris_optimal$score_medecins,
    iris_optimal$score_transport,
    iris_optimal$score_jeunes
  ),
  moyenne_rennes = c(
    mean(iris_scoring$score_distance,  na.rm = TRUE),
    mean(iris_scoring$score_densite,   na.rm = TRUE),
    mean(iris_scoring$score_seniors,   na.rm = TRUE),
    mean(iris_scoring$score_medecins,  na.rm = TRUE),
    mean(iris_scoring$score_transport, na.rm = TRUE),
    mean(iris_scoring$score_jeunes,    na.rm = TRUE)
  )
)

nom_optimal <- iris_optimal$NOM_IRIS

df_radar_long <- df_radar %>%
  pivot_longer(
    cols      = c(iris_optimal, moyenne_rennes),
    names_to  = "groupe",
    values_to = "score"
  ) %>%
  mutate(
    groupe = case_when(
      groupe == "iris_optimal"   ~ paste0("IRIS optimal (", nom_optimal, ")"),
      groupe == "moyenne_rennes" ~ "Moyenne Rennes (IRIS scorés)",
      TRUE                       ~ groupe
    )
  )

p_radar <- ggplot(df_radar_long,
                  aes(x = critere, y = score, fill = groupe)) +
  geom_col(position = "dodge", width = 0.6, alpha = 0.9) +
  geom_hline(yintercept = 50, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  scale_fill_manual(
    name   = "",
    values = c("#e63946", "#adb5bd")
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    labels = function(x) paste0(x, "/100")
  ) +
  labs(
    title    = "Profil de l'IRIS optimal vs moyenne de Rennes",
    subtitle = "Score normalisé 0-100 par critère",
    x        = NULL,
    y        = "Score (0-100)",
    caption  = "Ligne pointillée = score médian (50/100)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13, color = "#1a1a2e"),
    plot.subtitle   = element_text(size = 9, color = "#555555"),
    plot.caption    = element_text(size = 7, color = "#888888"),
    legend.position = "top",
    axis.text.x     = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("outputs/figures/18_profil_iris_optimal.png", p_radar,
       width = 10, height = 6, dpi = 200)
cat("  Sauvegardé : 18_profil_iris_optimal.png\n")

# =============================================================================
# 9. CARTE TOP 3 DÉTAILLÉE
# =============================================================================

cat("[3/3] Carte Top 3 zones candidates...\n")

top3_sf <- iris_rennes %>%
  filter(CODE_IRIS %in% top3_codes) %>%
  left_join(iris_ranked %>% select(CODE_IRIS, rang), by = "CODE_IRIS")

p_top3 <- ggplot() +
  # Fond : tous les IRIS gris clair
  geom_sf(
    data     = iris_rennes,
    fill     = "#f8f9fa",
    color    = "white",
    linewidth = 0.2
  ) +
  # Déserts pharmaceutiques en orange clair
  geom_sf(
    data     = iris_rennes %>% filter(desert_pharma == 1),
    fill     = "#ffd166",
    color    = "white",
    linewidth = 0.2,
    alpha    = 0.6
  ) +
  # Top 3 colorés par rang
  geom_sf(
    data     = top3_sf,
    aes(fill = factor(rang)),
    color    = "white",
    linewidth = 0.5
  ) +
  # Pharmacies
  geom_sf(
    data        = pharmacies,
    color       = "#1d3557",
    size        = 1.2,
    shape       = 3,
    inherit.aes = FALSE
  ) +
  # Hôpitaux
  geom_sf(
    data        = hopitaux,
    color       = "#2d6a4f",
    size        = 3,
    shape       = 15,
    inherit.aes = FALSE
  ) +
  # Labels top 3
  geom_sf_label(
    data          = top3_sf,
    aes(label     = paste0("#", rang, "\n", NOM_IRIS)),
    size          = 2.8,
    fontface      = "bold",
    color         = "white",
    fill          = "#1a1a2e",
    label.padding = unit(0.2, "cm"),
    nudge_y       = 200
  ) +
  scale_fill_manual(
    name   = "Rang",
    values = c("1" = "#e63946", "2" = "#f4a261", "3" = "#2a9d8f"),
    labels = c(
      "1" = "#1 — Priorité haute",
      "2" = "#2 — Priorité moyenne",
      "3" = "#3 — À considérer"
    )
  ) +
  labs(
    title    = "Top 3 zones d'implantation recommandées",
    subtitle = "Orange = déserts pharmaceutiques | Croix = pharmacies | Carrés verts = hôpitaux",
    caption  = "PharmOptim-Rennes | Sources : INSEE 2022, IGN 2024, OSM 2024"
  ) +
  theme_void() +
  theme(
    plot.title      = element_text(face = "bold", size = 13, color = "#1a1a2e"),
    plot.subtitle   = element_text(size = 9,  color = "#555555"),
    plot.caption    = element_text(size = 7,  color = "#888888"),
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("outputs/figures/19_top3_zones.png", p_top3,
       width = 9, height = 8, dpi = 200)
cat("  Sauvegardé : 19_top3_zones.png\n\n")

# =============================================================================
# 10. SAUVEGARDE FINALE
# =============================================================================

# Sauvegarder iris_rennes enrichi (avec scores) pour le rapport RMD
saveRDS(iris_rennes, "data/iris_rennes.rds")

# Sauvegarder tous les résultats du scoring
saveRDS(
  list(
    iris_optimal = iris_optimal,
    top10        = top10,
    tableau_comp = tableau_comp,
    df_radar     = df_radar,
    iris_ranked  = iris_ranked
  ),
  "outputs/stats/resultats_scoring.rds"
)

cat("=== SCORING TERMINÉ ===\n")
cat("Figures :\n")
cat("  outputs/figures/17_score_composite.png\n")
cat("  outputs/figures/18_profil_iris_optimal.png\n")
cat("  outputs/figures/19_top3_zones.png\n")
cat("Stats :\n")
cat("  outputs/stats/top10_iris.rds + .csv\n")
cat("  outputs/stats/tableau_comparatif.rds + .csv\n")
cat("  outputs/stats/resultats_scoring.rds\n")
cat("Données :\n")
cat("  data/iris_rennes.rds (enrichi avec scores)\n")
cat("\n")
cat("==============================================\n")
cat("  IRIS OPTIMAL RECOMMANDÉ :\n")
cat(sprintf("  %s — Score : %.1f/100\n",
            iris_optimal$NOM_IRIS,
            iris_optimal$SCORE_COMPOSITE))
cat(sprintf("  Population : %.0f hab. | Distance : %.0f m\n",
            iris_optimal$POP_TOTAL,
            iris_optimal$dist_min_pharmacie_m))
cat("==============================================\n")