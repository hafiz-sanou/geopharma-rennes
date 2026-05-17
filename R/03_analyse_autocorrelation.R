# =============================================================================
# PharmOptim-Rennes | Script 03 — Autocorrélation Spatiale
# Auteur : Hafiz | Statistique Spatiale — Magistère 2
# =============================================================================
# Méthodes (Chapitre 2 du cours) :
#   - Matrice de pondération spatiale (Queen, Rook, kNN)
#   - Indice de Moran global + test de significativité
#   - Diagramme de Moran (HH / LL / HL / LH)
#   - LISA — I de Moran local + carte des clusters significatifs
#   - Discussion MAUP (effet de l'échelle)
# Variables analysées :
#   - Densité de population
#   - Part des 65 ans et +
#   - Distance à la pharmacie la plus proche
#   - Nombre de pharmacies par IRIS
# =============================================================================

library(sf)
library(tidyverse)
library(spdep)
library(ggplot2)
library(viridis)
library(patchwork)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/stats",   recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 0. CHARGEMENT DES DONNÉES
# =============================================================================

cat("Chargement des données...\n")
iris_rennes <- readRDS("data/iris_rennes.rds")
cat(sprintf("  %d IRIS chargés\n\n", nrow(iris_rennes)))

# Thème graphique
theme_moran <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 13, color = "#1a1a2e"),
      plot.subtitle = element_text(size = 9,  color = "#555555"),
      plot.caption  = element_text(size = 7,  color = "#888888"),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# =============================================================================
# 1. MATRICES DE PONDÉRATION SPATIALE
# =============================================================================

cat("=== MATRICES DE PONDÉRATION SPATIALE ===\n\n")

# Voisinage Queen (partage d'un côté ou d'un sommet)
nb_queen <- poly2nb(iris_rennes, queen = TRUE)

# Voisinage Rook (partage d'un côté uniquement)
nb_rook  <- poly2nb(iris_rennes, queen = FALSE)

# Voisinage kNN (k plus proches voisins)
coords   <- st_coordinates(st_centroid(iris_rennes) %>% suppressWarnings())
nb_knn3  <- knearneigh(coords, k = 3) %>% knn2nb()
nb_knn5  <- knearneigh(coords, k = 5) %>% knn2nb()

# Matrices W normalisées (style W = row-standardized)
W_queen <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)
W_rook  <- nb2listw(nb_rook,  style = "W", zero.policy = TRUE)
W_knn3  <- nb2listw(nb_knn3,  style = "W", zero.policy = TRUE)
W_knn5  <- nb2listw(nb_knn5,  style = "W", zero.policy = TRUE)

cat("Résumé des voisinages :\n")
cat(sprintf("  Queen : %.1f voisins en moyenne\n", mean(card(nb_queen))))
cat(sprintf("  Rook  : %.1f voisins en moyenne\n", mean(card(nb_rook))))
cat(sprintf("  kNN-3 : %.1f voisins en moyenne\n", mean(card(nb_knn3))))
cat(sprintf("  kNN-5 : %.1f voisins en moyenne\n\n", mean(card(nb_knn5))))

# =============================================================================
# 2. INDICE DE MORAN GLOBAL
# =============================================================================

cat("=== INDICE DE MORAN GLOBAL ===\n\n")

# Variables à analyser
variables <- list(
  densite_pop          = iris_rennes$DENSITE_POP,
  part_65plus          = iris_rennes$PART_65P,
  dist_min_pharmacie   = iris_rennes$dist_min_pharmacie_m,
  nb_pharmacies        = as.numeric(iris_rennes$nb_pharmacies)
)

# Calcul du I de Moran pour chaque variable et chaque matrice W
resultats_moran <- list()

for (var_name in names(variables)) {
  y <- variables[[var_name]]

  # Remplacer les NA par la moyenne (nécessaire pour spdep)
  y[is.na(y)] <- mean(y, na.rm = TRUE)

  cat(sprintf("--- Variable : %s ---\n", var_name))

  for (w_name in c("Queen", "Rook", "kNN-3", "kNN-5")) {
    W <- switch(w_name,
      "Queen" = W_queen, "Rook" = W_rook,
      "kNN-3" = W_knn3,  "kNN-5" = W_knn5
    )

    # Test de Moran (hypothèse de randomisation)
    moran_test <- moran.test(y, W, randomisation = TRUE, zero.policy = TRUE)

    # Test de Moran par permutations Monte Carlo (999 simulations)
    moran_mc <- moran.mc(y, W, nsim = 999, zero.policy = TRUE)

    resultats_moran[[paste(var_name, w_name, sep = "_")]] <- list(
      variable  = var_name,
      voisinage = w_name,
      I         = moran_test$estimate["Moran I statistic"],
      E_I       = moran_test$estimate["Expectation"],
      Var_I     = moran_test$estimate["Variance"],
      p_value   = moran_test$p.value,
      p_mc      = moran_mc$p.value,
      significant = moran_test$p.value < 0.05
    )

    cat(sprintf("  %-8s | I = %+.4f | E[I] = %.4f | p = %.4f %s\n",
                w_name,
                moran_test$estimate["Moran I statistic"],
                moran_test$estimate["Expectation"],
                moran_test$p.value,
                ifelse(moran_test$p.value < 0.05, "***", "")))
  }
  cat("\n")
}

# Tableau récapitulatif
df_moran <- bind_rows(lapply(resultats_moran, as.data.frame)) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

saveRDS(df_moran,       "outputs/stats/moran_global.rds")
write.csv(df_moran,     "outputs/stats/moran_global.csv", row.names = FALSE)

# =============================================================================
# 3. DIAGRAMME DE MORAN — Variable principale : distance à la pharmacie
# =============================================================================

cat("=== DIAGRAMME DE MORAN ===\n\n")

# On travaille sur la distance min pharmacie (variable cible principale)
y_raw <- iris_rennes$dist_min_pharmacie_m
y_raw[is.na(y_raw)] <- mean(y_raw, na.rm = TRUE)

# Standardisation
y_std <- as.numeric(scale(y_raw))

# Lag spatial (moyenne des voisins)
y_lag <- lag.listw(W_queen, y_std, zero.policy = TRUE)

# Quadrants Moran
quadrant <- case_when(
  y_std >= 0 & y_lag >= 0 ~ "HH — Autocorrélation positive",
  y_std <  0 & y_lag <  0 ~ "LL — Autocorrélation positive",
  y_std >= 0 & y_lag <  0 ~ "HL — Autocorrélation négative",
  y_std <  0 & y_lag >= 0 ~ "LH — Autocorrélation négative",
  TRUE                     ~ "NS"
)

df_moran_plot <- data.frame(
  y_std    = y_std,
  y_lag    = y_lag,
  quadrant = quadrant,
  CODE_IRIS = iris_rennes$CODE_IRIS,
  NOM_IRIS  = iris_rennes$NOM_IRIS
)

# I de Moran global (pente de la droite)
moran_dist <- moran.test(y_raw, W_queen, randomisation = TRUE, zero.policy = TRUE)
I_val <- round(moran_dist$estimate["Moran I statistic"], 3)
p_val <- round(moran_dist$p.value, 4)

cat(sprintf("Distance pharmacie — I de Moran (Queen) : %.3f (p = %.4f)\n", I_val, p_val))

# Plot diagramme de Moran
p_moran_diag <- ggplot(df_moran_plot, aes(x = y_std, y = y_lag)) +
  # Lignes de quadrants
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5) +
  # Points colorés par quadrant
  geom_point(aes(color = quadrant), size = 2.5, alpha = 0.8) +
  # Droite de régression (pente = I de Moran)
  geom_smooth(method = "lm", se = TRUE, color = "#e63946",
              linewidth = 1, linetype = "solid") +
  # Labels quadrants
  annotate("text", x =  1.8, y =  2.2, label = "HH", fontface="bold",
           color = "#e63946", size = 5) +
  annotate("text", x = -1.8, y = -2.2, label = "LL", fontface="bold",
           color = "#457b9d", size = 5) +
  annotate("text", x =  1.8, y = -2.2, label = "HL", fontface="bold",
           color = "#f4a261", size = 5) +
  annotate("text", x = -1.8, y =  2.2, label = "LH", fontface="bold",
           color = "#2a9d8f", size = 5) +
  scale_color_manual(
    name = "Quadrant",
    values = c(
      "HH — Autocorrélation positive" = "#e63946",
      "LL — Autocorrélation positive" = "#457b9d",
      "HL — Autocorrélation négative" = "#f4a261",
      "LH — Autocorrélation négative" = "#2a9d8f"
    )
  ) +
  labs(
    title    = "Diagramme de Moran — Distance à la pharmacie la plus proche",
    subtitle = sprintf("I de Moran (Queen) = %.3f | p-value = %.4f | n = %d IRIS",
                       I_val, p_val, nrow(iris_rennes)),
    x        = "Distance standardisée (z-score)",
    y        = "Lag spatial (moyenne des voisins)",
    caption  = "Matrice de pondération : Queen | Hypothèse de randomisation"
  ) +
  theme_moran()

ggsave("outputs/figures/07_diagramme_moran.png", p_moran_diag,
       width = 9, height = 7, dpi = 200)

cat("  Diagramme de Moran sauvegardé\n\n")

# =============================================================================
# 4. DIAGRAMME DE MORAN — Densité de population
# =============================================================================

y2_raw <- iris_rennes$DENSITE_POP
y2_raw[is.na(y2_raw)] <- mean(y2_raw, na.rm = TRUE)
y2_std <- as.numeric(scale(y2_raw))
y2_lag <- lag.listw(W_queen, y2_std, zero.policy = TRUE)

moran_dens <- moran.test(y2_raw, W_queen, randomisation = TRUE, zero.policy = TRUE)
I2 <- round(moran_dens$estimate["Moran I statistic"], 3)
p2 <- round(moran_dens$p.value, 4)

cat(sprintf("Densité population — I de Moran (Queen) : %.3f (p = %.4f)\n\n", I2, p2))

df_moran2 <- data.frame(y_std = y2_std, y_lag = y2_lag)

p_moran_dens <- ggplot(df_moran2, aes(x = y_std, y = y_lag)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_point(color = "#1d3557", size = 2.5, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "#e63946", linewidth = 1) +
  labs(
    title    = "Diagramme de Moran — Densité de population",
    subtitle = sprintf("I = %.3f | p = %.4f | Matrice Queen", I2, p2),
    x = "Densité standardisée", y = "Lag spatial"
  ) +
  theme_moran()

ggsave("outputs/figures/08_diagramme_moran_densite.png", p_moran_dens,
       width = 8, height = 6, dpi = 200)

# =============================================================================
# 5. COMPARAISON DES MATRICES W (tableau du cours)
# =============================================================================

cat("=== COMPARAISON DES MATRICES W (effet du voisinage) ===\n\n")

y_dist <- iris_rennes$dist_min_pharmacie_m
y_dist[is.na(y_dist)] <- mean(y_dist, na.rm = TRUE)

comp_W <- data.frame(
  Voisinage  = c("QUEEN", "ROOK", "kNN-3", "kNN-5"),
  I_W        = c(
    moran.test(y_dist, W_queen, zero.policy=TRUE)$estimate[1],
    moran.test(y_dist, W_rook,  zero.policy=TRUE)$estimate[1],
    moran.test(y_dist, W_knn3,  zero.policy=TRUE)$estimate[1],
    moran.test(y_dist, W_knn5,  zero.policy=TRUE)$estimate[1]
  ),
  p_value = c(
    moran.test(y_dist, W_queen, zero.policy=TRUE)$p.value,
    moran.test(y_dist, W_rook,  zero.policy=TRUE)$p.value,
    moran.test(y_dist, W_knn3,  zero.policy=TRUE)$p.value,
    moran.test(y_dist, W_knn5,  zero.policy=TRUE)$p.value
  )
) %>%
  mutate(
    I_W     = round(I_W, 4),
    p_value = round(p_value, 4),
    H0      = ifelse(p_value < 0.05, "Rejetée", "Non rejetée")
  )

cat("Tableau comparatif (comme dans le cours) :\n")
print(comp_W)
saveRDS(comp_W, "outputs/stats/comparaison_matrices_W.rds")

# =============================================================================
# 6. LISA — I DE MORAN LOCAL
# =============================================================================

cat("\n=== LISA — INDICES LOCAUX DE MORAN ===\n\n")

# Calcul des LISA pour la distance à la pharmacie
lisa <- localmoran(y_dist, W_queen, zero.policy = TRUE,
                   alternative = "two.sided")

iris_rennes$lisa_I      <- lisa[, "Ii"]
iris_rennes$lisa_E_I    <- lisa[, "E.Ii"]
iris_rennes$lisa_p      <- lisa[, "Pr(z != E(Ii))"]
iris_rennes$lisa_sig    <- iris_rennes$lisa_p < 0.05

# Classification LISA (quadrants significatifs)
y_std_lisa <- as.numeric(scale(y_dist))
y_lag_lisa <- lag.listw(W_queen, y_std_lisa, zero.policy = TRUE)

iris_rennes$lisa_cluster <- case_when(
  iris_rennes$lisa_sig & y_std_lisa >  0 & y_lag_lisa >  0 ~ "HH",
  iris_rennes$lisa_sig & y_std_lisa <  0 & y_lag_lisa <  0 ~ "LL",
  iris_rennes$lisa_sig & y_std_lisa >  0 & y_lag_lisa <  0 ~ "HL",
  iris_rennes$lisa_sig & y_std_lisa <  0 & y_lag_lisa >  0 ~ "LH",
  TRUE                                                       ~ "NS"
)

# Comptes par cluster
cat("Distribution des clusters LISA :\n")
print(table(iris_rennes$lisa_cluster))

# =============================================================================
# 7. CARTE LISA
# =============================================================================

cat("\n[Carte LISA] Génération...\n")

iris_rennes$lisa_cluster <- factor(iris_rennes$lisa_cluster,
  levels = c("HH", "LL", "HL", "LH", "NS"))

p_lisa <- ggplot(iris_rennes) +
  geom_sf(aes(fill = lisa_cluster), color = "white", linewidth = 0.3) +
  scale_fill_manual(
    name   = "Cluster LISA",
    values = c(
      "HH" = "#e63946",   # Zone loin + voisins loin → priorité haute
      "LL" = "#457b9d",   # Zone proche + voisins proches → bien couvert
      "HL" = "#f4a261",   # Zone loin + voisins proches → outlier
      "LH" = "#2a9d8f",   # Zone proche + voisins loin → outlier
      "NS" = "#f0f0f0"    # Non significatif
    ),
    labels = c(
      "HH" = "HH — Désert entouré de déserts",
      "LL" = "LL — Bien couvert, voisins bien couverts",
      "HL" = "HL — Désert isolé",
      "LH" = "LH — Bien couvert, voisins éloignés",
      "NS" = "NS — Non significatif"
    ),
    drop = FALSE
  ) +
  labs(
    title    = "Carte LISA — Clusters spatiaux de distance aux pharmacies",
    subtitle = sprintf("I de Moran global = %.3f (p < 0.05) | %d clusters significatifs",
                       I_val, sum(iris_rennes$lisa_sig)),
    caption  = "LISA : Local Indicator of Spatial Autocorrelation | Matrice Queen | α = 0.05"
  ) +
  theme_void() +
  theme(
    plot.title      = element_text(face="bold", size=13, color="#1a1a2e"),
    plot.subtitle   = element_text(size=9, color="#555555"),
    plot.caption    = element_text(size=7, color="#888888"),
    legend.position = "bottom",
    legend.direction = "vertical",
    plot.background = element_rect(fill="white", color=NA)
  )

ggsave("outputs/figures/09_lisa_clusters.png", p_lisa,
       width = 8, height = 8, dpi = 200)

# =============================================================================
# 8. CARTE I DE MORAN LOCAUX (valeurs brutes)
# =============================================================================

p_lisa_vals <- ggplot(iris_rennes) +
  geom_sf(aes(fill = lisa_I), color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low      = "#457b9d",
    mid      = "#f0f0f0",
    high     = "#e63946",
    midpoint = 0,
    name     = "I local",
    na.value = "grey90"
  ) +
  labs(
    title    = "I de Moran locaux — Distance à la pharmacie",
    subtitle = "Valeurs positives = regroupement | Valeurs négatives = dispersion",
    caption  = "Matrice de pondération Queen | INSEE RP 2022 | OSM 2024"
  ) +
  theme_void() +
  theme(
    plot.title      = element_text(face="bold", size=13, color="#1a1a2e"),
    plot.subtitle   = element_text(size=9, color="#555555"),
    plot.caption    = element_text(size=7, color="#888888"),
    legend.position = "bottom",
    plot.background = element_rect(fill="white", color=NA)
  )

ggsave("outputs/figures/10_lisa_valeurs.png", p_lisa_vals,
       width = 8, height = 7, dpi = 200)

# =============================================================================
# 9. LISA SUR LA DENSITÉ DE POPULATION
# =============================================================================

lisa_dens <- localmoran(y2_raw, W_queen, zero.policy = TRUE,
                        alternative = "two.sided")

iris_rennes$lisa_dens_I   <- lisa_dens[, "Ii"]
iris_rennes$lisa_dens_p   <- lisa_dens[, "Pr(z != E(Ii))"]
iris_rennes$lisa_dens_sig <- iris_rennes$lisa_dens_p < 0.05

y_dens_lag <- lag.listw(W_queen, y2_std, zero.policy = TRUE)

iris_rennes$lisa_dens_cluster <- case_when(
  iris_rennes$lisa_dens_sig & y2_std >  0 & y_dens_lag >  0 ~ "HH",
  iris_rennes$lisa_dens_sig & y2_std <  0 & y_dens_lag <  0 ~ "LL",
  iris_rennes$lisa_dens_sig & y2_std >  0 & y_dens_lag <  0 ~ "HL",
  iris_rennes$lisa_dens_sig & y2_std <  0 & y_dens_lag >  0 ~ "LH",
  TRUE                                                         ~ "NS"
)

cat("\nDistribution clusters LISA (densité population) :\n")
print(table(iris_rennes$lisa_dens_cluster))

# =============================================================================
# 10. SAUVEGARDE FINALE
# =============================================================================

saveRDS(iris_rennes, "data/iris_rennes.rds")
saveRDS(list(
  W_queen = W_queen,
  W_rook  = W_rook,
  W_knn3  = W_knn3,
  W_knn5  = W_knn5,
  nb_queen = nb_queen
), "data/matrices_W.rds")

cat("\n=== AUTOCORRÉLATION SPATIALE TERMINÉE ===\n")
cat("Fichiers sauvegardés :\n")
cat("  outputs/figures/07_diagramme_moran.png\n")
cat("  outputs/figures/08_diagramme_moran_densite.png\n")
cat("  outputs/figures/09_lisa_clusters.png\n")
cat("  outputs/figures/10_lisa_valeurs.png\n")
cat("  outputs/stats/moran_global.csv\n")
cat("  outputs/stats/comparaison_matrices_W.rds\n")
cat("  data/iris_rennes.rds (enrichi LISA)\n")
cat("  data/matrices_W.rds\n")