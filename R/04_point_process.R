# =============================================================================
# PharmOptim-Rennes | Script 04 — Processus Ponctuels
# Auteur : Hafiz | Statistique Spatiale — Magistère 2
# =============================================================================
# Méthodes (Chapitre 3 du cours) :
#   - Estimation de l'intensité λ(x) — KDE (Kernel Density Estimation)
#   - Test CSR (Complete Spatial Randomness) — processus de Poisson homogène
#   - Fonction K de Ripley + L de Besag
#   - Fonction K inhomogène de Baddeley (Kinhom)
#   - Fonction D de Diggle (pharmacies vs médecins)
#   - Enveloppes Monte Carlo (999 simulations)
# =============================================================================

library(sf)
library(tidyverse)
library(spatstat)
library(spatstat.geom)
library(spatstat.explore)
library(ggplot2)
library(viridis)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/stats",   recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 0. CHARGEMENT & PRÉPARATION
# =============================================================================

cat("Chargement des données...\n")

iris_rennes <- readRDS("data/iris_rennes.rds")
pharmacies  <- readRDS("data/pharmacies.rds")
medecins    <- readRDS("data/medecins.rds")
hopitaux    <- readRDS("data/hopitaux.rds")
rennes_poly <- readRDS("data/rennes_contour.rds")

cat(sprintf("  %d pharmacies | %d médecins | %d hôpitaux\n\n",
            nrow(pharmacies), nrow(medecins), nrow(hopitaux)))

# -----------------------------------------------------------------------------
# Conversion en objets spatstat (ppp = planar point pattern)
# Fenêtre d'observation = contour de la commune de Rennes
# -----------------------------------------------------------------------------

# Contour Rennes → objet owin (spatstat)
rennes_owin <- as.owin(st_geometry(rennes_poly))

# Fonction de conversion sf → ppp
sf_to_ppp <- function(sf_points, win) {
  coords <- st_coordinates(sf_points)
  ppp(x = coords[, 1], y = coords[, 2], window = win,
      check = FALSE)
}

# Création des ppp
ppp_pharma   <- sf_to_ppp(pharmacies, rennes_owin)
ppp_medecins <- sf_to_ppp(medecins,   rennes_owin)

# Suppression des points dupliqués (spatstat l'exige)
ppp_pharma   <- unique(ppp_pharma)
ppp_medecins <- unique(ppp_medecins)

cat(sprintf("Processus ponctuels créés :\n"))
cat(sprintf("  Pharmacies : %d points dans la fenêtre\n", ppp_pharma$n))
cat(sprintf("  Médecins   : %d points dans la fenêtre\n\n", ppp_medecins$n))

# =============================================================================
# 1. INTENSITÉ — KDE (Kernel Density Estimation)
# =============================================================================

cat("=== ESTIMATION DE L'INTENSITÉ λ(x) — KDE ===\n\n")

# Estimation de l'intensité homogène λ̂
lambda_pharma   <- intensity(ppp_pharma)
lambda_medecins <- intensity(ppp_medecins)

cat(sprintf("Intensité homogène λ̂ :\n"))
cat(sprintf("  Pharmacies : %.4f points/m²  (= %.2f pharmacies/km²)\n",
            lambda_pharma, lambda_pharma * 1e6))
cat(sprintf("  Médecins   : %.4f points/m²  (= %.2f médecins/km²)\n\n",
            lambda_medecins, lambda_medecins * 1e6))

# KDE avec différentes bandes passantes (sigma)
# Sigma optimal par règle de Scott
sigma_scott <- bw.scott(ppp_pharma)
cat(sprintf("Bande passante optimale (Scott) : σ = %.0f m\n", sigma_scott[1]))

# KDE pharmacies
kde_pharma_300 <- density(ppp_pharma, sigma = 300, edge = TRUE)
kde_pharma_500 <- density(ppp_pharma, sigma = 500, edge = TRUE)
kde_pharma_opt <- density(ppp_pharma, sigma = sigma_scott[1], edge = TRUE)

# KDE médecins
kde_medecins   <- density(ppp_medecins, sigma = 500, edge = TRUE)

cat("KDE calculés (σ = 300m, 500m, optimal)\n\n")

# Conversion KDE → data.frame pour ggplot
kde_to_df <- function(kde_obj) {
  expand.grid(
    x = kde_obj$xcol,
    y = kde_obj$yrow
  ) %>%
  mutate(z = as.vector(t(kde_obj$v))) %>%
  filter(!is.na(z))
}

df_kde_pharma <- kde_to_df(kde_pharma_500)
df_kde_medecins <- kde_to_df(kde_medecins)

# Plot KDE pharmacies
p_kde_pharma <- ggplot() +
  geom_raster(data = df_kde_pharma,
              aes(x = x, y = y, fill = z), interpolate = TRUE) +
  geom_sf(data = rennes_poly, fill = NA, color = "white",
          linewidth = 0.8, inherit.aes = FALSE) +
  geom_sf(data = pharmacies, color = "white", size = 1.5,
          shape = 3, inherit.aes = FALSE) +
  scale_fill_viridis_c(
    option = "inferno", name = "Intensité λ(x)",
    na.value = "transparent"
  ) +
  coord_sf() +
  labs(
    title    = "Estimation de l'intensité — Pharmacies à Rennes",
    subtitle = sprintf("KDE (σ = 500m) | λ̂ homogène = %.2f pharmacies/km²",
                       lambda_pharma * 1e6),
    caption  = "Croix blanches = pharmacies | Source : OSM 2024"
  ) +
  theme_void() +
  theme(
    plot.title      = element_text(face="bold", size=13, color="#1a1a2e"),
    plot.subtitle   = element_text(size=9, color="#555555"),
    plot.caption    = element_text(size=7, color="#888888"),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.background = element_rect(fill="white", color=NA)
  )

ggsave("outputs/figures/11_kde_pharmacies.png", p_kde_pharma,
       width = 8, height = 7, dpi = 200)
cat("  KDE pharmacies sauvegardé\n")

# Plot KDE médecins
p_kde_medecins <- ggplot() +
  geom_raster(data = df_kde_medecins,
              aes(x = x, y = y, fill = z), interpolate = TRUE) +
  geom_sf(data = rennes_poly, fill = NA, color = "white",
          linewidth = 0.8, inherit.aes = FALSE) +
  geom_sf(data = medecins, color = "white", size = 1.5,
          shape = 3, inherit.aes = FALSE) +
  scale_fill_viridis_c(
    option = "mako", name = "Intensité λ(x)",
    na.value = "transparent"
  ) +
  coord_sf() +
  labs(
    title    = "Estimation de l'intensité — Médecins à Rennes",
    subtitle = sprintf("KDE (σ = 500m) | λ̂ homogène = %.2f médecins/km²",
                       lambda_medecins * 1e6),
    caption  = "Croix blanches = médecins | Source : OSM 2024"
  ) +
  theme_void() +
  theme(
    plot.title      = element_text(face="bold", size=13, color="#1a1a2e"),
    plot.subtitle   = element_text(size=9, color="#555555"),
    plot.caption    = element_text(size=7, color="#888888"),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.background = element_rect(fill="white", color=NA)
  )

ggsave("outputs/figures/12_kde_medecins.png", p_kde_medecins,
       width = 8, height = 7, dpi = 200)
cat("  KDE médecins sauvegardé\n\n")

# =============================================================================
# 2. TEST CSR — PROCESSUS DE POISSON HOMOGÈNE
# =============================================================================

cat("=== TEST CSR (Complete Spatial Randomness) ===\n\n")

# Quadrat count test (test du chi2 sur quadrillage)
qct <- quadrat.test(ppp_pharma, nx = 5, ny = 5)
cat("Quadrat test (pharmacies) :\n")
cat(sprintf("  Chi² = %.3f | df = %d | p = %.4f\n\n",
            qct$statistic, qct$parameter, qct$p.value))

if (qct$p.value < 0.05) {
  cat("  => H0 (CSR) REJETÉE — distribution non aléatoire\n\n")
} else {
  cat("  => H0 (CSR) non rejetée\n\n")
}

# =============================================================================
# 3. FONCTION K DE RIPLEY
# =============================================================================

cat("=== FONCTION K DE RIPLEY ===\n\n")

# Distances à tester (en mètres)
r_vals <- seq(0, 2000, by = 50)

# K de Ripley estimée avec correction des effets de bord (Ripley)
K_pharma <- Kest(ppp_pharma, r = r_vals,
                 correction = "Ripley")

cat("Fonction K de Ripley calculée\n")
cat(sprintf("  Distances : 0 à %.0f m (pas = 50m)\n", max(r_vals)))

# L de Besag (transformation de K)
L_pharma <- Lest(ppp_pharma, r = r_vals,
                 correction = "Ripley")

cat("Fonction L de Besag calculée\n\n")

# Enveloppes Monte Carlo (999 simulations CSR)
cat("Calcul des enveloppes Monte Carlo (999 simulations)...\n")
cat("  (peut prendre 1-2 minutes)\n")

set.seed(2024)
K_env <- envelope(ppp_pharma, Kest, r = r_vals,
                  nsim = 199,           # 199 sims pour niveau 5% bilatéral
                  correction = "Ripley",
                  verbose = FALSE)

L_env <- envelope(ppp_pharma, Lest, r = r_vals,
                  nsim = 199,
                  correction = "Ripley",
                  verbose = FALSE)

cat("  Enveloppes calculées\n\n")

# Conversion en data.frame pour ggplot
df_K <- data.frame(
  r     = K_env$r,
  obs   = K_env$obs,
  theo  = K_env$theo,
  lo    = K_env$lo,
  hi    = K_env$hi
)

df_L <- data.frame(
  r    = L_env$r,
  obs  = L_env$obs  - L_env$r,   # L(r) - r pour tester si = 0
  theo = L_env$theo - L_env$r,
  lo   = L_env$lo   - L_env$r,
  hi   = L_env$hi   - L_env$r
)

# Plot K de Ripley
p_K <- ggplot(df_K) +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi),
              fill = "#adb5bd", alpha = 0.4) +
  geom_line(aes(x = r, y = theo, linetype = "Théorique (CSR)"),
            color = "#e63946", linewidth = 0.8) +
  geom_line(aes(x = r, y = obs, linetype = "Observée"),
            color = "#1d3557", linewidth = 1) +
  scale_linetype_manual(
    name   = "",
    values = c("Théorique (CSR)" = "dashed", "Observée" = "solid")
  ) +
  scale_x_continuous(labels = scales::comma_format(suffix = "m")) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    title    = "Fonction K de Ripley — Pharmacies à Rennes",
    subtitle = "Zone grise = enveloppe Monte Carlo (199 sim.) | α = 5%",
    x        = "Distance r (m)",
    y        = "K̂(r)",
    caption  = "K̂(r) > K_Poisson(r) → agrégation | K̂(r) < K_Poisson(r) → répulsion"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face="bold", size=13, color="#1a1a2e"),
    plot.subtitle   = element_text(size=9, color="#555555"),
    legend.position = "top",
    plot.background = element_rect(fill="white", color=NA)
  )

ggsave("outputs/figures/13_K_ripley.png", p_K,
       width = 9, height = 6, dpi = 200)
cat("  Fonction K sauvegardée\n")

# Plot L de Besag (L(r) - r, attendu = 0 sous CSR)
p_L <- ggplot(df_L) +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi),
              fill = "#adb5bd", alpha = 0.4) +
  geom_hline(yintercept = 0, color = "#e63946",
             linewidth = 0.8, linetype = "dashed") +
  geom_line(aes(x = r, y = obs),
            color = "#1d3557", linewidth = 1) +
  scale_x_continuous(labels = scales::comma_format(suffix = "m")) +
  labs(
    title    = "Fonction L de Besag — L(r) - r",
    subtitle = "Attendu = 0 sous H0 (CSR) | Zone grise = enveloppe Monte Carlo",
    x        = "Distance r (m)",
    y        = "L̂(r) - r",
    caption  = "Au-dessus de 0 → agrégation | En-dessous de 0 → répulsion"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face="bold", size=13, color="#1a1a2e"),
    plot.subtitle   = element_text(size=9, color="#555555"),
    plot.background = element_rect(fill="white", color=NA)
  )

ggsave("outputs/figures/14_L_besag.png", p_L,
       width = 9, height = 6, dpi = 200)
cat("  Fonction L sauvegardée\n\n")

# =============================================================================
# 4. FONCTION K INHOMOGÈNE DE BADDELEY (Kinhom)
# =============================================================================

cat("=== FONCTION K INHOMOGÈNE DE BADDELEY ===\n\n")

# Estimation de l'intensité non-homogène par KDE
lambda_hat <- density(ppp_pharma, sigma = sigma_scott[1])

# Kinhom
K_inhom <- Kinhom(ppp_pharma, lambda = lambda_hat,
                  r = r_vals, correction = "Ripley")

# Enveloppe Kinhom
set.seed(2024)
K_inhom_env <- envelope(ppp_pharma, Kinhom,
                         lambda = lambda_hat,
                         r = r_vals, nsim = 199,
                         correction = "Ripley",
                         verbose = FALSE)

df_Kinhom <- data.frame(
  r    = K_inhom_env$r,
  obs  = K_inhom_env$obs,
  theo = K_inhom_env$theo,
  lo   = K_inhom_env$lo,
  hi   = K_inhom_env$hi
)

p_Kinhom <- ggplot(df_Kinhom) +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi),
              fill = "#adb5bd", alpha = 0.4) +
  geom_line(aes(x = r, y = theo, linetype = "Théorique"),
            color = "#e63946", linewidth = 0.8) +
  geom_line(aes(x = r, y = obs,  linetype = "Observée"),
            color = "#1d3557", linewidth = 1) +
  scale_linetype_manual(
    name   = "",
    values = c("Théorique" = "dashed", "Observée" = "solid")
  ) +
  scale_x_continuous(labels = scales::comma_format(suffix = "m")) +
  labs(
    title    = "Fonction K inhomogène de Baddeley — Pharmacies",
    subtitle = "Corrige l'hétérogénéité spatiale de l'intensité λ(x)",
    x        = "Distance r (m)",
    y        = "K̂_inhom(r)",
    caption  = "Baddeley et al. (2000) | Enveloppe Monte Carlo 199 sim."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face="bold", size=13, color="#1a1a2e"),
    plot.subtitle   = element_text(size=9, color="#555555"),
    legend.position = "top",
    plot.background = element_rect(fill="white", color=NA)
  )

ggsave("outputs/figures/15_K_inhom.png", p_Kinhom,
       width = 9, height = 6, dpi = 200)
cat("  Kinhom sauvegardé\n\n")

# =============================================================================
# 5. FONCTION D DE DIGGLE — Pharmacies vs Médecins
# =============================================================================

cat("=== FONCTION D DE DIGGLE — Pharmacies vs Médecins ===\n\n")
cat("D(r) = K_pharmacies(r) - K_medecins(r)\n")
cat("Si D(r) > 0 → pharmacies plus agrégées que médecins\n\n")

# K pour chaque type de point
K_pharma_d   <- Kest(ppp_pharma,   r = r_vals, correction = "Ripley")
K_medecins_d <- Kest(ppp_medecins, r = r_vals, correction = "Ripley")

# Fonction D de Diggle
df_D <- data.frame(
  r = K_pharma_d$r,
  D = K_pharma_d$iso - K_medecins_d$iso
)

# Enveloppe par permutation (mélange des labels)
ppp_combined <- superimpose(
  pharma   = ppp_pharma,
  medecins = ppp_medecins
)

set.seed(2024)
n_pharma <- ppp_pharma$n
n_total  <- ppp_combined$n

D_sims <- replicate(199, {
  idx   <- sample(n_total, n_pharma)
  ppp_s <- ppp_combined[idx]
  ppp_c <- ppp_combined[-idx]
  K_s   <- Kest(ppp_s, r = r_vals, correction = "Ripley")$iso
  K_c   <- Kest(ppp_c, r = r_vals, correction = "Ripley")$iso
  K_s - K_c
})

df_D$lo <- apply(D_sims, 1, min)
df_D$hi <- apply(D_sims, 1, max)

p_D <- ggplot(df_D) +
  geom_ribbon(aes(x = r, ymin = lo, ymax = hi),
              fill = "#adb5bd", alpha = 0.4) +
  geom_hline(yintercept = 0, color = "#e63946",
             linewidth = 0.8, linetype = "dashed") +
  geom_line(aes(x = r, y = D),
            color = "#1d3557", linewidth = 1) +
  scale_x_continuous(labels = scales::comma_format(suffix = "m")) +
  labs(
    title    = "Fonction D de Diggle — Pharmacies vs Médecins",
    subtitle = "D(r) = K_pharmacies(r) - K_médecins(r) | Zone grise = enveloppe permutations",
    x        = "Distance r (m)",
    y        = "D(r)",
    caption  = "D(r) > 0 → pharmacies plus agrégées | D(r) < 0 → médecins plus agrégés"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face="bold", size=13, color="#1a1a2e"),
    plot.subtitle   = element_text(size=9, color="#555555"),
    plot.background = element_rect(fill="white", color=NA)
  )

ggsave("outputs/figures/16_D_diggle.png", p_D,
       width = 9, height = 6, dpi = 200)
cat("  Fonction D de Diggle sauvegardée\n\n")

# =============================================================================
# 6. SAUVEGARDE DES RÉSULTATS
# =============================================================================

resultats_ripley <- list(
  lambda_pharma   = lambda_pharma,
  lambda_medecins = lambda_medecins,
  sigma_scott     = sigma_scott[1],
  quadrat_test    = list(
    chi2    = qct$statistic,
    df      = qct$parameter,
    p_value = qct$p.value
  ),
  K_pharma        = df_K,
  L_pharma        = df_L,
  K_inhom         = df_Kinhom,
  D_diggle        = df_D
)

saveRDS(resultats_ripley, "outputs/stats/resultats_ripley.rds")
saveRDS(kde_pharma_500,   "outputs/stats/kde_pharma_500.rds")
saveRDS(kde_medecins,     "outputs/stats/kde_medecins.rds")

cat("=== PROCESSUS PONCTUELS TERMINÉS ===\n")
cat("Fichiers sauvegardés :\n")
cat("  outputs/figures/11_kde_pharmacies.png\n")
cat("  outputs/figures/12_kde_medecins.png\n")
cat("  outputs/figures/13_K_ripley.png\n")
cat("  outputs/figures/14_L_besag.png\n")
cat("  outputs/figures/15_K_inhom.png\n")
cat("  outputs/figures/16_D_diggle.png\n")
cat("  outputs/stats/resultats_ripley.rds\n")
cat("  outputs/stats/kde_pharma_500.rds\n")
