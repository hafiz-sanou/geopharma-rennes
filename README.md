# GeoPharma Analytics 🗺️💊

> **Optimisation spatiale de l'implantation des pharmacies à Rennes**  
> Analyse géospatiale multicritère pour identifier les zones sous-desservies en pharmacies d'officine

---

## 📋 Présentation

Ce projet analyse la distribution spatiale des pharmacies sur le territoire de Rennes Métropole afin d'identifier les zones géographiques sous-desservies et de proposer des localisations optimales pour de nouvelles implantations. Il combine des méthodes d'analyse spatiale avancées (autocorrélation, statistiques de processus ponctuels) avec un scoring multicritère intégrant des données démographiques, socio-économiques et d'accessibilité.

L'analyse a été réalisée dans le cadre du cours *PharmOptim* du Master Mathématiques Appliquées & Statistiques — Magistère 2 (Université de Rennes, 2025).

---

## 🗂️ Structure du projet

```
geopharma-rennes/
├── data/                   # Données sources (non versionnées si sensibles)
│   ├── pharmacies_rennes.geojson
│   ├── iris_rennes.geojson
│   └── ...
├── R/                      # Scripts R d'analyse
├── outputs/                # Cartes et graphiques générés
├── rapport.Rmd             # Dashboard flexdashboard (source principale)
└── rapport.html            # Dashboard compilé
```

---

## 🛠️ Outils & Technologies

| Catégorie | Outils |
|-----------|--------|
| Langage | R |
| Dashboard | flexdashboard (storyboard) |
| Données spatiales | sf, OpenStreetMap, INSEE |
| Analyse spatiale | spdep, spatstat |
| Visualisation | tmap, leaflet, ggplot2, plotly |
| Statistiques | Moran's I, LISA, Ripley K, Besag L, KDE, Diggle D |

---

## 📊 Méthodologie

### 1. Collecte & préparation des données
- Extraction des pharmacies via OpenStreetMap
- Données socio-démographiques par IRIS (INSEE)
- Découpage spatial par quartiers de Rennes Métropole

### 2. Analyse d'autocorrélation spatiale
- **Indice de Moran global** : I = 0.278, p < 0.0001 → clustering spatial significatif
- **Statistiques LISA** : identification de 4 clusters HH (haute densité entourée de haute densité), révélant une concentration dans l'hypercentre

### 3. Statistiques de processus ponctuels
- **Fonction K de Ripley** & **L de Besag** : clustering net jusqu'à ~800m
- **KDE** (Kernel Density Estimation) : cartographie continue de la densité de pharmacies
- **Test de Diggle** : confirmation que la distribution n'est pas aléatoire (CSR rejetée, p < 0.001)

### 4. Scoring multicritère
Chaque IRIS reçoit un score composite (0–100) combinant :

| Critère | Poids |
|---------|-------|
| Densité de population | 25% |
| Part des 65 ans et plus | 20% |
| Taux de pauvreté | 20% |
| Distance à la pharmacie la plus proche | 35% |

Score final = `(densité × 0.25) + (seniors × 0.20) + (pauvreté × 0.20) + (distance × 0.35)`

---

## 📈 Résultats principaux

- **Autocorrélation spatiale significative** (Moran I = 0.278, p < 0.0001) → les pharmacies se regroupent spatialement et ne sont pas distribuées aléatoirement
- **4 clusters LISA HH** identifiés en hypercentre (centre-ville historique densément équipé)
- **Zone prioritaire** : IRIS Saint-Laurent avec le score composite le plus élevé (**41.4/100**) — population âgée, niveau de revenus faible, faible accessibilité aux pharmacies existantes
- Le clustering est net jusqu'à 800m de rayon, indiquant que les zones périphériques sont structurellement sous-équipées

---

## 🗺️ Dashboard interactif

Le fichier `rapport.html` contient un dashboard flexdashboard (format storyboard) avec :
- Carte Leaflet interactive des pharmacies et des IRIS
- Carte des scores multicritères par IRIS
- Visualisations des fonctions K/L de Ripley
- Correlograms de Moran

Pour le visualiser, ouvrir `rapport.html` directement dans un navigateur.

---

## ▶️ Reproduire l'analyse

### Prérequis

```r
install.packages(c(
  "flexdashboard", "sf", "spdep", "spatstat",
  "tmap", "leaflet", "plotly", "ggplot2",
  "osmdata", "tidyverse"
))
```

### Compilation du rapport

```r
rmarkdown::render("rapport.Rmd")
```

---

## 👤 Auteur

**Hafiz Sanou**  
Master Mathématiques Appliquées & Statistiques — Magistère 2  
Université de Rennes · 2025  

📧 [hafizsanou2004@gmail.com](mailto:hafizsanou2004@gmail.com)  
🔗 [Portfolio](https://portfolio-zeta-lemon-29.vercel.app/) · [GitHub](https://github.com/hafiz-sanou)

---

## 📄 Licence

Ce projet est à des fins académiques. Les données OpenStreetMap sont sous licence [ODbL](https://www.openstreetmap.org/copyright).
