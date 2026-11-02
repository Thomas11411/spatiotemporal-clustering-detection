# Spatial-Temporal Clustering Detection (`STclust`)

An R package for spatiotemporal cluster detection and epidemiological hotspot analysis using **Mixture Estimation Equations (MEE)**, specifically engineered to handle complex surveillance data with time-varying covariates, spatial heterogeneity, and excessive zeros (zero inflation).

## 📌 Background & Motivation

In the study of infectious diseases, cluster analysis is a useful statistical method to explore the relationship between disease hot spots and propagation patterns associated with risk factors. However, when covariates and excessive zeros also produce clustering effects for responses, traditional identification methods for geographic clusters may not address this confounding issue well.

Moreover, for complex surveillance systems, scrutinizing all risk factors that could simultaneously affect health outcomes and geographic clusters would be difficult, which could cause large spatial variations and thus create a convergency issue for the Bayesian algorithm.

To address these issues, a **mixture scan statistic** that is associated with two-part and cluster models is developed to identify geographic clusters for spatial data with covariates and excessive zeros. A model selection procedure is proposed to evaluate the significance of geographic clusters for the mixture scan statistic. Simulation studies and data analysis for associations between dengue infection and environmental risk factors are used to illustrate the proposed method.

## ✨ Key Features & Workflow

* **Spatiotemporal Window Scanning**: Automatically scans candidate spatiotemporal clusters across spatial polygons (`sf`) and extracts time-series trajectory steps directly from attributes.
* **Zero-Inflation & Covariate Adjustment**: Evaluates cluster candidates while controlling for zero inflation and time-varying risk factors using Mixture Estimation Equations (`MEE`).
* **Parallel Processing**: Speeds up candidate evaluations and spatial matrix calculations using multi-core processing (`parallel`).
* **Automated Hotspot Identification**: Evaluates candidate configurations using model-selection criteria to isolate primary transmission hot spots without over-clustering.
* **Spatial Visualization**: Directly accepts spatial objects and integrates with visualization utilities to map identified hotspots over target study regions.

## 📁 Package Structure

```text
spatial-temporal-clustering-detection/
├── R/             # Core functions (spatiotemporal_cluster_analysis, evaluate_candidate_cluster, MEE, etc.)
├── inst/          # Internal external data (inst/extdata/st_integrated.shp)
├── man/           # Documentation files (.Rd)
├── DESCRIPTION    # Package metadata and dependencies
└── NAMESPACE      # Exported functions and imported packages