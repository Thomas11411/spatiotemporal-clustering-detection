# Spatial-Temporal Clustering Detection (`STclust`)

An R package for spatiotemporal cluster detection and epidemiological hotspot analysis using **Mixture Estimation Equations (MEE)**, specifically engineered to handle complex surveillance data with time-varying covariates, spatial heterogeneity, and excessive zeros (zero inflation).

## 📌 Background & Motivation

In the study of infectious diseases, cluster analysis is a useful statistical method to explore the relationship between disease hot spots and propagation patterns associated with risk factors. However, when covariates and excessive zeros also produce clustering effects for responses, traditional identification methods for geographic clusters may not address this confounding issue well.

Moreover, for complex surveillance systems, scrutinizing all risk factors that could simultaneously affect health outcomes and geographic clusters would be difficult, which could cause large spatial variations and thus create a convergency issue for the Bayesian algorithm.

To address these issues, a **mixture scan statistic** that is associated with two-part and cluster models is developed to identify geographic clusters for spatial data with covariates and excessive zeros. A model selection procedure is proposed to evaluate the significance of geographic clusters for the mixture scan statistic. Simulation studies and data analysis for associations between dengue infection and environmental risk factors are used to illustrate the proposed method.

## ✨ Key Features & Workflow

* **Spatiotemporal Window Scanning**: Automatically scans candidate spatiotemporal clusters across spatial polygons (`sf`) over multiple time periods (`T.cnt`).
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
```

## 📦 Installation

### Install from GitHub

You can install the development version of `STclust` directly from GitHub using `remotes` or `devtools`:

```R
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("Thomas1141/spatial-temporal-clustering-detection")
```

## 🚀 Quick Start / Usage Example

Below is the standard execution workflow using the `STclust` (`spatial-temporal-clustering-detection`) package to identify spatiotemporal clusters:

```R
library(spatial.temporal.clustering.detection)
library(sf)

# 1. Access built-in spatial shapefile (or load your own sf object)
shp_path <- system.file("extdata", "st_integrated.shp", package = "spatial-temporal-clustering-detection")
my_map <- sf::st_read(shp_path)

# 2. Run Spatiotemporal Cluster Identification with Execution Timing
t0 <- Sys.time()

res <- spatiotemporal_cluster_analysis(
  shp.file = my_map,
  max.dist = 0.5,
  cutoff   = 1.5,
  T.cnt    = 3,
  nc       = 11
)

t1 <- Sys.time()

# Check computation duration
print(t1 - t0)

# 3. View Detection Results & Summary
print(res)
```

## ⚙️ Key Parameters

* `shp.file`: Spatial `sf` object containing geographic polygon boundaries and temporal attribute variables.
* `max.dist`: Maximum spatial distance threshold for candidate scanning windows.
* `cutoff`: Cutoff distance for spatial variogram/autocorrelation structure modeling.
* `T.cnt`: Number of temporal windows/periods for tracking dynamic trajectory variations.
* `nc`: Number of CPU cores allocated for parallel computing.

## 📄 Citation & License

If you use this package or methodology in your research, please cite the corresponding paper on mixture scan statistics for spatiotemporal cluster analysis.

**Author & Maintainer:** Li-Wei, Chen (`pig21411@gmail.com`)  
Institute of Population Health Sciences, National Health Research Institutes (NHRI), Taiwan.