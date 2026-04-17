# The Future of Conifer Trees in British Columbia's Changing Climate

Species distribution modelling study predicting how 17 Pinales (conifer) species in British Columbia will respond to climate change under four IPCC emissions scenarios (SSP1-2.6 through SSP5-8.5) for the period 2081–2100.

The full report is available in [`report_BC_pinales.pdf`](report_BC_pinales.pdf).

**Key findings:** The threatened whitebark pine (*Pinus albicaulis*) faces near-extinction under high-emissions scenarios. Dominant species like Douglas fir and western redcedar are projected to expand northward, while overall species richness shifts toward higher latitudes.

## Background

British Columbia's forests cover 57% of the province and span a wide range of bioclimatic zones — from coastal rainforest to boreal and subalpine. Increasing temperatures and more extreme weather are expected to substantially reshape tree species distributions. This project models those changes for all extant Pinales species with sufficient observational data in BC.

## Methods

Ensemble species distribution models were built with the **biomod2** package using three algorithms per species:
- Generalized Linear Model (GLM)
- Random Forest (RF)
- Maximum Entropy (MAXNET)

Each algorithm was trained twice with 80/20 random cross-validation. Individual models were combined into two ensemble models per species — committee averaging (EMca) and weighted mean (EMwmean) — retaining only models with TSS ≥ 0.7. Variable selection used the `select07` function from **mecofun** to remove highly correlated bioclimatic predictors.

**Occurrence data:** iNaturalist research-grade observations from GBIF (≥2000, ≥100 records per species) — ~57,948 observations across 17 species.

**Climate data:** WorldClim v2.1 at 2.5 arcmin resolution (current) and ACCESS-CM2 CMIP6 projections (SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5) for 2081–2100.

## Species Modelled

| Scientific Name | Common Name |
|---|---|
| *Pseudotsuga menziesii* | Douglas fir |
| *Thuja plicata* | Western redcedar |
| *Tsuga heterophylla* | Western hemlock |
| *Picea sitchensis* | Sitka spruce |
| *Abies grandis* | Grand fir |
| *Abies amabilis* | Pacific silver fir |
| *Picea engelmannii* | Engelmann spruce |
| *Pinus albicaulis* | Whitebark pine |
| *Pinus contorta* | Lodgepole pine |
| *Pinus monticola* | Western white pine |
| *Pinus ponderosa* | Ponderosa pine |
| *Larix occidentalis* | Western larch |
| *Larix lyallii* | Subalpine larch |
| *Tsuga mertensiana* | Mountain hemlock |
| *Taxus brevifolia* | Pacific yew |
| *Juniperus scopulorum* | Rocky Mountain juniper |
| *Xanthocyparis nootkatensis* | Alaska cedar |

## Repository Structure

```
├── pinales_BC.R          # Data preparation: climate grids, occurrence filtering, single-species prototype
├── model_loop.R          # Batch modelling loop: runs biomod2 pipeline for all 17 species
├── analysis.R            # Post-processing: projections, species richness change, climate variable plots
│
├── climate/              # WorldClim raster files (not tracked — download separately)
├── data/                 # Raw GBIF occurrence data (not tracked)
├── <Genus.species>/      # Per-species biomod2 model output directories
│
├── climate_grid_current.dat      # Preprocessed current bioclimate grid
├── climate_grid_future*.dat      # Preprocessed future bioclimate grids (one per SSP)
├── occurrence_iNat_100obs_y2000.Rds  # Filtered occurrence records
├── projections*.rds      # Saved species richness projections per SSP
│
└── results/              # Output figures (PDFs)
```

## Workflow

1. **`pinales_BC.R`** — Downloads and preprocesses WorldClim climate rasters, filters and formats GBIF occurrence data, and prototypes the full biomod2 pipeline for a single species.
2. **`model_loop.R`** — Iterates over all species: formats biomod2 data objects, selects uncorrelated bioclimatic variables, trains GLM/RF/MAXNET models, and builds ensemble models. Skips species whose ensemble models already exist on disk.
3. **`analysis.R`** — Loads saved models, runs ensemble forecasts for current and future climate scenarios, aggregates species richness across all species, and produces all result figures.

## Dependencies (R packages)

`terra`, `biomod2`, `mecofun`, `ggplot2`, `tidyverse`, `sf`, `bcmaps`, `rnaturalearth`, `geodata`, `gridExtra`

Install mecofun from the Macroecology Group, University of Potsdam; all others are available on CRAN.

## Data Sources

- Occurrence: [GBIF](https://www.gbif.org) — iNaturalist research-grade observations
- Current climate: [WorldClim v2.1](https://www.worldclim.org) (Fick & Hijmans 2017)
- Future climate: CMIP6 ACCESS-CM2 via WorldClim (O'Neill et al. 2017)
- Biogeoclimatic zones: [B.C. Ministry of Forests](https://www.for.gov.bc.ca/hre/becweb/)
