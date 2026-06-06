# Harbin manuscript reproducibility repository

This repository contains the analysis code used for the manuscript on forest functioning, hydro-climatic trends, hot-dry event resilience, and vulnerability across Northeast Asia.

The repository is intended for journal data availability and reproducibility. Large downloaded rasters, processed GeoTIFF outputs, shapefiles, QGIS projects, maps, and manuscript drafts are intentionally excluded from GitHub. The scripts document how the data were downloaded/exported, preprocessed, analyzed, and summarized.

## Repository contents

- `Scripts/rgee-01.R` - Google Earth Engine/rgee setup and example export workflows, including TerraClimate access and export to Google Drive.
- `Scripts/DataPrep.R` - preprocessing and alignment of NDVI, kNDVI, TerraClimate, and SPEI raster stacks.
- `Scripts/2-Trend-ndvikndvi.R` - pixel-wise Sen slope and Mann-Kendall trends for NDVI/kNDVI, forest masking, and tabulation.
- `Scripts/3- climate trends.R` - hydro-climatic trend analysis for precipitation, PET, water balance, temperature, and PDSI.
- `Scripts/3.1-Stage Hotdays.R` - hot-dry event, kNDVI functioning, resistance, and recovery calculations.
- `Scripts/4-correlation-Coupling.R` - coupling/correlation analysis between kNDVI and hydro-climatic variables.
- `Scripts/10_vulnerability_index.R` - composite forest functional vulnerability index.
- `Scripts/11_RF_vulnerability_drivers.R` - random forest analysis of vulnerability drivers.
- `Scripts/RevisionProfessor.R` - additional revision-stage analyses.

## Data sources

The workflow uses public geospatial datasets accessed through Google Earth Engine and local preprocessing:

- MODIS vegetation products used to derive annual NDVI/kNDVI.
- TerraClimate (`IDAHO_EPSCOR/TERRACLIMATE`) for precipitation, PET, PDSI, and temperature variables.
- SPEI monthly drought index data.
- MODIS land-cover forest-type masks.
- Public administrative/country boundary layers used only for spatial masking and country-level summaries.

## Reproducibility notes

1. Install R and the packages used in the scripts, including `terra`, `sf`, `dplyr`, `readr`, `tibble`, `Kendall`, `ranger`, `reticulate`, and `rgee`.
2. Authenticate Google Earth Engine before running the rgee workflow.
3. Export/download the raw raster stacks into the paths expected by `Scripts/DataPrep.R`.
4. Run the scripts in the numbered order shown above.
5. Generated rasters and tables will be written under `data/` and `outputs/`, which are ignored by Git because of file size and journal repository hygiene.

## Suggested journal data availability statement

The R scripts used for Google Earth Engine data access, preprocessing, trend analysis, coupling analysis, vulnerability-index construction, and random forest driver analysis are available in this GitHub repository. The underlying geospatial datasets are publicly available from their original providers, including Google Earth Engine-hosted MODIS and TerraClimate products, SPEI data, MODIS land-cover products, and public administrative boundary datasets. Large intermediate and derived raster outputs are not stored in the repository because of file-size constraints, but they can be reproduced by running the scripts in the documented order after downloading/exporting the source datasets.
