# # 1) Choose the env Python you want: candidate [2]
# py_path <- "C:/Users/Administrator/AppData/Local/r-miniconda/envs/r-reticulate/python.exe"
# file.exists(py_path)   # should be TRUE
# 
# # 2) Force reticulate to use this Python for this session
# Sys.setenv(RETICULATE_PYTHON = py_path)
# 
# library(reticulate)
# py_config()
# 
# 
# 
# 
# 
# 
# # Install needed Python packages into the 'r-reticulate' env
# py_install(
#   packages = c("earthengine-api", "ee_extra", "numpy"),
#   envname  = "r-reticulate",
#   method   = "conda"
# )
# 
# 
# library(rgee)
# 
# rgee::ee_install_set_pyenv(
#   py_path = py_path,         # same as above
#   py_env  = "r-reticulate"   # just a label; will be stored as EARTHENGINE_ENV
# )
# 
# 
# 
# 
# library(rgee)
# 
# # optional check
# Sys.getenv(c("EARTHENGINE_PYTHON", "EARTHENGINE_ENV", "RETICULATE_PYTHON"))
# reticulate::py_config()
# 
# # first-time initialization
# ee_Initialize()
# 
# 
# 

library(rgee)
ee_check()

ee_check_python()
ee_check_credentials()
ee_check_python_packages()



remotes::install_github("r-earthengine/rgeeExtra")


library(rgee)
ee_Initialize()

# Simple dataset (SRTM DEM)
dem <- ee$Image("CGIAR/SRTM90_V4")

# Print metadata
dem$bandNames()$getInfo()
dem$projection()$getInfo()



###############################################################
# Script: terrain_suitability_rgee.R
# Purpose: Identify and visualize terrain-suitable areas based on
#          elevation, slope, and aspect thresholds using 
#          Google Earth Engine (GEE) via rgee.
#
# Workflow:
#   1. Draw a region of interest (ROI) interactively
#   2. Convert the drawn ROI to Earth Engine FeatureCollection
#   3. Load SRTM (USGS SRTMGL1_003) elevation data
#   4. Derive slope and aspect layers using ee$Terrain
#   5. Create suitability masks:
#        - Elevation >= 2000 m
#        - Slope <= 30 degrees
#        - Aspect <= 90 degrees
#   6. Combine masks to identify suitable terrain pixels
#   7. Apply mask to elevation raster
#   8. Visualize suitable areas on the interactive map
#
# Author: Wyclife Agumba Oluoch
# YouTube: https://www.youtube.com/@wycology
# (Like, share, subscribe, comment)
# Date: 2025-11-14
###############################################################
# ---------------------------------------------------------------------------
# Load required libraries
# ---------------------------------------------------------------------------
library(rgee)
library(mapedit)  # no longer strictly needed if you don't draw, but harmless
library(dplyr)
library(sf)

ee_Initialize()

# ---------------------------------------------------------------------------
# 1. Use your own shapefile as ROI
# ---------------------------------------------------------------------------

# 1.1 Read your shapefile (adjust the path and filename)
roi_sf <- st_read("shp/NEAExtDissolve.shp")

# 1.2 Reproject to WGS84 (EPSG:4326) – GEE works in lat/long
roi_sf <- st_transform(roi_sf, 4326)

# (Optional) keep only the fields you care about and give a name field
roi_sf <- roi_sf |>
  mutate(name = "My_ROI") |>
  select(name)

# 1.3 Convert sf object to Earth Engine FeatureCollection
roi <- sf_as_ee(roi_sf)

# ---------------------------------------------------------------------------
# 2. Load Elevation Data (SRTM 30m)
# ---------------------------------------------------------------------------
elev <- ee$Image$Dataset$USGS_SRTMGL1_003$clip(roi)

slope  <- ee$Terrain$slope(elev)
aspect <- ee$Terrain$aspect(elev)

# ---------------------------------------------------------------------------
# 3. Create suitability masks
# ---------------------------------------------------------------------------
mask_elev   <- elev$gte(1000)    # Elevation >= 2000 m
mask_slope  <- slope$lte(30)     # Slope <= 30°
mask_aspect <- aspect$lte(180)   # Aspect <= 180° (e.g. E–S facing)

# ---------------------------------------------------------------------------
# 4. Combine masks and apply
# ---------------------------------------------------------------------------
suitable_mask <- mask_elev$
  And(mask_slope)$
  And(mask_aspect)

suitable_area <- elev$updateMask(suitable_mask)

# ---------------------------------------------------------------------------
# 5. Visualise
# ---------------------------------------------------------------------------
Map$centerObject(roi, zoom = 12)

Map$addLayer(
  suitable_area,
  list(palette = "green"),
  name = "Suitable terrain (elev>=2000, slope<=30, aspect<=180)"
) +
  Map$addLayer(
    eeObject = roi$style(
      color     = "red",
      width     = 4,
      fillColor = "#ffffff00"
    ),
    name = "ROI"
  )




























###############################################################
# Script: terraclimate_rgee_example.R
# Purpose: Example workflow for using TerraClimate in rgee
# Dataset: IDAHO_EPSCOR/TERRACLIMATE
#
# Steps:
#   1. Initialize Earth Engine (with Drive support for export)
#   2. Define region of interest (ROI) from a local shapefile
#   3. Load TerraClimate ImageCollection
#   4. Visualise a monthly TerraClimate layer (tmmx)
#   5. Derive a monthly PDSI time series over the ROI
#   6. Create annual precipitation (pr) composites
#   7. Export an annual composite to Google Drive
###############################################################

# -------------------------------------------------------------
# 0. Packages and initialization
# -------------------------------------------------------------
library(rgee)
library(sf)
library(dplyr)

# If you also want convenient extract/export helpers:
# install.packages("rgeeExtra")   # run once
# library(rgeeExtra)

# Initialize Earth Engine (enable Drive if you plan to export)
ee_Initialize(drive = FALSE)

# -------------------------------------------------------------
# 1. Region of interest (ROI) from local shapefile
# -------------------------------------------------------------
# Adjust this path to your own shapefile
roi_sf <- st_read("shp/NEAExtDissolve.shp")

# Ensure CRS is WGS84 (EPSG:4326) – required for GEE
roi_sf <- st_transform(roi_sf, 4326)

# Optionally keep only a simple attribute and rename
roi_sf <- roi_sf |>
  mutate(name = "ROI") |>
  select(name)

# Convert sf polygon to Earth Engine FeatureCollection
roi <- sf_as_ee(roi_sf)

# -------------------------------------------------------------
# 2. TerraClimate ImageCollection
# -------------------------------------------------------------
# TerraClimate ID from the Earth Engine catalogue
terraclim <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")

# We will often use this time window for examples
start_date <- "2001-01-01"
end_date   <- "2020-12-31"

# Optionally restrict by time
terraclim_sub <- terraclim$filterDate(start_date, end_date)

# -------------------------------------------------------------
# 3. Example: visualise monthly maximum temperature (tmmx)
# -------------------------------------------------------------
# TerraClimate monthly maximum air temperature (tmmx)
# Units: 0.1 °C (scale factor 0.1)
tmmx_example <- terraclim$
  filterDate("2017-07-01", "2017-08-01")$
  select("tmmx")$
  first()

tmmx_vis <- list(
  min     = -300,   # -30.0 °C (in 0.1 °C units)
  max     =  300,   #  30.0 °C
  palette = c(
    "1a3678", "2955bc", "5699ff", "8dbae9", "acd1ff", "caebff",
    "e5f9ff", "fdffb4", "ffe6a2", "ffc969", "ffa12d", "ff7c1f",
    "ca531a", "ff0000", "ab0000"
  )
)

Map$centerObject(roi, zoom = 5)

Map$addLayer(
  tmmx_example$clip(roi),
  tmmx_vis,
  "TerraClimate tmmx July 2017 (0.1 °C)"
) +
  Map$addLayer(
    roi$style(color = "black", width = 2, fillColor = "#ffffff00"),
    name = "ROI"
  )

# -------------------------------------------------------------
# 4. Monthly PDSI time series over the ROI
# -------------------------------------------------------------
# PDSI band: "pdsi" (scale factor 0.01)
tc_pdsi <- terraclim_sub$select("pdsi")

# We create a FeatureCollection where each feature holds:
#   - date (YYYY-MM-dd)
#   - mean PDSI over the ROI
pdsi_fc <- tc_pdsi$map(
  ee_utils_pyfunc(function(img) {
    # Mean over ROI polygon
    mean_dict <- img$reduceRegion(
      reducer  = ee$Reducer$mean(),
      geometry = roi$geometry(),
      scale    = 4638.3,   # native TerraClimate resolution
      maxPixels = 1e13
    )
    # Attach date as property
    img_date <- img$date()$format("YYYY-MM-dd")
    ee$Feature(NULL, mean_dict$set("date", img_date))
  })
)

# Bring the time series into R as a data.frame
pdsi_list <- pdsi_fc$getInfo()$features

pdsi_df <- lapply(pdsi_list, function(f) {
  props <- f$properties
  data.frame(
    date = as.Date(props$date),
    pdsi = as.numeric(props$pdsi),   # still in raw units (0.01 scale)
    stringsAsFactors = FALSE
  )
}) |>
  bind_rows() |>
  arrange(date) |>
  mutate(pdsi = pdsi * 0.01)  # convert to actual PDSI units

# Now pdsi_df is a monthly time series for the ROI:
head(pdsi_df)
plot(pdsi ~ date, data = pdsi_df, type = "l")

# -------------------------------------------------------------
# 5. Annual precipitation (pr) composites (2001–2020)
# -------------------------------------------------------------
# Precipitation band: "pr" (mm, monthly accumulation)
tc_pr <- terraclim_sub$select("pr")

years <- 2001:2020

annual_pr_ic <- ee$ImageCollection$fromImages(
  ee$List(years)$map(
    ee_utils_pyfunc(function(y) {
      y <- ee$Number(y)
      start <- ee$Date$fromYMD(y, 1, 1)
      end   <- start$advance(1, "year")
      
      # Sum precipitation over the year
      img <- tc_pr$
        filterDate(start, end)$
        sum()$
        set("year", y)
      
      img
    })
  )
)

# Example: pick year 2010 and visualise
pr_2010 <- annual_pr_ic$
  filter(ee$Filter$eq("year", 2010))$
  first()

pr_vis <- list(
  min     = 0,
  max     = 3000,    # adjust to your climate
  palette = c("f7fbff", "c6dbef", "6baed6", "2171b5", "08306b")
)

Map$addLayer(
  pr_2010$clip(roi),
  pr_vis,
  "Annual precipitation 2010 (mm)"
)

# -------------------------------------------------------------
# 6. Export an annual composite to Google Drive (optional)
# -------------------------------------------------------------
# Example: export pr_2010 to Drive as GeoTIFF
# # Make sure ee_Initialize(drive = TRUE) was used
# 
# task <- ee_image_to_drive(
#   image        = pr_2010,
#   description  = "TerraClimate_pr_2010",
#   folder       = "GEE_TerraClimate",   # Drive folder name (will be created if absent)
#   fileNamePrefix = "TerraClimate_pr_2010",
#   region       = roi$geometry(),
#   scale        = 4638.3,
#   maxPixels    = 1e13
# )
# 
# task$start()

# You can monitor exports with:
# ee_monitoring(task)

###############################################################
# End of script
###############################################################




























