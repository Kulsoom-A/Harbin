# -------------------------------------------------------------
# 01_GEE_data_export.R
# Google Earth Engine / rgee export workflow for the Harbin paper
# -------------------------------------------------------------
#
# Purpose:
#   Export the public remote-sensing and hydro-climate datasets used
#   in the manuscript workflow. Outputs are intended to be downloaded
#   locally and then processed by Scripts/DataPrep.R.
#
# Notes:
#   - Run this script after authenticating Google Earth Engine.
#   - Edit ROI path, export folder, years, and scale if needed.
#   - This script documents export logic for reproducibility. Some
#     collections/bands may require adjustment depending on the exact
#     product version used in the final manuscript.

library(rgee)
library(sf)
library(dplyr)

# -------------------------------------------------------------
# 0. Earth Engine initialization
# -------------------------------------------------------------
ee_Initialize(drive = TRUE)

# -------------------------------------------------------------
# 1. Region of interest
# -------------------------------------------------------------
roi_sf <- st_read("shp/NEAExtDissolve.shp")
roi_sf <- st_transform(roi_sf, 4326)
roi_sf <- roi_sf |>
  mutate(name = "Northeast_Asia") |>
  select(name)

roi <- sf_as_ee(roi_sf)

export_folder <- "Harbin_GEE_exports"
years <- 2001:2024

# -------------------------------------------------------------
# 2. TerraClimate annual exports
# -------------------------------------------------------------
terraclim <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")

make_annual_tc <- function(variable, year) {
  start <- ee$Date$fromYMD(year, 1, 1)
  end <- start$advance(1, "year")
  ic <- terraclim$filterDate(start, end)$select(variable)

  if (variable %in% c("pr", "aet", "def", "pet")) {
    img <- ic$sum()
  } else {
    img <- ic$mean()
  }

  img$set("year", year)$rename(paste0(variable, "_", year))
}

export_tc_variable <- function(variable) {
  annual_images <- lapply(years, function(y) make_annual_tc(variable, y))
  annual_stack <- ee$Image$cat(annual_images)$clip(roi)

  task <- ee_image_to_drive(
    image = annual_stack,
    description = paste0("TC_Annual_", variable, "_2001_2024_NEA"),
    folder = export_folder,
    fileNamePrefix = paste0("TC_Annual_", variable, "_2001_2024_NEA"),
    region = roi$geometry(),
    scale = 4638.3,
    maxPixels = 1e13
  )

  task$start()
  task
}

tc_variables <- c("pr", "pet", "aet", "def", "pdsi", "tmmn", "tmmx")
tc_tasks <- lapply(tc_variables, export_tc_variable)

# -------------------------------------------------------------
# 3. MODIS vegetation index annual exports
# -------------------------------------------------------------
# MOD13A2 contains 16-day NDVI at 1 km. The workflow later aligns
# exported vegetation rasters to the analysis grid in DataPrep.R.
modis_vi <- ee$ImageCollection("MODIS/061/MOD13A2")$select("NDVI")

make_annual_ndvi <- function(year) {
  start <- ee$Date$fromYMD(year, 1, 1)
  end <- start$advance(1, "year")
  modis_vi$
    filterDate(start, end)$
    mean()$
    multiply(0.0001)$
    rename(paste0("NDVI_", year))$
    set("year", year)
}

ndvi_images <- lapply(years, make_annual_ndvi)
ndvi_stack <- ee$Image$cat(ndvi_images)$clip(roi)

ndvi_task <- ee_image_to_drive(
  image = ndvi_stack,
  description = "NDVI_annual_2001_2024_NEA",
  folder = export_folder,
  fileNamePrefix = "NDVI_annual_2001_2024_NEA",
  region = roi$geometry(),
  scale = 1000,
  maxPixels = 1e13
)
ndvi_task$start()

# kNDVI is computed from NDVI as tanh(NDVI^2).
kndvi_images <- lapply(years, function(year) {
  ndvi <- make_annual_ndvi(year)
  ndvi$pow(2)$tanh()$rename(paste0("kNDVI_", year))$set("year", year)
})

kndvi_stack <- ee$Image$cat(kndvi_images)$clip(roi)

kndvi_task <- ee_image_to_drive(
  image = kndvi_stack,
  description = "kNDVI_annual_2001_2024_NEA",
  folder = export_folder,
  fileNamePrefix = "kNDVI_annual_2001_2024_NEA",
  region = roi$geometry(),
  scale = 1000,
  maxPixels = 1e13
)
kndvi_task$start()

# -------------------------------------------------------------
# 4. Monitor export tasks
# -------------------------------------------------------------
ee_monitoring()

