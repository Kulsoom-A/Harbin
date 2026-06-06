# -------------------------------------------------------------
# 01_prep_data.R
# Prepare NEA climate + vegetation stacks with proper time axis
# -------------------------------------------------------------
library(terra)

# -------------------------------------------------
# 0. Study area
# -------------------------------------------------
shape <- vect("shp/NEAFinal.shp") |> project("EPSG:4326")

# -------------------------------------------------
# 1. NDVI and kNDVI (annual)
# -------------------------------------------------
kndvi <- rast("data/kNDVI_NEA_annual_2001_2024_cf.nc")
ndvi  <- rast("data/NDVI_NEA_annual_2001_2024_cf.nc")

years_annual <- 2001:2024
dates_annual <- as.Date(paste0(years_annual, "-01-01"))

stopifnot(nlyr(kndvi) == length(years_annual))
stopifnot(nlyr(ndvi)  == length(years_annual))

time(kndvi) <- dates_annual
time(ndvi)  <- dates_annual

names(kndvi) <- paste0("kNDVI_", years_annual)
names(ndvi)  <- paste0("NDVI_",  years_annual)

# crop/mask NDVI grids to NEA final mask
kndvi <- crop(kndvi, shape) |> mask(shape)
ndvi  <- crop(ndvi,  shape) |> mask(shape)

# use first NDVI layer as the master grid
ref <- ndvi[[1]]

# -------------------------------------------------
# small helper: align any raster stack to NDVI grid
# -------------------------------------------------
align_to_ref <- function(r, ref, shp, method = "bilinear") {
  r |>
    # make sure extent/res match NDVI
    resample(ref, method = method) |>
    # mask to final shape
    mask(shp)
}

# -------------------------------------------------
# 2. TerraClimate annual stacks (2001–2024)
#    Variables: aet, def, pdsi, pet, pr, tmn, tmx
# -------------------------------------------------
tc_dir <- "data/TC_stacks_NEA"

tc_aet  <- rast(file.path(tc_dir, "TC_Annual_aet_2001_2024_0p1deg_NEA.tif"))
tc_def  <- rast(file.path(tc_dir, "TC_Annual_def_2001_2024_0p1deg_NEA.tif"))
tc_pdsi <- rast(file.path(tc_dir, "TC_Annual_pdsi_2001_2024_0p1deg_NEA.tif"))
tc_pet  <- rast(file.path(tc_dir, "TC_Annual_pet_2001_2024_0p1deg_NEA.tif"))
tc_pr   <- rast(file.path(tc_dir, "TC_Annual_pr_2001_2024_0p1deg_NEA.tif"))
tc_tmn  <- rast(file.path(tc_dir, "TC_Annual_tmmn_2001_2024_0p1deg_NEA.tif"))
tc_tmx  <- rast(file.path(tc_dir, "TC_Annual_tmmx_2001_2024_0p1deg_NEA.tif"))

set_time_names <- function(r, prefix, years) {
  stopifnot(nlyr(r) == length(years))
  dates <- as.Date(paste0(years, "-01-01"))
  time(r)  <- dates
  names(r) <- paste0(prefix, "_", years)
  r
}

tc_aet  <- set_time_names(tc_aet,  "aet",  years_annual)
tc_def  <- set_time_names(tc_def,  "def",  years_annual)
tc_pdsi <- set_time_names(tc_pdsi, "pdsi", years_annual)
tc_pet  <- set_time_names(tc_pet,  "pet",  years_annual)
tc_pr   <- set_time_names(tc_pr,   "pr",   years_annual)
tc_tmn  <- set_time_names(tc_tmn,  "tmn",  years_annual)
tc_tmx  <- set_time_names(tc_tmx,  "tmx",  years_annual)

# align all TerraClimate rasters to NDVI grid
tc_aet  <- align_to_ref(tc_aet,  ref, shape)
tc_def  <- align_to_ref(tc_def,  ref, shape)
tc_pdsi <- align_to_ref(tc_pdsi, ref, shape)
tc_pet  <- align_to_ref(tc_pet,  ref, shape)
tc_pr   <- align_to_ref(tc_pr,   ref, shape)
tc_tmn  <- align_to_ref(tc_tmn,  ref, shape)
tc_tmx  <- align_to_ref(tc_tmx,  ref, shape)

# -------------------------------------------------
# 3. SPEI12 monthly (2001–2022)
# -------------------------------------------------
spei12 <- rast("data/SPEI12_NEA/SPEI12_NEA_monthly_2001_2022_0p1deg.tif")

# expected: 22 years * 12 months = 264 layers
dates_spei <- seq(
  from = as.Date("2001-01-01"),
  to   = as.Date("2022-12-01"),
  by   = "1 month"
)
stopifnot(nlyr(spei12) == length(dates_spei))

time(spei12)  <- dates_spei
names(spei12) <- paste0("SPEI12_", format(dates_spei, "%Y_%m"))

# align SPEI to NDVI grid as well
spei12 <- align_to_ref(spei12, ref, shape)

# -------------------------------------------------
# 4. Check alignment of grids (geometry only)
# -------------------------------------------------
check_alignment <- function(r, ref) {
  ok <- compareGeom(r, ref, stopOnError = FALSE)
  if (!ok) stop("Geometry of raster does not match reference grid.")
}

check_alignment(kndvi, ref)
check_alignment(ndvi,  ref)
check_alignment(tc_pr, ref)
check_alignment(tc_pet, ref)
check_alignment(tc_aet, ref)
check_alignment(tc_def, ref)
check_alignment(tc_pdsi, ref)
check_alignment(tc_tmn, ref)
check_alignment(tc_tmx, ref)
check_alignment(spei12, ref)

# -------------------------------------------------
# 5. Save processed versions (optional but convenient)
# -------------------------------------------------
dir.create("data/processed", showWarnings = FALSE)

writeRaster(kndvi,  "data/processed/kNDVI_annual_2001_2024_NEA.tif", overwrite = TRUE)
writeRaster(ndvi,   "data/processed/NDVI_annual_2001_2024_NEA.tif",  overwrite = TRUE)
writeRaster(tc_pr,  "data/processed/TC_pr_annual_2001_2024_NEA.tif", overwrite = TRUE)
writeRaster(tc_pet, "data/processed/TC_pet_annual_2001_2024_NEA.tif",overwrite = TRUE)
writeRaster(tc_aet, "data/processed/TC_aet_annual_2001_2024_NEA.tif",overwrite = TRUE)
writeRaster(tc_def, "data/processed/TC_def_annual_2001_2024_NEA.tif",overwrite = TRUE)
writeRaster(tc_pdsi,"data/processed/TC_pdsi_annual_2001_2024_NEA.tif",overwrite = TRUE)
writeRaster(tc_tmn, "data/processed/TC_tmn_annual_2001_2024_NEA.tif",overwrite = TRUE)
writeRaster(tc_tmx, "data/processed/TC_tmx_annual_2001_2024_NEA.tif",overwrite = TRUE)
writeRaster(spei12, "data/processed/SPEI12_monthly_2001_2022_NEA.tif",overwrite = TRUE)

# You may also keep a list in memory for the next scripts
climate_list <- list(
  pr     = tc_pr,
  pet    = tc_pet,
  aet    = tc_aet,
  def    = tc_def,
  pdsi   = tc_pdsi,
  tmn    = tc_tmn,
  tmx    = tc_tmx,
  spei12 = spei12
)
