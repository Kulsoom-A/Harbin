# -------------------------------------------------------------
# 02_trend_veg.R
# Pixel-wise Sen's slope + Mann–Kendall trend for NDVI, kNDVI
# -------------------------------------------------------------
library(terra)
library(Kendall)

# -------------------------------------------------
# 1. Load processed annual NDVI / kNDVI
# -------------------------------------------------
kndvi <- rast("data/processed/kNDVI_annual_2001_2024_NEA.tif")
ndvi  <- rast("data/processed/NDVI_annual_2001_2024_NEA.tif")

years <- 2001:2024
stopifnot(nlyr(kndvi) == length(years))
stopifnot(nlyr(ndvi)  == length(years))

# -------------------------------------------------
# 2. Pre-compute index pairs for Theil–Sen slope
# -------------------------------------------------
n_years  <- length(years)
pair_idx <- combn(n_years, 2)  # 2 x 276 matrix

# -------------------------------------------------
# 3. Pixel function: Sen slope + MK tau + MK p
#    (note the extra argument `pair_idx`)
# -------------------------------------------------
trend_fun <- function(v, pair_idx) {
  v <- as.numeric(v)
  
  # if everything is NA, or too few data points
  if (all(is.na(v)) || sum(!is.na(v)) < 8) {
    return(c(NA, NA, NA))
  }
  
  ## --- Sen's slope (Theil–Sen) ---
  slopes <- (v[pair_idx[2, ]] - v[pair_idx[1, ]]) /
    (pair_idx[2, ] - pair_idx[1, ])
  
  sen_slope <- median(slopes, na.rm = TRUE)
  
  ## --- Mann–Kendall test ---
  mk  <- Kendall::MannKendall(v)
  tau <- as.numeric(mk$tau)
  p   <- as.numeric(mk$sl)
  
  c(sen_slope = sen_slope, mk_tau = tau, mk_p = p)
}

terraOptions(progress = 1, memfrac = 0.7)
dir.create("data/metrics_veg", showWarnings = FALSE)

# -------------------------------------------------
# 4. Apply to kNDVI
# -------------------------------------------------
kndvi_trend <- app(
  kndvi,
  fun      = trend_fun,
  pair_idx = pair_idx,          # <- passed here
  cores    = 4,
  filename = "data/metrics_veg/kNDVI_trend_Sen_MK_2001_2024_NEA1.tif",
  overwrite = TRUE
)
names(kndvi_trend) <- c("kndvi_senSlope", "kndvi_MK_tau", "kndvi_MK_p")
plot(kndvi_trend[[2]])
# -------------------------------------------------
# 5. Apply to NDVI
# -------------------------------------------------
ndvi_trend <- app(
  ndvi,
  fun      = trend_fun,
  pair_idx = pair_idx,          # <- passed here
  cores    = 4,
  filename = "data/metrics_veg/NDVI_trend_Sen_MK_2001_2024_NEA.tif",
  overwrite = TRUE
)
names(ndvi_trend) <- c("ndvi_senSlope", "ndvi_MK_tau", "ndvi_MK_p")

# -------------------------------------------------
# 6. Optional significance masks
# -------------------------------------------------
kndvi_sig <- kndvi_trend[[3]] < 0.05
names(kndvi_sig) <- "kndvi_trend_sig"

ndvi_sig  <- ndvi_trend[[3]] < 0.05
names(ndvi_sig) <- "ndvi_trend_sig"

writeRaster(
  c(kndvi_sig, ndvi_sig),
  "data/metrics_veg/veg_trend_significance_p05_2001_2024_NEA.tif",
  overwrite = TRUE
)

plot(kndvi_sig)






library(terra)

# kndvi_trend: layers = c("kndvi_senSlope", "kndvi_MK_tau", "kndvi_MK_p")
k_slope <- kndvi_trend[["kndvi_senSlope"]]
k_p     <- kndvi_trend[["kndvi_MK_p"]]

# start with all NA
k_class <- k_slope
k_class[] <- NA

# baseline: 0 = non-significant where we have a p-value
k_class <- ifel(!is.na(k_p), 0, NA)

# 1: weak increase (0.05 < p <= 0.10, slope > 0)
k_class <- ifel(k_p <= 0.10 & k_p > 0.05 & k_slope > 0, 1, k_class)

# 2: weak decrease (0.05 < p <= 0.10, slope < 0)
k_class <- ifel(k_p <= 0.10 & k_p > 0.05 & k_slope < 0, 2, k_class)

# 3: moderate increase (0.01 < p <= 0.05, slope > 0)
k_class <- ifel(k_p <= 0.05 & k_p > 0.01 & k_slope > 0, 3, k_class)

# 4: moderate decrease (0.01 < p <= 0.05, slope < 0)
k_class <- ifel(k_p <= 0.05 & k_p > 0.01 & k_slope < 0, 4, k_class)

# 5: strong increase (p <= 0.01, slope > 0)
k_class <- ifel(k_p <= 0.01 & k_slope > 0, 5, k_class)

# 6: strong decrease (p <= 0.01, slope < 0)
k_class <- ifel(k_p <= 0.01 & k_slope < 0, 6, k_class)

names(k_class) <- "kNDVI_trend_class"

# optional: raster attribute table for readability
rat_k <- data.frame(
  ID = 0:6,
  class = c(
    "NS",
    "Inc_p<=0.10",
    "Dec_p<=0.10",
    "Inc_p<=0.05",
    "Dec_p<=0.05",
    "Inc_p<=0.01",
    "Dec_p<=0.01"
  )
)
levels(k_class) <- rat_k

writeRaster(
  k_class,
  "data/metrics_veg/kNDVI_trend_class_7cat_2001_2024_NEA.tif",
  overwrite = TRUE
)

# quick visual check
plot(k_class)












n_slope <- ndvi_trend[["ndvi_senSlope"]]
n_p     <- ndvi_trend[["ndvi_MK_p"]]

n_class <- n_slope
n_class[] <- NA

n_class <- ifel(!is.na(n_p), 0, NA)

n_class <- ifel(n_p <= 0.10 & n_p > 0.05 & n_slope > 0, 1, n_class)
n_class <- ifel(n_p <= 0.10 & n_p > 0.05 & n_slope < 0, 2, n_class)

n_class <- ifel(n_p <= 0.05 & n_p > 0.01 & n_slope > 0, 3, n_class)
n_class <- ifel(n_p <= 0.05 & n_p > 0.01 & n_slope < 0, 4, n_class)

n_class <- ifel(n_p <= 0.01 & n_slope > 0, 5, n_class)
n_class <- ifel(n_p <= 0.01 & n_slope < 0, 6, n_class)

names(n_class) <- "NDVI_trend_class"

rat_n <- data.frame(
  ID = 0:6,
  class = c(
    "NS",
    "Inc_p<=0.10",
    "Dec_p<=0.10",
    "Inc_p<=0.05",
    "Dec_p<=0.05",
    "Inc_p<=0.01",
    "Dec_p<=0.01"
  )
)
levels(n_class) <- rat_n

writeRaster(
  n_class,
  "data/metrics_veg/NDVI_trend_class_7cat_2001_2024_NEA.tif",
  overwrite = TRUE
)

plot(n_class)





# -------------------------------------------------------------
# 03_forest_mask.R
# Align LC_Type5 forest layer and mask veg/trend rasters to forest
# -------------------------------------------------------------
library(terra)

# Study area and reference grid
shape <- vect("shp/NEAFinal.shp") |> project("EPSG:4326")

ndvi  <- rast("data/processed/NDVI_annual_2001_2024_NEA.tif")
kndvi <- rast("data/processed/kNDVI_annual_2001_2024_NEA.tif")
ref   <- ndvi[[1]]   # master grid (0.1° NEA)

# Raw forest type (LC_Type5 1–6, 0 = non-forest)
fType_raw <- rast("data/LC5_forest_1to6_2024_11km.tif")

# Check original classes (you already did this, but keep for reproducibility)
print(unique(fType_raw))

# Align to NDVI grid: nearest-neighbour for categorical data
fType <- fType_raw |>
  project(ref, method = "near") |>
  resample(ref, method = "near") |>
  crop(shape) |>
  mask(shape)

names(fType) <- "LC_Type5"

# Remove class 0 (non-forest): set to NA, keep 1–6
fType[fType == 0] <- NA

# Optional: attach labels to LC_Type5 = 1–6 if you want
# (edit these labels according to your LC_Type5 legend)
levels(fType) <- data.frame(
  ID    = 1:6,
  class = c(
    "ForestType1",
    "ForestType2",
    "ForestType3",
    "ForestType4",
    "ForestType5",
    "ForestType6"
  )
)

dir.create("data/processed", showWarnings = FALSE)

writeRaster(
  fType,
  "data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif",
  overwrite = TRUE
)

plot(fType, main = "LC_Type5 forest (1–6), non-forest removed")


















# -------------------------------------------------
# 2. Forest-masked NDVI / kNDVI and trend products
# -------------------------------------------------

# Reload if needed
ndvi        <- rast("data/processed/NDVI_annual_2001_2024_NEA.tif")
kndvi       <- rast("data/processed/kNDVI_annual_2001_2024_NEA.tif")
kndvi_trend <- rast("data/metrics_veg/kNDVI_trend_Sen_MK_2001_2024_NEA1.tif")
ndvi_trend  <- rast("data/metrics_veg/NDVI_trend_Sen_MK_2001_2024_NEA.tif")
kndvi_class <- rast("data/metrics_veg/kNDVI_trend_class_7cat_2001_2024_NEA.tif")
ndvi_class  <- rast("data/metrics_veg/NDVI_trend_class_7cat_2001_2024_NEA.tif")

fType <- rast("data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif")

# Mask everything to forest (LC_Type5 1–6)
ndvi_forest        <- mask(ndvi,        fType)
kndvi_forest       <- mask(kndvi,       fType)
kndvi_trend_forest <- mask(kndvi_trend, fType)
ndvi_trend_forest  <- mask(ndvi_trend,  fType)
kndvi_class_forest <- mask(kndvi_class, fType)
ndvi_class_forest  <- mask(ndvi_class,  fType)

dir.create("data/forest_only", showWarnings = FALSE)

writeRaster(
  ndvi_forest,
  "data/forest_only/NDVI_annual_2001_2024_NEA_forestOnly.tif",
  overwrite = TRUE
)
writeRaster(
  kndvi_forest,
  "data/forest_only/kNDVI_annual_2001_2024_NEA_forestOnly.tif",
  overwrite = TRUE
)

writeRaster(
  kndvi_trend_forest,
  "data/forest_only/kNDVI_trend_Sen_MK_2001_2024_NEA_forestOnly.tif",
  overwrite = TRUE
)
writeRaster(
  ndvi_trend_forest,
  "data/forest_only/NDVI_trend_Sen_MK_2001_2024_NEA_forestOnly.tif",
  overwrite = TRUE
)

writeRaster(
  kndvi_class_forest,
  "data/forest_only/kNDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif",
  overwrite = TRUE
)
writeRaster(
  ndvi_class_forest,
  "data/forest_only/NDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif",
  overwrite = TRUE
)

# quick visual sanity check
plot(ndvi_class_forest, main = "NDVI trend classes (forest only)")






########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################






library(terra)

# -------------------------------------------------------------
# Area (km²) and percentage of NDVI / kNDVI trend classes
# for forest-only rasters
# -------------------------------------------------------------

# 1. Load classified rasters (forest only)
ndvi_class_forest  <- rast("data/forest_only/NDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
kndvi_class_forest <- rast("data/forest_only/kNDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")

# 2. Descriptions for each class code
trend_desc_long <- c(
  "NS"          = "Non-significant trend (p > 0.10)",
  "Inc_p<=0.10" = "Positive trend, 0.05 < p \u2264 0.10",
  "Dec_p<=0.10" = "Negative trend, 0.05 < p \u2264 0.10",
  "Inc_p<=0.05" = "Positive trend, 0.01 < p \u2264 0.05",
  "Dec_p<=0.05" = "Negative trend, 0.01 < p \u2264 0.05",
  "Inc_p<=0.01" = "Positive trend, p \u2264 0.01",
  "Dec_p<=0.01" = "Negative trend, p \u2264 0.01"
)

# 3. Helper to compute area + percentage table for one raster
make_area_table <- function(r, index_name, desc_lookup) {
  ex <- expanse(r, unit = "km", byValue = TRUE)
  # 'value' may be factor; convert to character for safe matching
  ex$value <- as.character(ex$value)
  
  total_area <- sum(ex$area, na.rm = TRUE)
  ex$percent <- 100 * ex$area / total_area
  ex$description <- desc_lookup[ex$value]
  
  out <- ex[, c("value", "description", "area", "percent")]
  names(out) <- c("class_code", "class_description", "area_km2", "percent_forest")
  
  out$index <- index_name
  out$area_km2       <- round(out$area_km2, 1)
  out$percent_forest <- round(out$percent_forest, 2)
  
  out <- out[, c("index", "class_code", "class_description", "area_km2", "percent_forest")]
  out
}

# 4. Build tables for NDVI and kNDVI
ndvi_area_tbl  <- make_area_table(ndvi_class_forest,  "NDVI",  trend_desc_long)
kndvi_area_tbl <- make_area_table(kndvi_class_forest, "kNDVI", trend_desc_long)

# 5. Combined table
veg_trend_area_combined <- rbind(ndvi_area_tbl, kndvi_area_tbl)

# 6. Save to disk (ensure folders exist)
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

write.csv(
  ndvi_area_tbl,
  "outputs/tables/NDVI_trend_class_area_forestOnly_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  kndvi_area_tbl,
  "outputs/tables/kNDVI_trend_class_area_forestOnly_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  veg_trend_area_combined,
  "outputs/tables/NDVI_kNDVI_trend_class_area_forestOnly_2001_2024_combined.csv",
  row.names = FALSE
)

# 7. Inspect in console if you want
ndvi_area_tbl
kndvi_area_tbl
veg_trend_area_combined
















###############################################################
# 06_forestType_veg_trend_area.R
# Area (km²) and % of forest type in each NDVI/kNDVI trend class
###############################################################


###############################################################
# 06_forestType_veg_trend_area.R
# Area (km²) and % of forest type in each NDVI/kNDVI trend class
###############################################################





###############################################################
# 06_forestType_veg_trend_area.R
# Area (km²) and % of forest type in each NDVI/kNDVI trend class
###############################################################




library(terra)

# -------------------------------------------------------------
# 1. Load forest-type raster and trend-class rasters (forest only)
# -------------------------------------------------------------
# Forest type: LC_Type5 = 1–6, NA elsewhere
fType <- rast("data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif")

# NDVI / kNDVI trend classes, forest-only
ndvi_class_forest  <- rast("data/forest_only/NDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
kndvi_class_forest <- rast("data/forest_only/kNDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")

# -------------------------------------------------------------
# 2. Helper: build LC_Type5 × trend-class area table for one index
# -------------------------------------------------------------
make_ftype_trend_table <- function(class_raster, index_name, fType) {
  
  # encode unique combination: 10*forest_type + trend_code
  combo <- fType * 10 + class_raster
  
  # area per combination (km²)
  ex <- expanse(combo, unit = "km", byValue = TRUE)
  ex <- as.data.frame(ex)
  ex <- ex[!is.na(ex$value), ]  # remove NA combos if any
  
  # decode forest type and trend code
  ex$forest_type <- floor(ex$value / 10)
  ex$trend_code  <- ex$value %% 10
  
  # map trend codes to short and long labels
  trend_short <- c(
    "0" = "NS",
    "1" = "Inc_p<=0.10",
    "2" = "Dec_p<=0.10",
    "3" = "Inc_p<=0.05",
    "4" = "Dec_p<=0.05",
    "5" = "Inc_p<=0.01",
    "6" = "Dec_p<=0.01"
  )
  
  trend_long <- c(
    "NS"          = "Non-significant trend (p > 0.10)",
    "Inc_p<=0.10" = "Positive trend, 0.05 < p \u2264 0.10",
    "Dec_p<=0.10" = "Negative trend, 0.05 < p \u2264 0.10",
    "Inc_p<=0.05" = "Positive trend, 0.01 < p \u2264 0.05",
    "Dec_p<=0.05" = "Negative trend, 0.01 < p \u2264 0.05",
    "Inc_p<=0.01" = "Positive trend, p \u2264 0.01",
    "Dec_p<=0.01" = "Negative trend, p \u2264 0.01"
  )
  
  ex$trend_class <- trend_short[as.character(ex$trend_code)]
  ex$trend_desc  <- trend_long[ex$trend_class]
  ex$index       <- index_name
  
  # total forest area across all forest types for this index
  total_all <- sum(ex$area, na.rm = TRUE)
  
  # percentage of ALL forest area
  ex$perc_all_forest <- 100 * ex$area / total_all
  
  # percentage WITHIN each forest type (rows for a given type sum to 100)
  total_by_type <- tapply(ex$area, ex$forest_type, sum, na.rm = TRUE)
  ex$perc_within_type <- 100 * ex$area /
    total_by_type[as.character(ex$forest_type)]
  
  # rounding and tidy columns
  ex$area_km2          <- round(ex$area, 1)
  ex$perc_all_forest   <- round(ex$perc_all_forest, 2)
  ex$perc_within_type  <- round(ex$perc_within_type, 2)
  
  out <- ex[, c("index",
                "forest_type",
                "trend_code",
                "trend_class",
                "trend_desc",
                "area_km2",
                "perc_all_forest",
                "perc_within_type")]
  
  out[order(out$forest_type, out$trend_code), ]
}

# -------------------------------------------------------------
# 3. Build tables for NDVI and kNDVI
# -------------------------------------------------------------
ndvi_ftype_tbl  <- make_ftype_trend_table(ndvi_class_forest,  "NDVI",  fType)
kndvi_ftype_tbl <- make_ftype_trend_table(kndvi_class_forest, "kNDVI", fType)

combined_ftype_tbl <- rbind(ndvi_ftype_tbl, kndvi_ftype_tbl)

# -------------------------------------------------------------
# 4. Save to disk
# -------------------------------------------------------------
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

write.csv(
  ndvi_ftype_tbl,
  "outputs/tables/NDVI_trend_class_by_forestType_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  kndvi_ftype_tbl,
  "outputs/tables/kNDVI_trend_class_by_forestType_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  combined_ftype_tbl,
  "outputs/tables/NDVI_kNDVI_trend_class_by_forestType_2001_2024_combined.csv",
  row.names = FALSE
)

# Optional: inspect in console
ndvi_ftype_tbl
kndvi_ftype_tbl
combined_ftype_tbl
###############################################################



NAEcntry <- vect("shp/rough.shp")

str(NAEcntry)
NAEcntry
names(NAEcntry)








########################################################################
############### New analysis wid countries included





########################################################################
############### New analysis wid countries included





########################################################################
############### New analysis wid countries included














###############################################################
# 06_forestType_veg_trend_area.R (extended with labels)
###############################################################

library(terra)

# -------------------------------------------------------------
# 1. Load forest-type raster and trend-class rasters (forest only)
# -------------------------------------------------------------
# Forest type: LC_Type5 = 1–6, NA elsewhere
fType <- rast("data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif")
plot(fType)
# NDVI / kNDVI trend classes, forest-only
ndvi_class_forest  <- rast("data/forest_only/NDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
kndvi_class_forest <- rast("data/forest_only/kNDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")

# -------------------------------------------------------------
# 2. Helper: build LC_Type5 × trend-class area table for one index
# -------------------------------------------------------------
make_ftype_trend_table <- function(class_raster, index_name, fType) {
  
  # encode unique combination: 10*forest_type + trend_code
  combo <- fType * 10 + class_raster
  
  # area per combination (km²)
  ex <- expanse(combo, unit = "km", byValue = TRUE)
  ex <- as.data.frame(ex)
  ex <- ex[!is.na(ex$value), ]  # remove NA combos if any
  
  # decode forest type and trend code
  ex$forest_type <- floor(ex$value / 10)
  ex$trend_code  <- ex$value %% 10
  
  # map trend codes to short and long labels
  trend_short <- c(
    "0" = "NS",
    "1" = "Inc_p<=0.10",
    "2" = "Dec_p<=0.10",
    "3" = "Inc_p<=0.05",
    "4" = "Dec_p<=0.05",
    "5" = "Inc_p<=0.01",
    "6" = "Dec_p<=0.01"
  )
  
  trend_long <- c(
    "NS"          = "Non-significant trend (p > 0.10)",
    "Inc_p<=0.10" = "Positive trend, 0.05 < p \u2264 0.10",
    "Dec_p<=0.10" = "Negative trend, 0.05 < p \u2264 0.10",
    "Inc_p<=0.05" = "Positive trend, 0.01 < p \u2264 0.05",
    "Dec_p<=0.05" = "Negative trend, 0.01 < p \u2264 0.05",
    "Inc_p<=0.01" = "Positive trend, p \u2264 0.01",
    "Dec_p<=0.01" = "Negative trend, p \u2264 0.01"
  )
  
  ex$trend_class <- trend_short[as.character(ex$trend_code)]
  ex$trend_desc  <- trend_long[ex$trend_class]
  ex$index       <- index_name
  
  # forest-type labels for LC_Type5 = 1–6
  ftype_labels <- c(
    "1" = "Evergreen needleleaf forest",
    "2" = "Evergreen broadleaf forest",
    "3" = "Deciduous needleleaf forest",
    "4" = "Deciduous broadleaf forest",
    "5" = "Shrublands",
    "6" = "Grasslands"
  )
  ex$forest_type_label <- ftype_labels[as.character(ex$forest_type)]
  
  # total forest area across all forest types for this index
  total_all <- sum(ex$area, na.rm = TRUE)
  
  # percentage of ALL forest area
  ex$perc_all_forest <- 100 * ex$area / total_all
  
  # percentage WITHIN each forest type (rows for a given type sum to 100)
  total_by_type <- tapply(ex$area, ex$forest_type, sum, na.rm = TRUE)
  ex$perc_within_type <- 100 * ex$area /
    total_by_type[as.character(ex$forest_type)]
  
  # rounding and tidy columns
  ex$area_km2          <- round(ex$area, 1)
  ex$perc_all_forest   <- round(ex$perc_all_forest, 2)
  ex$perc_within_type  <- round(ex$perc_within_type, 2)
  
  out <- ex[, c("index",
                "forest_type",
                "forest_type_label",
                "trend_code",
                "trend_class",
                "trend_desc",
                "area_km2",
                "perc_all_forest",
                "perc_within_type")]
  
  out[order(out$forest_type, out$trend_code), ]
}

# -------------------------------------------------------------
# 3. Build tables for NDVI and kNDVI
# -------------------------------------------------------------
ndvi_ftype_tbl  <- make_ftype_trend_table(ndvi_class_forest,  "NDVI",  fType)
kndvi_ftype_tbl <- make_ftype_trend_table(kndvi_class_forest, "kNDVI", fType)

combined_ftype_tbl <- rbind(ndvi_ftype_tbl, kndvi_ftype_tbl)

# -------------------------------------------------------------
# 4. Save to disk (new folder tables2)
# -------------------------------------------------------------
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables2", showWarnings = FALSE, recursive = TRUE)

write.csv(
  ndvi_ftype_tbl,
  "outputs/tables2/NDVI_trend_class_by_forestType_2001_2024_labeled.csv",
  row.names = FALSE
)

write.csv(
  kndvi_ftype_tbl,
  "outputs/tables2/kNDVI_trend_class_by_forestType_2001_2024_labeled.csv",
  row.names = FALSE
)

write.csv(
  combined_ftype_tbl,
  "outputs/tables2/NDVI_kNDVI_trend_class_by_forestType_2001_2024_combined_labeled.csv",
  row.names = FALSE
)

# Optional check
ndvi_ftype_tbl
kndvi_ftype_tbl
combined_ftype_tbl







###############################################################
# 06_country_veg_trend_area.R
# Area (km²) and % of forest in each NDVI/kNDVI trend class
# stratified by country (China + Taiwan merged)
###############################################################

library(terra)

# -------------------------------------------------------------
# 1. Load trend-class rasters (forest only) and country polygons
# -------------------------------------------------------------
ndvi_class_forest  <- rast("data/forest_only/NDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
kndvi_class_forest <- rast("data/forest_only/kNDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")

# Country shapefile
NAEcntry <- vect("shp/rough.shp")

# Ensure same CRS as trend rasters
NAEcntry <- project(NAEcntry, crs(ndvi_class_forest))

# Inspect attribute names once in console to set this correctly:
# names(NAEcntry)
country_field <- "NAME"  # <-- change if your country name field is different

# -------------------------------------------------------------
# 2. Merge China + Taiwan into one category
# -------------------------------------------------------------
NAEcntry$COUNTRY_MRG <- NAEcntry[[country_field]]

# Values to be merged into "China" (adapt to your attribute values)
china_like <- c("China", "Mainland China", "Taiwan", "Taiwan Province of China")

NAEcntry$COUNTRY_MRG[NAEcntry[[country_field]] %in% china_like] <- "China"

# Create integer country IDs
cnames <- sort(unique(NAEcntry$COUNTRY_MRG))
NAEcntry$country_id <- match(NAEcntry$COUNTRY_MRG, cnames)

# -------------------------------------------------------------
# 3. Rasterise country_id onto forest trend grid
# -------------------------------------------------------------
ref <- ndvi_class_forest[[1]]

country_id_rast <- rasterize(NAEcntry, ref, field = "country_id")

# Mask to forest-only (trend rasters are already forest-only, this is just to be safe)
country_id_rast <- mask(country_id_rast, ndvi_class_forest)

# Optional: sanity check geometry
stopifnot(compareGeom(country_id_rast, ndvi_class_forest, stopOnError = FALSE))

# -------------------------------------------------------------
# 4. Helper: build country × trend-class area table
# -------------------------------------------------------------
make_country_trend_table <- function(class_raster, index_name, country_rast, country_names) {
  
  # encode unique combination: 10 * country_id + trend_code
  combo <- country_rast * 10 + class_raster
  
  ex <- expanse(combo, unit = "km", byValue = TRUE)
  ex <- as.data.frame(ex)
  ex <- ex[!is.na(ex$value), ]
  
  # decode ids
  ex$country_id <- floor(ex$value / 10)
  ex$trend_code <- ex$value %% 10
  
  # trend labels
  trend_short <- c(
    "0" = "NS",
    "1" = "Inc_p<=0.10",
    "2" = "Dec_p<=0.10",
    "3" = "Inc_p<=0.05",
    "4" = "Dec_p<=0.05",
    "5" = "Inc_p<=0.01",
    "6" = "Dec_p<=0.01"
  )
  
  trend_long <- c(
    "NS"          = "Non-significant trend (p > 0.10)",
    "Inc_p<=0.10" = "Positive trend, 0.05 < p \u2264 0.10",
    "Dec_p<=0.10" = "Negative trend, 0.05 < p \u2264 0.10",
    "Inc_p<=0.05" = "Positive trend, 0.01 < p \u2264 0.05",
    "Dec_p<=0.05" = "Negative trend, 0.01 < p \u2264 0.05",
    "Inc_p<=0.01" = "Positive trend, p \u2264 0.01",
    "Dec_p<=0.01" = "Negative trend, p \u2264 0.01"
  )
  
  ex$trend_class <- trend_short[as.character(ex$trend_code)]
  ex$trend_desc  <- trend_long[ex$trend_class]
  ex$country     <- country_names[ex$country_id]
  ex$index       <- index_name
  
  # total forest area across ALL countries
  total_all <- sum(ex$area, na.rm = TRUE)
  ex$perc_all_forest <- 100 * ex$area / total_all
  
  # percentage WITHIN each country (rows for a given country sum to 100)
  total_by_country <- tapply(ex$area, ex$country_id, sum, na.rm = TRUE)
  ex$perc_within_country <- 100 * ex$area /
    total_by_country[as.character(ex$country_id)]
  
  # rounding and tidy columns
  ex$area_km2            <- round(ex$area, 1)
  ex$perc_all_forest     <- round(ex$perc_all_forest, 2)
  ex$perc_within_country <- round(ex$perc_within_country, 2)
  
  out <- ex[, c("index",
                "country",
                "country_id",
                "trend_code",
                "trend_class",
                "trend_desc",
                "area_km2",
                "perc_all_forest",
                "perc_within_country")]
  
  out[order(out$country, out$trend_code), ]
}

# -------------------------------------------------------------
# 5. Build tables for NDVI and kNDVI
# -------------------------------------------------------------
ndvi_country_tbl  <- make_country_trend_table(ndvi_class_forest,  "NDVI",  country_id_rast, cnames)
kndvi_country_tbl <- make_country_trend_table(kndvi_class_forest, "kNDVI", country_id_rast, cnames)

combined_country_tbl <- rbind(ndvi_country_tbl, kndvi_country_tbl)

# -------------------------------------------------------------
# 6. Save to disk in outputs/tables2
# -------------------------------------------------------------
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables2", showWarnings = FALSE, recursive = TRUE)

write.csv(
  ndvi_country_tbl,
  "outputs/tables2/NDVI_trend_class_by_country_2001_2024_forestOnly.csv",
  row.names = FALSE
)

write.csv(
  kndvi_country_tbl,
  "outputs/tables2/kNDVI_trend_class_by_country_2001_2024_forestOnly.csv",
  row.names = FALSE
)

write.csv(
  combined_country_tbl,
  "outputs/tables2/NDVI_kNDVI_trend_class_by_country_2001_2024_forestOnly_combined.csv",
  row.names = FALSE
)

# Optional: inspect in console
ndvi_country_tbl
kndvi_country_tbl
combined_country_tbl













###############################################################
# 07_kNDVI_trend_stackedBars.R
# Publication-ready stacked bar plots for kNDVI trends
# - Forest types with abbreviated labels (ENT, EBT, etc.)
# - Countries with China and Taiwan merged
###############################################################

library(readr)
library(dplyr)
library(ggplot2)
library(forcats)
library(patchwork)

# -------------------------------------------------------------
# 1. Read combined tables (forest type + country)
# -------------------------------------------------------------
ftype_df <- read_csv(
  "outputs/tables2/NDVI_kNDVI_trend_class_by_forestType_2001_2024_combined_labeled.csv"
)

country_df <- read_csv(
  "outputs/tables2/NDVI_kNDVI_trend_class_by_country_2001_2024_forestOnly_combined.csv"
)

# -------------------------------------------------------------
# 2. Focus on kNDVI and set factor orders
# -------------------------------------------------------------
trend_levels <- c(
  "Dec_p<=0.01",
  "Dec_p<=0.05",
  "Dec_p<=0.10",
  "NS",
  "Inc_p<=0.10",
  "Inc_p<=0.05",
  "Inc_p<=0.01"
)

# Forest-type order (ecological sequence)
ftype_order <- c(
  "Evergreen needleleaf forest",
  "Evergreen broadleaf forest",
  "Deciduous needleleaf forest",
  "Deciduous broadleaf forest",
  "Shrublands",
  "Grasslands"
)

# Abbreviations for forest types
ftype_abbrev_map <- c(
  "Evergreen needleleaf forest" = "ENT",
  "Evergreen broadleaf forest"  = "EBT",
  "Deciduous needleleaf forest" = "DNT",
  "Deciduous broadleaf forest"  = "DBT",
  "Shrublands"                  = "SHB",
  "Grasslands"                  = "GRS"
)

abbr_order <- c("ENT", "EBT", "DNT", "DBT", "SHB", "GRS")

# ----------------- Forest types (kNDVI only) ------------------
ftype_k <- ftype_df |>
  filter(index == "kNDVI") |>
  mutate(
    trend_class       = factor(trend_class, levels = trend_levels),
    forest_type_label = factor(forest_type_label, levels = ftype_order),
    forest_type_abbr  = ftype_abbrev_map[as.character(forest_type_label)],
    forest_type_abbr  = factor(forest_type_abbr, levels = abbr_order)
  )

# ----------------- Countries (kNDVI only) ---------------------
# Merge Taiwan into China at this stage and recompute areas/percentages
country_k <- country_df |>
  filter(index == "kNDVI") |>
  mutate(
    # merge Taiwan into China (adjust pattern if your naming differs)
    country = case_when(
      grepl("Taiwan", country, ignore.case = TRUE) ~ "China",
      TRUE ~ country
    ),
    trend_class = factor(trend_class, levels = trend_levels)
  ) |>
  # recompute area per (country × trend_class)
  group_by(country, trend_class) |>
  summarise(area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop") |>
  # compute percentages within each country
  group_by(country) |>
  mutate(
    total_area_country   = sum(area_km2, na.rm = TRUE),
    perc_within_country  = 100 * area_km2 / total_area_country
  ) |>
  ungroup() |>
  mutate(
    country             = fct_reorder(country, total_area_country, .desc = TRUE),
    area_km2            = round(area_km2, 1),
    perc_within_country = round(perc_within_country, 2)
  )

# -------------------------------------------------------------
# 3. Colour palette for trend classes (diverging)
# -------------------------------------------------------------
trend_cols <- c(
  "Dec_p<=0.01" = "#67001F",
  "Dec_p<=0.05" = "#B2182B",
  "Dec_p<=0.10" = "#D6604D",
  "NS"          = "#F7F7F7",
  "Inc_p<=0.10" = "#92C5DE",
  "Inc_p<=0.05" = "#4393C3",
  "Inc_p<=0.01" = "#053061"
)

# -------------------------------------------------------------
# 4. Forest-type stacked bar plot (kNDVI, abbreviated labels)
#     - percent within forest type
#     - labels for segments >= 3%
# -------------------------------------------------------------
p_forest <- ggplot(
  ftype_k,
  aes(x = forest_type_abbr,
      y = perc_within_type,
      fill = trend_class)
) +
  geom_col(color = "grey20", width = 0.8) +
  geom_text(
    data = ftype_k |>
      filter(perc_within_type >= 3),
    aes(label = sprintf("%.0f", perc_within_type)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(
    values = trend_cols,
    name   = "kNDVI trend class"
  ) +
  labs(
    x     = "Forest type (ENT, EBT, DNT, DBT, SHB, GRS)",
    y     = "Share of forest area within type (%)",
    title = "A"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 25, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    plot.title         = element_text(hjust = 0, face = "bold")
  )

# -------------------------------------------------------------
# 5. Country-level stacked bar plot (kNDVI, China+Taiwan merged)
#     - percent within country
#     - labels for segments >= 3%
# -------------------------------------------------------------
p_country <- ggplot(
  country_k,
  aes(x = country,
      y = perc_within_country,
      fill = trend_class)
) +
  geom_col(color = "grey20", width = 0.8) +
  geom_text(
    data = country_k |>
      filter(perc_within_country >= 3),
    aes(label = sprintf("%.0f", perc_within_country)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(
    values = trend_cols,
    name   = "kNDVI trend class"
  ) +
  labs(
    x     = NULL,
    y     = "Share of forest area within country (%)",
    title = "B"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 25, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    plot.title         = element_text(hjust = 0, face = "bold")
  )

# -------------------------------------------------------------
# 6. Panel figure with shared legend (optional)
# -------------------------------------------------------------
p_panel <- p_forest + p_country +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# -------------------------------------------------------------
# 7. Save all figures
# -------------------------------------------------------------
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

ggsave(
  "outputs/figures/Fig_kNDVI_trend_forestType_NEA_2001_2024_abbr.png",
  p_forest, width = 7, height = 4, dpi = 600
)

ggsave(
  "outputs/figures/Fig_kNDVI_trend_country_NEA_2001_2024_ChinaMerged.png",
  p_country, width = 7, height = 4, dpi = 600
)

ggsave(
  "outputs/figures/Fig_kNDVI_trend_forestType_country_NEA_2001_2024_panel.png",
  p_panel, width = 8, height = 5, dpi = 600
)





























