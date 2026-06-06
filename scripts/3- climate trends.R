###############################################################
# 04_trend_climate.R
# Hydro-climatic trends (Sen + MK) and 7-class significance
# for PR, PET, water balance (PR-PET), Tmean, and PDSI
###############################################################

library(terra)
library(Kendall)

# -------------------------------------------------------------
# 1. Load processed annual TerraClimate stacks
#    (already aligned to NDVI grid in 01_prep_data.R)
# -------------------------------------------------------------
pr   <- rast("data/processed/TC_pr_annual_2001_2024_NEA.tif")
pet  <- rast("data/processed/TC_pet_annual_2001_2024_NEA.tif")
pdsi <- rast("data/processed/TC_pdsi_annual_2001_2024_NEA.tif")
tmn  <- rast("data/processed/TC_tmn_annual_2001_2024_NEA.tif")
tmx  <- rast("data/processed/TC_tmx_annual_2001_2024_NEA.tif")

years <- 2001:2024
stopifnot(
  nlyr(pr)   == length(years),
  nlyr(pet)  == length(years),
  nlyr(pdsi) == length(years),
  nlyr(tmn)  == length(years),
  nlyr(tmx)  == length(years)
)

# -------------------------------------------------------------
# 2. Derive water balance and mean annual temperature
# -------------------------------------------------------------
wb <- pr - pet
names(wb) <- paste0("wb_", years)

tmean <- (tmn + tmx) / 2
names(tmean) <- paste0("tmean_", years)

# Optional: rename originals for clarity
names(pr)   <- paste0("pr_", years)
names(pet)  <- paste0("pet_", years)
names(pdsi) <- paste0("pdsi_", years)

# -------------------------------------------------------------
# 3. Sen slope + Mann–Kendall trend function (pixel-wise)
# -------------------------------------------------------------
n_years  <- length(years)
pair_idx <- combn(n_years, 2)  # 2 x 276 for 24 years

trend_fun <- function(v, pair_idx) {
  v <- as.numeric(v)
  
  # if everything is NA, or too few data points
  if (all(is.na(v)) || sum(!is.na(v)) < 8) {
    return(c(NA, NA, NA))
  }
  
  # Theil–Sen slopes
  slopes <- (v[pair_idx[2, ]] - v[pair_idx[1, ]]) /
    (pair_idx[2, ] - pair_idx[1, ])
  
  sen_slope <- median(slopes, na.rm = TRUE)
  
  # Mann–Kendall
  mk  <- Kendall::MannKendall(v)
  tau <- as.numeric(mk$tau)
  p   <- as.numeric(mk$sl)
  
  c(sen_slope = sen_slope, mk_tau = tau, mk_p = p)
}

terraOptions(progress = 1, memfrac = 0.7)

dir.create("data/metrics_climate", showWarnings = FALSE)

# helper to run trend analysis and name layers
make_trend_stack <- function(r, prefix, out_path) {
  tr <- app(
    r,
    fun      = trend_fun,
    pair_idx = pair_idx,
    cores    = 4,
    filename = out_path,
    overwrite = TRUE
  )
  names(tr) <- c(
    paste0(prefix, "_senSlope"),
    paste0(prefix, "_MK_tau"),
    paste0(prefix, "_MK_p")
  )
  tr
}

# -------------------------------------------------------------
# 4. Compute trend stacks for each climate variable
# -------------------------------------------------------------
pr_trend   <- make_trend_stack(pr,    "pr",    "data/metrics_climate/PR_trend_Sen_MK_2001_2024_NEA.tif")
pet_trend  <- make_trend_stack(pet,   "pet",   "data/metrics_climate/PET_trend_Sen_MK_2001_2024_NEA.tif")
wb_trend   <- make_trend_stack(wb,    "wb",    "data/metrics_climate/WB_trend_Sen_MK_2001_2024_NEA.tif")
tmean_trend<- make_trend_stack(tmean, "tmean", "data/metrics_climate/Tmean_trend_Sen_MK_2001_2024_NEA.tif")
pdsi_trend <- make_trend_stack(pdsi,  "pdsi",  "data/metrics_climate/PDSI_trend_Sen_MK_2001_2024_NEA.tif")

# -------------------------------------------------------------
# 5. 7-class direction + significance classification
#    (same scheme as NDVI/kNDVI)
# -------------------------------------------------------------
# Classes:
# 0 = Non-significant (p > 0.10 or NA)
# 1 = Weak increase   (0.05 < p <= 0.10, slope > 0)
# 2 = Weak decrease   (0.05 < p <= 0.10, slope < 0)
# 3 = Moderate inc.   (0.01 < p <= 0.05, slope > 0)
# 4 = Moderate dec.   (0.01 < p <= 0.05, slope < 0)
# 5 = Strong increase (p <= 0.01, slope > 0)
# 6 = Strong decrease (p <= 0.01, slope < 0)

classify_trend <- function(tr_stack, prefix, out_path) {
  slope <- tr_stack[[1]]  # senSlope
  p     <- tr_stack[[3]]  # MK_p
  
  cls <- slope
  cls[] <- NA
  
  # baseline: 0 for all pixels with defined p
  cls <- ifel(!is.na(p), 0, NA)
  
  cls <- ifel(p <= 0.10 & p > 0.05 & slope > 0, 1, cls)
  cls <- ifel(p <= 0.10 & p > 0.05 & slope < 0, 2, cls)
  
  cls <- ifel(p <= 0.05 & p > 0.01 & slope > 0, 3, cls)
  cls <- ifel(p <= 0.05 & p > 0.01 & slope < 0, 4, cls)
  
  cls <- ifel(p <= 0.01 & slope > 0, 5, cls)
  cls <- ifel(p <= 0.01 & slope < 0, 6, cls)
  
  names(cls) <- paste0(prefix, "_trend_class")
  
  writeRaster(cls, out_path, overwrite = TRUE)
  cls
}

pr_class    <- classify_trend(pr_trend,    "PR",    "data/metrics_climate/PR_trend_class_7cat_2001_2024_NEA.tif")
pet_class   <- classify_trend(pet_trend,   "PET",   "data/metrics_climate/PET_trend_class_7cat_2001_2024_NEA.tif")
wb_class    <- classify_trend(wb_trend,    "WB",    "data/metrics_climate/WB_trend_class_7cat_2001_2024_NEA.tif")
tmean_class <- classify_trend(tmean_trend, "Tmean", "data/metrics_climate/Tmean_trend_class_7cat_2001_2024_NEA.tif")
pdsi_class  <- classify_trend(pdsi_trend,  "PDSI",  "data/metrics_climate/PDSI_trend_class_7cat_2001_2024_NEA.tif")

# -------------------------------------------------------------
# 6. Mask climate trend classes to forest area (LC_Type5 1–6)
# -------------------------------------------------------------
fType <- rast("data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif")

dir.create("data/forest_only", showWarnings = FALSE)

mask_to_forest <- function(r, prefix, out_path) {
  r_f <- mask(r, fType)
  writeRaster(r_f, out_path, overwrite = TRUE)
  r_f
}

pr_class_forest    <- mask_to_forest(pr_class,    "PR",    "data/forest_only/PR_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
pet_class_forest   <- mask_to_forest(pet_class,   "PET",   "data/forest_only/PET_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
wb_class_forest    <- mask_to_forest(wb_class,    "WB",    "data/forest_only/WB_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
tmean_class_forest <- mask_to_forest(tmean_class, "Tmean", "data/forest_only/Tmean_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
pdsi_class_forest  <- mask_to_forest(pdsi_class,  "PDSI",  "data/forest_only/PDSI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")

# quick sanity plots (optional)
plot(wb_class_forest, main = "Water balance trend classes (forest only)")
plot(tmean_class_forest, main = "Temperature trend classes (forest only)")
###############################################################












###############################################################
# 05_area_climate_trends.R
# Area (km²) and % of forest in hydro-climatic trend classes
###############################################################

library(terra)

# -------------------------------------------------------------
# 1. Load climate trend-class rasters (forest only)
#    (created in 04_trend_climate.R)
# -------------------------------------------------------------
pr_class_forest    <- rast("data/forest_only/PR_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
pet_class_forest   <- rast("data/forest_only/PET_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
wb_class_forest    <- rast("data/forest_only/WB_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
tmean_class_forest <- rast("data/forest_only/Tmean_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
pdsi_class_forest  <- rast("data/forest_only/PDSI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")

# Optional: also load NDVI/kNDVI trend classes for later comparison
ndvi_class_forest  <- rast("data/forest_only/NDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
kndvi_class_forest <- rast("data/forest_only/kNDVI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")

# -------------------------------------------------------------
# 2. Common descriptions for the 7 trend classes
#    (sign + significance; interpretation will differ by variable)
# -------------------------------------------------------------
trend_desc_long <- c(
  "NS"          = "Non-significant trend (p > 0.10)",
  "Inc_p<=0.10" = "Positive trend, 0.05 < p \u2264 0.10",
  "Dec_p<=0.10" = "Negative trend, 0.05 < p \u2264 0.10",
  "Inc_p<=0.05" = "Positive trend, 0.01 < p \u2264 0.05",
  "Dec_p<=0.05" = "Negative trend, 0.01 < p \u2264 0.05",
  "Inc_p<=0.01" = "Positive trend, p \u2264 0.01",
  "Dec_p<=0.01" = "Negative trend, p \u2264 0.01"
)

# -------------------------------------------------------------
# 3. Helper to compute area + percentage table for one raster
# -------------------------------------------------------------
make_area_table <- function(r, index_name, desc_lookup) {
  ex <- expanse(r, unit = "km", byValue = TRUE)
  # value may be a factor -> convert to character for matching
  ex$value <- as.character(ex$value)
  
  total_area <- sum(ex$area, na.rm = TRUE)
  
  ex$percent     <- 100 * ex$area / total_area
  ex$description <- desc_lookup[ex$value]
  
  out <- ex[, c("value", "description", "area", "percent")]
  names(out) <- c("class_code", "class_description", "area_km2", "percent_forest")
  
  out$index <- index_name
  out$area_km2       <- round(out$area_km2, 1)
  out$percent_forest <- round(out$percent_forest, 2)
  
  out <- out[, c("index", "class_code", "class_description", "area_km2", "percent_forest")]
  out
}

# -------------------------------------------------------------
# 4. Build tables for each climate variable
# -------------------------------------------------------------
pr_area_tbl    <- make_area_table(pr_class_forest,    "PR",    trend_desc_long)
pet_area_tbl   <- make_area_table(pet_class_forest,   "PET",   trend_desc_long)
wb_area_tbl    <- make_area_table(wb_class_forest,    "WB",    trend_desc_long)
tmean_area_tbl <- make_area_table(tmean_class_forest, "Tmean", trend_desc_long)
pdsi_area_tbl  <- make_area_table(pdsi_class_forest,  "PDSI",  trend_desc_long)

# Optional: area tables for NDVI/kNDVI in the same structure
ndvi_area_tbl  <- make_area_table(ndvi_class_forest,  "NDVI",  trend_desc_long)
kndvi_area_tbl <- make_area_table(kndvi_class_forest, "kNDVI", trend_desc_long)

# -------------------------------------------------------------
# 5. Combined table across all variables
# -------------------------------------------------------------
climate_trend_area_combined <- rbind(
  pr_area_tbl,
  pet_area_tbl,
  wb_area_tbl,
  tmean_area_tbl,
  pdsi_area_tbl
)

veg_trend_area_combined <- rbind(
  ndvi_area_tbl,
  kndvi_area_tbl
)

all_trend_area_combined <- rbind(
  climate_trend_area_combined,
  veg_trend_area_combined
)

# -------------------------------------------------------------
# 6. Save to disk
# -------------------------------------------------------------
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

write.csv(
  pr_area_tbl,
  "outputs/tables/PR_trend_class_area_forestOnly_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  pet_area_tbl,
  "outputs/tables/PET_trend_class_area_forestOnly_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  wb_area_tbl,
  "outputs/tables/WB_trend_class_area_forestOnly_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  tmean_area_tbl,
  "outputs/tables/Tmean_trend_class_area_forestOnly_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  pdsi_area_tbl,
  "outputs/tables/PDSI_trend_class_area_forestOnly_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  climate_trend_area_combined,
  "outputs/tables/Climate_trend_class_area_forestOnly_2001_2024_combined.csv",
  row.names = FALSE
)

write.csv(
  veg_trend_area_combined,
  "outputs/tables/NDVI_kNDVI_trend_class_area_forestOnly_2001_2024_combined.csv",
  row.names = FALSE
)

write.csv(
  all_trend_area_combined,
  "outputs/tables/ALL_trend_class_area_forestOnly_2001_2024_climate_veg.csv",
  row.names = FALSE
)

# -------------------------------------------------------------
# 7. Inspect in console (optional)
# -------------------------------------------------------------
pr_area_tbl
wb_area_tbl
tmean_area_tbl
ndvi_area_tbl
kndvi_area_tbl
###############################################################













###############################################################
# 07_forestType_climate_trend_area.R
# Area (km²) and % of forest type in each WB / PDSI trend class
###############################################################

library(terra)

# -------------------------------------------------------------
# 1. Load forest-type raster and climate trend-class rasters
# -------------------------------------------------------------
# Forest type: LC_Type5 = 1–6, NA elsewhere
fType <- rast("data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif")

# Water balance and PDSI trend classes, forest-only
wb_class_forest   <- rast("data/forest_only/WB_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
pdsi_class_forest <- rast("data/forest_only/PDSI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")

# -------------------------------------------------------------
# 2. Helper: LC_Type5 × trend-class area table for one index
# -------------------------------------------------------------
make_ftype_trend_table <- function(class_raster, index_name, fType) {
  
  # Encode unique combination: 10*forest_type + trend_code
  combo <- fType * 10 + class_raster
  
  # Area per combination (km²)
  ex <- expanse(combo, unit = "km", byValue = TRUE)
  ex <- as.data.frame(ex)
  ex <- ex[!is.na(ex$value), ]  # drop NA combos if any
  
  # Decode forest type and trend code
  ex$forest_type <- floor(ex$value / 10)
  ex$trend_code  <- ex$value %% 10
  
  # Map trend codes to short and long labels
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
  
  # Total forest area across all forest types (for this index)
  total_all <- sum(ex$area, na.rm = TRUE)
  
  # Percentage of ALL forest area
  ex$perc_all_forest <- 100 * ex$area / total_all
  
  # Percentage WITHIN each forest type (rows for a given type sum to 100)
  total_by_type <- tapply(ex$area, ex$forest_type, sum, na.rm = TRUE)
  ex$perc_within_type <- 100 * ex$area /
    total_by_type[as.character(ex$forest_type)]
  
  # Tidy and round
  ex$area_km2         <- round(ex$area, 1)
  ex$perc_all_forest  <- round(ex$perc_all_forest, 2)
  ex$perc_within_type <- round(ex$perc_within_type, 2)
  
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
# 3. Build tables for WB and PDSI
# -------------------------------------------------------------
wb_ftype_tbl   <- make_ftype_trend_table(wb_class_forest,   "WB",   fType)
pdsi_ftype_tbl <- make_ftype_trend_table(pdsi_class_forest, "PDSI", fType)

clim_ftype_tbl_combined <- rbind(wb_ftype_tbl, pdsi_ftype_tbl)

# -------------------------------------------------------------
# 4. Save to disk
# -------------------------------------------------------------
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

write.csv(
  wb_ftype_tbl,
  "outputs/tables/WB_trend_class_by_forestType_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  pdsi_ftype_tbl,
  "outputs/tables/PDSI_trend_class_by_forestType_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  clim_ftype_tbl_combined,
  "outputs/tables/WB_PDSI_trend_class_by_forestType_2001_2024_combined.csv",
  row.names = FALSE
)

# Optional: inspect in console
wb_ftype_tbl
pdsi_ftype_tbl
clim_ftype_tbl_combined
###############################################################






































################## STAGE -01



###############################################################
# 08_climate_trend_class_by_country.R
# Area (km²) and % of forest in climate trend classes by country
###############################################################

library(terra)
library(dplyr)
library(readr)

# -------------------------------------------------------------
# 1. Load climate trend-class rasters (forest only)
# -------------------------------------------------------------
wb_class_forest    <- rast("data/forest_only/WB_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
tmean_class_forest <- rast("data/forest_only/Tmean_trend_class_7cat_2001_2024_NEA_forestOnly.tif")
pdsi_class_forest  <- rast("data/forest_only/PDSI_trend_class_7cat_2001_2024_NEA_forestOnly.tif")

stopifnot(
  compareGeom(wb_class_forest, tmean_class_forest, stopOnError = FALSE),
  compareGeom(wb_class_forest, pdsi_class_forest,  stopOnError = FALSE)
)

ref <- wb_class_forest

# -------------------------------------------------------------
# 2. Countries (China + Taiwan merged)
# -------------------------------------------------------------
countries <- vect("shp/rough.shp")
countries <- project(countries, crs(ref))

countries$COUNTRY_MRG <- countries$NAME
countries$COUNTRY_MRG[
  grepl("Taiwan", countries$NAME, ignore.case = TRUE) |
    countries$NAME %in% c("China", "Mainland China")
] <- "China"

cnames <- sort(unique(countries$COUNTRY_MRG))
countries$country_id <- match(countries$COUNTRY_MRG, cnames)

country_id_rast <- rasterize(countries, ref, field = "country_id")
names(country_id_rast) <- "country_id"

# quick sanity check
tmp <- as.data.frame(country_id_rast, na.rm = TRUE)
print(table(cnames[tmp$country_id]))

# -------------------------------------------------------------
# 3. Trend labels
# -------------------------------------------------------------
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

# -------------------------------------------------------------
# 4. Helper: enforce plain numerics for area and % columns
# -------------------------------------------------------------
make_country_trend_table <- function(class_raster, index_name,
                                     country_rast, cnames,
                                     trend_short, trend_long) {
  combo <- country_rast * 10 + class_raster
  
  ex <- expanse(combo, unit = "km", byValue = TRUE)
  ex <- as.data.frame(ex)
  ex <- ex[!is.na(ex$value), ]
  
  ex$country_id <- floor(ex$value / 10)
  ex$trend_code <- ex$value %% 10
  ex$country    <- cnames[ex$country_id]
  
  ex$trend_class <- trend_short[as.character(ex$trend_code)]
  ex$trend_desc  <- trend_long[ex$trend_class]
  ex$index       <- index_name
  
  # total forest area per country (numeric vector)
  total_by_country <- tapply(as.numeric(ex$area), ex$country, sum, na.rm = TRUE)
  
  # ---- key lines: force to plain numeric before rounding ----
  ex$area_km2 <- as.numeric(ex$area)
  ex$area_km2 <- round(ex$area_km2, 1)
  
  ex$perc_within_country <- 100 * as.numeric(ex$area) /
    total_by_country[ex$country]
  ex$perc_within_country <- as.numeric(ex$perc_within_country)
  ex$perc_within_country <- round(ex$perc_within_country, 2)
  # -----------------------------------------------------------
  
  out <- ex[, c("index",
                "country",
                "trend_code",
                "trend_class",
                "trend_desc",
                "area_km2",
                "perc_within_country")]
  
  out <- as.data.frame(out)            # ensure base data.frame
  rownames(out) <- NULL
  out[order(out$country, out$trend_code), ]
}

# -------------------------------------------------------------
# 5. Build tables
# -------------------------------------------------------------
wb_country_tbl    <- make_country_trend_table(
  wb_class_forest, "WB", country_id_rast, cnames, trend_short, trend_long
)
tmean_country_tbl <- make_country_trend_table(
  tmean_class_forest, "Tmean", country_id_rast, cnames, trend_short, trend_long
)
pdsi_country_tbl  <- make_country_trend_table(
  pdsi_class_forest, "PDSI", country_id_rast, cnames, trend_short, trend_long
)

clim_country_tbl_combined <- rbind(
  wb_country_tbl,
  tmean_country_tbl,
  pdsi_country_tbl
)

# extra safety: drop any list columns if they slipped through
wb_country_tbl$perc_within_country    <- as.numeric(wb_country_tbl$perc_within_country)
tmean_country_tbl$perc_within_country <- as.numeric(tmean_country_tbl$perc_within_country)
pdsi_country_tbl$perc_within_country  <- as.numeric(pdsi_country_tbl$perc_within_country)
clim_country_tbl_combined$perc_within_country <- as.numeric(clim_country_tbl_combined$perc_within_country)

# -------------------------------------------------------------
# 6. Save tables
# -------------------------------------------------------------
dir.create("outputs/tables2", showWarnings = FALSE, recursive = TRUE)

write_csv(
  wb_country_tbl,
  "outputs/tables2/WB_trend_class_by_country_2001_2024_forestOnly.csv"
)

write_csv(
  tmean_country_tbl,
  "outputs/tables2/Tmean_trend_class_by_country_2001_2024_forestOnly.csv"
)

write_csv(
  pdsi_country_tbl,
  "outputs/tables2/PDSI_trend_class_by_country_2001_2024_forestOnly.csv"
)

write_csv(
  clim_country_tbl_combined,
  "outputs/tables2/Climate_trend_class_by_country_WB_Tmean_PDSI_2001_2024_forestOnly_combined.csv"
)

head(clim_country_tbl_combined, 20)
###############################################################




# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R
# 04b_climate_means_early_late.R


library(terra)

pr   <- rast("data/processed/TC_pr_annual_2001_2024_NEA.tif")
pet  <- rast("data/processed/TC_pet_annual_2001_2024_NEA.tif")
pdsi <- rast("data/processed/TC_pdsi_annual_2001_2024_NEA.tif")
tmn  <- rast("data/processed/TC_tmn_annual_2001_2024_NEA.tif")
tmx  <- rast("data/processed/TC_tmx_annual_2001_2024_NEA.tif")

years <- 2001:2024

# indices of periods
idx_early <- which(years %in% 2001:2010)
idx_mid   <- which(years %in% 2011:2020)
idx_late  <- which(years %in% 2021:2024)

wb    <- pr - pet
tmean <- (tmn + tmx) / 2

# early/mid/late means
pr_early    <- mean(pr[[idx_early]])
pr_mid      <- mean(pr[[idx_mid]])
pr_late     <- mean(pr[[idx_late]])

wb_early    <- mean(wb[[idx_early]])
wb_mid      <- mean(wb[[idx_mid]])
wb_late     <- mean(wb[[idx_late]])

tmean_early <- mean(tmean[[idx_early]])
tmean_mid   <- mean(tmean[[idx_mid]])
tmean_late  <- mean(tmean[[idx_late]])

pdsi_early  <- mean(pdsi[[idx_early]])
pdsi_mid    <- mean(pdsi[[idx_mid]])
pdsi_late   <- mean(pdsi[[idx_late]])

# save as multi-layer rasters if you like
dir.create("data/metrics_climate", showWarnings = FALSE)

wb_periods <- c(wb_early, wb_mid, wb_late)
names(wb_periods) <- c("wb_mean_2001_2010", "wb_mean_2011_2020", "wb_mean_2021_2024")
writeRaster(wb_periods, "data/metrics_climate/WB_early_mid_late_mean_NEA.tif", overwrite = TRUE)

tmean_periods <- c(tmean_early, tmean_mid, tmean_late)
names(tmean_periods) <- c("tmean_mean_2001_2010", "tmean_mean_2011_2020", "tmean_mean_2021_2024")
writeRaster(tmean_periods, "data/metrics_climate/Tmean_early_mid_late_mean_NEA.tif", overwrite = TRUE)

# similarly for pr, pdsi if you want those maps




























###############################################################
# 04c_WB_period_differences_by_country_forestType.R
# - Uses WB early/mid/late means
# - Computes ΔWB_mid and ΔWB_late (mid–early, late–early)
# - Summarises by country × forest type
###############################################################

library(terra)
library(dplyr)
library(readr)

# -------------------------------------------------------------
# 1. Load WB period means and compute differences
# -------------------------------------------------------------
wb_periods <- rast("data/metrics_climate/WB_early_mid_late_mean_NEA.tif")
names(wb_periods)
# "wb_mean_2001_2010" "wb_mean_2011_2020" "wb_mean_2021_2024"

wb_early <- wb_periods[[1]]
wb_mid   <- wb_periods[[2]]
wb_late  <- wb_periods[[3]]

# Differences (mid–early, late–early)
wb_dmid  <- wb_mid  - wb_early
wb_dlate <- wb_late - wb_early

names(wb_dmid)  <- "wb_diff_2011_2020_minus_2001_2010"
names(wb_dlate) <- "wb_diff_2021_2024_minus_2001_2010"

dir.create("data/metrics_climate", showWarnings = FALSE)

writeRaster(
  c(wb_dmid, wb_dlate),
  "data/metrics_climate/WB_diff_early_mid_late_NEA.tif",
  overwrite = TRUE
)

# -------------------------------------------------------------
# 2. Forest types (1–6) and country IDs on same grid
# -------------------------------------------------------------
# Forest type raster
fType <- rast("data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif")
levels(fType) <- NULL
names(fType) <- "forest_type"

# Country polygons
countries <- vect("shp/rough.shp")
countries <- project(countries, crs(wb_early))

# Build merged country name from NAME; merge Taiwan into China
countries$COUNTRY_MRG <- countries$NAME
countries$COUNTRY_MRG[
  grepl("Taiwan", countries$NAME, ignore.case = TRUE) |
    countries$NAME %in% c("China", "Mainland China")
] <- "China"

# (optional) keep only NEA countries you care about
nea_names <- c("China", "Japan", "Mongolia", "North Korea",
               "Russia", "South Korea")
countries <- countries[countries$COUNTRY_MRG %in% nea_names, ]

cnames <- sort(unique(countries$COUNTRY_MRG))
countries$country_id <- match(countries$COUNTRY_MRG, cnames)

# Rasterise country IDs to WB grid
ref <- wb_early
country_id_rast <- rasterize(countries, ref, field = "country_id")
names(country_id_rast) <- "country_id"

# Save for later reuse if you like
writeRaster(
  country_id_rast,
  "data/metrics_climate/country_id_0p1deg_NEA.tif",
  overwrite = TRUE
)

# quick sanity check: number of pixels per country
tmp <- as.data.frame(country_id_rast, na.rm = TRUE)
print(table(cnames[tmp$country_id]))

# -------------------------------------------------------------
# 3. Stack diffs + forest type + country, mask to forest
# -------------------------------------------------------------
stack_diff <- mask(
  c(wb_dmid, wb_dlate, fType, country_id_rast),
  fType
)

df <- as.data.frame(stack_diff, na.rm = TRUE)
names(df)
# "wb_diff_2011_2020_minus_2001_2010",
# "wb_diff_2021_2024_minus_2001_2010",
# "forest_type", "country_id"

# Country names
df$country <- cnames[df$country_id]

# Forest-type abbreviations
ftype_labels <- c(
  "1" = "ENT",
  "2" = "EBT",
  "3" = "DNT",
  "4" = "DBT",
  "5" = "SHB",
  "6" = "GRS"
)
df$forest_type_abbr <- ftype_labels[as.character(df$forest_type)]

df <- df[!is.na(df$forest_type_abbr), ]

# -------------------------------------------------------------
# 4. Summarise ΔWB by country × forest type
# -------------------------------------------------------------
summary_wb_diff <- df |>
  group_by(country, forest_type_abbr) |>
  summarise(
    n_cells        = n(),
    wb_dmid_med    = median(wb_diff_2011_2020_minus_2001_2010, na.rm = TRUE),
    wb_dmid_q25    = quantile(wb_diff_2011_2020_minus_2001_2010, 0.25, na.rm = TRUE),
    wb_dmid_q75    = quantile(wb_diff_2011_2020_minus_2001_2010, 0.75, na.rm = TRUE),
    wb_dlate_med   = median(wb_diff_2021_2024_minus_2001_2010, na.rm = TRUE),
    wb_dlate_q25   = quantile(wb_diff_2021_2024_minus_2001_2010, 0.25, na.rm = TRUE),
    wb_dlate_q75   = quantile(wb_diff_2021_2024_minus_2001_2010, 0.75, na.rm = TRUE),
    .groups        = "drop"
  )

dir.create("outputs/tables2", showWarnings = FALSE, recursive = TRUE)

write_csv(
  summary_wb_diff,
  "outputs/tables2/WB_period_differences_by_country_forestType_2001_2024.csv"
)

print(head(summary_wb_diff, 20))
###############################################################








library(dplyr)
library(tidyr)
library(ggplot2)

# starting from your summary_wb_diff object
# columns: country, forest_type_abbr, n_cells,
#          wb_dmid_med, wb_dmid_q25, wb_dmid_q75,
#          wb_dlate_med, wb_dlate_q25, wb_dlate_q75

# 1. Long format: one row per country × forest type × period
wb_diff_long <- bind_rows(
  summary_wb_diff %>%
    transmute(
      country,
      forest_type_abbr,
      period = "2011–2020 − 2001–2010",
      med    = wb_dmid_med,
      q25    = wb_dmid_q25,
      q75    = wb_dmid_q75
    ),
  summary_wb_diff %>%
    transmute(
      country,
      forest_type_abbr,
      period = "2021–2024 − 2001–2010",
      med    = wb_dlate_med,
      q25    = wb_dlate_q25,
      q75    = wb_dlate_q75
    )
)

# 2. Factor ordering (adjust if you prefer a different sequence)
wb_diff_long <- wb_diff_long %>%
  mutate(
    country = factor(
      country,
      levels = c("Russia", "Mongolia", "China", "North Korea", "South Korea", "Japan")
    ),
    forest_type_abbr = factor(
      forest_type_abbr,
      levels = c("ENT", "EBT", "DNT", "DBT", "SHB", "GRS")
    ),
    period = factor(
      period,
      levels = c("2011–2020 − 2001–2010", "2021–2024 − 2001–2010")
    )
  )

# 3. Plot: median ± IQR for each forest type within each country
position_dodge_width <- 0.6

p_wb_diff <- ggplot(wb_diff_long,
                    aes(x = country,
                        y = med,
                        colour = forest_type_abbr,
                        group = forest_type_abbr)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_linerange(
    aes(ymin = q25, ymax = q75),
    position = position_dodge(width = position_dodge_width),
    linewidth = 0.7
  ) +
  geom_point(
    position = position_dodge(width = position_dodge_width),
    size = 2.3
  ) +
  facet_wrap(~ period, ncol = 1, strip.position = "top") +
  scale_colour_brewer(
    palette = "Dark2",
    name = "Forest type"
  ) +
  labs(
    x = "",
    y = "Change in annual water balance (mm)",
    title = ""
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title.y    = element_text(face = "bold", size = 11),
    axis.text.x     = element_text(angle = 30, hjust = 1, face = "bold", colour = "black"),
    axis.text.y     = element_text(face = "bold", colour = "black"),
    strip.text      = element_text(face = "bold", size = 11),
    legend.title    = element_text(face = "bold"),
    legend.text     = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
#### captionChanges in water balance between early, mid, and late periods
# show in R
p_wb_diff

# save to file
ggsave(
  filename = "outputs/figures/Fig_WB_period_differences_country_forestType.png",
  plot     = p_wb_diff,
  width    = 8,
  height   = 6,
  dpi      = 300
)










###############################################################
# 04d_Tmean_period_differences_by_country_forestType.R
###############################################################

library(terra)
library(dplyr)
library(readr)
library(ggplot2)

# -------------------------------------------------------------
# 1. Load Tmean period means and compute differences
# -------------------------------------------------------------
tmean_periods <- rast("data/metrics_climate/Tmean_early_mid_late_mean_NEA.tif")
names(tmean_periods)
# "tmean_mean_2001_2010" "tmean_mean_2011_2020" "tmean_mean_2021_2024"

tmean_early <- tmean_periods[[1]]
tmean_mid   <- tmean_periods[[2]]
tmean_late  <- tmean_periods[[3]]

tmean_dmid  <- tmean_mid  - tmean_early
tmean_dlate <- tmean_late - tmean_early

names(tmean_dmid)  <- "tmean_diff_2011_2020_minus_2001_2010"
names(tmean_dlate) <- "tmean_diff_2021_2024_minus_2001_2010"

dir.create("data/metrics_climate", showWarnings = FALSE)

writeRaster(
  c(tmean_dmid, tmean_dlate),
  "data/metrics_climate/Tmean_diff_early_mid_late_NEA.tif",
  overwrite = TRUE
)

# -------------------------------------------------------------
# 2. Forest types and countries on same grid
# -------------------------------------------------------------
# Forest type 1–6
fType <- rast("data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif")
levels(fType) <- NULL
names(fType) <- "forest_type"

# Country polygons
countries <- vect("shp/rough.shp")
countries <- project(countries, crs(tmean_early))

# Merge Taiwan into China
countries$COUNTRY_MRG <- countries$NAME
countries$COUNTRY_MRG[
  grepl("Taiwan", countries$NAME, ignore.case = TRUE) |
    countries$NAME %in% c("China", "Mainland China")
] <- "China"

# Keep NEA countries of interest
nea_names <- c("China", "Japan", "Mongolia", "North Korea",
               "Russia", "South Korea")
countries <- countries[countries$COUNTRY_MRG %in% nea_names, ]

cnames <- sort(unique(countries$COUNTRY_MRG))
countries$country_id <- match(countries$COUNTRY_MRG, cnames)

# Rasterise to Tmean grid
ref <- tmean_early
country_id_rast <- rasterize(countries, ref, field = "country_id")
names(country_id_rast) <- "country_id"

# Optional: save for reuse
dir.create("data/metrics_climate", showWarnings = FALSE)
writeRaster(
  country_id_rast,
  "data/metrics_climate/country_id_0p1deg_NEA.tif",
  overwrite = TRUE
)

# -------------------------------------------------------------
# 3. Stack and mask to forest; convert to data frame
# -------------------------------------------------------------
stack_diff_t <- mask(
  c(tmean_dmid, tmean_dlate, fType, country_id_rast),
  fType
)

df_t <- as.data.frame(stack_diff_t, na.rm = TRUE)
names(df_t)

df_t$country <- cnames[df_t$country_id]

ftype_labels <- c(
  "1" = "ENT",
  "2" = "EBT",
  "3" = "DNT",
  "4" = "DBT",
  "5" = "SHB",
  "6" = "GRS"
)
df_t$forest_type_abbr <- ftype_labels[as.character(df_t$forest_type)]
df_t <- df_t[!is.na(df_t$forest_type_abbr), ]

# -------------------------------------------------------------
# 4. Summarise ΔTmean by country × forest type
# -------------------------------------------------------------
summary_tmean_diff <- df_t |>
  group_by(country, forest_type_abbr) |>
  summarise(
    n_cells         = n(),
    tmean_dmid_med  = median(tmean_diff_2011_2020_minus_2001_2010, na.rm = TRUE),
    tmean_dmid_q25  = quantile(tmean_diff_2011_2020_minus_2001_2010, 0.25, na.rm = TRUE),
    tmean_dmid_q75  = quantile(tmean_diff_2011_2020_minus_2001_2010, 0.75, na.rm = TRUE),
    tmean_dlate_med = median(tmean_diff_2021_2024_minus_2001_2010, na.rm = TRUE),
    tmean_dlate_q25 = quantile(tmean_diff_2021_2024_minus_2001_2010, 0.25, na.rm = TRUE),
    tmean_dlate_q75 = quantile(tmean_diff_2021_2024_minus_2001_2010, 0.75, na.rm = TRUE),
    .groups         = "drop"
  )

dir.create("outputs/tables2", showWarnings = FALSE, recursive = TRUE)

write_csv(
  summary_tmean_diff,
  "outputs/tables2/Tmean_period_differences_by_country_forestType_2001_2024.csv"
)

print(head(summary_tmean_diff, 20))

# -------------------------------------------------------------
# 5. Publication-ready plot (mirror of WB figure)
# -------------------------------------------------------------
tmean_diff_long <- bind_rows(
  summary_tmean_diff %>%
    transmute(
      country,
      forest_type_abbr,
      period = "2011–2020 − 2001–2010",
      med    = tmean_dmid_med,
      q25    = tmean_dmid_q25,
      q75    = tmean_dmid_q75
    ),
  summary_tmean_diff %>%
    transmute(
      country,
      forest_type_abbr,
      period = "2021–2024 − 2001–2010",
      med    = tmean_dlate_med,
      q25    = tmean_dlate_q25,
      q75    = tmean_dlate_q75
    )
)

tmean_diff_long <- tmean_diff_long %>%
  mutate(
    country = factor(
      country,
      levels = c("Russia", "Mongolia", "China", "North Korea", "South Korea", "Japan")
    ),
    forest_type_abbr = factor(
      forest_type_abbr,
      levels = c("ENT", "EBT", "DNT", "DBT", "SHB", "GRS")
    ),
    period = factor(
      period,
      levels = c("2011–2020 − 2001–2010", "2021–2024 − 2001–2010")
    )
  )

pos_dodge <- 0.6

p_tmean_diff <- ggplot(tmean_diff_long,
                       aes(x = country,
                           y = med,
                           colour = forest_type_abbr,
                           group = forest_type_abbr)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_linerange(
    aes(ymin = q25, ymax = q75),
    position = position_dodge(width = pos_dodge),
    linewidth = 0.7
  ) +
  geom_point(
    position = position_dodge(width = pos_dodge),
    size = 2.3
  ) +
  facet_wrap(~ period, ncol = 1, strip.position = "top") +
  scale_colour_brewer(
    palette = "Dark2",
    name = "Forest type"
  ) +
  labs(
    x = "",
    y = "Change in mean annual temperature (°C)",
    title = ""
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    strip.text   = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_tmean_diff

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = "outputs/figures/Fig_Tmean_period_differences_country_forestType.png",
  plot     = p_tmean_diff,
  width    = 8,
  height   = 6,
  dpi      = 300
)
###############################################################








###############################################################
# 04e_PDSI_trendClass_by_country_plot.R
###############################################################

library(dplyr)
library(readr)
library(ggplot2)

# 1. Read combined climate trend-class table, keep only PDSI
clim_country_tbl_combined <- read_csv(
  "outputs/tables2/Climate_trend_class_by_country_WB_Tmean_PDSI_2001_2024_forestOnly_combined.csv"
)

pds_country <- clim_country_tbl_combined %>%
  filter(index == "PDSI")

# 2. Factor ordering for countries and classes
pds_country <- pds_country %>%
  mutate(
    country = factor(
      country,
      levels = c("Russia", "Mongolia", "China", "North Korea", "South Korea", "Japan")
    ),
    trend_class = factor(
      trend_class,
      levels = c("Dec_p<=0.01", "Dec_p<=0.05", "Dec_p<=0.10",
                 "NS",
                 "Inc_p<=0.10", "Inc_p<=0.05", "Inc_p<=0.01")
    )
  )

# 3. Colour palette consistent with your 7-class scheme
trend_cols <- c(
  "Dec_p<=0.01" = "#b2182b",  # strong negative
  "Dec_p<=0.05" = "#ef8a62",
  "Dec_p<=0.10" = "#fddbc7",
  "NS"          = "#d9d9d9",
  "Inc_p<=0.10" = "#c7eae5",
  "Inc_p<=0.05" = "#5ab4ac",
  "Inc_p<=0.01" = "#01665e"   # strong positive
)

# 4. Stacked bar plot: % of forest in each PDSI trend class, by country
p_pds_bar <- ggplot(pds_country,
                    aes(x = country,
                        y = perc_within_country,
                        fill = trend_class)) +
  geom_col(colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = trend_cols,
                    name = "PDSI trend class") +
  labs(
    x = "",
    y = "Forest area (%)",
    title = "PDSI trend classes in forested areas (2001–2024)"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold", colour = "black"),
    axis.text.y  = element_text(face = "bold", colour = "black"),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_pds_bar

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = "outputs/figures/Fig_PDSI_trend_classes_country_forestOnly.png",
  plot     = p_pds_bar,
  width    = 7,
  height   = 4.5,
  dpi      = 300
)
###############################################################



















