


























# -------------------------------------------------------------
# 1. Load annual forest-only kNDVI / NDVI and climate stacks
# -------------------------------------------------------------
kndvi_forest <- rast("data/forest_only/kNDVI_annual_2001_2024_NEA_forestOnly.tif")
ndvi_forest  <- rast("data/forest_only/NDVI_annual_2001_2024_NEA_forestOnly.tif") # optional, not used yet

pr   <- rast("data/processed/TC_pr_annual_2001_2024_NEA.tif")
pet  <- rast("data/processed/TC_pet_annual_2001_2024_NEA.tif")
pdsi <- rast("data/processed/TC_pdsi_annual_2001_2024_NEA.tif")

fType <- rast("data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif")

years <- 2001:2024
stopifnot(
  nlyr(kndvi_forest) == length(years),
  nlyr(pr)           == length(years),
  nlyr(pet)          == length(years),
  nlyr(pdsi)         == length(years)
)

# -------------------------------------------------------------
# 2. Derive water balance and mask climate stacks to forests
# -------------------------------------------------------------
wb <- pr - pet
names(wb)   <- paste0("wb_", years)
names(pdsi) <- paste0("pdsi_", years)

wb_forest   <- mask(wb,   fType)
pdsi_forest <- mask(pdsi, fType)

# optional sanity checks
ref <- kndvi_forest[[1]]
stopifnot(all(res(wb_forest)   == res(ref)))
stopifnot(all(ext(wb_forest)   == ext(ref)))
stopifnot(all(crs(wb_forest)   == crs(ref)))
stopifnot(all(res(pdsi_forest) == res(ref)))
stopifnot(all(ext(pdsi_forest) == ext(ref)))
stopifnot(all(crs(pdsi_forest) == crs(ref)))

# -------------------------------------------------------------
# 3. Demean each pixel's time series (anomalies)
# -------------------------------------------------------------
demean_fun <- function(v) {
  if (all(is.na(v))) return(v)
  m <- mean(v, na.rm = TRUE)
  v - m
}

kndvi_anom <- app(kndvi_forest, demean_fun)
wb_anom    <- app(wb_forest,    demean_fun)
pdsi_anom  <- app(pdsi_forest,  demean_fun)

n_years <- length(years)

# -------------------------------------------------------------
# 4. Correlation function: returns r and p-value
# -------------------------------------------------------------
cor_fun <- function(v, n_years) {
  kn   <- v[1:n_years]
  clim <- v[(n_years + 1):(2 * n_years)]
  
  good <- !is.na(kn) & !is.na(clim)
  n <- sum(good)
  if (n < 8) {
    return(c(r = NA, p = NA))
  }
  
  r <- cor(kn[good], clim[good])
  
  # two-sided p-value from t-statistic
  tval <- r * sqrt((n - 2) / (1 - r^2))
  p <- 2 * pt(-abs(tval), df = n - 2)
  
  c(r = r, p = p)
}

terraOptions(progress = 1, memfrac = 0.7)
dir.create("data/metrics_coupling", showWarnings = FALSE)

# -------------------------------------------------------------
# 5. kNDVI–WB correlation (forest only)
# -------------------------------------------------------------
stack_kndvi_wb <- c(kndvi_anom, wb_anom)

kndvi_wb_cor <- app(
  stack_kndvi_wb,
  fun      = cor_fun,
  n_years  = n_years,
  cores    = 4,
  filename = "data/metrics_coupling/kNDVI_WB_cor_p_2001_2024_NEA_forestOnly.tif",
  overwrite = TRUE
)
names(kndvi_wb_cor) <- c("r_kNDVI_WB", "p_kNDVI_WB")

# -------------------------------------------------------------
# 6. kNDVI–PDSI correlation (forest only)
# -------------------------------------------------------------
stack_kndvi_pdsi <- c(kndvi_anom, pdsi_anom)

kndvi_pdsi_cor <- app(
  stack_kndvi_pdsi,
  fun      = cor_fun,
  n_years  = n_years,
  cores    = 4,
  filename = "data/metrics_coupling/kNDVI_PDSI_cor_p_2001_2024_NEA_forestOnly.tif",
  overwrite = TRUE
)
names(kndvi_pdsi_cor) <- c("r_kNDVI_PDSI", "p_kNDVI_PDSI")

# -------------------------------------------------------------
# 7. Optional: quick look at correlation distributions
# -------------------------------------------------------------
hist(values(kndvi_wb_cor[["r_kNDVI_WB"]]),  breaks = 40, main = "r(kNDVI, WB)")
hist(values(kndvi_pdsi_cor[["r_kNDVI_PDSI"]]), breaks = 40, main = "r(kNDVI, PDSI)")
###############################################################





















###############################################################
# 09_coupling_class_area.R
# Classify kNDVI–WB and kNDVI–PDSI correlations and
# compute area statistics (overall + by forest type)
###############################################################

library(terra)

# -------------------------------------------------------------
# 1. Load correlation rasters and forest-type raster
# -------------------------------------------------------------
kndvi_wb_cor   <- rast("data/metrics_coupling/kNDVI_WB_cor_p_2001_2024_NEA_forestOnly.tif")
kndvi_pdsi_cor <- rast("data/metrics_coupling/kNDVI_PDSI_cor_p_2001_2024_NEA_forestOnly.tif")

names(kndvi_wb_cor)   <- c("r_kNDVI_WB",   "p_kNDVI_WB")
names(kndvi_pdsi_cor) <- c("r_kNDVI_PDSI", "p_kNDVI_PDSI")

fType <- rast("data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif")

# -------------------------------------------------------------
# 2. Classify correlation into 7 categories (sign + p-value)
#    0 = non-significant (p > 0.10)
#    1 = weak positive     (0.05 < p <= 0.10, r > 0)
#    2 = weak negative     (0.05 < p <= 0.10, r < 0)
#    3 = moderate positive (0.01 < p <= 0.05, r > 0)
#    4 = moderate negative (0.01 < p <= 0.05, r < 0)
#    5 = strong positive   (p <= 0.01, r > 0)
#    6 = strong negative   (p <= 0.01, r < 0)
# -------------------------------------------------------------
classify_corr <- function(r_stack, prefix, out_path) {
  r <- r_stack[[1]]
  p <- r_stack[[2]]
  
  cls <- r
  cls[] <- NA
  
  # 0: non-significant (p > 0.10)
  cls <- ifel(!is.na(p) & p > 0.10, 0, cls)
  
  # 1–2: weak
  cls <- ifel(p <= 0.10 & p > 0.05 & r > 0, 1, cls)
  cls <- ifel(p <= 0.10 & p > 0.05 & r < 0, 2, cls)
  
  # 3–4: moderate
  cls <- ifel(p <= 0.05 & p > 0.01 & r > 0, 3, cls)
  cls <- ifel(p <= 0.05 & p > 0.01 & r < 0, 4, cls)
  
  # 5–6: strong
  cls <- ifel(p <= 0.01 & r > 0, 5, cls)
  cls <- ifel(p <= 0.01 & r < 0, 6, cls)
  
  names(cls) <- paste0(prefix, "_corr_class")
  
  writeRaster(cls, out_path, overwrite = TRUE)
  cls
}

dir.create("data/metrics_coupling", showWarnings = FALSE)

kndvi_wb_class <- classify_corr(
  kndvi_wb_cor,
  prefix   = "kNDVI_WB",
  out_path = "data/metrics_coupling/kNDVI_WB_corr_class_7cat_2001_2024_NEA_forestOnly.tif"
)

kndvi_pdsi_class <- classify_corr(
  kndvi_pdsi_cor,
  prefix   = "kNDVI_PDSI",
  out_path = "data/metrics_coupling/kNDVI_PDSI_corr_class_7cat_2001_2024_NEA_forestOnly.tif"
)

# -------------------------------------------------------------
# 3. Descriptions for correlation classes
# -------------------------------------------------------------
corr_desc_long <- c(
  "NS"          = "Non-significant correlation (p > 0.10)",
  "Inc_p<=0.10" = "Positive correlation, 0.05 < p \u2264 0.10",
  "Dec_p<=0.10" = "Negative correlation, 0.05 < p \u2264 0.10",
  "Inc_p<=0.05" = "Positive correlation, 0.01 < p \u2264 0.05",
  "Dec_p<=0.05" = "Negative correlation, 0.01 < p \u2264 0.05",
  "Inc_p<=0.01" = "Positive correlation, p \u2264 0.01",
  "Dec_p<=0.01" = "Negative correlation, p \u2264 0.01"
)

# mapping from numeric code to short label (same as trend tables)
code_to_short <- c(
  "0" = "NS",
  "1" = "Inc_p<=0.10",
  "2" = "Dec_p<=0.10",
  "3" = "Inc_p<=0.05",
  "4" = "Dec_p<=0.05",
  "5" = "Inc_p<=0.01",
  "6" = "Dec_p<=0.01"
)

# -------------------------------------------------------------
# 4. Helper: overall area table for one correlation-class raster
# -------------------------------------------------------------
make_overall_area_table <- function(r, index_name, code_to_short, desc_long) {
  ex <- expanse(r, unit = "km", byValue = TRUE)
  ex <- as.data.frame(ex)
  ex$value <- as.character(ex$value)
  
  ex$class_code   <- ex$value
  ex$class_short  <- code_to_short[ex$class_code]
  ex$description  <- desc_long[ex$class_short]
  
  total_area <- sum(ex$area, na.rm = TRUE)
  ex$percent <- 100 * ex$area / total_area
  
  out <- ex[, c("class_code", "class_short", "description", "area", "percent")]
  names(out) <- c("class_code", "class_label", "class_description",
                  "area_km2", "percent_forest")
  out$index           <- index_name
  out$area_km2        <- round(out$area_km2, 1)
  out$percent_forest  <- round(out$percent_forest, 2)
  out <- out[, c("index", "class_code", "class_label",
                 "class_description", "area_km2", "percent_forest")]
  out[order(out$class_code), ]
}

wb_corr_area    <- make_overall_area_table(kndvi_wb_class,   "kNDVI_WB",   code_to_short, corr_desc_long)
pdsi_corr_area  <- make_overall_area_table(kndvi_pdsi_class, "kNDVI_PDSI", code_to_short, corr_desc_long)

overall_corr_area <- rbind(wb_corr_area, pdsi_corr_area)

# -------------------------------------------------------------
# 5. Helper: forest-type × correlation-class table
# -------------------------------------------------------------
make_ftype_corr_table <- function(class_raster, index_name,
                                  fType, code_to_short, desc_long) {
  
  combo <- fType * 10 + class_raster  # encode type + class
  ex <- expanse(combo, unit = "km", byValue = TRUE)
  ex <- as.data.frame(ex)
  ex <- ex[!is.na(ex$value), ]
  
  ex$forest_type <- floor(ex$value / 10)
  ex$class_code  <- ex$value %% 10
  ex$class_code  <- as.character(ex$class_code)
  
  ex$class_label <- code_to_short[ex$class_code]
  ex$class_desc  <- desc_long[ex$class_label]
  ex$index       <- index_name
  
  total_all <- sum(ex$area, na.rm = TRUE)
  ex$perc_all_forest <- 100 * ex$area / total_all
  
  total_by_type <- tapply(ex$area, ex$forest_type, sum, na.rm = TRUE)
  ex$perc_within_type <- 100 * ex$area /
    total_by_type[as.character(ex$forest_type)]
  
  ex$area_km2         <- round(ex$area, 1)
  ex$perc_all_forest  <- round(ex$perc_all_forest, 2)
  ex$perc_within_type <- round(ex$perc_within_type, 2)
  
  out <- ex[, c("index", "forest_type", "class_code",
                "class_label", "class_desc",
                "area_km2", "perc_all_forest", "perc_within_type")]
  
  out[order(out$forest_type, as.numeric(out$class_code)), ]
}

wb_corr_ftype   <- make_ftype_corr_table(kndvi_wb_class,   "kNDVI_WB",   fType,
                                         code_to_short, corr_desc_long)
pdsi_corr_ftype <- make_ftype_corr_table(kndvi_pdsi_class, "kNDVI_PDSI", fType,
                                         code_to_short, corr_desc_long)

corr_ftype_combined <- rbind(wb_corr_ftype, pdsi_corr_ftype)

# -------------------------------------------------------------
# 6. Save tables
# -------------------------------------------------------------
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)

write.csv(
  wb_corr_area,
  "outputs/tables/kNDVI_WB_corr_class_area_forestOnly_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  pdsi_corr_area,
  "outputs/tables/kNDVI_PDSI_corr_class_area_forestOnly_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  overall_corr_area,
  "outputs/tables/kNDVI_WB_PDSI_corr_class_area_forestOnly_2001_2024_combined.csv",
  row.names = FALSE
)

write.csv(
  wb_corr_ftype,
  "outputs/tables/kNDVI_WB_corr_class_by_forestType_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  pdsi_corr_ftype,
  "outputs/tables/kNDVI_PDSI_corr_class_by_forestType_2001_2024.csv",
  row.names = FALSE
)

write.csv(
  corr_ftype_combined,
  "outputs/tables/kNDVI_WB_PDSI_corr_class_by_forestType_2001_2024_combined.csv",
  row.names = FALSE
)

# Optional: quick look in console
wb_corr_area
pdsi_corr_area
wb_corr_ftype[wb_corr_ftype$forest_type == 2, ]   # example: forest type 2
###############################################################





























