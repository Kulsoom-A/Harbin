###############################################################
# 08_kNDVI_functioning_resilience.R
# kNDVI mean, trend, stability (CV) + resistance & recovery
# by country × forest type, NE Asia
###############################################################

library(terra)
library(dplyr)
library(readr)

terraOptions(progress = 1, memfrac = 0.7)

# -------------------------------------------------------------
# 0. Base rasters and forest mask
# -------------------------------------------------------------
# Annual kNDVI, 2001–2024 (already aligned)
kndvi <- rast("data/processed/kNDVI_annual_2001_2024_NEA.tif")
years_all <- 2001:2024
stopifnot(nlyr(kndvi) == length(years_all))

# Forest mask: LC_Type5 = 1–6
fType <- rast("data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif")
levels(fType) <- NULL          # ensure numeric, not factor
names(fType) <- "forest_type"

# Reference grid
ref <- kndvi[[1]]

dir.create("data/metrics_climate", showWarnings = FALSE, recursive = TRUE)
dir.create("data/metrics_veg",     showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/tables2",      showWarnings = FALSE, recursive = TRUE)

country_id_path <- "data/metrics_climate/country_id_0p1deg_NEA.tif"
cnames_path     <- "data/metrics_climate/country_names.rds"

# -------------------------------------------------------------
# 1. Country raster and names (from iso_a2, Taiwan → China)
# -------------------------------------------------------------
# Always rebuild here so we are not stuck with the broken NAME field
NAEcntry <- vect("shp/rough.shp")
NAEcntry <- project(NAEcntry, crs(ref))

# Inspect useful fields (one-off check if needed)
# print(names(NAEcntry))
# head(as.data.frame(NAEcntry)[, c("fid", "iso_a2", "NAME")])

codes <- as.character(NAEcntry$iso_a2)
# Check what we have
print(unique(codes))

# Map ISO A2 codes to country names
country_raw <- dplyr::case_when(
  codes %in% c("CN")             ~ "China",
  codes %in% c("TW")             ~ "China",        # merge Taiwan into China
  codes %in% c("JP")             ~ "Japan",
  codes %in% c("MN")             ~ "Mongolia",
  codes %in% c("KP")             ~ "North Korea",
  codes %in% c("KR")             ~ "South Korea",
  codes %in% c("RU")             ~ "Russia",
  TRUE                           ~ NA_character_
)

if (any(is.na(country_raw))) {
  warning("Some polygons have unknown iso_a2 codes; please check:")
  print(as.data.frame(NAEcntry)[is.na(country_raw), c("fid", "iso_a2", "NAME")])
}

cnames <- sort(unique(na.omit(country_raw)))
print(cnames)
# e.g. "China" "Japan" "Mongolia" "North Korea" "Russia" "South Korea"

NAEcntry$country_id <- match(country_raw, cnames)

country_id_rast <- rasterize(NAEcntry, ref, field = "country_id")
names(country_id_rast) <- "country_id"

writeRaster(country_id_rast, country_id_path, overwrite = TRUE)
saveRDS(cnames, cnames_path)

# Quick diagnostic: should now have several IDs, not just 1
print(freq(country_id_rast))

# -------------------------------------------------------------
# 2. Forest-type abbreviations
# -------------------------------------------------------------
ftype_labels <- c(
  "1" = "ENT",  # Evergreen needleleaf
  "2" = "EBT",  # Evergreen broadleaf
  "3" = "DNT",  # Deciduous needleleaf
  "4" = "DBT",  # Deciduous broadleaf
  "5" = "SHB",  # Shrublands
  "6" = "GRS"   # Grasslands
)

# =============================================================
# 3. BASIC kNDVI FUNCTIONING METRICS (mean, Sen slope, CV)
# =============================================================

# Mask kNDVI to forest
kndvi_forest <- mask(kndvi, fType)

# --- 3.1 Mean and CV -----------------------------------------
kndvi_mean <- app(kndvi_forest, fun = mean, na.rm = TRUE)
names(kndvi_mean) <- "kNDVI_mean"

kndvi_sd <- app(kndvi_forest, fun = sd, na.rm = TRUE)
names(kndvi_sd) <- "kNDVI_sd"

kndvi_cv <- kndvi_sd / kndvi_mean
names(kndvi_cv) <- "kNDVI_CV"

# Optional: avoid extreme CV in very low mean pixels
kndvi_cv[kndvi_mean < 0.05] <- NA

# --- 3.2 Sen slope of kNDVI (trend) --------------------------
n_years <- length(years_all)
pair_idx <- combn(n_years, 2)  # index pairs for Sen slope

sen_slope_fun <- function(v, pair_idx) {
  v <- as.numeric(v)
  if (all(is.na(v)) || sum(!is.na(v)) < 8) {
    return(NA_real_)
  }
  slopes <- (v[pair_idx[2, ]] - v[pair_idx[1, ]]) /
    (pair_idx[2, ] - pair_idx[1, ])
  median(slopes, na.rm = TRUE)
}

kndvi_slope <- app(
  kndvi_forest,
  fun = sen_slope_fun,
  pair_idx = pair_idx,
  cores = 4
)
names(kndvi_slope) <- "kNDVI_senSlope"

# --- 3.3 Save kNDVI functioning rasters ----------------------
writeRaster(
  c(kndvi_mean, kndvi_slope, kndvi_cv),
  "data/metrics_veg/kNDVI_mean_slope_CV_2001_2024_NEA_forestOnly.tif",
  overwrite = TRUE
)

# --- 3.4 Summarise by country × forest type ------------------
cnames <- readRDS(cnames_path)              # ensure consistent order
country_id_rast <- rast(country_id_path)
names(country_id_rast) <- "country_id"

stack_fun <- c(kndvi_mean, kndvi_slope, kndvi_cv, fType, country_id_rast)
fun_df <- as.data.frame(stack_fun, xy = FALSE, na.rm = TRUE)

fun_df$forest_type_int  <- as.integer(fun_df$forest_type)
fun_df$forest_type_abbr <- ftype_labels[as.character(fun_df$forest_type_int)]
fun_df$country          <- cnames[fun_df$country_id]

fun_df <- fun_df[!is.na(fun_df$forest_type_abbr) & !is.na(fun_df$country), ]

summary_fun <- fun_df |>
  group_by(country, forest_type_abbr) |>
  summarise(
    n_cells       = n(),
    kmean_med     = median(kNDVI_mean,     na.rm = TRUE),
    kmean_q25     = quantile(kNDVI_mean,   0.25, na.rm = TRUE),
    kmean_q75     = quantile(kNDVI_mean,   0.75, na.rm = TRUE),
    slope_med     = median(kNDVI_senSlope, na.rm = TRUE),
    slope_q25     = quantile(kNDVI_senSlope, 0.25, na.rm = TRUE),
    slope_q75     = quantile(kNDVI_senSlope, 0.75, na.rm = TRUE),
    cv_med        = median(kNDVI_CV,       na.rm = TRUE),
    cv_q25        = quantile(kNDVI_CV,     0.25, na.rm = TRUE),
    cv_q75        = quantile(kNDVI_CV,     0.75, na.rm = TRUE),
    .groups       = "drop"
  )

write_csv(
  summary_fun,
  "outputs/tables2/kNDVI_functioning_summary_by_country_forestType_2001_2024.csv"
)

# =============================================================
# 4. HOT–DRY YEARS & RESISTANCE / RECOVERY (2001–2022)
# =============================================================

# --- 4.1 SPEI12 growing-season minima ------------------------
spei12 <- rast("data/processed/SPEI12_monthly_2001_2022_NEA.tif")

dates_spei <- seq(
  from = as.Date("2001-01-01"),
  to   = as.Date("2022-12-01"),
  by   = "1 month"
)
stopifnot(nlyr(spei12) == length(dates_spei))

years_event <- 2001:2022
grow_months <- 5:9  # May–September

spei_min_list <- vector("list", length(years_event))
names(spei_min_list) <- as.character(years_event)

for (i in seq_along(years_event)) {
  yr <- years_event[i]
  idx <- which(
    format(dates_spei, "%Y") == yr &
      as.integer(format(dates_spei, "%m")) %in% grow_months
  )
  spei_min_list[[i]] <- app(spei12[[idx]], fun = min, na.rm = TRUE)
}

spei_gs_min <- rast(spei_min_list)
names(spei_gs_min) <- paste0("SPEImin_", years_event)

# Align to kNDVI grid if needed
spei_gs_min <- resample(spei_gs_min, ref, method = "bilinear")

# --- 4.2 Tmax anomalies vs 2001–2010 baseline ----------------
tmx <- rast("data/processed/TC_tmx_annual_2001_2024_NEA.tif")
tmx <- tmx[[1:length(years_event)]]  # 2001–2022 only
stopifnot(nlyr(tmx) == length(years_event))

tmx_base_idx <- which(years_event %in% 2001:2010)
tmx_base <- mean(tmx[[tmx_base_idx]])

tmx_anom <- tmx - tmx_base
names(tmx_anom) <- paste0("tmx_anom_", years_event)

# --- 4.3 Resistance and recovery function --------------------
resilience_fun <- function(v, years, base_years,
                           drought_thr = -1.0,
                           heat_thr    = 0.75) {
  # v = c(kNDVI[years], SPEImin[years], tmx_anom[years])
  n <- length(years)
  k    <- v[1:n]
  spei <- v[(n + 1):(2 * n)]
  tmax <- v[(2 * n + 1):(3 * n)]
  
  if (all(is.na(k))) {
    return(c(R_med = NA_real_, Rec_med = NA_real_, n_events = 0))
  }
  
  base_idx <- which(years %in% base_years)
  mu  <- mean(k[base_idx], na.rm = TRUE)
  sig <- sd(k[base_idx],   na.rm = TRUE)
  
  if (is.na(mu) || is.na(sig) || sig == 0) {
    return(c(R_med = NA_real_, Rec_med = NA_real_, n_events = 0))
  }
  
  A <- (k - mu) / sig  # standardized kNDVI anomalies
  
  hotdry <- which(
    !is.na(A) &
      !is.na(spei) & !is.na(tmax) &
      spei <= drought_thr &
      tmax >= heat_thr
  )
  
  if (length(hotdry) == 0) {
    return(c(R_med = NA_real_, Rec_med = NA_real_, n_events = 0))
  }
  
  # Resistance: anomaly during hot–dry year
  R_vals <- A[hotdry]
  
  # Recovery: mean anomaly in t+1 and t+2
  Rec_vals <- rep(NA_real_, length(hotdry))
  for (i in seq_along(hotdry)) {
    t <- hotdry[i]
    nxt <- t + 1:2
    nxt <- nxt[nxt <= n]
    if (length(nxt) > 0) {
      Rec_vals[i] <- mean(A[nxt], na.rm = TRUE)
    }
  }
  
  R_med   <- median(R_vals,  na.rm = TRUE)
  Rec_med <- median(Rec_vals, na.rm = TRUE)
  n_ev    <- sum(!is.na(R_vals))
  
  c(R_med = R_med, Rec_med = Rec_med, n_events = n_ev)
}

# --- 4.4 Apply resilience_fun pixel-wise ----------------------
# Restrict kNDVI to 2001–2022 and forest
kndvi_2001_2022 <- kndvi_forest[[1:length(years_event)]]

big_stack <- c(kndvi_2001_2022, spei_gs_min, tmx_anom)

res_stack <- app(
  big_stack,
  fun = resilience_fun,
  years      = years_event,
  base_years = 2001:2010,
  drought_thr = -1.0,
  heat_thr    = 0.75,
  cores = 4
)

names(res_stack) <- c("R_hotdry", "Rec_hotdry", "n_hotdry")

# Mask to forest explicitly
res_stack_forest <- mask(res_stack, fType)

writeRaster(
  res_stack_forest,
  "data/metrics_veg/kNDVI_resistance_recovery_hotdry_2001_2022_NEA_forestOnly.tif",
  overwrite = TRUE
)

# --- 4.5 Summarise resistance & recovery by country × forest type ----
stack_res <- c(res_stack_forest, fType, country_id_rast)
res_df <- as.data.frame(stack_res, xy = FALSE, na.rm = TRUE)

res_df$forest_type_int  <- as.integer(res_df$forest_type)
res_df$forest_type_abbr <- ftype_labels[as.character(res_df$forest_type_int)]
res_df$country          <- cnames[res_df$country_id]

res_df <- res_df[
  !is.na(res_df$forest_type_abbr) &
    !is.na(res_df$country) &
    !is.na(res_df$R_hotdry),
]

summary_res <- res_df |>
  group_by(country, forest_type_abbr) |>
  summarise(
    n_cells       = n(),
    R_med         = median(R_hotdry,   na.rm = TRUE),
    R_q25         = quantile(R_hotdry, 0.25, na.rm = TRUE),
    R_q75         = quantile(R_hotdry, 0.75, na.rm = TRUE),
    Rec_med       = median(Rec_hotdry,   na.rm = TRUE),
    Rec_q25       = quantile(Rec_hotdry, 0.25, na.rm = TRUE),
    Rec_q75       = quantile(Rec_hotdry, 0.75, na.rm = TRUE),
    n_events_med  = median(n_hotdry,   na.rm = TRUE),
    n_events_mean = mean(n_hotdry,     na.rm = TRUE),
    .groups       = "drop"
  )

write_csv(
  summary_res,
  "outputs/tables2/kNDVI_resistance_recovery_hotdry_summary_by_country_forestType_2001_2022.csv"
)

###############################################################
# End of script
###############################################################













###############################################################
# 09_kNDVI_classify_functioning_resilience.R
# Classify kNDVI mean, trend, CV, resistance & recovery
###############################################################

library(terra)

terraOptions(progress = 1, memfrac = 0.7)

# -------------------------------------------------------------
# 1. Load existing metric rasters
# -------------------------------------------------------------
fun_path <- "data/metrics_veg/kNDVI_mean_slope_CV_2001_2024_NEA_forestOnly.tif"
res_path <- "data/metrics_veg/kNDVI_resistance_recovery_hotdry_2001_2022_NEA_forestOnly.tif"

fun_rast <- rast(fun_path)
res_rast <- rast(res_path)

# Make sure names are as expected (adjust if needed)
names(fun_rast) <- c("kNDVI_mean", "kNDVI_senSlope", "kNDVI_CV")
names(res_rast) <- c("R_hotdry", "Rec_hotdry", "n_hotdry")

kndvi_mean  <- fun_rast[["kNDVI_mean"]]
kndvi_slope <- fun_rast[["kNDVI_senSlope"]]
kndvi_cv    <- fun_rast[["kNDVI_CV"]]

R_hotdry    <- res_rast[["R_hotdry"]]
Rec_hotdry  <- res_rast[["Rec_hotdry"]]
n_hotdry    <- res_rast[["n_hotdry"]]

# -------------------------------------------------------------
# 2. kNDVI mean – productivity classes
# -------------------------------------------------------------
# breaks in kNDVI units
br_mean <- c(0, 0.3, 0.5, 0.7, 0.85, 1.1)

kmean_class <- app(kndvi_mean, fun = function(x) {
  as.integer(cut(
    x,
    breaks = br_mean,
    labels = FALSE,
    include.lowest = TRUE
  ))
})

names(kmean_class) <- "kmean_class"

lev_mean <- data.frame(
  value = 1:5,
  class = c(
    "very low (≤ 0.30)",
    "low (0.30–0.50)",
    "moderate (0.50–0.70)",
    "high (0.70–0.85)",
    "very high (> 0.85)"
  )
)
levels(kmean_class)[[1]] <- lev_mean

# -------------------------------------------------------------
# 3. kNDVI Sen slope – greening / browning classes
# -------------------------------------------------------------
# per-year Sen slope classes
br_slope <- c(-1, -0.004, -0.001, 0.001, 0.004, 1)

ktrend_class <- app(kndvi_slope, fun = function(x) {
  as.integer(cut(
    x,
    breaks = br_slope,
    labels = FALSE,
    include.lowest = TRUE
  ))
})

names(ktrend_class) <- "ktrend_class"

lev_trend <- data.frame(
  value = 1:5,
  class = c(
    "strong decline (≤ −0.004 yr⁻¹)",
    "weak decline (−0.004 to −0.001)",
    "stable (−0.001 to 0.001)",
    "weak greening (0.001 to 0.004)",
    "strong greening (> 0.004 yr⁻¹)"
  )
)
levels(ktrend_class)[[1]] <- lev_trend

# -------------------------------------------------------------
# 4. kNDVI CV – functional stability classes
# -------------------------------------------------------------
br_cv <- c(0, 0.10, 0.20, 0.30, 0.50, 10)

kcv_class <- app(kndvi_cv, fun = function(x) {
  as.integer(cut(
    x,
    breaks = br_cv,
    labels = FALSE,
    include.lowest = TRUE
  ))
})

names(kcv_class) <- "kcv_class"

lev_cv <- data.frame(
  value = 1:5,
  class = c(
    "very stable (CV ≤ 0.10)",
    "stable (0.10–0.20)",
    "moderate variability (0.20–0.30)",
    "high variability (0.30–0.50)",
    "very high variability (> 0.50)"
  )
)
levels(kcv_class)[[1]] <- lev_cv

# -------------------------------------------------------------
# 5. Resistance & recovery – σ-based classes
# -------------------------------------------------------------
# Clamp extreme anomalies before classification
R_clamp   <- clamp(R_hotdry,   lower = -3, upper = 3, values = TRUE)
Rec_clamp <- clamp(Rec_hotdry, lower = -3, upper = 3, values = TRUE)

# Breaks in SD units
br_R <- c(-3, -1.0, -0.5, 0.5, 1.5, 3)

R_class <- app(R_clamp, fun = function(x) {
  as.integer(cut(
    x,
    breaks = br_R,
    labels = FALSE,
    include.lowest = TRUE
  ))
})
names(R_class) <- "R_class"

lev_R <- data.frame(
  value = 1:5,
  class = c(
    "strong loss during hot–dry (≤ −1.0 σ)",
    "moderate loss (−1.0 to −0.5 σ)",
    "near-neutral (−0.5 to 0.5 σ)",
    "buffered (>0.5 to 1.5 σ)",
    "over-compensating (> 1.5 σ)"
  )
)
levels(R_class)[[1]] <- lev_R

Rec_class <- app(Rec_clamp, fun = function(x) {
  as.integer(cut(
    x,
    breaks = br_R,
    labels = FALSE,
    include.lowest = TRUE
  ))
})
names(Rec_class) <- "Rec_class"

lev_Rec <- data.frame(
  value = 1:5,
  class = c(
    "strongly incomplete recovery (≤ −1.0 σ)",
    "partial recovery (−1.0 to −0.5 σ)",
    "full recovery (−0.5 to 0.5 σ)",
    "good recovery (>0.5 to 1.5 σ)",
    "over-recovery (> 1.5 σ)"
  )
)
levels(Rec_class)[[1]] <- lev_Rec

# -------------------------------------------------------------
# 6. n_hotdry – exposure to compound extremes
# -------------------------------------------------------------
nhot_class <- app(n_hotdry, fun = function(x) {
  as.integer(cut(
    x,
    breaks = c(-0.5, 0.5, 3.5, 6.5, 50),
    labels = FALSE,
    include.lowest = TRUE
  ))
})

names(nhot_class) <- "nhot_class"

lev_nhot <- data.frame(
  value = 1:4,
  class = c(
    "no hot–dry events",
    "occasional (1–3)",
    "frequent (4–6)",
    "very frequent (≥7)"
  )
)
levels(nhot_class)[[1]] <- lev_nhot

# -------------------------------------------------------------
# 7. Save all classified layers
# -------------------------------------------------------------
dir.create("data/metrics_veg/classes", showWarnings = FALSE, recursive = TRUE)

class_stack <- c(
  kmean_class,
  ktrend_class,
  kcv_class,
  R_class,
  Rec_class,
  nhot_class
)

writeRaster(
  class_stack,
  "data/metrics_veg/classes/kNDVI_functioning_resilience_classes_NEA_forestOnly.tif",
  overwrite = TRUE
)

###############################################################
# End of script
###############################################################








###############################################################
# Save kNDVI functioning & resilience classes for QGIS
###############################################################

library(terra)

# If you still have the individual class rasters in memory:
# kmean_class, ktrend_class, kcv_class, R_class, Rec_class, nhot_class

# 1. Make sure each layer has a clear name
names(kmean_class) <- "kmean_class"
names(ktrend_class) <- "ktrend_class"
names(kcv_class)    <- "kcv_class"
names(R_class)      <- "R_class"
names(Rec_class)    <- "Rec_class"
names(nhot_class)   <- "nhot_class"

# (Optional) check levels are still attached
# levels(kmean_class)
# levels(ktrend_class)
# ...

# 2. Build a stack with proper band names (optional but useful)
class_stack <- c(
  kmean_class,
  ktrend_class,
  kcv_class,
  R_class,
  Rec_class,
  nhot_class
)

# Check
class_stack
names(class_stack)
levels(class_stack[[1]])  # e.g. productivity classes

# 3. Create folder for class rasters
dir.create("data/metrics_veg/classes", showWarnings = FALSE, recursive = TRUE)

# 4. Save the full multi-band stack (one file, six bands)
writeRaster(
  class_stack,
  "data/metrics_veg/classes/kNDVI_functioning_resilience_classes_NEA_forestOnly.tif",
  overwrite = TRUE,
  datatype  = "INT1U"  # small integer codes, good for classes
)

# 5. Save each layer separately – easiest for QGIS styling
writeRaster(
  kmean_class,
  "data/metrics_veg/classes/kNDVI_mean_class.tif",
  overwrite = TRUE,
  datatype  = "INT1U"
)

writeRaster(
  ktrend_class,
  "data/metrics_veg/classes/kNDVI_trend_class.tif",
  overwrite = TRUE,
  datatype  = "INT1U"
)

writeRaster(
  kcv_class,
  "data/metrics_veg/classes/kNDVI_CV_class.tif",
  overwrite = TRUE,
  datatype  = "INT1U"
)

writeRaster(
  R_class,
  "data/metrics_veg/classes/kNDVI_R_hotdry_class.tif",
  overwrite = TRUE,
  datatype  = "INT1U"
)

writeRaster(
  Rec_class,
  "data/metrics_veg/classes/kNDVI_Rec_hotdry_class.tif",
  overwrite = TRUE,
  datatype  = "INT1U"
)

writeRaster(
  nhot_class,
  "data/metrics_veg/classes/kNDVI_n_hotdry_class.tif",
  overwrite = TRUE,
  datatype  = "INT1U"
)

###############################################################
# End of saving block
###############################################################


















