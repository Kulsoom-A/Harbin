###############################################################
# Replot Figures 7 and 8 by vegetation type only
# Supervisor comment: do not compare on national level;
# compare different vegetation types and use 9-point filtering.
###############################################################

library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# -------------------------------------------------------------
# 1. Inputs already saved from previous processing
# -------------------------------------------------------------

wb_diff_path <- "data/metrics_climate/WB_diff_early_mid_late_NEA.tif"
tmean_diff_path <- "data/metrics_climate/Tmean_diff_early_mid_late_NEA.tif"

fType_path <- "data/processed/LC5_forest_1to6_2024_NEA_0p1deg.tif"

wb_diff <- rast(wb_diff_path)
tmean_diff <- rast(tmean_diff_path)
fType <- rast(fType_path)

levels(fType) <- NULL
names(fType) <- "forest_type"

# Align forest type raster if needed
if (!compareGeom(wb_diff, fType, stopOnError = FALSE)) {
  fType_wb <- resample(fType, wb_diff[[1]], method = "near")
} else {
  fType_wb <- fType
}

if (!compareGeom(tmean_diff, fType, stopOnError = FALSE)) {
  fType_tmean <- resample(fType, tmean_diff[[1]], method = "near")
} else {
  fType_tmean <- fType
}

# -------------------------------------------------------------
# 2. Vegetation labels and filtering threshold
# -------------------------------------------------------------

ftype_labels <- c(
  "1" = "ENT",
  "2" = "EBT",
  "3" = "DNT",
  "4" = "DBT",
  "5" = "SHB",
  "6" = "GRS"
)

min_cells <- 9   # 9-point/cell filtering as suggested

# -------------------------------------------------------------
# 3. Function to summarise period differences by vegetation type
# -------------------------------------------------------------

summarise_by_vegetation_type <- function(diff_stack, ftype_raster, variable_name) {
  
  names(diff_stack) <- c(
    "diff_2011_2020_minus_2001_2010",
    "diff_2021_2024_minus_2001_2010"
  )
  
  s <- c(diff_stack, ftype_raster)
  names(s)[3] <- "forest_type"
  
  df <- as.data.frame(s, na.rm = TRUE)
  
  df$forest_type <- as.integer(df$forest_type)
  df$forest_type_abbr <- ftype_labels[as.character(df$forest_type)]
  
  df <- df |>
    filter(!is.na(forest_type_abbr))
  
  df_long <- df |>
    pivot_longer(
      cols = c(
        diff_2011_2020_minus_2001_2010,
        diff_2021_2024_minus_2001_2010
      ),
      names_to = "period",
      values_to = "value"
    ) |>
    mutate(
      period = recode(
        period,
        "diff_2011_2020_minus_2001_2010" = "2011–2020 − 2001–2010",
        "diff_2021_2024_minus_2001_2010" = "2021–2024 − 2001–2010"
      ),
      period = factor(
        period,
        levels = c(
          "2011–2020 − 2001–2010",
          "2021–2024 − 2001–2010"
        )
      ),
      forest_type_abbr = factor(
        forest_type_abbr,
        levels = c("ENT", "EBT", "DNT", "DBT", "SHB", "GRS")
      )
    )
  
  summary_tbl <- df_long |>
    group_by(forest_type_abbr, period) |>
    summarise(
      variable = variable_name,
      n_cells = sum(!is.na(value)),
      med = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(n_cells >= min_cells)
  
  summary_tbl
}

# -------------------------------------------------------------
# 4. Summarise WB and Tmean by vegetation type only
# -------------------------------------------------------------

wb_veg_summary <- summarise_by_vegetation_type(
  diff_stack = wb_diff,
  ftype_raster = fType_wb,
  variable_name = "Water balance"
)

tmean_veg_summary <- summarise_by_vegetation_type(
  diff_stack = tmean_diff,
  ftype_raster = fType_tmean,
  variable_name = "Mean annual temperature"
)

dir.create("outputs/tables2", showWarnings = FALSE, recursive = TRUE)

write_csv(
  wb_veg_summary,
  "outputs/tables2/WB_period_differences_by_vegetationType_2001_2024.csv"
)

write_csv(
  tmean_veg_summary,
  "outputs/tables2/Tmean_period_differences_by_vegetationType_2001_2024.csv"
)

# -------------------------------------------------------------
# 5. Common plotting function
# -------------------------------------------------------------

plot_vegetation_change <- function(summary_tbl, y_lab, out_file) {
  
  p <- ggplot(
    summary_tbl,
    aes(
      x = forest_type_abbr,
      y = med,
      colour = forest_type_abbr
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_linerange(
      aes(ymin = q25, ymax = q75),
      linewidth = 0.8
    ) +
    geom_point(size = 2.8) +
    facet_wrap(~ period, ncol = 1, strip.position = "top") +
    scale_colour_brewer(
      palette = "Dark2",
      name = "Vegetation type"
    ) +
    labs(
      x = "Vegetation type",
      y = y_lab
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(face = "bold", size = 11),
      axis.title.y = element_text(face = "bold", size = 11),
      axis.text.x  = element_text(face = "bold", colour = "black", size = 10),
      axis.text.y  = element_text(face = "bold", colour = "black", size = 10),
      strip.text   = element_text(face = "bold", size = 11),
      legend.title = element_text(face = "bold"),
      legend.text  = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 6.5,
    height = 5.2,
    dpi = 300
  )
  
  return(p)
}

# -------------------------------------------------------------
# 6. Figure 7: Water balance by vegetation type
# -------------------------------------------------------------

p_wb_veg <- plot_vegetation_change(
  summary_tbl = wb_veg_summary,
  y_lab = "Change in annual water balance (mm)",
  out_file = "outputs/figures/Figure7_WB_period_differences_by_vegetationType.png"
)

p_wb_veg

# -------------------------------------------------------------
# 7. Figure 8: Mean annual temperature by vegetation type
# -------------------------------------------------------------

p_tmean_veg <- plot_vegetation_change(
  summary_tbl = tmean_veg_summary,
  y_lab = "Change in mean annual temperature (°C)",
  out_file = "outputs/figures/Figure8_Tmean_period_differences_by_vegetationType.png"
)

p_tmean_veg

