# ---------------------------------------------------------------
# 01_belgium_spatial_baseline.R
#
# Spatial baseline for the Fulbright Belgium / Ghent University
# small-area hypertension proposal.
#
# What this script does, in order:
#   1. Loads the Statbel statistical-sector shapefile.
#   2. Cleans and validates the geometry.
#   3. Simulates a spatially-structured hypertension prevalence
#      outcome at the sector level.
#   4. Runs Global Moran's I (analytic + 999-permutation MC).
#   5. Fits a baseline SAR lag model as a smoke test.
#
# Author: Romario Joseph (rjoseph3@bu.edu)
# Stack:  R 4.x, sf, spdep, dplyr, ggplot2
# Data:   Public Statbel sector layer. Outcome is simulated.
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(spdep)
  library(dplyr)
  library(ggplot2)
})

set.seed(42)

# ---------------------------------------------------------------
# 0. Paths and config
# ---------------------------------------------------------------
DATA_DIR    <- "data"
OUTPUTS_DIR <- "outputs"
dir.create(DATA_DIR,    showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Statbel 2023 statistical-sector layer (open data).
# If the URL ever changes, update it in data/README.md too.
STATBEL_URL <- "https://statbel.fgov.be/sites/default/files/files/opendata/Statistische%20sectoren/sh_statbel_statistical_sectors_3812_20230101.shp.zip"
LOCAL_ZIP   <- file.path(DATA_DIR, "statbel_sectors_2023.zip")
LOCAL_SHP   <- file.path(DATA_DIR, "statbel_sectors_2023.shp")

# Belgian Lambert 2008, the official planar CRS for Belgium.
TARGET_CRS  <- 3812

# ---------------------------------------------------------------
# 1. Geometry ingest
# ---------------------------------------------------------------
if (!file.exists(LOCAL_SHP)) {
  if (!file.exists(LOCAL_ZIP)) {
    message("Downloading Statbel sector layer ...")
    download.file(STATBEL_URL, LOCAL_ZIP, mode = "wb")
  }
  unzip(LOCAL_ZIP, exdir = DATA_DIR)
}

sectors_raw <- sf::st_read(LOCAL_SHP, quiet = TRUE)

# Pull a stable sector ID. Statbel uses CS01012023 / CD_REFNIS_S
# depending on release; we pick whichever is present.
id_col <- intersect(
  c("CS01012023", "CD_SECTOR", "cs012023", "CNIS5_2023"),
  names(sectors_raw)
)[1]
if (is.na(id_col)) stop("Could not locate a sector ID column.")

sectors_raw$sector_id <- as.character(sectors_raw[[id_col]])

# ---------------------------------------------------------------
# 2. Geometry cleaning + reprojection
# ---------------------------------------------------------------
n_before <- nrow(sectors_raw)

sectors <- sectors_raw |>
  sf::st_make_valid() |>
  dplyr::filter(!sf::st_is_empty(geometry)) |>
  dplyr::mutate(area_m2 = as.numeric(sf::st_area(geometry))) |>
  dplyr::filter(area_m2 > 0) |>
  sf::st_transform(TARGET_CRS)

n_after <- nrow(sectors)
message(sprintf(
  "Sector cleaning: kept %d / %d polygons (%d dropped).",
  n_after, n_before, n_before - n_after
))

# ---------------------------------------------------------------
# 3. Simulate a spatially-structured outcome
# ---------------------------------------------------------------
# Strategy: draw a smooth spatial trend by averaging each polygon's
# neighbors' centroids' coordinates (cheap proxy for a Gaussian
# random field), pass through a logit so values live in (0, 1),
# and calibrate the national mean to ~0.28 (recent Belgian HIS).

centroids <- suppressWarnings(sf::st_centroid(sectors))
coords    <- sf::st_coordinates(centroids)

# Standardize coordinates so the latent surface is unitless.
xs <- scale(coords[, 1])[, 1]
ys <- scale(coords[, 2])[, 1]

latent_trend <- 0.8 * xs - 0.4 * ys + 0.3 * xs * ys
noise        <- rnorm(nrow(sectors), mean = 0, sd = 0.4)

# Calibrate intercept so mean prevalence is roughly 0.28.
linpred  <- latent_trend + noise
intercept <- log(0.28 / (1 - 0.28)) - mean(linpred)
sectors$htn_prev <- plogis(intercept + linpred)

# Sanity checks before we touch Moran.
stopifnot(all(!is.na(sectors$htn_prev)))
stopifnot(all(sectors$htn_prev > 0 & sectors$htn_prev < 1))

# Also simulate a "deprivation index" covariate for the SAR demo.
sectors$deprivation <- 0.6 * xs + rnorm(nrow(sectors), 0, 0.5)

# ---------------------------------------------------------------
# 4. Spatial weights + Global Moran's I
# ---------------------------------------------------------------
message("Building queen-contiguity neighbors ...")
nb <- spdep::poly2nb(sectors, queen = TRUE)

# Some Belgian sectors are coastal islands and have zero neighbors.
# We log them and drop them from the Moran test rather than silently
# zero-weighting.
n_islands <- sum(spdep::card(nb) == 0)
message(sprintf("Detected %d sector(s) with zero neighbors.", n_islands))

keep <- spdep::card(nb) > 0
sectors_keep <- sectors[keep, ]
nb_keep      <- spdep::poly2nb(sectors_keep, queen = TRUE)
lw           <- spdep::nb2listw(nb_keep, style = "W", zero.policy = FALSE)

moran_analytic <- spdep::moran.test(
  sectors_keep$htn_prev, lw,
  randomisation = TRUE,
  alternative   = "greater"
)

moran_mc <- spdep::moran.mc(
  sectors_keep$htn_prev, lw,
  nsim = 999,
  alternative = "greater"
)

message("Global Moran's I (analytic):  I = ",
        round(moran_analytic$estimate[1], 4),
        ", p = ", signif(moran_analytic$p.value, 3))
message("Global Moran's I (999 perms): I = ",
        round(moran_mc$statistic, 4),
        ", p = ", signif(moran_mc$p.value, 3))

# ---------------------------------------------------------------
# 5. Moran scatterplot
# ---------------------------------------------------------------
lag_htn <- spdep::lag.listw(lw, sectors_keep$htn_prev)

scatter_df <- data.frame(
  htn     = sectors_keep$htn_prev,
  lag_htn = lag_htn
)

p <- ggplot2::ggplot(scatter_df, ggplot2::aes(htn, lag_htn)) +
  ggplot2::geom_point(alpha = 0.4, size = 0.7) +
  ggplot2::geom_smooth(method = "lm", se = FALSE) +
  ggplot2::geom_hline(yintercept = mean(scatter_df$lag_htn),
                      linetype = "dashed") +
  ggplot2::geom_vline(xintercept = mean(scatter_df$htn),
                      linetype = "dashed") +
  ggplot2::labs(
    title = "Moran scatterplot, simulated Belgian sector HTN prevalence",
    x     = "Sector prevalence",
    y     = "Spatial lag of prevalence"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(
  file.path(OUTPUTS_DIR, "moran_scatterplot.png"),
  p, width = 6, height = 5, dpi = 150
)

# ---------------------------------------------------------------
# 6. Baseline SAR lag model (smoke test)
# ---------------------------------------------------------------
sar_fit <- tryCatch(
  spdep::spautolm(
    htn_prev ~ deprivation,
    data    = sectors_keep,
    listw   = lw,
    family  = "SAR"
  ),
  error = function(e) {
    message("SAR fit failed: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(sar_fit)) {
  message("SAR lag model fit OK. lambda = ",
          round(sar_fit$lambda, 4))
}

# ---------------------------------------------------------------
# 7. Write results table
# ---------------------------------------------------------------
results <- data.frame(
  metric = c("Moran I (analytic)", "Moran I (MC, 999)",
             "p (analytic)", "p (MC)", "n sectors used", "n islands dropped"),
  value  = c(
    round(moran_analytic$estimate[1], 4),
    round(moran_mc$statistic, 4),
    signif(moran_analytic$p.value, 3),
    signif(moran_mc$p.value, 3),
    nrow(sectors_keep),
    n_islands
  )
)

write.csv(results,
          file.path(OUTPUTS_DIR, "moran_results.csv"),
          row.names = FALSE)

message("Done. See outputs/ for the scatterplot and the results CSV.")
