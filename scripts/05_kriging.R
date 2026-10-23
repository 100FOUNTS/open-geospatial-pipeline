# Step 5. Geostatistics in R: variogram, kriging, IDW and Thiessen, with validation.
#
# Idea: take a sample of summer NDMI points (as if they were field sensors), hold out 20 % for
# validation and compare three ways of reconstructing the full map.
# If you have real measurements (for example soil-moisture sensors), replace "sample" with your
# points: a file with an "ndmi" column (or change the variable name).
#
# Since the full NDMI already exists, each method can be graded: its estimate is compared with
# the true value at the held-out points (which the method never saw).
#
# The three methods:
#   Thiessen : copies the value of the nearest sensor (a patchwork map).
#   IDW      : averages nearby sensors, weighting the closer ones more (a fixed formula).
#   Kriging  : averages with weights learned from the data (the variogram) and gives the uncertainty.
#
# Run from the project root:   Rscript scripts/05_kriging.R

suppressPackageStartupMessages({
  library(terra)   # rasters
  library(sf)      # points and polygons
  library(gstat)   # variogram, kriging, IDW
  library(spdep)   # Moran's I
  library(yaml)    # read config.yml
})

cfg <- read_yaml("config.yml")
k   <- cfg$kriging
set.seed(k$seed)               # same seed = same random points every time
dir.create("outputs", showWarnings = FALSE)

# ---- 1. Data ---------------------------------------------------------------
# The summer NDMI is clipped to the basin(s) and n random points are drawn inside.
r <- rast("data/processed/ndmi_summer.tif")
names(r) <- "ndmi"
basins <- st_read("data/processed/basins.gpkg", quiet = TRUE)
basins <- st_transform(basins, crs(r))
r <- mask(crop(r, vect(basins)), vect(basins))

n_cells <- global(!is.na(r), "sum")$sum
n <- min(k$n_points, n_cells)
sample_pts <- spatSample(r, size = n, method = "random", na.rm = TRUE,
                         as.points = TRUE, values = TRUE)
sample_pts <- st_as_sf(sample_pts)

# 80 % to fit the methods (training), 20 % held out to grade them (test).
idx   <- sample(seq_len(nrow(sample_pts)), size = round(k$train_fraction * nrow(sample_pts)))
train <- sample_pts[idx, ]
test  <- sample_pts[-idx, ]
cat(sprintf("Points: %d training, %d test\n", nrow(train), nrow(test)))

# ---- 2. Spatial autocorrelation (Moran's I) --------------------------------
# Measures whether nearby points have similar values (-1 to 1). If it is high and its p-value low,
# there is spatial structure to exploit. "Neighbours" = the 8 nearest points of each point.
coords    <- st_coordinates(sample_pts)
neighbours <- knn2nb(knearneigh(coords, k = 8))
weights   <- nb2listw(neighbours, style = "W")
moran     <- moran.test(sample_pts$ndmi, weights)
cat(sprintf("Moran's I = %.3f (p = %.2g)\n", moran$estimate[["Moran I statistic"]], moran$p.value))

# ---- 3. Variogram -----------------------------------------------------------
# Empirical variogram: for pairs of points at each distance, how much their values differ.
# A model is then fitted (spherical, exponential or Gaussian) and the one with the lowest error (SSErr) is kept.
#   nugget = difference between almost touching points (noise)
#   range  = distance beyond which points no longer resemble each other
#   sill   = level where the curve flattens
v <- variogram(ndmi ~ 1, train)
models <- c("Sph", "Exp", "Gau")
fits <- lapply(models, function(m) {
  tryCatch(
    fit.variogram(v, vgm(psill = max(v$gamma), model = m,
                         range = max(v$dist) / 3, nugget = min(v$gamma))),
    error = function(e) NULL, warning = function(w) NULL)
})
sse <- sapply(fits, function(a) if (is.null(a)) Inf else attr(a, "SSErr"))
if (all(is.infinite(sse))) stop("No variogram model could be fitted. Check the data.")
best <- fits[[which.min(sse)]]
cat(sprintf("Chosen model: %s | range = %.0f m | nugget = %.4f | partial sill = %.4f\n",
            as.character(best$model[2]), best$range[2], best$psill[1], best$psill[2]))

png("outputs/variogram.png", width = 900, height = 600, res = 120)
print(plot(v, best, main = "Summer NDMI variogram"))
dev.off()

# ---- 4. Validation on held-out points ---------------------------------------
# Each method estimates the NDMI at the test points using only the training set.
# nmax = how many neighbouring sensors it uses per estimate.
pred_kri <- krige(ndmi ~ 1, train, test, model = best, nmax = k$max_neighbors, debug.level = 0)$var1.pred
pred_idw <- idw(ndmi ~ 1, train, test, idp = 2, nmax = k$max_neighbors, debug.level = 0)$var1.pred
pred_thi <- idw(ndmi ~ 1, train, test, idp = 2, nmax = 1, debug.level = 0)$var1.pred   # Thiessen = nearest neighbour

# RMSE and MAE: mean error (lower = better). Bias: >0 overestimates, <0 underestimates.
metrics <- function(obs, pred) c(RMSE = sqrt(mean((obs - pred)^2)),
                                 MAE  = mean(abs(obs - pred)),
                                 Bias = mean(pred - obs))
tab <- rbind(Kriging  = metrics(test$ndmi, pred_kri),
             IDW      = metrics(test$ndmi, pred_idw),
             Thiessen = metrics(test$ndmi, pred_thi))
tab <- round(tab, 4)
cat("\nError at the test points:\n")
print(tab)
write.csv(tab, "outputs/interpolation_metrics.csv")

# Cross-validation: drop each training point, estimate it from the rest and compare.
cv <- krige.cv(ndmi ~ 1, train, model = best, nmax = k$max_neighbors)
cat(sprintf("\nKriging cross-validation (training): RMSE = %.4f\n", sqrt(mean(cv$residual^2))))
if (tab["Kriging", "RMSE"] >= tab["IDW", "RMSE"]) {
  cat("Warning: kriging does not beat IDW. The variogram adds little on this data.\n")
}

# ---- 5. Full maps -----------------------------------------------------------
# NDMI is estimated at the centre of each basin cell (a grid) and saved as a raster.
# Kriging returns two things: the prediction (var1.pred) and its variance (var1.var), which is
# the uncertainty: higher far from the sensors.
cells <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
grid  <- st_as_sf(cells, coords = c("x", "y"), crs = st_crs(train))

kri   <- krige(ndmi ~ 1, train, grid, model = best, nmax = k$max_neighbors, debug.level = 0)
idw_m <- idw(ndmi ~ 1, train, grid, idp = 2, nmax = k$max_neighbors, debug.level = 0)

to_raster <- function(pts, field, name) {
  rr <- rasterize(vect(pts), r, field = field)
  names(rr) <- name
  rr
}
writeRaster(to_raster(kri,   "var1.pred", "ndmi_kriging"),  "data/processed/ndmi_kriging.tif",     overwrite = TRUE)
writeRaster(to_raster(kri,   "var1.var",  "ndmi_variance"), "data/processed/ndmi_kriging_var.tif", overwrite = TRUE)
writeRaster(to_raster(idw_m, "var1.pred", "ndmi_idw"),      "data/processed/ndmi_idw.tif",         overwrite = TRUE)

# The points are saved so they can be viewed in QGIS (where the "sensors" were).
st_write(train, "data/processed/training_points.gpkg", delete_dsn = TRUE, quiet = TRUE)
st_write(test,  "data/processed/test_points.gpkg",     delete_dsn = TRUE, quiet = TRUE)

cat("\nDone. New files in data/processed/ and outputs/\n")
