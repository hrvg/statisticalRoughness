library("statisticalRoughness")
library("tictoc")

# PARAMETERS
tic()
print("Saving parameters...")
out_dir <- "out10"
if(!dir.exists(out_dir)) dir.create(out_dir)

.selected <- c("beta2", "alpha1.y", "alpha1.x", "alpha2.y", "alpha2.x", "w.y", "w.x", "xi.y", "xi.x", "zeta1", "zeta2", "median_Pe", "H", "q", "geol_diversity", "inv.fc")
cparams <- data.frame(
  selected = .selected,
  n_sigma = rep(3, length(.selected)),
  lower = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  upper = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  lower_clamp = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 1E-16, -Inf, -Inf, -Inf),
  upper_clamp = c(0, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf)
)
print(cparams)
saveRDS(cparams, file.path(out_dir, "cparams.Rds"))
print("Parameters saved.")
toc()

# LOADING
tic()
print("Loading zeta results...")
results_directory <- file.path("F:/hguillon/research/exploitation/out/run140/")
zeta_results <- read_zeta_raster(out_dir = results_directory, Lmax = 1E4, raster_resolution = 10, .len = 32, proxy = FALSE)
zeta_results$spatial_scales <- zeta_results$spatial_scales[-1]
band_names <- read.csv(file.path(results_directory, "band_names.csv")) %>% names()
if ("theta" %in% .selected) zeta_results$raster_list <- zeta_results$raster_list %>% make_angular(match("theta", band_names))
print("Zeta results loaded.")
toc()

# CLAMPING
tic()
print("Clamping rasters...")
sliced_clamped_raster <- slice_clamp(
  raster_list = zeta_results$raster_list,
  att_names = band_names,
  selected = .selected,
  clamp_raster = TRUE, 
  clamp_params = cparams)
print("Raster clamped.")
toc()

# MASK q from H
tic()
print("Masking rasters...")
sliced_clamped_raster <- mask_layer_from_layer(sliced_clamped_raster, match("q", .selected), match("H", .selected), option = "NA")
print("Rasters masked.")
toc()

# # CLIPING CENTRAL VALLEY
tic()
print("Clipping Central Valley...")
shp <- raster::shapefile("F:/hguillon/research/data/safe_water/central_valley_alluvial_boundary/Alluvial_Bnd.shp")
shp <- sp::spTransform(shp, raster::crs(sliced_clamped_raster[[31]] %>% as("Raster")))
sliced_clamped_raster <- lapply(seq_along(sliced_clamped_raster), function(n){
	s <- sliced_clamped_raster[[n]] %>% as("Raster")
	s <- raster::mask(s, shp, inverse = TRUE)
	return(stars::st_as_stars(s, proxy = FALSE))
})  
print("Clipped Central Valley.")
toc()

# PAR PERP
tic()
print("Deriving parallel and perpendicular directions...")
sliced_clamped_raster <- find_par_perp(sliced_clamped_raster, match("alpha1.x", .selected), match("alpha1.y", .selected))
.selected <- c(.selected, "alpha_perp1", "alpha_par1")
sliced_clamped_raster <- find_par_perp(sliced_clamped_raster, match("alpha2.x", .selected), match("alpha2.y", .selected))
.selected <- c(.selected, "alpha_perp2", "alpha_par2")
print("Parallel and perpendicular directions derived.")
toc()

# CORRECT ZETA
tic()
print("Correcting zeta...")
sliced_clamped_raster <- correct_zeta(sliced_clamped_raster, match("zeta1", .selected))
.selected <- c(.selected, "zeta_perp1")
sliced_clamped_raster <- correct_zeta(sliced_clamped_raster, match("zeta2", .selected))
.selected <- c(.selected, "zeta_perp2")
print("Corrected zeta.")
toc()

# LOG TRANSFORM
tic()
print("Applying log-transforms...")
.sliced_clamped_raster <- sliced_clamped_raster
for (x in c("median_Pe", "w.x", "w.y", "xi.x", "xi.y", "inv.fc")){
    .sliced_clamped_raster <- logtransform_layer(.sliced_clamped_raster, match(x, .selected))
}
print("Applied log-transforms.")
toc()

# REDUCE NOISE
# tic()
# print("Reducing noise...")
# sliced_clamped_raster <- reduce_spatial_noise(sliced_clamped_raster, .NAonly = TRUE)
# print("Reduced noise...")
# toc()

# FOUR VALUES CHECK
tic()
print("Performing four-values check...")
ind <- four_values_check(.sliced_clamped_raster)
.sliced_clamped_raster <- .sliced_clamped_raster[ind]
spatial_scales <- zeta_results$spatial_scales[ind]
print("Performed four-values check.")
toc()

# OUTPUT: MAPS
tic()
print("Outputting maps...")
for(x in .selected){
  print(x)
  raster_select(sliced_clamped_raster, band_id = match(x, .selected)) %>% 
  save_map(ttl = x, groups = zeta_results$spatial_scales, begin = 0.1, end = 0.95, direction = -1, option = "viridis", col_style = "cont", out_path = out_dir)
}
print("Outputted maps.")
toc()

# OUTPUT: DISTRIBUTION
tic()
print("Outputting distributions...")
pdf(file.path(out_dir, "distributions.pdf"))
for(x in .selected){
  gb <- make_all_plots(.sliced_clamped_raster, spatial_scales, match(x, .selected), x, begin = 0.1, end = 0.85, direction = 1, option = "viridis")
  layout_mat <- matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3), byrow = FALSE, nrow = 5, ncol = 2)
  gridExtra::grid.arrange(grobs = gb, layout_matrix = layout_mat)
}
dev.off()
print("Outputted distributions.")
toc()

# OUTPUT: PEARSON CORRELATIONS
tic()
print("Outputting correlations...")
crossscale_corr <- crossscale_correlations(
  raster_list = .sliced_clamped_raster,
  selected = .selected,
  spatial_scales = spatial_scales,
  corr_type = "pearson")
pdf(file.path(out_dir, "pearson.pdf"))
crossscale_corr$p
dev.off()
saveRDS(crossscale_corr$graphics_df, file.path(out_dir, "graphics_df.Rds"))
print("Outputted correlations.")
toc()