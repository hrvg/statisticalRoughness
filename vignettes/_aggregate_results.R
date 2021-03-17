library("statisticalRoughness")

# PARAMETERS
out_dir <- "out8"
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

# LOADING
results_directory <- file.path("F:/hguillon/research/exploitation/out/run140/")
zeta_results <- read_zeta_raster(out_dir = results_directory, Lmax = 1E4, raster_resolution = 10, .len = 32, proxy = FALSE)
zeta_results$spatial_scales <- zeta_results$spatial_scales[-1]
band_names <- read.csv(file.path(results_directory, "band_names.csv")) %>% names()
if ("theta" %in% .selected) zeta_results$raster_list <- zeta_results$raster_list %>% make_angular(match("theta", band_names))

# CLAMPING
sliced_clamped_raster <- slice_clamp(
  raster_list = zeta_results$raster_list,
  att_names = band_names,
  selected = .selected,
  clamp_raster = TRUE, 
  clamp_params = cparams)

# MASK q from H
sliced_clamped_raster <- mask_layer_from_layer(sliced_clamped_raster, match("q", .selected), match("H", .selected), option = "NA")

# CLIPING CENTRAL VALLEY
shp <- raster::shapefile("F:/hguillon/research/data/safe_water/central_valley_alluvial_boundary/Alluvial_Bnd.shp")
shp <- sp::spTransform(shp, raster::crs(sliced_clamped_raster[[31]] %>% as("Raster")))
sliced_clamped_raster <- lapply(seq_along(sliced_clamped_raster), function(n){
	s <- sliced_clamped_raster[[n]] %>% as("Raster")
	s <- raster::mask(s, shp, inverse = TRUE)
	return(stars::st_as_stars(s, proxy = FALSE))

# PAR PERP
sliced_clamped_raster <- find_par_perp(sliced_clamped_raster, match("alpha1.x", .selected), match("alpha1.y", .selected))
.selected <- c(.selected, "alpha_perp1", "alpha_par1")
sliced_clamped_raster <- find_par_perp(sliced_clamped_raster, match("alpha2.x", .selected), match("alpha2.y", .selected))
.selected <- c(.selected, "alpha_perp2", "alpha_par2")

# CORRECT ZETA
sliced_clamped_raster <- correct_zeta(sliced_clamped_raster, match("zeta1", .selected))
.selected <- c(.selected, "zeta_perp1")
sliced_clamped_raster <- correct_zeta(sliced_clamped_raster, match("zeta2", .selected))
.selected <- c(.selected, "zeta_perp2")


# LOG TRANSFORM
.sliced_clamped_raster <- sliced_clamped_raster
for (x in c("median_Pe", "w.x", "w.y", "xi.x", "xi.y", "inv.fc")){
    .sliced_clamped_raster <- logtransform_layer(.sliced_clamped_raster, match(x, .selected))
}

# REDUCE NOISE
sliced_clamped_raster <- reduce_spatial_noise(sliced_clamped_raster)

# FOUR VALUES CHECK
ind <- four_values_check(.sliced_clamped_raster)
.sliced_clamped_raster <- .sliced_clamped_raster[ind]
spatial_scales <- zeta_results$spatial_scales[ind]

# OUTPUT: MAPS
for(x in .selected){
  print(x)
  raster_select(sliced_clamped_raster, band_id = match(x, .selected)) %>% 
  save_map(ttl = x, groups = zeta_results$spatial_scales, begin = 0.1, end = 0.95, direction = -1, option = "viridis", col_style = "cont", out_path = out_dir)
}

# OUTPUT: DISTRIBUTION
pdf(file.path(out_dir, "distributions.pdf"))
for(x in .selected){
  gb <- make_all_plots(.sliced_clamped_raster, spatial_scales, match(x, .selected), x, begin = 0.1, end = 0.85, direction = 1, option = "viridis")
  layout_mat <- matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3), byrow = FALSE, nrow = 5, ncol = 2)
  gridExtra::grid.arrange(grobs = gb, layout_matrix = layout_mat)
}
dev.off()

# OUTPUT: PEARSON CORRELATIONS
pdf(file.path(out_dir, "pearson.pdf"))
crossscale_correlations(
  raster_list = .sliced_clamped_raster,
  selected = .selected,
  spatial_scales = spatial_scales,
  corr_type = "pearson")
dev.off()