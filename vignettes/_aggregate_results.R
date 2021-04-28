library("statisticalRoughness")
library("tictoc")
library("foreach")

# NEW FUNCTIONS
is_band1_gtlt_band2 <- function(band1, band2, rasters, .selected, spatial_scales, option = "gt", knit = TRUE){
  allx <- raster_select(rasters, band_id = match(band1, .selected))$rasters
  ally <- raster_select(rasters, band_id = match(band2, .selected))$rasters
  xy <- foreach(x = allx, y = ally, .combine = rbind) %do% {
    x <- as(x, "Raster")
    y <- as(y, "Raster")
    s <- raster::stack(x, y)
    if (option == "gt") res <- raster::calc(s, fun = function(x){x[1] > x[2]}) %>% raster::getValues() %>% rstatix::freq_table() %>% dplyr::filter(group == "TRUE") %>% dplyr::select(-group)
    if (option == "lt") res <- raster::calc(s, fun = function(x){x[1] < x[2]}) %>% raster::getValues() %>% rstatix::freq_table() %>% dplyr::filter(group == "TRUE") %>% dplyr::select(-group)
    if(nrow(res) == 0) res <- data.frame(n = 0, prop = 0)
    res
  }
  mean_xy <- xy %>% dplyr::summarize(n = mean(n), prop = mean(prop))
  xy <- rbind(xy, mean_xy)
  rownames(xy) <- c(spatial_scales, "mean")
  if (knit){
    xy  %>% 
    knitr::kable("html", caption = paste0("Is ", band1, " ", ifelse(option == "gt", "greater", "lesser"), " than ", band2, "?")) %>% 
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover"),
      position = "center",
      full_width = TRUE)
  } else {
    xy
  }
}

# PARAMETERS
tic()
print("Saving parameters...")
ZETA_MASKING <- FALSE
REDUCE_NOISE <- FALSE
NOISE_FUN <- mean
PAR_PERP <- FALSE
MAPS <- TRUE
DIST <- TRUE
CORR <- FALSE
SUB <- TRUE
COMPARE <- FALSE
AVERAGE <- TRUE

out_dir <- paste("out", ifelse(ZETA_MASKING, "zeta-masking", "no-masking"), ifelse(REDUCE_NOISE, "noise-reduction", "no-noise-reduction"), sep = "_")
if(!dir.exists(out_dir)) dir.create(out_dir)

.selected <- c("beta2", "alpha1.y", "alpha1.x", "alpha2.y", "alpha2.x", "w.y", "w.x", "xi.y", "xi.x", "zeta1", "zeta2", "theta.x", "theta.y", "q", "H", "inv.fc", "median_Pe")
cparams <- data.frame(
  selected = .selected,
  n_sigma = rep(3, length(.selected)),
  lower = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
  upper = c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
  lower_clamp = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 1, 1, -Inf, -Inf, -Inf, 1e-16, -Inf, -Inf),
  upper_clamp = c(0, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf)
)
print(cparams)
saveRDS(cparams, file.path(out_dir, "cparams.Rds"))
print("Parameters saved.")
toc()

# LOADING
tic()
print("Loading zeta results...")
results_directory <- file.path("F:/hguillon/research/exploitation/out/run149/")
zeta_results <- read_zeta_raster(out_dir = results_directory, Lmax = 1E4, raster_resolution = 10, .len = 32, proxy = FALSE)
zeta_results$spatial_scales <- zeta_results$spatial_scales[3:32]
band_names <- read.csv(file.path(results_directory, "band_names.csv")) %>% names()
if ("theta" %in% .selected) zeta_results$raster_list <- zeta_results$raster_list %>% make_angular(match("theta", band_names))
print("Zeta results loaded.")
toc()

# AVERAGING
#' Average layer1 and layer2
#' @param raster_list a `list` of `stars` objects
#' @param xtarget_id `numeric`, optional, the id of the target layer with the x direction
#' @param ytarget_id `numeric`, optional, the id of the target layer with the y direction
#' @return a `list` of `stars`
#' @export
create_average_layer <- function(raster_list, xtarget_id, ytarget_id){
  transformed_raster_list <- lapply(seq_along(raster_list), function(n){
    s <- raster_list[[n]] %>% as("Raster")
    r_target_x <- s[[xtarget_id]]
    r_target_y <- s[[ytarget_id]]
    tmp_s <- raster::stack(r_target_x, r_target_y)
    res <- raster::calc(tmp_s, fun = mean, na.rm = TRUE)
    s <- raster::stack(s, res)
    return(stars::st_as_stars(s))
  })
  return(transformed_raster_list)
}

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

# AVERAGING
if(AVERAGE){
  tic()
  print("Averaging...")
  sliced_clamped_raster <- create_average_layer(sliced_clamped_raster, match("alpha1.x", .selected), match("alpha1.y", .selected))
  .selected <- c(.selected, "alpha1")
  sliced_clamped_raster <- create_average_layer(sliced_clamped_raster, match("alpha2.x", .selected), match("alpha2.y", .selected))
  .selected <- c(.selected, "alpha2")
  sliced_clamped_raster <- create_average_layer(sliced_clamped_raster, match("xi.x", .selected), match("xi.y", .selected))
  .selected <- c(.selected, "xi")
  sliced_clamped_raster <- create_average_layer(sliced_clamped_raster, match("w.x", .selected), match("w.y", .selected))
  .selected <- c(.selected, "w")
  print("Averaging done.")
  toc()
}


# MASK q from H
tic()
print("Masking q rasters from H...")
sliced_clamped_raster <- mask_layer_from_layer(sliced_clamped_raster, match("q", .selected), match("H", .selected), option = "NA")
print("Rasters masked.")
toc()

# MASK ZETA
if(ZETA_MASKING){
  tic()
  print("Masking all rasters zeta1...")
  for (x in .selected[!grepl("zeta1", .selected)]){
    sliced_clamped_raster <- mask_layer_from_layer(sliced_clamped_raster, match(x, .selected), match("zeta1", .selected), option = "NA")
  }
  print("Rasters masked.")
  toc() 
}

# CLIPING CENTRAL VALLEY
tic()
print("Clipping Central Valley...")
shp <- raster::shapefile("F:/hguillon/research/data/safe_water/central_valley_alluvial_boundary/Alluvial_Bnd.shp")
shp <- sp::spTransform(shp, raster::crs(sliced_clamped_raster[[length(sliced_clamped_raster)]] %>% as("Raster")))
sliced_clamped_raster <- lapply(seq_along(sliced_clamped_raster), function(n){
	s <- sliced_clamped_raster[[n]] %>% as("Raster")
	s <- raster::mask(s, shp, inverse = TRUE)
	return(stars::st_as_stars(s, proxy = FALSE))
})  
print("Clipped Central Valley.")
toc()

# LOG TRANSFORM
tic()
print("Applying log-transforms...")
var_to_transforms <- c("median_Pe", "w.x", "w.y", "xi.x", "xi.y", "inv.fc", "w", "xi")
for (x in var_to_transforms[var_to_transforms %in% .selected]){
    sliced_clamped_raster <- logtransform_layer(sliced_clamped_raster, match(x, .selected))
}
print("Applied log-transforms.")
toc()

# THETA
tic()
print("Getting theta values...")
#' Get theta_y from theta_x
#' @param raster_list a `list` of `stars` objects
#' @param target_id the id of the layer to be masked
#' @return a `list` of `stars`
#' @export
get_theta_from_theta <- function(raster_list, theta_x, theta_y){
  transformed_raster_list <- lapply(seq_along(raster_list), function(n){
    s <- raster_list[[n]] %>% as("Raster")
    theta_y_values <- raster::getValues(s[[theta_x]]) + 90
    s <- raster::setValues(s, theta_y_values, layer = theta_y)
    return(stars::st_as_stars(s))
  })
  return(transformed_raster_list)
}

if ("theta.x" %in% .selected) sliced_clamped_raster  <- sliced_clamped_raster  %>% make_angular(match("theta.x", .selected))
sliced_clamped_raster <- get_theta_from_theta(sliced_clamped_raster, match("theta.x", .selected), match("theta.y", .selected))
print("Got theta values.")
toc()


# # REDUCE NOISE BEFORE PAR_PERP STEP
if (REDUCE_NOISE && !PAR_PERP){
  tic()
  print("Addressing noise...")
  sliced_clamped_raster <- reduce_spatial_noise(sliced_clamped_raster, .FUN = NOISE_FUN)
  print("Addressed spatial noise.")
  toc() 
}

# COMPARING X-Y DIRECTIONS
if (COMPARE){
  tic()
  print("Comparing x-y directions...")
  vars <- .selected[grepl(".x", .selected)]
  opts <- sapply(grepl("theta", vars), function(x) ifelse(x, "lt", "gt"))
  res_xy <- foreach(var = vars, opt = opts) %do% {
    is_band1_gtlt_band2(var, gsub(".x", ".y", var), sliced_clamped_raster, .selected, zeta_results$spatial_scales, option = opt)
  }
  print("Compared x-y directions.")
  toc()
}

# PAR PERP
if (PAR_PERP){
  tic()
  print("Deriving parallel and perpendicular directions...")

  res1 <- find_par_perp(sliced_clamped_raster, match("alpha1.x", .selected), match("alpha1.y", .selected))
  res2 <- find_par_perp(sliced_clamped_raster, match("alpha2.x", .selected), match("alpha2.y", .selected))
  sliced_clamped_raster <- res1$raster_list
  .selected <- c(.selected, "alpha_perp1", "alpha_par1")

  sliced_clamped_raster <- find_par_perp(sliced_clamped_raster, match("alpha1.x", .selected), match("alpha1.y", .selected), match("alpha2.x", .selected), match("alpha2.y", .selected), index_perp = res1$index_perp, index_par = res1$index_par)
  .selected <- c(.selected, "alpha_perp2", "alpha_par2")

  sliced_clamped_raster <- find_par_perp(sliced_clamped_raster, match("alpha1.x", .selected), match("alpha1.y", .selected), match("w.x", .selected), match("w.y", .selected), index_perp = res1$index_perp, index_par = res1$index_par)
  .selected <- c(.selected, "w_perp", "w_par")

  sliced_clamped_raster <- find_par_perp(sliced_clamped_raster, match("alpha1.x", .selected), match("alpha1.y", .selected), match("xi.x", .selected), match("xi.y", .selected), index_perp = res1$index_perp, index_par = res1$index_par)
  .selected <- c(.selected, "xi_perp", "xi_par")

  sliced_clamped_raster <- find_par_perp(sliced_clamped_raster, match("alpha1.x", .selected), match("alpha1.y", .selected), match("theta.x", .selected), match("theta.y", .selected), index_perp = res1$index_perp, index_par = res1$index_par)
  .selected <- c(.selected, "theta_perp", "theta_par")

  sliced_clamped_raster <- get_theta_from_theta(sliced_clamped_raster, match("theta_perp", .selected), match("theta_par", .selected))
  sliced_clamped_raster  <- sliced_clamped_raster %>% make_angular(match("theta_par", .selected), mode = "full")


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

  # COMPARING PAR PERP DIRECTIONS
  tic()
  print("Comparing par-perp directions...")
  vars <- .selected[grepl("_perp", .selected)]
  vars <- vars[!grepl("zeta", vars)]
  opts <- sapply(grepl("theta", vars), function(x) ifelse(x, "lt", "gt"))
  res_perp <- foreach(var = vars, opt = opts) %do% {
  is_band1_gtlt_band2(var, gsub("_perp", "_par", var), sliced_clamped_raster, .selected, zeta_results$spatial_scales, option = opt)
  }
  print("Compared par-perp directions.")
  toc()
}

# # REDUCE NOISE AFTER PAR_PERP STEP
if (REDUCE_NOISE && PAR_PERP){
  tic()
  print("Addressing noise...")
  sliced_clamped_raster <- reduce_spatial_noise(sliced_clamped_raster, .FUN = NOISE_FUN)
  print("Addressed spatial noise.")
  toc() 
}

# FOUR VALUES CHECK
tic()
print("Performing four-values check...")
ind <- four_values_check(sliced_clamped_raster)
sliced_clamped_raster <- sliced_clamped_raster[ind]
spatial_scales <- zeta_results$spatial_scales[ind]
print("Performed four-values check.")
toc()

# OUTPUT: DISTRIBUTION
if (DIST){
  tic()
  print("Outputting distributions...")
  pdf(file.path(out_dir, "distributions.pdf"))
  for(x in .selected){
    gb <- make_all_plots(sliced_clamped_raster, spatial_scales, match(x, .selected), x, begin = 0.1, end = 0.85, direction = 1, option = "viridis")
    layout_mat <- matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3), byrow = FALSE, nrow = 5, ncol = 2)
    gridExtra::grid.arrange(grobs = gb, layout_matrix = layout_mat)
  }
  dev.off()
  print("Outputted distributions.")
  toc()
}

# OUTPUT: MAPS
if (MAPS){
  tic()
  print("Outputting maps...")
  for(x in .selected){
    print(x)
    raster_select(sliced_clamped_raster, band_id = match(x, .selected)) %>% 
    save_map(ttl = x, groups = zeta_results$spatial_scales, begin = 0.1, end = 0.95, direction = -1, option = "viridis", col_style = "cont", out_path = out_dir)
  }
  print("Outputted maps.")
  toc()
}

# OUTPUT: PEARSON CORRELATIONS
if (CORR){
  tic()
  print("Outputting correlations...")
  crossscale_corr <- crossscale_correlations(
    raster_list = sliced_clamped_raster,
    selected = .selected,
    spatial_scales = spatial_scales,
    corr_type = "pearson")
  saveRDS(crossscale_corr$graphics_df, file.path(out_dir, "graphics_df.Rds"))
  # pdf(file.path(out_dir, "pearson.pdf"))
  # crossscale_corr$p
  # dev.off()
  print("Outputted correlations.")
  toc()  
}



# SUBSELECTING (to reduce the output from the correlations)
if (SUB){
  # .subselected <- .selected[which(!(grepl("\\.x", .selected) | grepl("\\.y", .selected) | .selected == "zeta1" | .selected == "zeta2"))]
  .subselected <- c("H", "q", "w", "xi", "alpha1", "alpha2", "theta.x", "zeta1")
  tic()
  print("Subselecting rasters...")
  sub_sliced_clamped_raster <- slice_clamp(
    raster_list = sliced_clamped_raster,
    att_names = .selected,
    selected = .subselected,
    clamp_raster = FALSE, 
    clamp_params = NULL)
  print("Raster subselected.")
  toc()
}

# OUTPUT: MAPS
if (SUB && MAPS){
  tic()
  print("Outputting maps...")
  for(x in .subselected){
    print(x)
    raster_select(sub_sliced_clamped_raster, band_id = match(x, .subselected)) %>% 
    save_map(ttl = x, groups = zeta_results$spatial_scales, begin = 0.1, end = 0.95, direction = -1, option = "viridis", col_style = "cont", out_path = out_dir)
  }
  print("Outputted maps.")
  toc()
}

# OUTPUT: DISTRIBUTION
if (SUB && DIST){
  tic()
  print("Outputting distributions...")
  pdf(file.path(out_dir, "distributions_sub.pdf"))1
  for(x in .subselected){
    gb <- make_all_plots(sub_sliced_clamped_raster, spatial_scales, match(x, .subselected), x, begin = 0.1, end = 0.85, direction = 1, option = "viridis")
    layout_mat <- matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3), byrow = FALSE, nrow = 5, ncol = 2)
    gridExtra::grid.arrange(grobs = gb, layout_matrix = layout_mat)
  }
  dev.off()
  print("Outputted distributions.")
  toc()
}

# OUTPUT: PEARSON CORRELATIONS
if (SUB && CORR){
  tic()
  print("Outputting correlations...")
  crossscale_corr <- crossscale_correlations(
    raster_list = sub_sliced_clamped_raster,
    selected = .subselected,
    spatial_scales = spatial_scales,
    corr_type = "pearson")
  pdf(file.path(out_dir, "pearson_sub.pdf"))
  crossscale_corr$p
  dev.off()
  print("Outputted correlations.")
  toc()
}

# # OUTPUT: MAPS
# tic()
# print("Outputting maps...")
# for(x in .subselected){
#   print(x)
#   raster_select(sub_sliced_clamped_raster, band_id = match(x, .subselected)) %>% 
#   save_map(ttl = x, groups = zeta_results$spatial_scales, begin = 0.1, end = 0.95, direction = -1, option = "viridis", col_style = "cont", out_path = out_dir)
# }
# print("Outputted maps.")
# toc()




