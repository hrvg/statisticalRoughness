devtools::load_all()

ANGLE_STEP <- 5
NITERS <- 2^seq(7)
N_SAMPLE <- 30
raster_resolution <- 10

gabilan_mesa <- raster::raster(file.path(system.file("extdata/rasters/", package = "statisticalRoughness"), "gabilan_mesa.tif"))

best_iter_results <- foreach(NITER = NITERS, .combine = rbind) %do% {
	res <- foreach(k = seq(N_SAMPLE), .combine = rbind) %do% {
		get_zeta(gabilan_mesa, raster_resolution, .mode = "radial", angle_step = ANGLE_STEP, niter = NITER)		
	} %>%
	dplyr::mutate(NITER = NITER)
}

best_iter_results <- best_iter_results %>% dplyr::mutate(NITER = as.factor(NITER))

p <- ggplot2::ggplot(best_iter_results, ggplot2::aes(x = theta.x, group = NITER, color = NITER)) +
ggplot2::geom_density() +
ggpubr::theme_classic2()
plotly::ggplotly(p)