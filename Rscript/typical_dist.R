library(tictoc)
library(statisticalRoughness)

results_directory <- file.path("F:/hguillon/research/exploitation/out/run140/")
zeta_results <- read_zeta_raster(out_dir = results_directory, Lmax = 1E4, raster_resolution = 10, .len = 32, proxy = FALSE)
zeta_results$spatial_scales <- zeta_results$spatial_scales[-1]
band_names <- read.csv(file.path(results_directory, "band_names.csv")) %>% names()

.ssl <- sf::read_sf("F:/hguillon/research/data/california-rivers/gis-files/Flowline_CA18_SO/NHDFlowline_SO.shp") %>% sf::st_zm()
len = sapply(seq(7), function(x) .ssl %>% dplyr::filter(StreamOrde >= x) %>% sf::st_length() %>% sum())

rl <- raster_select(zeta_results$raster_list, match("H", band_names))
rl <- rl$rasters

tic()
shp <- sf::read_sf("F:/hguillon/research/data/safe_water/central_valley_alluvial_boundary/Alluvial_Bnd.shp") %>% sf::st_transform(sf::st_crs(3310))
bb <- sf::st_bbox(rl[[31]]) %>% sf::st_as_sfc() %>% sf::st_transform(sf::st_crs(3310))
shp <- sf::st_difference(bb, shp)
toc()

tic()
rl <- lapply(seq_along(rl), function(n){
	s <- rl[[n]] %>% sf::st_transform(sf::st_crs(3310)) %>% sf::st_crop(shp, crop = FALSE)
	return(s)
})
toc()

pct_covering_streamlines = sapply(seq_along(rl), function(i){
	print(i)
	r <- rl[[i]]
	pct_non_na <- sum(!is.na(r[[1]])) / prod(dim(r))
	r <- sf::st_crop(rl[[i]], ssl, crop = FALSE)
	pct_non_na_cropped <- sum(!is.na(r[[1]])) / prod(dim(r))
	ratio <- pct_non_na_cropped / pct_non_na
toc	return(ratio)
})
plot(zeta_results$spatial_scales, pct_covering_streamlines, log = "xy")

ssl <- .ssl %>% 
	dplyr::filter(StreamOrde >= 6) %>% 
	sf::st_geometry() %>% 
	sf::st_transform(sf::st_crs(3310))
shp_CA_dir <- "F:/hguillon/research/data/california-rivers/geomorph_input_data_(Elaheh-White)/DWR_CA_Boundary"
shp_CA <- sf::st_read(file.path(shp_CA_dir, "California_Boundary_TA.shp")) %>% 
	sf::st_transform(sf::st_crs(3310))
CA_area <-  sf::st_intersection(shp, shp_CA)
ssl <- sf::st_union(ssl, CA_area)

# scales <- pracma::logspace(log10(250), 5, n = 20)
scales <- zeta_results$spatial_scales
areas <- rep(NA, length(scales))
for(i in seq_along(scales)){
	tictoc::tic()
	print(i)
	print(scales[i])
	bf <- ssl %>% 
		sf::st_buffer(scales[i], endCapStyle="SQUARE", nQuadSegs=1) %>%
		sf::st_union() %>%
		sf::st_intersection(CAnoCV)
	areas[i] <- bf %>% sf::st_area()
	tictoc::toc()
}

.areas <- areas

Carea <- CAnoCV %>% sf::st_area() %>% as.numeric()
areas <- (Carea - areas) / Carea

df = data.frame(scales = (scales), areas = (areas))
.df = data.frame(scales = log10(scales), areas = (areas))
obj <- lm(areas ~ scales, data = .df)

n <- 1
r2 <- 0
while(signif(r2, 3) < 0.999){
	print(n)
	segmented_fit <- segmented::segmented(obj, npsi = n, control = segmented::seg.control(n.boot = 100))
	s <- summary(segmented_fit)
	r2 <- s$r.squared
	print(r2) 
	n <- n + 1
}

confint <- segmented::confint.segmented(segmented_fit, level = 0.99)

p <- ggplot2::ggplot(df, ggplot2::aes(x = scales, y = areas)) +
ggplot2::geom_point() +
ggplot2::geom_line() +
ggplot2::scale_x_log10(
   breaks = scales::trans_breaks("log10", function(x) 10^x),
   labels = scales::trans_format("log10", scales::math_format(10^.x))
 )
for (i in seq(nrow(confint))){
	p <- p +
		ggplot2::geom_vline(xintercept = 10^(confint[i,2]), lty = 2) +
		ggplot2::geom_vline(xintercept = 10^(confint[i,1]), lty = 1) +
		ggplot2::geom_vline(xintercept = 10^(confint[i,3]), lty = 2)
}
p <- p +
	ggpubr::theme_pubclean() +
	ggplot2::labs(x = "Buffer size (m)", y = "Proportion of non-fluvial landscape (%)")
p