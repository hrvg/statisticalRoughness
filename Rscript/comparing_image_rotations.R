devtools::load_all()

rstr <- raster::raster(file.path(system.file("extdata/rasters/", package = "statisticalRoughness"), "gabilan_mesa.tif"))
rstr <- rstr %>% raster::trim()
ext <- raster::extent(rstr)
circle <- sf::st_sfc(sf::st_buffer(sf::st_point(c(ext[1]+(ext[2]-ext[1])/2, ext[3]+(ext[4]-ext[3])/2)), mean(dim(rstr)[1:2]) * mean(raster::res(rstr)) / 2), crs = sf::st_crs(rstr)) %>% as("Spatial")
rstr <- raster::mask(rstr, circle)

rotate_raster_imager <- function(rstr, ang_fourier){
	rstr <- rstr %>% 
		raster::calc(fun = function(x) x + 1e4) %>% 
		imager::as.cimg() %>% 
		imager::imrotate(ang_fourier, interpolation = 0) %>%  # default to clockwise rotation
		as.matrix() %>% 
		t() %>% 
		raster::raster() %>% 
		raster::reclassify(matrix(c(-Inf, 0, NA))) %>%
		raster::calc(fun = function(x) x - 1e4) %>% 
		raster::as.matrix()
	return(rstr)
}

rotate_raster_openimager <- function(rstr, ang_fourier){
	ang_fourier <- - ang_fourier %% 360
	rstr <- rstr %>%
		raster::as.matrix() %>%
		 OpenImageR::rotateImage(ang_fourier, method = "bilinear")  # default to counter-clockwise rotation
	rstr[rstr == 0] <- NA	
	return(rstr)
}

ang_fourier <- - 45

p1_imager <- rstr %>% rotate_raster_imager(ang_fourier) %>% view_matrix(ply = FALSE)
p1_openimager <- rstr %>% rotate_raster_openimager(ang_fourier) %>% view_matrix(ply = FALSE)

gridExtra::grid.arrange(
	grobs = list(
		rstr %>% raster::as.matrix() %>% view_matrix(ply = FALSE),
		p1_imager,
		p1_openimager
	),
	ncol = 3
)

library(microbenchmark)
mbm = microbenchmark(
  imager = rstr %>% rotate_raster_imager(ang_fourier),
  openimager = rstr %>% rotate_raster_openimager(ang_fourier),
  times = 100
)
mbm
ggplot2::autoplot(mbm)
