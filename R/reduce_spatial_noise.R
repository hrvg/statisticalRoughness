#' @param raster_list a `list` of `stars` or `RasterStack` objects
#' @return  a `list` of `stars` or `RasterStack` objects
#' @export
reduce_spatial_noise <- function(raster_list){
	raster_list <- lapply(raster_list, function(s){
		s <- as(s, "Raster")
		lr <- lapply(seq(raster::nlayers(s)), function(i) raster::focal(s[[i]], w = matrix(1, nrow = 3, ncol = 3), na.rm = TRUE, pad = TRUE, fun = mean))
		s <- do.call(raster::stack, lr)
		s %>% stars::st_as_stars()
	})
	return(raster_list)
}