#' Convolves with a 3x3 mean kernel with weights of 1
#' @param raster_list a `list` of `stars` or `RasterStack` objects
#' @param .NAonly `logical`. If `TRUE`, only cell values that are NA are replaced with the computed focal values
#' @return  a `list` of `stars` or `RasterStack` objects
#' @export
#' @keywords postprocessing
reduce_spatial_noise <- function(raster_list, .NAonly = FALSE){
	raster_list <- lapply(raster_list, function(s){
		s <- as(s, "Raster")
		lr <- lapply(seq(raster::nlayers(s)), function(i) raster::focal(s[[i]], w = matrix(1, nrow = 3, ncol = 3), na.rm = TRUE, pad = TRUE, fun = mean, NAonly = .NAonly))
		s <- do.call(raster::stack, lr)
		s %>% stars::st_as_stars()
	})
	return(raster_list)
}