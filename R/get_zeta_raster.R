#' Get the value of the anisotropy exponent for a tile of elevation
#' @param DEM an elevation raster
#' @param tiles a `RasterLayer` or a `SpatialPolygons` to tiles `DEM`
#' @param .zeta_df null or a data.frame with column zeta and a number of rows equals to the number of cells of tiles
#' @param ... to pass the `crs_ref` and `vertical_accuracy` argument to `get_zeta_df()`
#' @return stack of raster of size tiles with anisotropy exponent values and a number of layers equals to the number of column in .zeta_df
#' @export
get_zeta_raster <- function(DEM, tiles, .zeta_df = NULL, ...){
	# class check
	if(!class(tiles) %in% c("SpatRaster")) stop("invalid class: tiles is not of class 'SpatRaster'")
	if (is.null(.zeta_df)){
		.zeta_df <- get_zeta_df(DEM, tiles, ...)
	}
	zeta_raster <- lapply(seq(ncol(.zeta_df)), function(i){
		.zeta_raster <- tiles
		vals <- rep(NA, terra::ncell(.zeta_raster))
		vals[which(!is.na(terra::values(.zeta_raster)))] <- unlist(.zeta_df[[i]])
		.zeta_raster <- terra::setValues(.zeta_raster, vals)
	})
	zeta_raster <- do.call(c, zeta_raster)
	names(zeta_raster) <- colnames(.zeta_df)
	return(zeta_raster)
}
