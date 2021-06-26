#' Get the value of the anisotropy exponent for a tile of elevation
#' @param DEM an elevation raster
#' @param tiles a `RasterLayer` or a `SpatialPolygons` to tiles `DEM`
#' @param raster_resolution `numeric` the resolution of `rstr` in meters
#' @param .zeta_df null or a data.frame with column zeta and a number of rows equals to the number of cells of tiles
#' @param ... to pass the `vertical_accuracy` argument to `get_zeta_df()`
#' @return stack of raster of size tiles with anisotropy exponent values and a number of layers equals to the number of column in .zeta_df
#' @export
#' @keywords zeta
get_zeta_raster <- function(DEM, tiles, raster_resolution, .zeta_df = NULL, ...) {
  # class check
  if (!class(tiles) %in% c("RasterLayer", "SpatialPolygonsDataFrame")) stop("invalid class: tiles is not of class 'RasterLayer' or 'SpatialPolygonsDataFrame'")
  if (is.null(.zeta_df)) {
    .zeta_df <- get_zeta_df(DEM, tiles, raster_resolution, ...)
  }
  zeta_raster <- lapply(seq(ncol(.zeta_df)), function(i) {
    .zeta_raster <- tiles
    vals <- rep(NA, terra::ncell(.zeta_raster))
    vals[which(!is.na(terra::values(.zeta_raster)))] <- unlist(.zeta_df[[i]])
    .zeta_raster <- terra::setValues(.zeta_raster, vals)
  })
  zeta_raster <- do.call(raster::stack, zeta_raster)
  names(zeta_raster) <- colnames(.zeta_df)
  return(zeta_raster)
}
