#' Get the value of the anisotropy exponent and associated metrics for a tile of elevation
#' @param DEM an elevation raster
#' @param tiles a `RasterLayer` or a `SpatialPolygons` to tiles `DEM`
#' @param crs_ref a `raster::crs`, a coordinate reference system, default to `terra::crs("EPSG:3310")`
#' @param vertical_accuracy `numeric`, the vertical accuracy of `DEM`, default to 1.87 for the 10-m National Elevation Dataset
#' @return a data.frame with a number of rows equals to the length of `tiles`
#' @import doFuture
#' @import doRNG
#' @import future
#' @import foreach
#' @export
get_zeta_df <- function(DEM, tiles, crs_ref = terra::crs("EPSG:3310"), vertical_accuracy = 1.87){
	# class check
	if(!class(DEM) == "SpatRaster") stop("invalid class: DEM is not of class 'SpatRaster'")
	if(!class(tiles) %in% c("SpatRaster", "SpatVector")) stop("invalid class: tiles is not of class 'SpatRaster' or 'SpatVector'")
	if(!class(vertical_accuracy) == 'numeric') stop("invalid class: vertical_accuracy is not of class 'numeric'")
	if(!class(crs_ref) == 'CRS') stop("invalid class: crs_ref is not of class 'CRS'")
	
	# coercion
	if (class(tiles) == "SpatRaster") tiles <- terra::as.polygons(tiles, dissolve = FALSE)
	tiles <- terra::crop(tiles, DEM)
	
	# spatial checks
	if (length(tiles) == 0) stop("There are no tiles covering your raster.")
	if (length(tiles) > terra::ncell(DEM)) stop("There are more tiles than cells in your raster.")

	# setting up parallelization
	registerDoFuture()
	if (length(tiles) < availableCores() %/% 2){
		plan(sequential)
	} else {
		plan(multicore, workers = availableCores() %/% 2)
	}

	# main
	i <- NULL
	zeta_dfs <-	foreach(i = seq_along(tiles), .combine = rbind, .inorder = TRUE) %dorng% {
		cropped_DEM <- terra::crop(DEM, tiles[i, ])
		cropped_DEM_values <- terra::values(cropped_DEM)
		pct_non_na <- sum(is.finite(cropped_DEM_values)) / terra::ncell(cropped_DEM)
		if ((pct_non_na > 1/3) & !(stats::IQR(stats::na.omit(cropped_DEM_values)) < vertical_accuracy)){ 
			mean_res <- mean(terra::res(terra::project(cropped_DEM, as.character(crs_ref))))
			cropped_DEM <- raster::raster(cropped_DEM)
			zeta_df <- get_zeta(cropped_DEM, raster_resolution = mean_res)
		}
		zeta_df
	}
	rownames(zeta_dfs) <- NULL
	# stopping parallelization
	plan(sequential)
	return(zeta_dfs)
}