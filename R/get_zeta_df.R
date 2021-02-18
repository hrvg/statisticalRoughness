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
	if(!class(DEM) == "RasterLayer") stop("invalid class: DEM is not of class 'RasterLayer'")
	if(!class(tiles) %in% c("RasterLayer", "SpatialPolygonsDataFrame")) stop("invalid class: tiles is not of class 'RasterLayer' or 'SpatialPolygonsDataFrame'")
	if(!class(vertical_accuracy) == 'numeric') stop("invalid class: vertical_accuracy is not of class 'numeric'")
	if(!class(crs_ref) == 'CRS') stop("invalid class: crs_ref is not of class 'CRS'")
	
	# coercion
	if (class(tiles) == "RasterLayer") tiles <- raster::rasterToPolygons(tiles, dissolve = FALSE)
	tiles <- raster::crop(tiles, DEM)
	
	# spatial checks
	if (length(tiles) == 0) stop("There are no tiles covering your raster.")
	if (length(tiles) > raster::ncell(DEM)) stop("There are more tiles than cells in your raster.")

	# setting up parallelization
	registerDoFuture()
	if (length(tiles) < availableCores() %/% 2){
		plan(sequential)
	} else {
		if(.Platform$OS.type == "unix"){
			plan(multicore, workers = availableCores() - 1)
		} else {
			plan(multisession, workers = availableCores() - 1)
		}
	}

	# main
	i <- NULL
	zeta_dfs <-	foreach(i = seq_along(tiles), .combine = rbind, .inorder = TRUE) %dorng% {
		cropped_DEM <- raster::crop(DEM, tiles[i, ])
		cropped_DEM_values <- raster::getValues(cropped_DEM)
		pct_non_na <- sum(is.finite(cropped_DEM_values)) / raster::ncell(cropped_DEM)
		if ((pct_non_na > 1/3) & !(stats::IQR(stats::na.omit(cropped_DEM_values)) < vertical_accuracy)){ 
			mean_res <- mean(raster::res(raster::projectRaster(cropped_DEM, crs = crs_ref)))
			zeta_df <- get_zeta(cropped_DEM, raster_resolution = mean_res)
		} else {
			zeta_df <- matrix(NA, nrow = 1, ncol = 19) %>% as.data.frame()
			colnames(zeta_df) <- c("beta1", "beta2", "alpha1", "alpha1.x", "alpha1.y", "zeta1", "alpha2", "alpha2.x", "alpha2.y", "zeta2", "theta", "inv.fc", "rc", "xi", "xi.x", "xi.y", "w", "w.x", "w.y")
		}
		zeta_df
	}
	# rownames(zeta_dfs) <- NULL
	# stopping parallelization
	plan(sequential)
	return(zeta_dfs)
}
