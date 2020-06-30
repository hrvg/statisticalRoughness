#' Compute the mean resolution of lat long raster in meters
#' @param r aggegated raster
#' @param L aggregation factor
#' @return numeric, mean resolution in m
#' @export
get_mean_res <- function(r, L, crs.ref = crs_ref){
	.r <- raster::projectRaster(raster::raster(r), crs = crs.ref)
	mean_res <- mean(raster::res(.r)) / L
	return(mean_res)
}
