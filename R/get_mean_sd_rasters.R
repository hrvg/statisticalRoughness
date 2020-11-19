#' Compute box-counting mean standard deviation rasters
#' @param dem initial DEM
#' @param L integer
#' @param R vector of values
#' @return list of rasters
#' @importFrom magrittr %>%
#' @importFrom stats sd
#' @export
#' @keywords Hurst
get_mean_sd_rasters <- function(dem, L, R){
	mean_sd_rasters <- lapply(R, function(r) {
		terra::aggregate(dem, fact = r, fun = sd, na.rm = TRUE) %>% 
		terra::aggregate(fact = L / r, fun = mean, na.rm = TRUE)
	})
	return(mean_sd_rasters)
}