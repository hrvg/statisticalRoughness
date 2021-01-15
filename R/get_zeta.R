#' Internal function to derive the anisotropy exponent from two roughness exponents
#' @param alpha_1 `numeric`, a roughness coefficient
#' @param alpha_2 `numeric`, a roughness coefficient
#' @return a `numeric` anisotropy exponent
#' @keywords zeta
#' @export
get_zeta_ <- function(alpha_1, alpha_2){
	zeta <- max(c(alpha_1, alpha_2), na.rm = TRUE) / min(c(alpha_1, alpha_2), na.rm = TRUE)
	return(zeta)
}

#' Wrapper function to derive the anisotropy exponent from a non-rotated raster
#' @param rstr `Raster` object
#' @param raster_resolution `numeric` the resolution of `rstr` in meters
#' @param nbin `numeric`, the number of bins
#' @param .Hann `logical`, if `TRUE` performs a Hann windowing, default to `TRUE`, passed to `fft2D()`
#' @param .quantile_prob vector of quantile probabilities, default to `c(0.9999)`, passed to `filter_spectral_power_matrix()`
#' @param .prob, `numeric` the value of the probabilities passed to `quantile()`, passed to `filter_alpha()`
#' @importFrom magrittr %>%
#' @return a `data.frame`
#' @export
#' @keywords zeta
get_zeta <- function(rstr, raster_resolution = 9.015, nbin = 50, .Hann = TRUE, .quantile_prob = c(0.9999), .prob = .999){
	rstr <- rstr %>% detrend_dem() 
	FT2D <- fft2D(raster::as.matrix(rstr), dx = raster_resolution, dy = raster_resolution, Hann = .Hann)
	binned_power_spectrum <- bin(log10(FT2D$radial_frequency_vector), log10(FT2D$spectral_power_vector), nbin) %>% na.omit()
	normalized_spectral_power_matrix <- get_normalized_spectral_power_matrix(binned_power_spectrum, FT2D)
	filtered_spectral_power_matrix <- filter_spectral_power_matrix(normalized_spectral_power_matrix, FT2D, quantile_prob = .quantile_prob)
	ang_fourier <- get_fourier_angle(filtered_spectral_power_matrix, FT2D)
	rotated_raster <- rotate_raster(raster::as.matrix(rstr), ang_fourier)
	hhcf_x <- get_hhcf(rotated_raster, margin = 1)
	hhcf_y <- get_hhcf(rotated_raster, margin = 2)
	alpha_x <- get_all_alpha(hhcf_x, raster_resolution) %>% summarise_alpha()
	colnames(alpha_x) <- paste0(colnames(alpha_x), ".x")
	alpha_y <- get_all_alpha(hhcf_y, raster_resolution) %>% summarise_alpha()
	colnames(alpha_y) <- paste0(colnames(alpha_y), ".y")
	res <- dplyr::bind_cols(alpha_x, alpha_y) %>% dplyr::mutate(
		zeta1 = get_zeta_(alpha_x$slope1_mean, alpha_y$slope1_mean),
		zeta2 = get_zeta_(alpha_x$slope2_mean, alpha_y$slope2_mean),
		theta = ang_fourier
		)
	return(res)
}