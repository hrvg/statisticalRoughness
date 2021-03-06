#' Internal function to derive the anisotropy exponent from two roughness exponents
#' 
#' If any IQR is higher than 0.225, this function returns `NA`. The threshold value of 0.225 is the IQR of a normal distribution with a standard deviation equals to 0.5/3, a three-sigma range between 0 and 1.
#' 
#' @param alpha_1 `numeric`, a roughness coefficient
#' @param alpha_2 `numeric`, a roughness coefficient
#' @param IQR_1 `numeric`, an IQR
#' @param IQR_2 `numeric`, an IQR
#' @param kruskal_flag `logical` flag returned by `get_kruskal_flag()`
#' @return a `numeric` anisotropy exponent
#' @keywords zeta
#' @export
get_zeta_ <- function(alpha_1, alpha_2, IQR_1, IQR_2, kruskal_flag){
	if(any(is.na(c(IQR_1, IQR_2)))) return(NA)
	if(IQR_1 <= 0.225 & IQR_2 <= 0.225 & kruskal_flag){
		# zeta <- max(c(alpha_1, alpha_2), na.rm = TRUE) / min(c(alpha_1, alpha_2), na.rm = TRUE)
		zeta <- max(c(alpha_1), na.rm = TRUE) / min(c(alpha_2), na.rm = TRUE)
	} else {
		zeta <- NA
	}
	return(zeta)
}

#' Computes a Kruskal-Wallis ranksum test for significant differences between the distributions of alpha
#' @param alpha_x `data.frame` from get_all_alpha_()`
#' @param alpha_y `data.frame` from get_all_alpha_()`
#' @param var `character` the variable to compare
#' @importFrom rlang .data
#' @return a `logical` value: `TRUE` if the distribution are statistically different, `FALSE` otherwise
#' @keywords zeta
#' @export
get_kruskal_flag <- function(alpha_x, alpha_y, var){
	alpha_x <- alpha_x %>% 
		dplyr::select(var) %>% 
		dplyr::mutate(type = "x") %>% 
		stats::na.omit()
	alpha_y <- alpha_y %>%
		dplyr::select(var) %>%
		dplyr::mutate(type = "y") %>% 
		stats::na.omit()
	if(nrow(alpha_x) > 1 & nrow(alpha_y) > 1){
		kruskal <- rbind(alpha_x, alpha_y) %>% dplyr::mutate(type = as.factor(.data$type))
		colnames(kruskal) <- c("var", "type")
		kruskal <- kruskal %>% rstatix::kruskal_test(var ~ type) 
		return(kruskal$p < 0.05)
	} else {
		return(FALSE)		
	}
}

#' Wrapper function to derive the anisotropy exponent from a non-rotated raster
#' @param rstr `Raster` object
#' @param raster_resolution `numeric` the resolution of `rstr` in meters
#' @param nbin `numeric`, the number of bins
#' @param .Hann `logical`, if `TRUE` performs a Hann windowing, default to `TRUE`, passed to `fft2D()`
#' @param .quantile_prob vector of quantile probabilities, default to `c(0.9999)`, passed to `filter_spectral_power_matrix()`
#' @param .prob, `numeric` the value of the probabilities passed to `quantile()`, passed to `filter_alpha()`
#' @param full, `logical` if `TRUE` all results are returned, else, the default, results of interest are returned
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return a `data.frame`
#' @export
#' @keywords zeta
get_zeta <- function(rstr, raster_resolution, nbin = 20, .Hann = TRUE, .quantile_prob = c(0.9999), .prob = .999, full = FALSE){
	FT2D <- rstr %>% detrend_dem() %>% raster::as.matrix() %>% 
		fft2D(dx = raster_resolution, dy = raster_resolution, Hann = .Hann)
	binned_power_spectrum <- bin(log10(FT2D$radial_frequency_vector), log10(FT2D$spectral_power_vector), nbin) %>% stats::na.omit()
	beta <- get_beta(binned_power_spectrum, FT2D)
	normalized_spectral_power_matrix <- get_normalized_spectral_power_matrix(binned_power_spectrum, FT2D)
	filtered_spectral_power_matrix <- filter_spectral_power_matrix(normalized_spectral_power_matrix, FT2D, quantile_prob = .quantile_prob)
	ang_fourier <- get_fourier_angle(filtered_spectral_power_matrix, FT2D)
	rotated_raster <- rotate_raster(rstr, ang_fourier)
	hhcf_x <- get_hhcf_(rotated_raster, raster_resolution, margin = 1)
	hhcf_y <- get_hhcf_(rotated_raster, raster_resolution, margin = 2)
	w_x <- mean(hhcf_x$rms, na.rm = TRUE)
	w_y <- mean(hhcf_y$rms, na.rm = TRUE)
	rms_flag <- w_x > w_y
	if(is.na(rms_flag)) rms_flag <- TRUE
	if(rms_flag){
		alpha_x <- get_all_alpha_(hhcf_x, raster_resolution) 
		alpha_y <- get_all_alpha_(hhcf_y, raster_resolution)
		xi_x <- mean(hhcf_x$autocorr_len, na.rm = TRUE)
		xi_y <- mean(hhcf_y$autocorr_len, na.rm = TRUE)
	} else { # invert so that x has the highest rms
		alpha_x <- get_all_alpha_(hhcf_y, raster_resolution)
		alpha_y <- get_all_alpha_(hhcf_x, raster_resolution)
		xi_x <- mean(hhcf_y$autocorr_len, na.rm = TRUE)
		xi_y <- mean(hhcf_x$autocorr_len, na.rm = TRUE)
		w_x <- mean(hhcf_y$rms, na.rm = TRUE)
		w_y <- mean(hhcf_x$rms, na.rm = TRUE)
	}
	kruskal_flag1 <- get_kruskal_flag(alpha_x, alpha_y, "alpha1")
	kruskal_flag2 <- get_kruskal_flag(alpha_x, alpha_y, "alpha2")
	alpha_x <- alpha_x %>% summarise_alpha()
	alpha_y <- alpha_y %>% summarise_alpha()
	colnames(alpha_x) <- paste0(colnames(alpha_x), ".x")
	colnames(alpha_y) <- paste0(colnames(alpha_y), ".y")
	res <- dplyr::bind_cols(beta, alpha_x, alpha_y) %>% 
		dplyr::mutate(
			zeta1 = get_zeta_(alpha_x$alpha1_mean.x, alpha_y$alpha1_mean.y, alpha_x$alpha1_IQR.x, alpha_y$alpha1_IQR.y, kruskal_flag1),
			zeta2 = get_zeta_(alpha_x$alpha2_mean.x, alpha_y$alpha2_mean.y, alpha_x$alpha2_IQR.x, alpha_y$alpha2_IQR.y, kruskal_flag2),
			theta = ang_fourier,
			rc = mean(c(alpha_x$rc_mean.x, alpha_y$rc_mean.y)),
			w.x = w_x,
			w.y = w_y,
			w = mean(c(.data$w.x, .data$w.y)),
			xi.x = xi_x,
			xi.y = xi_y,
			xi = mean(c(.data$xi.x, .data$xi.y))
		)
	if (full) return(res)
	res <- res %>%
		dplyr::rename(
			alpha1.x = .data$alpha1_mean.x,
			alpha2.x = .data$alpha2_mean.x,
			alpha1.y = .data$alpha1_mean.y,
			alpha2.y = .data$alpha2_mean.y
		) %>%
		dplyr::mutate(
			alpha1 = mean(c(.data$alpha1.x, .data$alpha1.y)),
			alpha2 = mean(c(.data$alpha2.x, .data$alpha2.y)),
			inv.fc = 1/.data$fc
		) %>%
		dplyr::select(
			.data$beta1,
			.data$beta2,
			.data$alpha1,
			.data$alpha1.x,
			.data$alpha1.y,
			.data$zeta1,
			.data$alpha2,
			.data$alpha2.x,
			.data$alpha2.y,
			.data$zeta2,
			.data$theta,
			.data$inv.fc,
			.data$rc,
			.data$xi,
			.data$xi.x,
			.data$xi.y,
			.data$w,
			.data$w.x,
			.data$w.y
		)
	return(res)
}