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
	# if(IQR_1 <= 0.225 & IQR_2 <= 0.225 & kruskal_flag){
	if(kruskal_flag){
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
get_zeta <- function(rstr, raster_resolution, .mode = "radial", angle_step = 5, nbin = 20, .Hann = TRUE, .quantile_prob = c(0.99), .prob = .999, full = FALSE){
	if(.mode == "fourier"){
		res <- get_zeta_fourier(rstr, raster_resolution, nbin, .Hann, .quantile_prob, .prob, full)
	} else if(.mode == "radial"){
		res <- get_zeta_radial(rstr, raster_resolution, angle_step, full)
	}
	return(res)
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
get_zeta_fourier <- function(rstr, raster_resolution, nbin, .Hann, .quantile_prob, .prob, full){
	FT2D <- rstr %>% raster::trim() %>% detrend_dem() %>% raster::as.matrix() %>% 
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
	theta_x <- ang_fourier %% 360
	theta_y <- (theta_x + 90) %% 360
	alpha_x <- get_all_alpha_(hhcf_x, raster_resolution) 
	alpha_y <- get_all_alpha_(hhcf_y, raster_resolution)
	xi_x <- mean(hhcf_x$autocorr_len, na.rm = TRUE)
	xi_y <- mean(hhcf_y$autocorr_len, na.rm = TRUE)
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
			xi = mean(c(.data$xi.x, .data$xi.y)),
			theta.x = theta_x,
			theta.y = theta_y
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
			.data$theta.x,
			.data$theta.y,
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
get_zeta_radial <- function(rstr, raster_resolution, angle_step, full){
	ext <- raster::extent(rstr)
	circle <- sf::st_sfc(sf::st_buffer(sf::st_point(c(ext[1]+(ext[2]-ext[1])/2, ext[3]+(ext[4]-ext[3])/2)), mean(dim(rstr)[1:2]) * mean(raster::res(rstr)) / 2), crs = sf::st_crs(rstr)) %>% as("Spatial")
	rstr <- raster::mask(rstr %>% raster::trim(), circle)
	radial_res <- get_radial_angle(rstr, raster_resolution, angle_step)
	rotation_angle <- radial_res$theta_perp
	rotated_raster <- rotate_raster(rstr, rotation_angle)
	hhcf_x <- get_hhcf_(rotated_raster, raster_resolution, margin = 1)
	hhcf_y <- get_hhcf_(rotated_raster, raster_resolution, margin = 2)
	w_x <- mean(hhcf_x$rms, na.rm = TRUE)
	w_y <- mean(hhcf_y$rms, na.rm = TRUE)
	theta_x <- rotation_angle %% 360
	theta_y <- (theta_x + 90) %% 360
	alpha_x <- get_all_alpha_(hhcf_x, raster_resolution) 
	alpha_y <- get_all_alpha_(hhcf_y, raster_resolution)
	xi_x <- mean(hhcf_x$autocorr_len, na.rm = TRUE)
	xi_y <- mean(hhcf_y$autocorr_len, na.rm = TRUE)
	kruskal_flag1 <- get_kruskal_flag(alpha_x, alpha_y, "alpha1")
	kruskal_flag2 <- get_kruskal_flag(alpha_x, alpha_y, "alpha2")
	alpha_x <- alpha_x %>% summarise_alpha()
	alpha_y <- alpha_y %>% summarise_alpha()
	colnames(alpha_x) <- paste0(colnames(alpha_x), ".x")
	colnames(alpha_y) <- paste0(colnames(alpha_y), ".y")
	res <- dplyr::bind_cols(alpha_x, alpha_y) %>% 
		dplyr::mutate(
			zeta1 = get_zeta_(alpha_x$alpha1_mean.x, alpha_y$alpha1_mean.y, alpha_x$alpha1_IQR.x, alpha_y$alpha1_IQR.y, kruskal_flag1),
			zeta2 = get_zeta_(alpha_x$alpha2_mean.x, alpha_y$alpha2_mean.y, alpha_x$alpha2_IQR.x, alpha_y$alpha2_IQR.y, kruskal_flag2),
			theta = rotation_angle,
			rc = mean(c(alpha_x$rc_mean.x, alpha_y$rc_mean.y)),
			w.x = w_x,
			w.y = w_y,
			w = mean(c(.data$w.x, .data$w.y)),
			xi.x = xi_x,
			xi.y = xi_y,
			xi = mean(c(.data$xi.x, .data$xi.y)),
			theta.x = theta_x,
			theta.y = theta_y
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
			alpha1_median = radial_res$alpha1_median,
			alpha1_mad = radial_res$alpha1_mad
		) %>%
		dplyr::select(
			.data$alpha1,
			.data$alpha1.x,
			.data$alpha1.y,
			.data$zeta1,
			.data$alpha2,
			.data$alpha2.x,
			.data$alpha2.y,
			.data$zeta2,
			.data$theta.x,
			.data$theta.y,
			.data$alpha1_median,
			.data$alpha1_mad,
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

#' Wrapper function to derive the anisotropy exponent from a non-rotated raster
#' @param rstr `Raster` object
#' @param raster_resolution `numeric` the resolution of `rstr` in meters
#' @param angle_step `numeric`, angular step
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @return a `data.frame`
#' @export
#' @keywords zeta
get_radial_angle <- function(rstr, raster_resolution, angle_step){
	angles <- seq(0, 90 - angle_step, angle_step)
	res <- foreach(rotation_angle = angles, .combine = cbind) %do% {
		rotated_raster <- rotate_raster(rstr, rotation_angle)
		hhcf_x <- get_hhcf_(rotated_raster, raster_resolution, margin = 1, average = TRUE)
		hhcf_y <- get_hhcf_(rotated_raster, raster_resolution, margin = 2, average = TRUE)
		theta_x <- rotation_angle %% 360
		theta_y <- (theta_x + 90) %% 360
		alpha_x <- get_all_alpha_(hhcf_x, raster_resolution) 
		alpha_y <- get_all_alpha_(hhcf_y, raster_resolution)
		colnames(alpha_x) <- paste0(colnames(alpha_x), ".x")
		colnames(alpha_y) <- paste0(colnames(alpha_y), ".y")
		res <- dplyr::bind_cols(alpha_x, alpha_y)
		res <- res %>%
			dplyr::select(
				.data$alpha1.x,
				.data$alpha1.y
			)
		colnames(res) <- gsub(".x", paste0(".", theta_x), colnames(res))
		colnames(res) <- gsub(".y", paste0(".", theta_y), colnames(res))
		return(res)
	}
	alpha1 <- res %>% dplyr::select(dplyr::contains("alpha1.")) %>% as.matrix() %>% c()
	alpha1_median <- stats::median(alpha1, na.rm = TRUE)
	alpha1_mad <- stats::mad(alpha1, constant = 1, na.rm = TRUE)
	theta <- res %>% dplyr::select(dplyr::contains("alpha1.")) %>% colnames() %>% gsub("alpha1.", "", .) %>% as.numeric()
	theta <- theta[!is.na(alpha1)]
	alpha1 <- alpha1[!is.na(alpha1)]
	alpha1 <- c(alpha1, alpha1)
	theta <- c(theta, theta %>% +(180))
	ww <- scales::rescale(alpha1)
	ww <- ww / sum(ww)
	circular_density <- spatstat.core::circdensity(theta, weights = ww, bw = angle_step)
	theta_perp <- circular_density$x[which.max(circular_density$y)] %% 180
	return(list(alpha1_median = alpha1_median, alpha1_mad = alpha1_mad, theta_perp = theta_perp))
}