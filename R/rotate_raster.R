#' Bins a power spectrum
#' 
#' This code is (largely) a port from Matlab to `R` of the code from [2DSpecTools](https://web.mit.edu/perron/www/files/2DSpecTools-v1.1.zip): A 2D spectral analysis toolkit for Matlab by Taylor Perron et al. 
#' 
#' @param x `numeric`
#' @param y `numeric`
#' @param nbin `numeric`, the number of bins
#' @return a matrix with binned values and summary statistics
#' @export
#' @keywords rotate_raster
bin <- function(x, y, nbin){
	sorted <- pracma::sortrows(cbind(c(x), c(y)), 1)
	x <- sorted[, 1]
	y <- sorted[, 2]
	xmin <- head(x ,1)
	xmax <- tail(x, 1)
	xrange <- xmax - xmin
	w <- xrange / nbin
	s <- pracma::zeros(nbin, 8)
	for (i in (1:nbin)){
		xlo <- xmin + (i-1) * w
		xhi <- xlo + w
		x.window <- which((x >= xlo) & (x <= xhi))
		if (length(x.window) > 1){
			mini <- min(x.window)
			maxi <- max(x.window)
			s[i, ] <- c(mean(xlo, xhi), 
				mean(y[mini:maxi]),
				sd(y[mini:maxi]),
				sd(y[mini:maxi])/sqrt(maxi-mini+1),
				maxi-mini+1,
				max(y[mini:maxi]),
				min(y[mini:maxi]),
				median(y[mini:maxi])
				)
		} else {
			s[i, ] <- c(mean(xlo, xhi), NA, NA, NA, NA, NA, NA, NA)
		}
	}
	return(s)
}

#' Returns a quantile-filtered raster of the logarithm of the real-part of the normalized spectral power matrix
#' @param normalized_spectral_power_matrix `matrix`, normalized spectral power matrix
#' @param quantile_prob vector of quantile probabilities, default to `c(0.99)`
#' @return a `RasterLayer` of the normalized spectral power matrix filtered for components above `quantile_prob`
#' @export
#' @keywords rotate_raster
filter_spectral_power_matrix <- function(normalized_spectral_power_matrix, FT2D, quantile_prob = c(0.99)){
	log_Pmn <- log10(Re(normalized_spectral_power_matrix))
	log_Pmn[!is.finite(log_Pmn)] <- NA

	nfx <- ncol(FT2D$radial_frequency_matrix)
	nfy <- nrow(FT2D$radial_frequency_matrix)
	nyq <- FT2D$radial_frequency_matrix[(nfy / 2 + 1), 1]
	x <- pracma::linspace(-nyq, nyq, nfx)
	y <- pracma::linspace(-nyq, nyq, nfy)

	ras <- raster::raster(nrow = nfx, ncol = nfy)
	raster::extent(ras) <- c(min(x), max(x), min(y), max(y))
	ras_Pmn <- raster::setValues(ras, c(t(log_Pmn)))

	q_log_Pmn <- quantile(log_Pmn, probs = c(0.99), na.rm = TRUE)
	ras_Pmn_clamped <- raster::clamp(ras_Pmn, lower = q_log_Pmn[1], useValues = FALSE)
	return(ras_Pmn_clamped)
}

#' Normalizes the spectral power matrix with the background spectrum
#' @param B `matrix` from `bin()`, binned radial frequency and spectral power vectors
#' @param FT2D `list`, results from `fft2d()` with objects radial frequency and spectral power matrices
#' @return `matrix`, normalized spectral power matrix
#' @export
#' @keywords rotate_raster
get_normalized_spectral_power_matrix <- function(B, FT2D){
	rlm_fit <- MASS::rlm(B[, 2] ~ B[, 1], maxit = 100)
	background_sprectrum <- (10^rlm_fit$coefficients[1] * FT2D$radial_frequency_matrix^rlm_fit$coefficients[2])
	normalized_spectral_power_matrix <- FT2D$spectral_power_matrix / background_sprectrum
	return(normalized_spectral_power_matrix)
}


#' Get the angle of the 2D Fourier spectrum
#' @param Pmn_truncated truncated normalized power spectrum
#' @return angle of the main component of the Fourier spectrum
#' @export
#' @keywords rotate_raster
get_fourier_angle <- function(Pmn_truncated){
	points <- raster::rasterToPoints(Pmn_truncated, spatial = TRUE)
	coord <- raster::coordinates(points)
	coord_polar <- useful::cart2pol(coord[, 1], coord[, 2], degree = TRUE)

	power_ww <- na.omit(raster::getValues(Pmn_truncated))
	power_ww <- log10(power_ww)^2
	power_ww <- power_ww / sum(power_ww)
	cDens_fourier <- spatstat::circdensity(coord_polar$theta, weights = power_ww, bw = 5)
	ang_fourier <- cDens_fourier$x[which.max(cDens_fourier$y)]
	return(ang_fourier)
}

#' Rotates the detrended DEM according to the main direction of the Fourier spectrum
#' @param mm `matrix`, matrix corresponding to a DEM
#' @param ang_fourier `numeric`, angle of the main component of the 2D Fourier spectrum
#' @return a `matrix` corresponding to the DEM rotated in the main direction of the Fourier spectrum
#' @export
#' @keywords rotate_raster
rotate_raster <- function(mm, ang_fourier){
	mm.max <- max(mm, na.rm = TRUE)
	mm.min <- .9 * min(mm, na.rm = TRUE)
	mm <- (mm - mm.min) / (mm.max - mm.min)
	im <- imager::as.cimg(mm)
	im.r <- imager::imrotate(im, -ang_fourier, interpolation = 0)
	mm <- as.matrix(im.r)
	dim(mm) <- NULL
	mm[mm == 0] <- NA
	mm <- mm * (mm.max - mm.min) + mm.min
	mm <- matrix(mm, nrow = dim(im.r)[1], ncol = dim(im.r)[2])
	return(mm)
}