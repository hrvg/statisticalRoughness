#' Detrends a DEM with `pracma::odregress()`
#' @param dem a `RasterLayer`
#' @return the detrended raster as a `RasterLayer`
#' @export
#' @keywords fft2d
detrend_dem <- function(dem){
	coord <- raster::coordinates(dem)
	x <- coord[, 1]
	y <- coord[, 2]
	z <- raster::getValues(dem)
	ind <- which(!is.na(z))

	odr <- pracma::odregress(coord[ind, ], z[ind])
	z_detrend <- z
	z_detrend[ind] <- z[ind] - odr$fitted

	detrended_dem <- dem
	detrended_dem <- raster::setValues(detrended_dem, z_detrend)
	detrended_dem[is.na(detrended_dem)] <- 0
	return(detrended_dem)
}

#' Pads a matrix with zeros so it has a given dimension
#' 
#' This function is useful to pad a matrix so that the dimension is a factor of 2.
#' 
#' @param x a `matrix`
#' @param n the desired number of rows of `x`
#' @param m the desired number of columns of `x`
#' @return a padded version of `x` with dimension `n,m`
#' @export
#' @keywords fft2d
pad_mat <- function(x, n, m){
	padx <- n - nrow(x)  
	pady <- m - ncol(x) 
	padx <- padx / 2
	pady <- pady / 2
	x <- matlab::padarray(x, c(floor(padx), floor(pady)), direction = "pre")
	x <- matlab::padarray(x, c(ceiling(padx), ceiling(pady)), direction = "post")
	return(x)
}

#' A two-dimensional Hann window
#' 
#' This code is (largely) a port from Matlab to `R` of the code from [2DSpecTools](https://web.mit.edu/perron/www/files/2DSpecTools-v1.1.zip): A 2D spectral analysis toolkit for Matlab by Taylor Perron et al. 
#' 
#' @param M a `matrix`
#' @return a list with two elements the windowed matrix `H` and the sum of the coefficients `Wss`
#' @export
#' @keywords fft2d
Hann2D <- function(M){
	nx <- ncol(M)
	ny <- nrow(M)
	# matrix coordinates of centroid of M
	a <- (nx + 1) / 2 
	b <- (ny + 1) / 2 
	XY <- pracma::meshgrid((1:nx), (1:ny))
	X <- XY$X
	Y <- XY$Y
	theta <- (X == a) * (pi / 2) + (X != a) * atan2((Y - b), (X - a))
	r <- sqrt((Y - b)^2 + (X - a)^2)
	rprime <- sqrt((a^2) * (b^2) * (b^2 * (cos(theta))^2 + a^2 * (sin(theta))^2)^(-1))
	hanncoeff <- (r < rprime) * (0.5 * (1 + cos(pi * r / rprime)))
	H <- M * hanncoeff
	Wss <- sum(sum(hanncoeff^2))
	return(list(H = H, Wss = Wss))
}

#' Perfoms a two-dimensional Fourier transform
#' 
#' This code is (largely) a port from Matlab to `R` of the code from [2DSpecTools](https://web.mit.edu/perron/www/files/2DSpecTools-v1.1.zip): A 2D spectral analysis toolkit for Matlab by Taylor Perron et al. 
#' As `dx` and `dy` refers to spacing in map units, it is recommended, if applicable, that the original raster is projected from latitude, longitude to a planar projection.
#' 
#' @param M a `matrix`
#' @param dx the spacing in the x direction in map units (usually meters) 
#' @param dy the spacing in the y direction in map units (usually meters)
#' @param Hann `logical`, if `TRUE` performs a Hann windowing, default to `FALSE`
#' @return a list with 4 elements: `spectral_power_matrix`, `radial_frequency_matrix`, `spectral_power_vector` and `radial_frequency_vector`.
#' @export
#' @keywords fft2d
fft2D <- function(M, dx, dy, Hann = TRUE){
	nx <- ncol(M)
	ny <- nrow(M)
	if (Hann){
		l <- Hann2D(M)
		M <- l$H
		Wss <- l$Wss
	} else {
		Wss <- sum(sum(pracma::ones(ny, nx)))
	}
	Lx <- 2^(ceiling(log(max(c(nx, ny)))/log(2)))
	Ly <- Lx
	dfx <-  1 / (dx * Lx)
	dfy <-  1 / (dy * Ly)
	M <- pad_mat(M, Lx, Ly)

	M <- fft(M)
	M <-  mrbsizeR::fftshift(M)
	M[(Ly/2 + 1), (Lx/2 + 1)] <- 0
	M <- M * Conj(M) / (Lx * Ly * Wss)
	spectral_power_matrix <- M 
	xc <- Lx / 2 + 1 # matrix indices of zero frequency
	yc <- Ly / 2 + 1
	XY <- pracma::meshgrid((1:Lx), (1:Ly)) # matrices of column and row indices
	cols <- XY$X
	rows <- XY$Y
	
	radial_frequency_matrix <- sqrt((dfy * (rows - yc))^2 + (dfx*(cols - xc)) ^2)
	# Creating sorted, non-redundant vectors of frequency and power 
	M <- M[ , (1:(Lx/2+1))]
	radial_frequency_vector <- radial_frequency_matrix[, (1:(Lx/2+1))]
	radial_frequency_vector[((yc+1):Ly), xc] <- -1 # This half-column is redundant_ Set the frequencies to negative values so they will be clipped out below
	radial_frequency_vector <- pracma::sortrows(cbind(c(radial_frequency_vector), c(Re(M))), 1) # concatenate column vectors of frequency and PSD and sort by frequency
	radial_frequency_vector <- radial_frequency_vector[radial_frequency_vector[ , 1] > 0 , ] # Take only positive frequencies_ This gets rid of DC (zero freq) as well as the redundant frequencies we set to -1 above
	spectral_power_vector <-  2 * radial_frequency_vector[ , 2] # Separate into power and frequency vectors and assign to output arguments; # the factor of 2 corrects for the fact that we have taken only half of the 2D spectrum_ sum(spectral_power_vectorec) should now equal sum(Pmat(:))_
	radial_frequency_vector <- radial_frequency_vector[, 1]
	l <- list(
		spectral_power_matrix = spectral_power_matrix,
		radial_frequency_matrix = radial_frequency_matrix,
		spectral_power_vector = spectral_power_vector,
		radial_frequency_vector = radial_frequency_vector)
	return(l)
}

#' Plot a radial Fourier spectrum
#' @param binned_power_spectrum `matrix`, a binned spectrum from `bin()`
#' @param FT2D a `list` from `fft2d()`
#' @return a `ggplot` object
#' @export
#' @keywords fft2d
spectrum_plot <- function(binned_power_spectrum, FT2D, xdecades = 3, ydecades = 3){
		df_bin <- data.frame(frequency = 10^binned_power_spectrum[, 1], normalized_power = 10^binned_power_spectrum[, 2])
		df_complete <- data.frame(frequency = FT2D$radial_frequency_vector, normalized_power = FT2D$spectral_power_vector)
		p1 <- ggplot2::ggplot() +
			ggplot2::geom_point(data = df_complete, ggplot2::aes(x = frequency , y = normalized_power, alpha = 0.1)) +
			ggplot2::geom_point(data = df_bin, ggplot2::aes(x = frequency , y = normalized_power, color = "red")) +
			ggplot2::labs(title = "Topography power spectrum") +
			ggplot2::labs(y= bquote("DFT mean squared amplitude "~(m^2)), x = bquote("Wavenumber "~(m^{-1}))) +
			ggplot2::scale_x_log10(
				breaks = scales::trans_breaks(n = xdecades, 'log10', function(x) 10^x),
	            labels = scales::trans_format('log10', scales::math_format(10^.x)),
	            sec.axis = ggplot2::sec_axis(~1/., 
	            	breaks = scales::trans_breaks(n = xdecades, 'log10', function(x) 10^x), 
	            	labels = scales::trans_format('log10', scales::math_format(10^.x)),
	            	name = "Wavelength (m)")
	            	) +
			ggplot2::scale_y_log10(
				breaks = scales::trans_breaks(n = ydecades, 'log10', function(x) 10^x),
	            labels = scales::trans_format('log10', scales::math_format(10^.x))
	            ) + 
			ggpubr::theme_pubr() +
			ggplot2::guides(color = FALSE, alpha = FALSE)
	    return(p1)
}