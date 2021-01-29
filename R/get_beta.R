#' Derive the slope of the binned power spectrum
#' @param binned_power_spectrum `matrix`, a binned spectrum from `bin()`
#' @param FT2D a `list` from `fft2d()`
#' @param do_plot `logical`, plot the fit, default to `FALSE`
#' @return a `data.frame`
#' @export
#' @keywords fft2d
get_beta <- function(binned_power_spectrum, FT2D, do_plot = FALSE){
	x <- binned_power_spectrum[, 1]
	y <- binned_power_spectrum[, 2]
	segmented.fit <- segmented::segmented(lm(y ~ x, weights = abs(x)))
	beta <- data.frame(
		fc = 10^segmented.fit$psi[1,2],
		beta1 = segmented.fit$coefficients[2],
		beta2 = segmented.fit$coefficients[2:3] %>% sum(),
		beta.r2 = summary(segmented.fit)$adj.r.squared
	)
	# if(any(summary(segmented.fit)$Ttable[, 4] %>% stats::na.omit() > 0.05)){
	# 	beta <- data.frame(
	# 		change_point = NA,
	# 		slope1 = NA,
	# 		slope2 = NA,
	# 		adj.r.squared = NA
	# 	)
	# }
	if (do_plot){
		spectrum_plot(binned_power_spectrum, FT2D) + ggplot2::geom_vline(xintercept = beta$fc, lty = 2)
	}
	return(beta)
}