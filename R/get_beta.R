#' Derive the slope of the binned power spectrum
#' @param 
#' 
#' 
get_beta <- function(binned_power_spectrum, change_point, max_point){
	change_frequency <- 1 / log10(change_point)
	maxfrequency <- 1 / log10(max_point)
	df_bin <- data.frame(frequency = binned_power_spectrum[, 1], normalized_power = binned_power_spectrum[, 2])
	pw_fun1 <- splinefun(x = df_bin$frequency, y = df_bin$normalized_power)
	d_P_df <- pw_fun1(df_bin$frequency, deriv = 1)

	segmented.fit <- 	
}

#' Internal companion function to `get_alpha`.
#' @param df a `data.frame` with the non-binned data
#' @param df_bin a `data.frame` with the binned data
#' @param change_point `numeric`, the change point identified by `get_alpha`
#' @param xdecades `numeric`, how decades should be plotted in the x-direction
#' @param ydecades `numeric`, how decades should be plotted in the y-direction
#' @return a `ggplot` object
#' @export
#' @keywords zeta
beta_plot <- function(df, change_point, xdecades= 3,  ydecades = 3){
	p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$dr, y = .data$hhcf)) +
		ggplot2::geom_point(ggplot2::aes(alpha = 0.1)) +
		ggplot2::scale_x_log10(
			breaks = scales::trans_breaks(n = xdecades, 'log10', function(x) 10^x),
			labels = scales::trans_format('log10', scales::math_format(10^.x))
		) +
		ggplot2::scale_y_log10(
			breaks = scales::trans_breaks(n = ydecades, 'log10', function(x) 10^x),
			labels = scales::trans_format('log10', scales::math_format(10^.x))
		) +
		ggplot2::geom_vline(xintercept = change_point, linetype = 2) +
		ggpubr::theme_pubr() +
		ggplot2::coord_equal() +
		ggplot2::guides(color = FALSE, alpha = FALSE)
  return(p)
}
