#' Wrapper to extract all roughness exponents from the HHCFs
#' @param hhcf `list` from `get_hhcf()`
#' @param dr `numeric`, spacing of the values along the axis
#' @return a `data.frame`
#' @export
#' @keywords zeta
get_all_alpha <- function(hhcf, dr){
	hhcf <- hhcf$hhcf
	all_alpha <- lapply(seq_along(hhcf), function(k) get_alpha(hhcf[k, ], dr))
	all_alpha <- do.call(rbind, all_alpha)
	return(all_alpha)
}

#' Estimates the roughness exponent from the HHCF from spline fit.
#' 
#' This function will return `NA` if the HHFC is monotically increasing.
#' 
#' @param row array of hhcf values
#' @param dr `numeric`, spacing of the values along the axis
#' @param do_plot `logical`, plot the fit, default to `FALSE`
#' @export
#' @keywords zeta
get_alpha <- function(row, dr, do_plot = FALSE){
	if (length(na.omit(row)) < 30){
		return(data.frame(slope1 = NA, slope2 = NA, break1 = NA))
	}
	df <- data.frame(dr = seq_along(row) * dr, hhcf = unlist(row))
	binned_hhcf <- bin(log10(df$dr), log10(df$hhcf), 60)
	if (nrow(na.omit(binned_hhcf)) < 30){
		return(data.frame(slope1 = NA, slope2 = NA, break1 = NA))
	}
	df_bin <- data.frame(dr = 10^binned_hhcf[, 1], hhcf = 10^binned_hhcf[, 2])

	hhcf_fun <- splinefun(x = log10(df_bin$dr), y = log10(df_bin$hhcf))
	d_hhcf_dr <- hhcf_fun(log10(df_bin$dr), deriv = 1)

	ind_neg <- which(d_hhcf_dr < 0)
	if (length(ind_neg) < 1){
		return(data.frame(slope1 = NA, slope2 = NA, break1 = NA))
	} else {
		df_bin <- df_bin[1:(ind_neg %>% head(1)), ]
		min_dr <- stats::na.omit(df_bin) %>% dplyr::pull(.data$dr) %>% min()
		df_bin <- dplyr::bind_rows(df %>% dplyr::filter(.data$dr < min_dr), na.omit(df_bin))

		hhcf_fun <- splinefun(x = log10(df_bin$dr), y = log10(df_bin$hhcf))
		dd_hhcf_dr <- hhcf_fun(log10(df_bin$dr), deriv = 2)

		change_point <- df_bin$dr[which.min(dd_hhcf_dr)]
		alpha_1 <- hhcf_fun(log10(df_bin$dr)[which(df_bin$dr <= change_point)], deriv = 1) %>% mean()
		alpha_2 <- hhcf_fun(log10(df_bin$dr)[which(df_bin$dr >= change_point)], deriv = 1) %>% mean()
	if(do_plot){
		p <- alpha_plot(df, df_bin, change_point)
		print(p)
	} 
	return(data.frame(slope1 = alpha_1, slope2 = alpha_2, break1 = change_point))
	}
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
alpha_plot <- function(df, df_bin, change_point, xdecades= 3,  ydecades = 3){
	p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$dr, y = .data$hhcf)) +
		ggplot2::geom_point(ggplot2::aes(alpha = 0.1)) +
		ggplot2::geom_point(data = df_bin, ggplot2::aes(x = .data$dr , y = .data$hhcf, color = "red")) +
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