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
	if (length(stats::na.omit(log10(unlist(row)))) < 30){
		return(data.frame(slope1 = NA, slope2 = NA, change_point = NA, max_point = NA, adj.r.squared = NA))
	}
	df <- data.frame(dr = seq_along(row) * dr, hhcf = unlist(row))
	# binned_hhcf <- bin(log10(df$dr), log10(df$hhcf), 60)
	# if (nrow(stats::na.omit(binned_hhcf)) < 30){
	# 	return(data.frame(slope1 = NA, slope2 = NA, change_point = NA, max_point = NA))
	# }
	# df_bin <- data.frame(dr = 10^binned_hhcf[, 1], hhcf = 10^binned_hhcf[, 2])

	hhcf_fun1 <- stats::splinefun(x = log10(df$dr), y = log10(df$hhcf))
	d_hhcf_dr <- hhcf_fun1(log10(df$dr), deriv = 1)

	ind_neg <- which(d_hhcf_dr < 0)
	if (length(ind_neg) < 1){
		return(data.frame(slope1 = NA, slope2 = NA, change_point = NA, max_point = NA, adj.r.squared = NA))
	} else {
		df_filtered <- df[1:(ind_neg %>% utils::head(1)), ]
		if (nrow(df_filtered) < 5){
			return(data.frame(slope1 = NA, slope2 = NA, change_point = NA, max_point = NA, adj.r.squared = NA))
		}

		x <- log10(df_filtered$dr)
		y <- log10(df_filtered$hhcf)
		segmented.fit <- segmented::segmented(lm(y ~ x))
		alpha <- data.frame(
			change_point = 10^segmented.fit$psi[1,2],
			slope1 = segmented.fit$coefficients[2],
			slope2 = segmented.fit$coefficients[2:3] %>% sum(),
			max_point = max(df_filtered$dr, na.rm = TRUE),
			adj.r.squared = summary(segmented.fit)$adj.r.squared
		)
		if(any(summary(segmented.fit)$Ttable[, 4] %>% stats::na.omit() > 0.05)){
			alpha <- data.frame(
				change_point = NA,
				slope1 = NA,
				slope2 = NA,
				max_point = NA,
				adj.r.squared = NA
			)
		}

		# min_dr <- stats::stats::na.omit(df_bin) %>% dplyr::pull(.data$dr) %>% min()
		# df_bin <- df %>% dplyr::filter(.data$dr < min_dr)
		# df_bin <- dplyr::bind_rows(df %>% dplyr::filter(.data$dr < min_dr), stats::na.omit(df_bin))

		# binned_hhcf <- bin(log10(df_filtered$dr), log10(df_filtered$hhcf), 30)
		# if (nrow(stats::na.omit(binned_hhcf)) < 5){
		# 	return(data.frame(slope1 = NA, slope2 = NA, change_point = NA, max_point = NA, adj.r.squared = NA))
		# }
		# df_bin <- data.frame(dr = 10^binned_hhcf[, 1], hhcf = 10^binned_hhcf[, 2]) %>% stats::na.omit()

		# hhcf_fun2 <- stats::splinefun(x = log10(df_bin$dr), y = log10(df_bin$hhcf))
		# dd_hhcf_dr <- hhcf_fun2(log10(df_bin$dr), deriv = 2)

		# ind_neg <- which(dd_hhcf_dr < 0)
		# if (length(ind_neg) < 1){
		# 	return(data.frame(slope1 = NA, slope2 = NA, change_point = NA, max_point = NA))
		# } else {
		# 	change_point <- df_bin$dr[ind_neg %>% utils::head(1)]
		# 	alpha_1 <- hhcf_fun1(log10(df_filtered$dr)[which(df_filtered$dr <= change_point)], deriv = 1) %>% median()
		# 	alpha_2 <- hhcf_fun1(log10(df_filtered$dr)[which(df_filtered$dr >= change_point)], deriv = 1) %>% median()
		if(do_plot){
			p <- alpha_plot(df, df_filtered, alpha$change_point)
			print(p)
		} 
			# return(data.frame(slope1 = alpha_1, slope2 = alpha_2, change_point = change_point, max_point = max(df_filtered$dr, na.rm = TRUE)))
		return(alpha)
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
	.x <- NULL
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

#' Filters the slope values of alpha returned by `get_all_alpha()` between zero and the 99.9% quantile of the initial slopes
#' @param alpha, a `data.frame` returned by `get_all_alpha()`
#' @param prob, `numeric` the value of the probabilities passed to `quantile()`
#' @return a `data.frame`
#' @export
#' @keywords zeta
#' @importFrom rlang .data
filter_alpha <- function(alpha, prob = .999){
	alpha <- stats::na.omit(alpha)
	threshold1 <- stats::quantile(alpha$slope1, probs = prob[1])
	threshold2 <- stats::quantile(alpha$slope2, probs = prob[1])
	alpha <- alpha %>% dplyr::filter(.data$slope1 > 0, .data$slope2 > 0, .data$slope1 <= threshold1, .data$slope2 <= threshold2)
	return(alpha)
}

#' Wrapper to summarize `alpha` values
#' @param alpha, a `data.frame` returned by `get_all_alpha()` or `filter_alpha()`
#' @return a summarized `data.frame`
#' @export
#' @keywords zeta
#' @importFrom rlang .data
summarise_alpha <- function(alpha){
	alpha %>% stats::na.omit(alpha) %>% dplyr::summarise(dplyr::across(dplyr::everything(), list(min = min, mean = mean, median = stats::median, max = max, sd = sd, mode = modeest::meanshift)))
}
