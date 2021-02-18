#' Wrapper to extract all roughness exponents from the HHCFs
#' @param hhcf `list` from `get_hhcf()`
#' @param dr `numeric`, spacing of the values along the axis
#' @return a `data.frame`
#' @export
#' @keywords zeta
get_all_alpha <- function(hhcf, dr){
	hhcf <- hhcf$hhcf
	all_alpha <- lapply(seq_along(hhcf), function(k){
		get_alpha(hhcf[k, ], dr)
	})
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
#' @importFrom methods is
#' @export
#' @keywords zeta
get_alpha <- function(row, dr, do_plot = FALSE){
	if (length(stats::na.omit(log10(unlist(row)))) < 30){
		return(data.frame(alpha1 = NA, alpha2 = NA, rc = NA, rmax = NA, alpha.r2 = NA))
	}
	df <- data.frame(dr = seq_along(row) * dr, hhcf = unlist(row))
	# binned_hhcf <- bin(log10(df$dr), log10(df$hhcf), 60)
	# if (nrow(stats::na.omit(binned_hhcf)) < 30){
	# 	return(data.frame(alpha1 = NA, alpha2 = NA, rc = NA, rmax = NA))
	# }
	# df_bin <- data.frame(dr = 10^binned_hhcf[, 1], hhcf = 10^binned_hhcf[, 2])

	hhcf_fun1 <- stats::splinefun(x = log10(df$dr), y = log10(df$hhcf))
	d_hhcf_dr <- hhcf_fun1(log10(df$dr), deriv = 1)

	ind_neg <- which(d_hhcf_dr < 0)
	if (length(ind_neg) < 1){
		return(data.frame(alpha1 = NA, alpha2 = NA, rc = NA, rmax = NA, alpha.r2 = NA))
	} else {
		df_filtered <- df[1:(ind_neg %>% utils::head(1)), ]
		if (nrow(df_filtered) < 10){
			return(data.frame(alpha1 = NA, alpha2 = NA, rc = NA, rmax = NA, alpha.r2 = NA))
		}

		x <- log10(df_filtered$dr)
		y <- log10(df_filtered$hhcf)
		segmented.fit <- tryCatch(segmented::segmented(lm(y ~ x, weights = 1 / x)), error = function(e) e, warning = function(w) w)
		if(is(segmented.fit, "warning") | is(segmented.fit, "error")){
			return(data.frame(alpha1 = NA, alpha2 = NA, rc = NA, rmax = NA, alpha.r2 = NA))
		}
		alpha <- data.frame(
			rc = 10^segmented.fit$psi[1,2],
			alpha1 = segmented.fit$coefficients[2],
			alpha2 = segmented.fit$coefficients[2:3] %>% sum(),
			rmax = max(df_filtered$dr, na.rm = TRUE),
			alpha.r2 = summary(segmented.fit)$adj.r.squared
		)
		if(any(summary(segmented.fit)$Ttable[, 4] %>% stats::na.omit() > 0.05)){
			alpha <- data.frame(
				rc = NA,
				alpha1 = NA,
				alpha2 = NA,
				rmax = NA,
				alpha.r2 = NA
			)
		}

		# min_dr <- stats::stats::na.omit(df_bin) %>% dplyr::pull(.data$dr) %>% min()
		# df_bin <- df %>% dplyr::filter(.data$dr < min_dr)
		# df_bin <- dplyr::bind_rows(df %>% dplyr::filter(.data$dr < min_dr), stats::na.omit(df_bin))

		# binned_hhcf <- bin(log10(df_filtered$dr), log10(df_filtered$hhcf), 30)
		# if (nrow(stats::na.omit(binned_hhcf)) < 5){
		# 	return(data.frame(alpha1 = NA, alpha2 = NA, rc = NA, rmax = NA, alpha.r2 = NA))
		# }
		# df_bin <- data.frame(dr = 10^binned_hhcf[, 1], hhcf = 10^binned_hhcf[, 2]) %>% stats::na.omit()

		# hhcf_fun2 <- stats::splinefun(x = log10(df_bin$dr), y = log10(df_bin$hhcf))
		# dd_hhcf_dr <- hhcf_fun2(log10(df_bin$dr), deriv = 2)

		# ind_neg <- which(dd_hhcf_dr < 0)
		# if (length(ind_neg) < 1){
		# 	return(data.frame(alpha1 = NA, alpha2 = NA, rc = NA, rmax = NA))
		# } else {
		# 	rc <- df_bin$dr[ind_neg %>% utils::head(1)]
		# 	alpha_1 <- hhcf_fun1(log10(df_filtered$dr)[which(df_filtered$dr <= rc)], deriv = 1) %>% median()
		# 	alpha_2 <- hhcf_fun1(log10(df_filtered$dr)[which(df_filtered$dr >= rc)], deriv = 1) %>% median()
		if(do_plot){
			p <- alpha_plot(df, df_filtered, alpha$rc)
			print(p)
		} 
			# return(data.frame(alpha1 = alpha_1, alpha2 = alpha_2, rc = rc, rmax = max(df_filtered$dr, na.rm = TRUE)))
		return(alpha)
	}
}

#' Internal companion function to `get_alpha`.
#' @param df a `data.frame` with the non-binned data
#' @param df_bin a `data.frame` with the binned data
#' @param rc `numeric`, the change point identified by `get_alpha`
#' @param xdecades `numeric`, how decades should be plotted in the x-direction
#' @param ydecades `numeric`, how decades should be plotted in the y-direction
#' @return a `ggplot` object
#' @export
#' @keywords zeta
alpha_plot <- function(df, df_bin, rc, xdecades= 3,  ydecades = 3){
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
		ggplot2::geom_vline(xintercept = rc, linetype = 2) +
		ggpubr::theme_pubr() +
		ggplot2::coord_equal() +
		ggplot2::guides(color = FALSE, alpha = FALSE)
  return(p)
}

#' Filters the alpha values of alpha returned by `get_all_alpha()` between zero and the 99.9% quantile of the initial alphas
#' @param alpha, a `data.frame` returned by `get_all_alpha()`
#' @param prob, `numeric` the value of the probabilities passed to `quantile()`
#' @param rsquared_filter, `logical`, if `TRUE`, filters out the `alpha.r2` values lower than 0.99
#' @return a `data.frame`
#' @export
#' @keywords zeta
#' @importFrom rlang .data
filter_alpha <- function(alpha, prob = .999, rsquared_filter = TRUE){
	alpha <- stats::na.omit(alpha)
	threshold1 <- stats::quantile(alpha$alpha1, probs = prob[1])
	threshold2 <- stats::quantile(alpha$alpha2, probs = prob[1])
	alpha <- alpha %>% dplyr::filter(
		.data$alpha1 > 0,
		.data$alpha2 > 0,
		.data$alpha1 <= threshold1,
		.data$alpha2 <= threshold2
	)
	if (rsquared_filter){
		alpha <- alpha %>% dplyr::filter(
			.data$alpha.r2 >= 0.99
		)
	}
	return(alpha)
}

#' Wrapper to summarize `alpha` values
#' @param alpha, a `data.frame` returned by `get_all_alpha()` or `filter_alpha()`
#' @return a summarized `data.frame`
#' @export
#' @keywords zeta
#' @importFrom rlang .data
summarise_alpha <- function(alpha){
	alpha <- alpha %>% dplyr::filter_all(dplyr::all_vars(is.finite(.)))
	.list = list(min = min, mean = mean, max = max, sd = sd, IQR = stats::IQR)
	if (nrow(alpha) > 1){
		alpha <- alpha %>% dplyr::summarise(dplyr::across(dplyr::everything(), .list))
	} else {
		cnames <- outer(colnames(alpha), names(.list), paste, sep = "_") %>% t() %>% c()
		alpha <- matrix(NA, ncol = length(cnames), nrow = 1) %>% as.data.frame()
		colnames(alpha) <- cnames
	}
	return(alpha)
}
