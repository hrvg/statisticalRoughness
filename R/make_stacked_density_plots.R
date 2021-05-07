#' Makes stacked density plot
#' @param raster_list a `list` of `stars` or `RasterStack` objects
#' @param band_id `numeric`, identifier of the band of the `stars` object
#' @param var_name `character`, used for the title of the graph
#' @param limx `numeric` limits of the x-axis, default to `c(0,1)`
#' @param logscale `logical` if `TRUE` the x-axis is log-transformed
#' @param spatial_scales `numeric`, the spatial scales corresponding to each raster in `raster_list`
#' @param ... passed to `ggplot2::scale_color_viridis_c()`
#' @importFrom methods as
#' @importFrom rlang .data
#' @return a `ggplot()` object
#' @export
#' @keywords postprocessing
make_stacked_density_plot <- function(raster_list, band_id = 1, var_name = "variable", limx = c(0, 1), logscale = FALSE, spatial_scales = NULL, ...){
	..density.. <- .x <- x <- NULL
	if(length(raster_list) != length(spatial_scales)) stop("Length of `spatial_scales` is not compatible with `raster_list`")
	l.df <- lapply(seq_along(raster_list), function(i){
		df <- data.frame(values = raster_list[[i]][[1]][,,band_id] %>% c())
		# df <- data.frame(values = raster::getValues(as(raster_list[[i]], "Raster")[[band_id]]))
		p <- ggplot2::ggplot(df, ggplot2::aes(x = values)) + 
			ggplot2::geom_density(ggplot2::aes(y = ..density..), bw = "SJ")
		if (logscale){
			p <- p + ggplot2::scale_x_log10(
				breaks = scales::trans_breaks(n = 3, 'log10', function(x) 10^x),
	            labels = scales::trans_format('log10', scales::math_format(10^.x)),
	            limits = limx
	        )
		}	
		p_df <- ggplot2::ggplot_build(p)$data[[1]]	
		df <- data.frame(x = p_df$x, y = p_df$y)
		df$scale <- rep(spatial_scales[i], nrow(df))
		return(df)
	})
	df <- do.call(rbind, l.df)
	p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, group = .data$scale, colour = .data$scale)) + 
			ggplot2::geom_line() +
			ggplot2::scale_color_viridis_c(...) +
			ggpubr::theme_pubr() +
			ggplot2::labs(x = paste(var_name, "values"), y = "density")
	return(p)
}

#' Extract the modal values of the stacked density plot
#' @param plot_object a `ggplot` object
#' @param spatial_scales `numeric`, the spatial scales corresponding to each raster in `raster_list`
#' @param var_name `character`, used for the title of the graph
#' @param ... passed to `ggplot2::scale_color_viridis_c()`
#' @return a `list` with three `ggplot` elements: a plot in the variable - density space, a plot in the variable mode - scale space, a density plot with the modes added.
#' @importFrom rlang .data
#' @export
#' @keywords postprocessing
modes_from_stacked_density <- function(plot_object, spatial_scales = NULL, var_name = "variable", ...){
	x <- .x <- NULL
	df <- ggplot2::ggplot_build(plot_object)$data[[1]]
	if(length(unique(df$group)) != length(spatial_scales)) stop("Length of `spatial_scales` is not compatible with the plot")
	df <- df %>% 
		dplyr::mutate(scale = spatial_scales[.data$group]) %>%
		dplyr::select(c("x", "y", "scale")) %>%
		dplyr::group_by(.data$scale) %>%
		dplyr::group_modify( ~ data.frame(y = max(.$y), x = .$x[which.max(.$y)]))
	p_mode_density <-  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, group = .data$scale, colour = .data$scale, label = .data$scale)) +
		ggplot2::geom_point() +
		ggplot2::scale_color_viridis_c(...) +
		ggpubr::theme_pubr() +
		ggplot2::labs(x = paste(var_name, "values"), y = "density") +
		ggrepel::geom_text_repel()
	p_mode_scale <- ggplot2::ggplot(df, ggplot2::aes(x = .data$scale, y = .data$x)) +
		ggplot2::geom_point() +
		ggpubr::theme_pubr() +
		ggplot2::labs(x = "spatial scale (m)", y = paste(var_name, "mode")) + 
		ggplot2::scale_x_log10(
				breaks = scales::trans_breaks(n = 5, 'log10', function(x) 10^x),
	            labels = scales::trans_format('log10', scales::math_format(10^.x))
	    )
	p_density_with_mode <- plot_object + ggplot2::geom_point(data = df, ggplot2::aes(x = .data$x, y = .data$y, group = .data$scale, colour = .data$scale))
	return(list(p_mode_density = p_mode_density, p_mode_scale = p_mode_scale, p_density_with_mode = p_density_with_mode))
}

#' Wrapper function around `make_stacked_density_plot()` and `modes_from_stacked_density()`
#' @param raster_list a `list` of `stars` or `RasterStack` objects
#' @param .spatial_scales `numeric`, the spatial scales corresponding to each raster in `raster_list`
#' @param .band_id `numeric`, identifier of the band of the `stars` object
#' @param .var_name `character`, used for the title of the graph
#' @param .base_size `numeric`, base size for the font, passed to `ggpubr::theme_pubr()`
#' @param ... passed to `make_stacked_density_plot()` and  `modes_from_stacked_density()` 
#' @return a list of grobs to be used by `gridExtra::grid.arrange()`
#' @export
#' @keywords postprocessing
make_all_plots <- function(raster_list, .spatial_scales, .band_id, .var_name, .base_size = 12, ...){
	stacked_density_plot <- make_stacked_density_plot(raster_list, spatial_scales = .spatial_scales, band_id = .band_id, var_name = .var_name, ...)
	mode_plots <- modes_from_stacked_density(stacked_density_plot, spatial_scales = .spatial_scales, var_name = .var_name, ...)
	gb <- list(
	  mode_plots$p_density_with_mode +
		  ggpubr::theme_pubr(base_size = .base_size) +
		  ggplot2::theme(legend.key.width = ggplot2::unit(.5, "in")),
	  mode_plots$p_mode_density  + 
	  	ggpubr::theme_pubr(base_size = .base_size) +
	  	ggplot2::theme(legend.key.width = ggplot2::unit(.5, "in")),
	  mode_plots$p_mode_scale + ggpubr::theme_pubr(base_size = .base_size)
	)
	return(gb)
}