#' Creates a facet plot of correlation between `selected` attributes from a `raster_list` representing different `spatial_scales`
#' @param raster_list a `list` of `stars` or `RasterStack` objects
#' @param att_names `character`, the attributes corresponding to each band
#' @param selected `character`, the selected attributes
#' @param spatial_scales `numeric`, the spatial scales corresponding to each raster in `raster_list`
#' @param corr_type `character`, the type of correlation to compute, 'pearson' or 'spearman'
#' @param clamp_raster `logical`, if `TRUE` clamps the rasters using `clamp_params`
#' @return if `clamp` is `FALSE` a `ggplot` object, otherwise a named `list` with two elements, a `ggplot` object and a list of clamped rasters.
#' @importFrom rlang .data
#' @export
#' @keywords postprocessing
crossscale_correlations <- function(raster_list = NULL, att_names = NULL, selected = NULL, spatial_scales = NULL, corr_type = "spearman", clamp_raster = FALSE){
	# checks and tests
	.x <- NULL
	if(class(raster_list) != "list") stop("`raster_list` is not a `list`.")
	if(!length(raster_list) >= 1) stop("`raster_list` is empty.")
	if(!all(sapply(raster_list, function(obj) any(class(obj) %in% c("RasterStack", "stars", "stars_proxy"))))) stop("`raster_list` does not contain only `stars` or `Raster` objects.")
	if(class(selected) != "character") stop("`selected` is not a `character`.")
	if(!class(spatial_scales) %in% c("integer", "numeric")) stop("`spatial_scales` is not `integer` or `numeric`.")
	if(length(unique(sapply(raster_list, function(obj) dim(obj)[3]))) != 1) stop("Not all rasters have the same number of bands.")
	if(length(raster_list) != length(spatial_scales)) stop("`raster_list` and `spatial_scales` have different lengths.")
	if(!all(selected %in% att_names)) stop("Some or all `selected` attributes are not contained in `att_names`.")
	if(!length(selected) >= 2) stop("At least two attributes should be selected.")
	if(! corr_type %in% c("spearman", "pearson")) stop("`corr_type should be either 'spearman' or 'pearson'.")

	list_values <- lapply(raster_list, function(s){
		df <- as.data.frame(s) %>% dplyr::select(-c("x", "y")) %>% dplyr::group_by(.data$band) %>% dplyr::group_map(~.)
		df <- do.call(cbind, df)
		names(df) <- selected
		return(df)
	})
	list_correlations <- lapply(list_values, function(df){
		df <- stats::na.omit(df)
		corr <- Hmisc::rcorr(as.matrix(df), type = corr_type)$r 
		corr[lower.tri(corr, diag = TRUE)] <- NA
		as.data.frame(corr)
	})
	list_graphics <- lapply(seq_along(list_correlations), function(i){
		list_correlations[[i]] %>%
		dplyr::add_rownames() %>% 
		reshape2::melt(id.names = .data$rowname) %>% 
		stats::na.omit() %>%
		dplyr::rename(var1 = .data$rowname, var2 = .data$variable) %>%
		dplyr::mutate(scale = spatial_scales[[i]])
	})
	graphics_df <- do.call(rbind, list_graphics) %>% dplyr::arrange(.data$var1, .data$var2)

	p <- ggplot2::ggplot(graphics_df, ggplot2::aes(x = .data$scale, y = .data$value)) +
		ggplot2::scale_x_log10(
				breaks = scales::trans_breaks(n=3, 'log10', function(x) 10^x),
	            labels = scales::trans_format('log10', scales::math_format(10^.x))
	            	) +
		ggplot2::stat_smooth(fill = NA, alpha = .3) +
		ggplot2::guides(colour = FALSE) +
	    ggplot2::stat_smooth(lwd = 1, colour = 'black', alpha = .9) +
		ggplot2::facet_wrap(.data$var1 ~ .data$var2) +
		ggplot2::geom_point(alpha = 1) +
		ggpubr::theme_pubr() +
		ggplot2::ylim(c(min(graphics_df$value), max(graphics_df$value))) +
		ggplot2::labs(x = "spatial scale (m)", y = paste(corr_type, "correlation"), title = "cross-scale correlations")
	p <- p + ggforce::facet_wrap_paginate(.data$var1 ~ .data$var2, nrow = 2, ncol = 2)
		p <- lapply(seq(ggforce::n_pages(p)), function(i) p + ggforce::facet_wrap_paginate(.data$var1 ~ .data$var2, nrow = 2, ncol = 2, scales = "free", page = i))
	if(clamp_raster){
		return(list(p = p, rasters = raster_list))
	} else {
		return(p)
	}
}