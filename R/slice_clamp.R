#' Slicing and clamping of rasters
#' @param raster_list a `list` of `stars` or `RasterStack` objects
#' @param att_names `character`, the attributes corresponding to each band
#' @param selected `character`, the selected attributes
#' @param clamp_raster `logical`, if `TRUE` clamps the rasters using `clamp_params`
#' @param clamp_params `data.frame` of parameters to clamp rasters using `clamp_raster_sigmas()`
#' @return if `clamp` is `FALSE` a `ggplot` object, otherwise a named `list` with two elements, a `ggplot` object and a list of clamped rasters.
#' @importFrom methods as
#' @export
slice_clamp <- function(raster_list = NULL, att_names = NULL, selected = NULL, clamp_raster = FALSE, clamp_params = NULL){
	# checks and tests
	if(class(raster_list) != "list") stop("`raster_list` is not a `list`.")
	if(!length(raster_list) >= 1) stop("`raster_list` is empty.")
	if(!all(sapply(raster_list, function(obj) any(class(obj) %in% c("RasterStack", "stars", "stars_proxy"))))) stop("`raster_list` does not contain only `stars` or `Raster` objects.")
	if(class(att_names) != "character") stop("`att_names` is not a `character`.")
	if(class(selected) != "character") stop("`selected` is not a `character`.")
	if(length(unique(sapply(raster_list, function(obj) dim(obj)[3]))) != 1) stop("Not all rasters have the same number of bands.")
	if(dim(raster_list[[1]])[3] != length(att_names)) stop("`raster_list` and `att_names` have different lengths.")
	if(!all(selected %in% att_names)) stop("Some or all `selected` attributes are not contained in `att_names`.")
	if(!length(selected) >= 2) stop("At least two attributes should be selected.")
	if(class(clamp_raster) != "logical") stop("`clamp_raster` is not `logical`")
	if(clamp_raster){
		if(class(clamp_params) != "data.frame") stop("`clamp_params` is not a `data.frame`")
		if(nrow(clamp_params) != length(selected)) stop("Some or all `clamp_params` are missing.")
		if(!(all(colnames(clamp_params) %in% c("selected", "n_sigma", "lower", "upper", "lower_clamp", "upper_clamp")))) stop("`clamp_params` does not have the correct format.")
	}

	sliced_rasters <- lapply(raster_list, function(s) dplyr::slice(s, along = "band", index = match(selected, att_names)))

	if(clamp_raster){
		clamped_rasters <- lapply(seq_along(selected), function(i){
			pars <- clamp_params %>% dplyr::filter(selected == selected[i])
			lr <- clamp_raster_sigmas(sliced_rasters, n_sigma = pars$n_sigma, band_id = i, lower = pars$lower, upper = pars$upper, lower_clamp = pars$lower_clamp, upper_clamp = pars$upper_clamp)
			return(lr$rasters)
		})
		clamped_rasters <- do.call(function(...) Map(list, ...), clamped_rasters)
		clamped_rasters <- lapply(clamped_rasters, function(lr){
			lr <- lapply(lr, function(r) as(r, "Raster"))
			s <- do.call(raster::stack, lr)
			stars::st_as_stars(s) %>% stars::st_set_dimensions(3, values = NULL)
		})
		sliced_rasters <- clamped_rasters
	}
	return(sliced_rasters)
}


#' Filter outliers from a list of raster data using data from all rasters
#' 
#' This function uses a `list` as input so that the object can have vastly different resolutions.
#' 
#' @param raster_list a `list` of `stars` objects
#' @param n_sigma `numeric`, number of standard deviations to consider to calculate the percentile corresponding to values to filter, default to 3
#' @param band_id `numeric`, identifier of the band of the `stars` object
#' @param lower `logical`, if `TRUE` a lower clamp value is computed using the specified number of standard deviations
#' @param upper `logical`, if `TRUE` an upper clamp value is computed using the specified number of standard deviations
#' @param lower_clamp `numeric`, a value for the lower clamp, default to `-Inf`
#' @param upper_clamp `numeric`, a value for the upper clamp, default to `Inf`
#' @param use_values logical. If `FALSE` values outside the clamping range become `NA`, if `TRUE`, they get the extreme values
#' @return a named `list` of with the clamped `stars` objects and the clamped values of all the rasters (from ploting purposes)
#' @importFrom methods as
#' @export
clamp_raster_sigmas <- function(raster_list, n_sigma = 3, band_id = 1, lower = TRUE, upper = TRUE, lower_clamp = -Inf, upper_clamp = Inf, use_values = FALSE){
	assign("band_id", band_id, envir = .GlobalEnv) # this is not best practice but there is an issue with slice I cannot solve
	sigma_level <- stats::pnorm(n_sigma * c(-1,1), mean = 0, sd = 1, log = FALSE) %>% diff()
	all_values <- sapply(seq_along(raster_list), function(i){
		v <- dplyr::slice(raster_list[[i]], along = "band", index = band_id) %>% as.data.frame()
		v[[ncol(v)]]
	})
	all_values <- unlist(all_values) %>% as.numeric() %>% stats::na.omit()
	if (lower) lower_clamp <- stats::quantile(all_values, probs = 1 - sigma_level)
	if (upper) upper_clamp <- stats::quantile(all_values, probs = sigma_level)
	clamped_rasters <- lapply(seq_along(raster_list), function(i){
		band_index <- band_id
		r <- dplyr::slice(raster_list[[i]], along = "band", index = band_id)
		as(r, "Raster") %>% 
		raster::clamp(lower = lower_clamp, upper = upper_clamp, useValues = use_values) %>% 
		stars::st_as_stars()
	})
	clamped_values <- all_values
	clamped_values[clamped_values < lower_clamp] <- NA
	clamped_values[clamped_values > upper_clamp] <- NA
	clamped_values <- stats::na.omit(clamped_values)
	return(list(rasters = clamped_rasters, values = clamped_values))
}