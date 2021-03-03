#' Check that at least four values of the raster stack are not `NA` for each element of `raster_list`
#' @param raster_list a `list` of `stars` objects
#' @export
#' @return `logical`
#' @keywords postprocessing
four_values_check <- function(raster_list){
	sapply(raster_list, function(s){
		# notNA <- sapply(seq(dim(s)["band"]), function(i) sum(!is.na(c(s[,,,i][[1]]))))
	  	# all(notNA > 4)
		df <- as.data.frame(s) %>% dplyr::select(-c("x", "y")) %>% dplyr::group_by(.data$band) %>% dplyr::group_map(~.)
		df <- do.call(cbind, df)
		nrow(stats::na.omit(df)) > 4
	})
}

