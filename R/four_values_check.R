#' Check that at least four values of the raster stack are not `NA` for each element of `raster_list`
#' @param raster_list a `list` of `RasterStack` objects
#' @export
#' @return `logical`
four_values_check <- function(raster_list){
	sapply(raster_list, function(s){
	  notNA <- sapply(seq(dim(s)["band"]), function(i) sum(!is.na(c(s[,,,i][[1]]))))
	  all(notNA > 4)
	})
}

