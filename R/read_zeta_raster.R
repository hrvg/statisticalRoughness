#' Read zeta raster from destination
#' @param out_dir `character` file path to folder
#' @param Lmax `numeric` maximum scale
#' @param raster_resolution `numeric` raster resolution, to determine the scales
#' @param .len numbers of scales, passed to `get_all_R_L()` with a default of 5 factors and only with factors of 12
#' @param ... passed to `stars::read_stars()` to allow `proxy = FALSE`
#' @return a `list` with two elements containing the rasters as a `list` of `stars` objects and the spatial scales
#' @export
read_zeta_raster <- function(out_dir = "/out/run128/out", Lmax = 1E4, raster_resolution = 1, .len = 32, ...){
	lf <- list.files(path = file.path(out_dir), pattern = ".tif")
	lf_ext <- unname(sapply(lf, tools::file_ext))
	lf <- lf[which(lf_ext == "tif")]
	spatial_scales <- unlist(lapply(lf, function(f) utils::tail(unlist(strsplit(tools::file_path_sans_ext(f), "_")), 1)))
	spatial_scales <- unname(sapply(spatial_scales, as.numeric))
	lf <- lf[order(spatial_scales)]
	raster_list <- lapply(lf, function(f) stars::read_stars(file.path(out_dir,f), ...))
	allLR <- statisticalRoughness::get_all_R_L(Lmax, 5, only = 12, len = .len)
	spatial_scales <- allLR$allL * raster_resolution 
	return(list(raster_list = raster_list, spatial_scales = spatial_scales))
}