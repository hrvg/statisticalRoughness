#' This function transforms all rasters found into a directory into SpatialPolygons and save them to shapefiles.
#' 
#' The main use of this function is to parallelize processing for each pixels of each rasters.
#' 
#' @param input_dir `character`, input directory
#' @param output_dir `character`, output directory
#' @param reference_path optional, `character`, path to a reference raster used to crop the rasters found in `input_dir`
#' @return `list` of `SpatialPolygons
#' @import foreach
#' @import doFuture
#' @import raster
#' @import sp
#' @import rgdal
all_rasters_to_polygons <- function(input_dir, output_dir, reference_path = NULL){
	supported_ext <- c("tif", "grd")
	raster_file <- NULL
	if(!class(input_dir) == "character") stop("invalid class: input_dir is not of class 'character'")
	if(!class(output_dir) == "character") stop("invalid class: output_dir is not of class 'character'")
	if(!dir.exists(input_dir)) stop(paste(input_dir, "does not exist"))
	if(!dir.exists(output_dir)) stop(paste(output_dir, "does not exist"))
	if(!is.null(reference_path)){
		if(!class(reference_path) == "character") stop("invalid class: reference_path is not of class 'character'")
		if(!tools::file_ext(reference_path) %in% supported_ext){
			stop(paste0("invalid path: reference_path is not a path to a valid raster file; suported extensions: ", paste(supported_ext, collapse = ", ")))
		} else {
			if(!file.exists(reference_path)) stop(paste("invalid path: I cannot find a reference raster at:", reference_path))
		}
	}
	file_list <- list.files(input_dir, full.names = TRUE)
	file_list_ext <- sapply(file_list, tools::file_ext)
	file_list <- file_list[which(file_list_ext %in% supported_ext)]
	if(length(file_list) == 0) stop(paste("No supported file found at:", input_dir))
	if(!is.null(reference_path)) reference_raster <- raster::raster(reference_path)
	all_polygons <- foreach(raster_file = file_list, .packages = c("raster", "sp", "rgdal")) %do% {
		r <- raster::raster(raster_file)
		if(!is.null(reference_path)) r <- raster::crop(r, reference_raster)
		pol <- raster::rasterToPolygons(r, na.rm = FALSE)
		raster::shapefile(pol, filename = file.path(output_dir, paste0(tools::file_path_sans_ext(basename(raster_file)), ".shp")), overwrite = TRUE)
		return(pol)
	}
	return(all_polygons)
}