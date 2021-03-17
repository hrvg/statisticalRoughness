#' Add parallel and perpendicular directions
#' @param raster_list a `list` of `stars` objects
#' @param xdir_id the id of the layer with the x direction
#' @param ydir_id the id of the layer with the y direction
#' @return a `list` of `stars` 
#' @export
find_par_perp <- function(raster_list, xdir_id, ydir_id){
	par_perp_raster_list <- lapply(seq_along(raster_list), function(n){
		s <- raster_list[[n]] %>% as("Raster")
		r_perp <- s[[xdir_id]]
		r_par <- s[[ydir_id]]
		s <- raster::stack(s, r_perp, r_par)

		index_perp <- which(raster::getValues(s[[xdir_id]]) > raster::getValues(s[[ydir_id]]))
		perp_values <- raster::getValues(s[[ydir_id]])
		perp_values[index_perp] <- raster::getValues(s[[xdir_id]])[index_perp]
		s <- raster::setValues(s, perp_values, layer = raster::nlayers(s)-1)

		index_par <- which(raster::getValues(s[[ydir_id]]) < raster::getValues(s[[xdir_id]]))
		par_values <- raster::getValues(s[[xdir_id]])
		par_values[index_par] <- raster::getValues(s[[ydir_id]])[index_par]
		s <- raster::setValues(s, par_values, layer = raster::nlayers(s))
		return(stars::st_as_stars(s, proxy = FALSE))
	})
	return(par_perp_raster_list)
}