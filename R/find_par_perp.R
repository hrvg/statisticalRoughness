#' Add parallel and perpendicular directions and modifies the values of the x,y raster so that the raster in the perpendicular direction has the highest value for the defining direction.
#' 
#' Here the perpendicular direction is defined so that alpha_1_perp > alpha_1_par.
#' 
#' @param raster_list a `list` of `stars` objects
#' @param xdir_id `numeric`, the id of the layer with the x direction
#' @param ydir_id `numeric`, the id of the layer with the y direction
#' @param xtarget_id `numeric`, optional, the id of the target layer with the x direction
#' @param ytarget_id `numeric`, optional, the id of the target layer with the y direction
#' @param index_perp `numeric`, optional, the index of values in the perpendicular direction
#' @param index_par `numeric`, optional, the index of values in the parallel direction
#' @return a `list` of `stars` 
#' @export
find_par_perp <- function(raster_list, xdir_id, ydir_id, xtarget_id = NULL, ytarget_id = NULL, index_perp = NULL, index_par = NULL){
	res <- lapply(seq_along(raster_list), function(n){
		s <- raster_list[[n]] %>% as("Raster")

		if (is.null(xtarget_id)) xtarget_id <- xdir_id
		r_target_perp <- s[[xtarget_id]]

		if (is.null(ytarget_id)) ytarget_id <- ydir_id
		r_target_par <- s[[ytarget_id]]

		s <- raster::stack(s, r_target_perp, r_target_par)

		if(is.null(index_perp)){
			i_perp <- which(raster::getValues(s[[xdir_id]]) > raster::getValues(s[[ydir_id]]))
		} else {
			i_perp <- index_perp[[n]]
		}
		perp_values <- raster::getValues(s[[ytarget_id]]) # pre-assign when predicate is FALSE
		perp_values[i_perp] <- raster::getValues(s[[xtarget_id]])[i_perp] # assign when predicate is TRUE
		s <- raster::setValues(s, perp_values, layer = raster::nlayers(s)-1)

		if(is.null(index_par)){
			i_par <- which(raster::getValues(s[[ydir_id]]) < raster::getValues(s[[xdir_id]]))
		} else {
			i_par <- index_par[[n]]
		}
		par_values <- raster::getValues(s[[xtarget_id]]) # pre-assign when predicate is FALSE
		par_values[i_par] <- raster::getValues(s[[ytarget_id]])[i_par] # assign when predicate is TRUE
		s <- raster::setValues(s, par_values, layer = raster::nlayers(s))
		
		return(list(r = stars::st_as_stars(s, proxy = FALSE), index_perp = i_perp, ratio_perp = length(i_perp) / length(stats::na.omit(perp_values)), index_par = i_par))
	})
	
	par_perp_raster_list <- lapply(res, function(l) l$r)

	if(!is.null(index_perp) & !is.null(index_par)){
		return(par_perp_raster_list)
	} else {
		ratio_perp <- lapply(res, function(l) l$ratio_perp)
		index_perp <- lapply(res, function(l) l$index_perp)
		index_par <- lapply(res, function(l) l$index_par)
		return(list(raster_list = par_perp_raster_list, index_perp = index_perp, index_par = index_par, ratio_perp = ratio_perp))
	}
}