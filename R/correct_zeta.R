#' Correct zeta
#' @param raster_list a `list` of `stars` objects
#' @param target_id the id of the layer to be masked
#' @return a `list` of `stars`
#' @export
correct_zeta <- function(raster_list, target_id) {
  corrected_raster_list <- lapply(seq_along(raster_list), function(n) {
    s <- raster_list[[n]] %>% as("Raster")
    r_zeta <- s[[target_id]]
    s <- raster::stack(s, r_zeta)
    index_values <- which(raster::getValues(s[[target_id]]) < 1)
    zeta_values <- raster::getValues(s[[target_id]])
    zeta_values[index_values] <- 1 / raster::getValues(s[[target_id]])[index_values]
    s <- raster::setValues(s, zeta_values, layer = raster::nlayers(s))
    return(stars::st_as_stars(s, proxy = FALSE))
  })
  return(corrected_raster_list)
}
