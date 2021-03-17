#' Apply a log transform to a given layer
#' @param raster_list a `list` of `stars` objects
#' @param target_id the id of the layer to be masked
#' @return a `list` of `stars`
#' @export
logtransform_layer <- function(raster_list, target_id){
  transformed_raster_list <- lapply(seq_along(raster_list), function(n){
    s <- raster_list[[n]] %>% as("Raster")
    transformed_values <- raster::getValues(s[[target_id]])
    transformed_values <- log10(transformed_values)
    transformed_values[!is.finite(transformed_values)] <- NA
    s <- raster::setValues(s, transformed_values, layer = target_id)
    return(stars::st_as_stars(s))
  })
  return(transformed_raster_list)
}
