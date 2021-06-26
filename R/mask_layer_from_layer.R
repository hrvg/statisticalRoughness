#' Mask a layer with another
#' @param raster_list a `list` of `stars` objects
#' @param target_id the id of the layer to be masked
#' @param mask_id the id of the layer to use as a mask, currently default to pvalue mask with a threshold at 0.05
#' @param option `character` one of `pval`, `r2`, `NA`
#' @param pval_threshold `numeric` threshold of p-values, default to 0.05
#' @param r2_threshold `numeric` threshold of r squared, default to 0.95
#' @return a `list` of `stars`
#' @export
mask_layer_from_layer <- function(raster_list, target_id, mask_id, option = "pval", pval_threshold = 0.05, r2_threshold = 0.95) {
  masked_raster_list <- lapply(seq_along(raster_list), function(n) {
    s <- raster_list[[n]] %>% as("Raster")
    if (option == "pval") {
      index_values <- which(raster::getValues(s[[mask_id]]) > pval_threshold)
    } else if (option == "r2") {
      index_values <- which(raster::getValues(s[[mask_id]]) < r2_threshold)
    } else if (option == "NA") {
      index_values <- which(is.na(raster::getValues(s[[mask_id]])))
    }
    masked_values <- raster::getValues(s[[target_id]])
    masked_values[index_values] <- NA
    s <- raster::setValues(s, masked_values, layer = target_id)
    return(stars::st_as_stars(s))
  })
  return(masked_raster_list)
}
