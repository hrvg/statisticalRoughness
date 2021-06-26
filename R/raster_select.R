#' Selecting raster bands
#' @param raster_list a `list` of `stars` or `RasterStack` objects
#' @param band_id `numeric`, identifier of the band of the `stars` object
#' @export
#' @keywords postprocessing
raster_select <- function(raster_list, band_id) {
  selected_rasters <- lapply(seq_along(raster_list), function(i) {
    dplyr::slice(raster_list[[i]], along = "band", index = band_id)
  })
  all_values <- sapply(seq_along(raster_list), function(i) {
    v <- dplyr::slice(raster_list[[i]], along = "band", index = band_id) %>% as.data.frame()
    v[[ncol(v)]]
  })
  all_values <- unlist(all_values) %>%
    as.numeric() %>%
    stats::na.omit()
  return(list(rasters = selected_rasters, values = all_values))
}
