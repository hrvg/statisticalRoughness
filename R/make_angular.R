#' Transform a 360 degrees variable to a 180 degrees
#' @param raster_list a `list` of `stars` objects
#' @param target_id the id of the layer to be masked
#' @param mode `character`, should the function consider a full or a half-circle
#' @return a `list` of `stars`
#' @export
#' @keywords postprocessing
make_angular <- function(raster_list, target_id, mode = "half") {
  if (!mode %in% c("half", "full", "sin2", "cos", "sin")) stop("make_angular: `mode` should be `half`, `full`, `sin2`, `cos`, or `sin`.")
  transformed_raster_list <- lapply(seq_along(raster_list), function(n, .mode = mode) {
    s <- raster_list[[n]] %>% as("Raster")
    transformed_values <- raster::getValues(s[[target_id]])
    if (mode == "half") transformed_values <- sapply(transformed_values, function(x) x %% 180)
    if (mode == "full") transformed_values <- sapply(transformed_values, function(x) x %% 360)
    if (mode == "sin2") transformed_values <- sapply(transformed_values, function(x) pracma::sind(x)^2)
    if (mode == "sin") transformed_values <- sapply(transformed_values, function(x) pracma::sind(x))
    if (mode == "cos") transformed_values <- sapply(transformed_values, function(x) pracma::cosd(x))
    # transformed_values <-  sapply(transformed_values, function(x) pracma::sind((x %% 180)))
    # transformed_values <-  sapply(transformed_values, function(x) ifelse(x > 180, x - 180, x))
    # transformed_values[!is.finite(transformed_values)] <- NA
    s <- raster::setValues(s, transformed_values, layer = target_id)
    return(stars::st_as_stars(s))
  })
  return(transformed_raster_list)
}
