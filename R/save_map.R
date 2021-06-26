#' Creates and save a map of a given variablw
#' @param raster_list a named list with two elements `rasters` and `values` supplied from `raster_select()
#' @param groups the groups to use in the maps
#' @param n_class `numeric` number of classes
#' @param ttl `character` the title of the map
#' @param begin `numeric` passed to `tmap::tm_raster()`
#' @param end `numeric` passed to `tmap::tm_raster()`
#' @param direction `numeric` passed to `tmap::tm_raster()`
#' @param option `character` passed to `tmap::tm_raster()`
#' @param out_path `character` path to save the images
#' @param col_style `character`, default to `NULL`, passed as `style` to `tmap::tm_raster()`
#' @param cstyle `character`, default to `equal`, passed as `style` to `classInt::classIntervals()`
#' @export
save_map <- function(raster_list, groups, n_class = 10, ttl = "zeta", begin = 0.1, end = 0.95, direction = 1, option = "inferno", out_path = "out", col_style = NULL, cstyle = "equal") {
  groups <- sapply(groups, as.character)
  if (!dir.exists(out_path)) dir.create(out_path)
  out_path <- file.path(out_path, ttl)
  if (!dir.exists(out_path)) dir.create(out_path)
  if (direction == -1) option <- paste0("-", option)
  cInt <- classInt::classIntervals(raster_list$values, n_class, style = cstyle)
  brk <- cInt$brks
  for (ind_H in seq_along(raster_list$rasters)) {
    main_title <- paste0(ttl, " at scale: ", groups[ind_H], "-m")
    r <- raster_list$rasters[[ind_H]]
    if (is.null(col_style)) {
      ma <- tmap::tm_shape(r, name = ttl) + tmap::tm_raster(palette = option, n = n_class, breaks = brk, contrast = c(begin, end), alpha = 1, legend.format = list(format = "f", digits = 1), title = ttl) + tmap::tm_layout(legend.position = c("right", "top"), main.title = main_title)
    } else {
      ma <- tmap::tm_shape(r, name = ttl) + tmap::tm_raster(palette = option, breaks = brk, style = col_style, contrast = c(begin, end), alpha = 1, legend.format = list(format = "f", digits = 2), title = ttl) + tmap::tm_layout(legend.position = c("right", "top"), main.title = main_title)
    }
    tmap::tmap_save(ma, file.path(out_path, paste0(ind_H, ".jpg")), dpi = 600, units = "px", width = 3840 * 2, height = 2160 * 2)
  }
}
