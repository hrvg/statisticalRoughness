#' merge all Hurst raster
#' @param Hurst_dir a file.path to a directory containing tiled Hurst rasters
#' @importFrom magrittr %>%
#' @export
#' @keywords Hurst
merge_Hurst_rasters <- function(Hurst_dir = file.path("F:/hguillon/research/exploitation/out/run127")) {
  path <- xmin <- xmax <- ymin <- ymax <- NULL
  out_dir <- file.path(Hurst_dir, "out")
  if (!dir.exists(out_dir)) dir.create(out_dir)
  lf <- list.files(Hurst_dir, full.names = TRUE, pattern = ".tif")
  .lf <- list.files(Hurst_dir, full.names = FALSE, pattern = ".tif")
  tile_id <- sapply(.lf, function(fn) unlist(strsplit(fn, "_"))[[3]]) %>% as.numeric()
  L_id <- sapply(.lf, function(fn) {
    unlist(strsplit(fn, "_"))[[4]] %>%
      tools::file_path_sans_ext()
  }) %>% as.numeric()
  tiled_file_df <- data.frame(path = lf, tile_id = tile_id, L_id = L_id) %>%
    dplyr::arrange(tile_id, L_id) %>%
    dplyr::group_by(L_id)
  for (l in unique(tiled_file_df$L_id)) {
    raster_paths <- tiled_file_df %>%
      dplyr::filter(L_id == l) %>%
      dplyr::pull(path) %>%
      as.character()
    rasters <- lapply(raster_paths, raster::raster)
    xt <- lapply(rasters, function(r) {
      xt <- raster::extent(r)
      data.frame(xmin = xt@xmin, xmax = xt@xmax, ymin = xt@ymin, ymax = xt@ymax)
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::summarize(xmin = min(xmin), xmax = max(xmax), ymin = min(ymin), ymax = max(ymax)) %>%
      as.numeric()
    template <- raster::raster(raster::extent(xt))
    raster::projection(template) <- as.character(raster::crs(rasters[[1]]))
    out_path <- file.path(out_dir, paste0("Hurst_raster_", l, ".tif"))
    raster::writeRaster(template, file = out_path, format = "GTiff", overwrite = TRUE)
    gdalUtils::mosaic_rasters(gdalfile = raster_paths, dst_dataset = out_path, of = "GTiff")
  }
}
