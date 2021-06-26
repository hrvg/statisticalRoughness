#' Produces a leaflet map
#' @param clamped named `list` from `clamp_raster_sigmas()`
#' @param ttl `character`, title
#' @param n_class `numeric`, number of bins for the color scheme
#' @param circular `logical`, if `TRUE` the value is transformed so be between 0 and 180, default to `FALSE`
#' @param style `character`, changes between continuous and binned values for the colorscale, default to `continuous`
#' @param groups `character`, overlay groups to adjust map visuals
#' @importFrom methods as
#' @return a `leaflet` object
#' @export
#' @keywords postprocessing
make_leaflet_map <- function(clamped, ttl = "values", n_class = 10, circular = FALSE, style = "continuous", groups) {
  groups <- sapply(groups, as.character)
  ma <- NULL
  ma <- leaflet::leaflet()
  cInt <- classInt::classIntervals(clamped$values, n_class, style = "fisher")
  brk <- cInt$brks
  if (circular) {
    clamped$values <- sapply(clamped$values, function(x) ifelse(x > 180, x - 180, x))
    if (style == "continuous") {
      cols <- leaflet::colorNumeric("Spectral", clamped$values, na.color = NA)
    } else {
      cols <- leaflet::colorBin("Spectral", clamped$values, bins = brk, na.color = NA)
    }
  } else {
    if (style == "continuous") {
      cols <- leaflet::colorNumeric("viridis", clamped$values, na.color = NA)
    } else {
      cols <- leaflet::colorBin("viridis", clamped$values, bins = brk, na.color = NA)
    }
  }
  for (ind_H in seq_along(clamped$rasters)) {
    r <- clamped$rasters[[ind_H]]
    if (circular) {
      r <- as(r, "Raster") %>%
        raster::calc(function(x) ifelse(x > 180, x - 180, x)) %>%
        stars::st_as_stars()
    }
    ma <- ma %>% leafem::addStarsImage(r, colors = cols, opacity = 1, project = TRUE, maxBytes = Inf, group = groups[ind_H])
  }
  ma <- ma %>%
    leaflet::addLegend(pal = cols, values = clamped$values, title = ttl) %>%
    leaflet::addProviderTiles(
      leaflet::providers$Stamen.TerrainBackground,
      group = "Stamen.TerrainBackground",
      options = leaflet::providerTileOptions(maxZoom = 22)
    ) %>%
    leaflet::addProviderTiles(
      leaflet::providers$Esri.WorldImagery,
      group = "Esri.WorldImagery",
      options = leaflet::providerTileOptions(maxZoom = 22)
    ) %>%
    leaflet.extras::addResetMapButton() %>%
    leaflet.extras::addFullscreenControl(position = "topleft", pseudoFullscreen = FALSE) %>%
    leaflet::addLayersControl(
      baseGroups = c("Stamen.TerrainBackground", "Esri.WorldImagery"),
      overlayGroups = groups,
      options = leaflet::layersControlOptions(collapsed = FALSE)
    ) %>%
    leaflet::hideGroup(utils::tail(groups, -1))
  return(ma)
}
