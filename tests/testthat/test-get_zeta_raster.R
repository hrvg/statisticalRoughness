test_that("get_zeta_raster works", {
	test_path <- system.file("extdata/rasters/", package = "statisticalRoughness")
	test_raster <- raster::raster(file.path(test_path, "gabilan_mesa.tif"))
	raster_resolution <- 10
	test_tiles <- raster::aggregate(test_raster, fact = 100)
	test_tiles_vect <- raster::rasterToPolygons(test_tiles, dissolve = FALSE)
	results <- get_zeta_df(test_raster, test_tiles, raster_resolution)
	expect_is(results, 'data.frame')
	expect_error(get_zeta_raster(test_raster, "test_tiles_vect", raster_resolution, .zeta_df = results), "invalid class: tiles is not of class 'RasterLayer' or 'SpatialPolygonsDataFrame'")
	expect_is(get_zeta_raster(test_raster, test_tiles, raster_resolution, .zeta_df = results), 'RasterStack')
	expect_is(get_zeta_raster(test_raster, test_tiles, raster_resolution), 'RasterStack')
})
