test_that("get_zeta_df works", {
	test_path <- system.file("extdata/rasters/", package = "statisticalRoughness")
	test_raster <- raster::raster(file.path(test_path, "gabilan_mesa.tif"))
	raster_resolution <- 10
	test_raster2 <- raster::raster(file.path(test_path, "yosemite.tif"))
	test_tiles <- raster::aggregate(test_raster, fact = 100)
	results <- get_zeta_df(test_raster, test_tiles, raster_resolution)
	test_tiles_vect <- raster::rasterToPolygons(test_tiles, dissolve = FALSE)
	results <- get_zeta_df(test_raster, test_tiles_vect, raster_resolution)
	expect_is(results, 'data.frame')
	expect_equal(nrow(results), raster::ncell(test_tiles))
	expect_error(
		get_zeta_df("test_raster", test_tiles, raster_resolution),
		"invalid class: DEM is not of class 'RasterLayer'"
	)
	expect_error(
		get_zeta_df(test_raster, "test_tiles", raster_resolution),
		"invalid class: tiles is not of class 'RasterLayer' or 'SpatialPolygonsDataFrame'"
	)
	expect_error(
		get_zeta_df(test_raster, test_tiles, "raster_resolution"),
		"invalid class: raster_resolution is not of class 'numeric'"
	)
	expect_error(
		get_zeta_df(test_raster, test_tiles, raster_resolution, vertical_accuracy = "abc"),
		"invalid class: vertical_accuracy is not of class 'numeric'"
	)
	expect_error(
		get_zeta_df(test_raster, raster::disaggregate(test_raster, fact = 2), raster_resolution),
		"There are more tiles than cells in your raster."
	)
})
