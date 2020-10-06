test_that("all_rasters_to_polygons throws an error when classes are wrong", {
  expect_error(
  	all_rasters_to_polygons(TRUE, get_wd()),
  	"invalid class: input_dir is not of class 'character'")
  expect_error(
  	all_rasters_to_polygons(getwd(), TRUE),
  	"invalid class: output_dir is not of class 'character'")
  expect_error(
  	all_rasters_to_polygons("xyz", getwd()),
	paste("xyz", "does not exist"))
    expect_error(
  	all_rasters_to_polygons(getwd(), "xyz"),
	paste("xyz", "does not exist"))
  expect_error(
  	all_rasters_to_polygons(getwd(), getwd(), reference_path = TRUE),
  	"invalid class: reference_path is not of class 'character'")
  supported_ext <- c("tif", "grd")
  expect_error(
  	all_rasters_to_polygons(getwd(), getwd(), reference_path = "raster.netcdf"),
  	paste0("invalid path: reference_path is not a path to a valid raster file; suported extensions: ", paste(supported_ext, collapse = ", ")))
  expect_error(
  	all_rasters_to_polygons(getwd(), getwd(), reference_path = "raster.tif"),
  	paste("invalid path: I cannot find a reference raster at:", "raster.tif"))
  expect_error(
  	all_rasters_to_polygons(getwd(), getwd()),
  	paste("No supported file found at:", getwd())
  	)
})

test_that("all_rasters_to_polygons works", {
	test_path <- system.file("extdata/run127b/", package = "statisticalRoughness")
	test_raster <- raster::raster(file.path(test_path, "Hurst_raster_23.tif"))
	polygons <- all_rasters_to_polygons(test_path, test_path)
	expect_is(polygons, "list")
	expect_is(polygons[[1]], "SpatialPolygons")
	expect_equal(raster::crs(polygons[[1]]), raster::crs(test_raster))
	expect_equal(length(polygons[[1]]), raster::ncell(test_raster))
})