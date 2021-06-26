test_that("all_rasters_to_polygons throws an error when classes are wrong", {
  expect_error(
    all_rasters_to_polygons(TRUE, get_wd()),
    "invalid class: input_dir is not of class 'character'"
  )
  expect_error(
    all_rasters_to_polygons(getwd(), TRUE),
    "invalid class: output_dir is not of class 'character'"
  )
  expect_error(
    all_rasters_to_polygons("xyz", getwd()),
    paste("xyz", "does not exist")
  )
  expect_error(
    all_rasters_to_polygons(getwd(), "xyz"),
    paste("xyz", "does not exist")
  )
  expect_error(
    all_rasters_to_polygons(getwd(), getwd(), reference_path = TRUE),
    "invalid class: reference_path is not of class 'character'"
  )
  supported_ext <- c("tif", "grd")
  expect_error(
    all_rasters_to_polygons(getwd(), getwd(), reference_path = "raster.netcdf"),
    paste0("invalid path: reference_path is not a path to a valid raster file; suported extensions: ", paste(supported_ext, collapse = ", "))
  )
  expect_error(
    all_rasters_to_polygons(getwd(), getwd(), reference_path = "raster.tif"),
    paste("invalid path: I cannot find a reference raster at:", "raster.tif")
  )
  expect_error(
    all_rasters_to_polygons(getwd(), getwd()),
    paste("No supported file found at:", getwd())
  )
})

test_that("all_rasters_to_polygons works", {
  test_path <- system.file("extdata/run127b/", package = "statisticalRoughness")
  test_raster1 <- terra::rast(file.path(test_path, "Hurst_raster_23.tif"))
  test_raster2 <- terra::rast(file.path(test_path, "Hurst_raster_24.tif"))
  pols <- all_rasters_to_polygons(test_path, test_path)
  expect_is(pols, "list")
  expect_is(pols[[1]], "SpatVector")
  expect_equal(terra::crs(pols[[1]]), terra::crs(test_raster1))
  expect_equal(terra::crs(pols[[2]]), terra::crs(test_raster2))
  pols <- all_rasters_to_polygons(test_path, test_path, reversed = TRUE)
  expect_equal(length(pols[[2]]), sum(!is.na(terra::values(test_raster1[[2]]))))
  expect_equal(length(pols[[1]]), sum(!is.na(terra::values(test_raster2[[2]]))))
})
