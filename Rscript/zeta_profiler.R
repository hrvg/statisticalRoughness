# libraries
devtools::load_all()
tictoc::tic()

# init
root_dir <- '/home/hguillon/research'
crs_ref <- raster::crs("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 + +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")
options(future.globals.maxSize= 2 * 1024^3) # 2 GiB

l <- 28
n <- 30

# dem reading
fname <- "CA_DEM.grd"
fname_sans_ext <- tools::file_path_sans_ext(fname)
dem_dir <- 'data/california-rivers/gis-files/NED/'
dem_dir <- file.path(dem_dir, paste0(fname_sans_ext, "_tiles/"))
rstr <- raster::raster(file.path(root_dir, dem_dir, paste0(fname_sans_ext, "_tile_", n, ".tif")))
raster_resolution <- 10

vertical_accuracy = 1.8
mode = "radial"

Hurst_dir <- "out/run127"
if(root_dir != "/home/hguillon") Hurst_dir <- paste0("exploitation/", Hurst_dir) 
Hurst_dir <- file.path(root_dir, Hurst_dir)

run_dir <- "out/run152"
if(root_dir != "/home/hguillon") run_dir <- paste0("exploitation/", run_dir) 
out_dir <- file.path(root_dir, run_dir)
if(!dir.exists(out_dir)) dir.create(out_dir)

# main
Hurst_raster <- raster::stack(file.path(Hurst_dir, paste0("Hurst_raster_", n, "_",l , ".tif")))
tiles <- Hurst_raster[[2]]

# 

# coercion
DEM <- rstr
if (class(tiles) == "RasterLayer") tiles <- raster::rasterToPolygons(tiles, dissolve = FALSE)
DEM <- raster::extend(DEM, tiles)

# spatial check
if (length(tiles) > raster::ncell(DEM)) stop("There are more tiles than cells in your raster.")

# i <- seq_along(tiles) %>% sample(1)
i <- length(tiles) %/% 2

registerDoFuture()
	# if (length(tiles) < availableCores() - 2){
	# 	plan(sequential)
	# } else {
if(.Platform$OS.type == "unix"){
	plan(multicore, workers = availableCores())
} else {
	plan(multisession, workers = availableCores())
}

cropped_DEM <- raster::crop(DEM, tiles[i, ])
cropped_DEM_values <- raster::getValues(cropped_DEM)
pct_non_na <- sum(is.finite(cropped_DEM_values)) / raster::ncell(cropped_DEM)
if ((pct_non_na > 1/3) & !(stats::IQR(stats::na.omit(cropped_DEM_values)) < vertical_accuracy)){ 
	zeta_df <- get_zeta(cropped_DEM, raster_resolution, .mode = mode)
}

print(zeta_df)
tictoc::toc()