# libraries
devtools::load_all()

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

fpath <- file.path(out_dir, paste0("zeta_raster_", n, "_", l, ".tif"))
if (!file.exists(fpath)){
	zeta_df <- get_zeta_df(
	  rstr, tiles, 
	  raster_resolution, vertical_accuracy = 1.87
	)
	zeta_raster <- get_zeta_raster(
	  rstr, tiles, 
	  .zeta_df = zeta_df, raster_resolution, vertical_accuracy = 1.87
	)
	raster::writeRaster(zeta_raster, fpath)
}