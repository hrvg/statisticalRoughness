# libraries
devtools::load_all()

# init
root_dir <- '/home/hguillon/research'
crs_ref <- raster::crs("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 + +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")
options(future.globals.maxSize= 2 * 1024^3) # 2 GiB

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


df <- foreach(l = seq(1,30), .combine = rbind) %do% {
	Hurst_raster <- raster::stack(file.path(Hurst_dir, paste0("Hurst_raster_", n, "_",l , ".tif")))
	tiles <- Hurst_raster[[2]]
	i <- sample(length(tiles), 1)

	# coercion
	DEM <- rstr
	if (class(tiles) == "RasterLayer") tiles <- raster::rasterToPolygons(tiles, dissolve = FALSE)
	DEM <- raster::extend(DEM, tiles)

	# spatial check
	if (length(tiles) > raster::ncell(DEM)) stop("There are more tiles than cells in your raster.")


	cropped_DEM <- raster::crop(DEM, tiles[i, ])
	cropped_DEM_values <- raster::getValues(cropped_DEM)
	pct_non_na <- sum(is.finite(cropped_DEM_values)) / raster::ncell(cropped_DEM)
	etime <- system.time({
		if ((pct_non_na > 1/3) & !(stats::IQR(stats::na.omit(cropped_DEM_values)) < vertical_accuracy)){ 
			zeta_df <- get_zeta(cropped_DEM, raster_resolution, .mode = mode)
		}
	})
	cbind(
		data.frame(l = l, ncell = raster::ncell(tiles)),
		rbind(etime)
	)
}

df <- df %>% dplyr::mutate(child = user.child + sys.child)

p_elapsed <- ggplot2::ggplot(df, ggplot2::aes(x = l, y = elapsed)) +
ggplot2::geom_point() +
ggpubr::theme_classic2()

p_child <- ggplot2::ggplot(df, ggplot2::aes(x = l, y = child)) +
ggplot2::geom_point() +
ggpubr::theme_classic2()

CORES <- 24
CORES <- availableCores()
time_unit <- 3600

df <- df %>%
	dplyr::mutate(
		total_elapsed = ncell * elapsed / time_unit,
		total_child = ncell * child / CORES / time_unit
	)

p_total <- ggplot2::ggplot(
	df %>% reshape2::melt(measure.vars = c("total_elapsed", "total_child")), 
	ggplot2::aes(x = l, y = value, color = variable, group = variable)) +
ggplot2::geom_point() +
ggpubr::theme_classic2()
plotly::ggplotly(p_total)
