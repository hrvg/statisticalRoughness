devtools::load_all()
library(tictoc)

# region_names <- c("fort_bragg", "yosemite", "gabilan_mesa", "modoc")
region_names <- c("fort_bragg", "yosemite")

for (.name in region_names){
	print(.name)
	rstr <- raster::raster(file.path("./inst/extdata/big_rasters/", paste0(.name, ".tif")))
	out_dir <- file.path(file.path("./inst/extdata/zeta_results"), .name)

	if(!dir.exists(file.path(out_dir))) dir.create(out_dir)

	Lmax <- min(head(dim(rstr), 2))
	spatial_scales <- get_all_R_L(Lmax, 5, len = 13)$allL %>% head(-1)
	spatial_scales <- spatial_scales[spatial_scales >= 72]

	for (n in rev(spatial_scales)){
		tic()
		print(n)
		if(!file.exists(file.path(out_dir, paste0(.name, "_", n, ".tif")))){
			tiles <- raster::aggregate(rstr, fact = n)
			raster_resolution <- 10
			zeta_df <- get_zeta_df(
			  rstr, tiles, 
			  raster_resolution, vertical_accuracy = 1.87
			)
			zeta_raster <- get_zeta_raster(
			  rstr, tiles, 
			  .zeta_df = zeta_df, raster_resolution, vertical_accuracy = 1.87
			)
			raster::writeRaster(zeta_raster, file.path(out_dir, paste0(.name, "_", n, ".tif")), overwrite = TRUE)	
			write.csv(head(zeta_df, 1), file.path(out_dir, paste0("band_names.csv")), row.names = FALSE)
		}
		toc()
	}

}
