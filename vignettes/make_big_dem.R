library(raster)

root_dir <- 'F:/hguillon/research/'
dem_dir <- 'data/california-rivers/gis-files/NED' 
DEM <- raster::raster(file.path(root_dir, dem_dir, "CA_DEM.grd"))

getpol <- function(i, pts, dl = 1000,.crs = crs(DEM)){
	x_min <- pts@coords[i,1] - dl
	x_max <- pts@coords[i,1] + dl
	y_min <- pts@coords[i,2] - dl
	y_max <- pts@coords[i,2] + dl
	coords = matrix(c(x_min, y_min,
	               x_min, y_max,
	               x_max, y_max,
	               x_max, y_min,
	               x_min, y_min), 
	             ncol = 2, byrow = TRUE)
	p <-  Polygon(coords)
	sp1 <-  SpatialPolygons(list(Polygons(list(p), ID = "a")), proj4string=.crs)
	return(sp1)
}

get_dem <- function(region, f){
	if (region == "gabilan_mesa"){
		lat <- 35.9287904
		lon <- -120.8124603
	} else if (region == "yosemite"){
		lat <- 37.7541922
		lon <- -119.5437924
	} else if (region == "modoc"){
		lat <- 41.6503415
		lon <- -121.0992454
	} else if (region == "oregon"){
		lat <- 47.40
		lon <- -125.45
	} else if (region == "marble"){
		lat <- 36.7788543
		lon <- -110.767441
	}

	lonlat <- cbind(lon,lat)
	crdref <- CRS('+proj=longlat +datum=WGS84')
	pts <- SpatialPoints(lonlat, proj4string = crdref)
	pts <- spTransform(pts, crs(DEM))
	r <- res(DEM)[1]
	.dl <- f * r 
	sf_pol <- getpol(1, pts, dl = .dl)
	dem.c <- raster::crop(DEM, sf_pol)
	plot(dem.c)
	writeRaster(dem.c, paste0(region, ".tif"), overwrite = TRUE)
}

# for (region in c("gabilan_mesa", "yosemite", "modoc")) get_dem(region, 1000)


dem_dir <- "data/california-rivers/gis-files/nw_pacific_crm-v1"
DEM <- raster::raster(file.path(root_dir, dem_dir, "nw_pacific_crm_v1.asc"))
get_dem("oregon", 300)


dlx = 10e3/10/2 * res(DEM)[1] 
dly = 25e3/10/2 * res(DEM)[2]
i = 1
x_min <- pts@coords[i,1] - dlx
x_max <- pts@coords[i,1] + dlx
y_min <- pts@coords[i,2] - dly
y_max <- pts@coords[i,2] + dly
coords = matrix(c(x_min, y_min,
               x_min, y_max,
               x_max, y_max,
               x_max, y_min,
               x_min, y_min), 
             ncol = 2, byrow = TRUE)
p <-  Polygon(coords)
sp1 <-  SpatialPolygons(list(Polygons(list(p), ID = "a")), proj4string=.crs)