#devtools::install_version("rgrass7", version = "0.2.6")
#remotes::install_github("rsbivand/rgrass", ref="init_no_SG_53")

library("gtbox")
library("terra")
library("rgrass7")
library("sp")
library("raster")

packageVersion("rgrass7")

data("rivers_cevennes",package = "gtbox")

data("dem_cevennes",package = "gtbox")
dem = rast(dem_cevennes)
rm(dem_cevennes)

gisbase = "/usr/lib/grass78/"
#dem1 = rast("/home/vincent/jobs/gabilan_mesa3/LIDAR_GROUND.tif")
start_grass(dem,"dem",gisbase)

th_px = 2000 # stream initiation threshold in pixels
st = process_dem(dem,th_px,gisbase=gisbase)


gisbase = "/usr/lib/grass78/"
rgrass7::initGRASS(gisbase,home=tempdir(),mapset="PERMANENT",override = TRUE)


dem0 <- as(raster(dem), "SpatialGridDataFrame")
dem0 <- as(raster(dem), "SpatialGrid")

rgrass7::initGRASS(gisbase,SG=dem,home=tempdir(),mapset="PERMANENT",override = TRUE)

rgrass7::initGRASS(gisbase,SG=dem0,home=tempdir(),mapset="PERMANENT",override = TRUE)

writeRaster(dem,"/tmp/dem.tif")
