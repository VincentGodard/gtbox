library("raster")
library("rgdal")

dem_cevennes = raster::raster("dev/data/srtm_cevennes.tif")
outlets_cevennes = rgdal::readOGR("dev/data/outlets_cevennes.gpkg")
profile_cevennes = rgdal::readOGR("dev/data/profile_cevennes.gpkg")


#export
usethis::use_data(dem_cevennes,profile_cevennes,outlets_cevennes, internal = FALSE, overwrite = TRUE)

