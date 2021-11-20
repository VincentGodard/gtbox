library("raster")
library("rgdal")

dem_cevennes = raster("devs/data/srtm_cevennes.tif")
outlets_cevennes = readOGR("devs/data/outlets_cevennes.gpkg")
profile_cevennes = readOGR("devs/data/profile_cevennes.gpkg")


#export
usethis::use_data(dem_cevennes,profile_cevennes,outlets_cevennes, internal = FALSE, overwrite = TRUE)

