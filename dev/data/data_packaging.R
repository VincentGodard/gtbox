library("terra")
library("gtbox")

tmp = terra::rast("dev/data/srtm_cevennes_utm.tif")
dem_cevennes = trim_na(tmp)

#outlets_cevennes = rgdal::readOGR("dev/data/outlets_cevennes.gpkg")
#profile_cevennes = rgdal::readOGR("dev/data/profile_cevennes.gpkg")


#export
#usethis::use_data(dem_cevennes,profile_cevennes,outlets_cevennes, internal = FALSE, overwrite = TRUE)
usethis::use_data(dem_cevennes, internal = FALSE, overwrite = TRUE)

