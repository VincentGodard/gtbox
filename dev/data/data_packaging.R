library("terra")
library("gtbox")

tmp = terra::rast("dev/data/srtm_cevennes_utm.tif")
dem_cevennes = trim_na(tmp)

outlets_cevennes = terra::vect("dev/data/outlets_cevennes.gpkg")
profile_cevennes = terra::vect("dev/data/profile_cevennes.gpkg")


#export
usethis::use_data(dem_cevennes,profile_cevennes,outlets_cevennes, internal = FALSE, overwrite = TRUE)

