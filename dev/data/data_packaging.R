library("terra")
library("gtbox")

dem_cevennes= wrap(trim_na(terra::rast("dev/data/srtm_cevennes_utm.tif")))

outlets_cevennes = wrap(terra::vect("dev/data/outlets_cevennes.gpkg"))
profile_cevennes = wrap(terra::vect("dev/data/profile_cevennes.gpkg"))


#export
usethis::use_data(dem_cevennes,profile_cevennes,outlets_cevennes, internal = FALSE, overwrite = TRUE)

