library("gtbox")
library("terra")

isbase = "/usr/lib/grass78/"
th_px = 2000 # stream initiation threshold in pixels

data("dem_cevennes",package = "gtbox")
dem0 = dem_cevennes
rm(dem_cevennes)

dem <- terra::project(dem0, crs="EPSG:32631",res=30,method="bilinear")
dem = trim_na(dem) # we trim to get a nice rectangular DEM with no NA on the sides
plot(dem)
rm(dem0)

# compute stack and network
st = process_dem(dem,th_px,gisBase=gisbase)
