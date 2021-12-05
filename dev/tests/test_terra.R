library("gtbox")

gisbase = "/usr/lib/grass78/"
th_px = 2000 # stream initiation threshold in pixels

data("dem_cevennes",package = "gtbox")
dem = dem_cevennes
rm(dem_cevennes)



# compute stack and network
st = process_dem(dem,th_px,gisBase=gisbase)
