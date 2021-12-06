library("gtbox")

gisbase = "/usr/lib/grass78/"
th_px = 2000 # stream initiation threshold in pixels

data("dem_cevennes",package = "gtbox")
dem = rast(dem_cevennes)
rm(dem_cevennes)

precip = dem/dem


# compute stack and network (,precip=precip)
st = process_dem(dem,th_px,gisBase=gisbase,to_net=FALSE)

tmp = get_next(st$st_id,st$dir)



