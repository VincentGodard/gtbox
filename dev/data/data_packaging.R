library("terra")
library("gtbox")

# dem cevennes ----
dem_cevennes= wrap(trim_na(terra::rast("dev/data/srtm_cevennes_utm.tif")))

# outlets cevennes -----
outlets_cevennes = wrap(terra::vect("dev/data/outlets_cevennes.gpkg"))

# extract streams profiles ----
gisbase = "/usr/lib/grass78/"
# get Digital Elevation Model
data("dem_cevennes",package = "gtbox")
dem = rast(dem_cevennes)
rm(dem_cevennes)

# compute accumulation, flow direction, etc ...
st = process_dem(dem,th_px=2000,gisbase=gisbase)
# get network
net = get_network(st,gisbase)
# get largest basin
out = get_outlets(st,large = 10e4)
bas = get_basins(st,out,gisbase)
#
for (i in 1:length(bas)){
  bas0 = bas[i,] #
  st0 = terra::mask(terra::crop(st,bas0),bas0) # clip and crop the stack
  ids = unique(st0$st_id)$st_id # get the stream ids located in the clipped stack
  net0 = net[net$stream%in%ids,] # subset the network according to that
  #
  id0 = net0[which.max(net0$flow_accum)]$stream # get id of stream from which we start
  ids0 = get_streams(net0,id0,mode="up") # get upstream ids
  #
  df0 = as.data.frame(st0,xy=TRUE)
  df0 = df0[df0$st_id%in%ids0,]
  df0$river = i
  if(i==1){rivers_cevennes=df0}else{rivers_cevennes=rbind(rivers_cevennes,df0)}
}


#export
usethis::use_data(dem_cevennes,outlets_cevennes,rivers_cevennes,internal = FALSE, overwrite = TRUE)

