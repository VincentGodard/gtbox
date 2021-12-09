library("gtbox")
library("terra")
library("sf")


gisbase = "/usr/lib/grass78/"
th_px = 2000 # stream initiation threshold in pixels

data("dem_cevennes",package = "gtbox")
dem = rast(dem_cevennes)
rm(dem_cevennes)
writeRaster(dem,"/tmp/dem.tif")

# Swath profiles ----
x1 = 583695
y1 = 4875921
x2 = 550939
y2 = 4907255

w_buf = 5e3

l1 = terra::vect(matrix(c(x1,y1,x2,y2),ncol =2,byrow=T),
                 type="lines",crs=crs(dem))

bf = terra::vect(st_buffer(st_as_sf(l1),dist=w_buf,
                    endCapStyle = "FLAT",joinStyle = "ROUND"))
#bf$test = 0
plot(dem)
lines(l1)
lines(bf,lty=2)

data = swath_profile(dem,x1,y1,x2,y2,w_buf,100) # computation of the swath profile

# compute curvature ----
C0 =  compute_curvature(dem,win=5,gisbase)
S0 =  compute_slope(dem,win=5,gisbase)

# compute stack and network (,precip=precip) ----
#precip = dem/dem
st = process_dem(dem,th_px,gisbase=gisbase,to_net=FALSE)

net = get_network(st,gisbase,clip=TRUE)

writeRaster(st$st_id,"/tmp/st.tif")

#tmp = get_next(st$st_id,st$dir)


slope <- terrain(st$z, "slope", unit="radians")
aspect <- terrain(st$z, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)

plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(st$z, col=rainbow(25, alpha=0.35), add=TRUE)
lines(net,lwd=net$strahler,col="blue")

# outlets strahler
out = get_outlets(st,strahler = 2)
bas = get_basins(st,out,gisbase)
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(st$z, col=rainbow(25, alpha=0.35), add=TRUE)
lines(net,lwd=net$strahler,col="blue")
lines(bas,lwd=2)
points(out)

# zonal stats
library(tictoc)
tic()
tmp = compute_zonal_stats(bas,st$z,"alti",gisbase)
toc()
tic()
tmp2 =extract(st$z,bas,fun="mean")
toc()

# outlets large
out = get_outlets(st,large = 10e4)
bas = get_basins(st,out,gisbase)
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(st$z, col=rainbow(25, alpha=0.35), add=TRUE)
lines(bas)
points(out)

# crop the stack
bas = bas[1,] # largest basin
st0 = terra::mask(terra::crop(st,bas),bas) # clip and crop the stack
ids = unique(st0$st_id)$st_id # get the stream ids located in the clipped stack
net0 = net[net$stream%in%ids,] # subset the network according to that
#
id0 = net0[which.max(net0$flow_accum)]$stream # get id of stream from which we start
ids0 = get_streams(net0,722,mode="down") # get upstream ids




# outlets elevation
out = get_outlets(st,elevation = 600)
bas = get_basins(st,out,gisbase)
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(st$z, col=rainbow(25, alpha=0.35), add=TRUE)
lines(bas)
points(out)

# outlet custom
data("outlets_cevennes",package = "gtbox")
out = vect(outlets_cevennes)
rm(outlets_cevennes)
bas = get_basins(st,out,gisbase)
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(st$z, col=rainbow(25, alpha=0.35), add=TRUE)
lines(bas)
points(out)
writeVector(bas,"/tmp/bas.gpkg")
