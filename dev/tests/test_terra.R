library("gtbox")

gisbase = "/usr/lib/grass78/"
th_px = 2000 # stream initiation threshold in pixels

data("dem_cevennes",package = "gtbox")
dem = rast(dem_cevennes)
rm(dem_cevennes)

precip = dem/dem


# compute stack and network (,precip=precip)
st = process_dem(dem,th_px,gisBase=gisbase,to_net=FALSE)

net = get_network(st,gisbase,clip=TRUE)


#tmp = get_next(st$st_id,st$dir)

slope <- terrain(st$z, "slope", unit="radians")
aspect <- terrain(st$z, "aspect", unit="radians")
hill <- shade(slope, aspect, 40, 270)

plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(st$z, col=rainbow(25, alpha=0.35), add=TRUE)
lines(net,lwd=net$strahler,col="blue")

# outlets
out = get_outlets(st,strahler = 3)
bas = get_basins(st,out,gisbase)
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(st$z, col=rainbow(25, alpha=0.35), add=TRUE)
lines(net,lwd=net$strahler,col="blue")
lines(bas,lwd=2)
points(out)

out = get_outlets(st,large = 10e4)
bas = get_basins(st,out,gisbase)
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(st$z, col=rainbow(25, alpha=0.35), add=TRUE)
lines(bas)
points(out)

out = get_outlets(st,elevation = 600)
bas = get_basins(st,out,gisbase)
plot(hill, col=grey(0:100/100), legend=FALSE, mar=c(2,2,1,4))
plot(st$z, col=rainbow(25, alpha=0.35), add=TRUE)
lines(bas)
points(out)


