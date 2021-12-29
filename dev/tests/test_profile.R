library("gtbox")
library("terra")

# get river data
# see DEM processing vignette for example of the extraction of such data from a DEM
data("rivers_cevennes",package = "gtbox")

# get DEM (for plotting)
data("dem_cevennes",package = "gtbox")
dem = rast(dem_cevennes)
rm(dem_cevennes)

df = rivers_cevennes[rivers_cevennes$river==1,] # select a river
df = df[order(df$dist),] # order data according to distance along river
df$area = df$acc*30^2 # compute drainage area
df$chi = profile_chi(df$area,df$dist,aref=1000,mn=0.5) # compute chi
df$z2 = profile_smooth(df$dist,df$z,span=0.1,radius=10e3) # smooth profile
df$z3 = profile_gaussian(df$dist,df$z2,span=0.1,sigma=500) # filter profile
df$gradient = profile_gradient(df$dist,df$z3,window=200) # compute gradient

par(mfrow=c(2,2))
plot(dem)
lines(df$x,df$y,lwd=3,col="blue")
par(mar = c(4,4,0.5,0.5))
plot(df$dist,df$z,type="l",xlab="Distance along river (m)",ylab="Elevation (m)")
lines(df$dist,df$z2,col="firebrick")
lines(df$dist,df$z3,col="dodgerblue")
legend("topleft",c("Raw","Smoothed","Filtered"),lty=1,col=c("black","firebrick","dodgerblue"))
plot(df$dist,df$gradient,type="l",xlab="Distance along river (m)",ylab="Gradient (m/m)")
plot(df$chi,df$z,type="l",xlab="Chi (m)",ylab="Elevation (m)")

# # get profile
# bas0 = bas[1,]
# st0 = terra::mask(terra::crop(st,bas0),bas0) # clip and crop the SpatRasters
# ids = unique(st0$st_id)$st_id # get the stream ids located in the basin
# net0 = net[net$stream%in%ids,] # subset the network SpatVector according to that
# #
# id0 = net0[which.max(net0$flow_accum)]$stream # get id of stream segment from which we start (largest accumulation)
# ids0 = get_streams(net0,id0,mode="up") # get upstream ids
# #
# df0 = as.data.frame(st0,xy=TRUE) # convert the SpatRaster to data frame
# df0 = df0[df0$st_id%in%ids0,] # select the relevant stream segments according to their id
# df0 = df0[order(df0$dist),] # order according to distance
# df0$area = df0$acc*res(st0)[1]*res(st0)[2] # compute drainage area
# df0$chi = profile_chi(df0$area,df0$dist,aref=1000,mn=0.5) # compute chi
# #
# plot(df0$dist,df0$z,type="l",xlab="Distance (m)",ylab="Elevation (m)")
# plot(df0$chi,df0$z,type="l",xlab="Chi (m)",ylab="Elevation (m)")


#
for (i in 1:length(bas)){
bas0 = bas[i,] # largest basin
st0 = terra::mask(terra::crop(st,bas0),bas0) # clip and crop the stack
ids = unique(st0$st_id)$st_id # get the stream ids located in the clipped stack
net0 = net[net$stream%in%ids,] # subset the network according to that
#
id0 = net0[which.max(net0$flow_accum)]$stream # get id of stream from which we start
ids0 = get_streams(net0,id0,mode="up") # get upstream ids
#
df0 = as.data.frame(st0,xy=TRUE)
df0 = df0[df0$st_id%in%ids0,]
df0$drain = i
if(i==1){df=df0}else{df=rbind(df,df0)}
df0 = df0[order(df0$dist),]
df0$area = df0$acc*res(st0)[1]*res(st0)[2]
df0$chi = profile_chi(df0$area,df0$dist,aref=1000,mn=0.5)
plot(df0$dist,df0$z,type="l",main=i)
plot(df0$chi,df0$z,type="l",main=i)
}





