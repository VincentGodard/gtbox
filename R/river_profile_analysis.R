
#' Smooth channel profile by removing artifacts
#'
#' @param distance distance along channel network (m)
#' @param elevation elevation along channel network (m)
#' @param span span parameter for the loess fit used to detrend the profile (fraction of the total points used for the local fits)
#' @param radius radius of the rolling circle used to smooth the profile (m)
#'
#' @return a vector (same length as the inputs) with corrected elevation values along the river
#' @export
#'
#' @examples
#' # get river data
#' # see DEM processing vignette for example of the extraction of such data from a DEM
#' data("rivers_cevennes",package = "gtbox")
#'
#' # get DEM (for plotting)
#' data("dem_cevennes",package = "gtbox")
#' dem = rast(dem_cevennes)
#' rm(dem_cevennes)
#'
#' df = rivers_cevennes[rivers_cevennes$river==1,] # select a river
#' df = df[order(df$dist),] # order data according to distance along river
#' df$area = df$acc*30^2 # compute drainage area
#' df$chi = profile_chi(df$area,df$dist,aref=1000,mn=0.5) # compute chi
#' df$z2 = profile_smooth(df$dist,df$z,span=0.1,radius=10e3) # smooth profile
#' df$z3 = profile_gaussian(df$dist,df$z2,span=0.1,sigma=500) # filter profile
#' df$gradient = profile_gradient(df$dist,df$z3,window=200) # compute gradient
#'
#' par(mfrow=c(2,2))
#' plot(dem)
#' lines(df$x,df$y,lwd=3,col="blue")
#' par(mar = c(4,4,0.5,0.5))
#' plot(df$dist,df$z,type="l",xlab="Distance along river (m)",ylab="Elevation (m)")
#' lines(df$dist,df$z2,col="firebrick")
#' lines(df$dist,df$z3,col="dodgerblue")
#' legend("topleft",c("Raw","Smoothed","Filtered"),lty=1,col=c("black","firebrick","dodgerblue"))
#' plot(df$dist,df$gradient,type="l",xlab="Distance along river (m)",ylab="Gradient (m/m)")
#' plot(df$chi,df$z,type="l",xlab="Chi (m)",ylab="Elevation (m)")
profile_smooth<-function(distance,elevation,span,radius){
  tmp=as.data.frame(cbind(distance,elevation))
  colnames(tmp)<-c("distance","elevation")
  rk=rank(tmp$distance) # retain memory of the initial order in case the profile is not ordered
  tmp=tmp[with(tmp, order(distance)), ] #order data
  fit1=loess(tmp$elevation ~ tmp$distance,span=span)
  tmp$residu=tmp$elevation - predict(fit1)
  tmp$flag=0
  #for(i in 1:nrow(tmp)) {
  for(i in 2:(nrow(tmp)-1)) {
    xc = tmp$distance[i]
    yc = tmp$residu[i]-radius
    tmp$dist = sqrt((tmp$distance-xc)^2+(tmp$residu-yc)^2)-radius
    tmp2 = subset(tmp,dist<=0)
    if (dim(tmp2)[1]>1) {tmp$flag[i]=1}
  }
  tmp3 = subset(tmp,flag==0)
  tmp$alt2 = approx(tmp3$distance,tmp3$residu,xout=tmp$distance,rule=2)$y + predict(fit1)
  return(tmp$alt2[rk])
}
# other approach for detrending
#x = seq(min(tmp$distance),max(tmp$distance),length.out = npt)
#   y = approx(tmp$distance,tmp$elevation,x)$y
#   y2 = approx(x,y,tmp$distance)$y
#   tmp$residu=tmp$elevation - y2






#' profile_gaussian
#'
#' Perform Gaussian filtering of river profile
#'
#' @param distance distance along channel network (m)
#' @param elevation elevation along channel network (m)
#' @param span span parameter for the loess fit used to detrend the profile
#' @param sigma standard deviation of the Gaussian window (m)
#'
#' @return a vector (same length as the inputs) with filtered elevation values along the river
#' @export
#'
#' @examples
#' # get river data
#' # see DEM processing vignette for example of the extraction of such data from a DEM
#' data("rivers_cevennes",package = "gtbox")
#'
#' # get DEM (for plotting)
#' data("dem_cevennes",package = "gtbox")
#' dem = rast(dem_cevennes)
#' rm(dem_cevennes)
#'
#' df = rivers_cevennes[rivers_cevennes$river==1,] # select a river
#' df = df[order(df$dist),] # order data according to distance along river
#' df$area = df$acc*30^2 # compute drainage area
#' df$chi = profile_chi(df$area,df$dist,aref=1000,mn=0.5) # compute chi
#' df$z2 = profile_smooth(df$dist,df$z,span=0.1,radius=10e3) # smooth profile
#' df$z3 = profile_gaussian(df$dist,df$z2,span=0.1,sigma=500) # filter profile
#' df$gradient = profile_gradient(df$dist,df$z3,window=200) # compute gradient
#'
#' par(mfrow=c(2,2))
#' plot(dem)
#' lines(df$x,df$y,lwd=3,col="blue")
#' par(mar = c(4,4,0.5,0.5))
#' plot(df$dist,df$z,type="l",xlab="Distance along river (m)",ylab="Elevation (m)")
#' lines(df$dist,df$z2,col="firebrick")
#' lines(df$dist,df$z3,col="dodgerblue")
#' legend("topleft",c("Raw","Smoothed","Filtered"),lty=1,col=c("black","firebrick","dodgerblue"))
#' plot(df$dist,df$gradient,type="l",xlab="Distance along river (m)",ylab="Gradient (m/m)")
#' plot(df$chi,df$z,type="l",xlab="Chi (m)",ylab="Elevation (m)")
profile_gaussian<-function(distance,elevation,span,sigma){
  tmp=as.data.frame(cbind(distance,elevation))
  colnames(tmp)<-c("distance","elevation")
  rk=rank(tmp$distance) # retain memory of the initial order
  tmp=tmp[with(tmp, order(distance)), ] #order data
  #  fit1=lm(tmp$elevation ~ poly(tmp$distance,order))
  fit1=loess(tmp$elevation ~ tmp$distance,span=span)
  tmp$residu=tmp$elevation - predict(fit1)
  for(i in 1:nrow(tmp)) {
    abs0=tmp$distance[i]
    tmp$gaussian=1/(sqrt(2*pi)*sigma)*exp(-1*(tmp$distance-abs0)^2/(2*sigma^2))
    tmp$filt[i]=sum(tmp$residu*tmp$gaussian)/sum(tmp$gaussian)
  }
  tmp$alt2=tmp$filt + predict(fit1)
  return(tmp$alt2[rk])
}




#' profile_gradient
#'
#' Compute river profile gradient over a given window
#'
#' @param distance distance along channel network (m)
#' @param elevation elevation along channel network (m)
#' @param window size of the window for the computation of the gradient
#'
#' @return a vector (same length as the inputs) with gradient values along the river
#' @export
#'
#' @examples
#' # get river data
#' # see DEM processing vignette for example of the extraction of such data from a DEM
#' data("rivers_cevennes",package = "gtbox")
#'
#' # get DEM (for plotting)
#' data("dem_cevennes",package = "gtbox")
#' dem = rast(dem_cevennes)
#' rm(dem_cevennes)
#'
#' df = rivers_cevennes[rivers_cevennes$river==1,] # select a river
#' df = df[order(df$dist),] # order data according to distance along river
#' df$area = df$acc*30^2 # compute drainage area
#' df$chi = profile_chi(df$area,df$dist,aref=1000,mn=0.5) # compute chi
#' df$z2 = profile_smooth(df$dist,df$z,span=0.1,radius=10e3) # smooth profile
#' df$z3 = profile_gaussian(df$dist,df$z2,span=0.1,sigma=500) # filter profile
#' df$gradient = profile_gradient(df$dist,df$z3,window=200) # compute gradient
#'
#' par(mfrow=c(2,2))
#' plot(dem)
#' lines(df$x,df$y,lwd=3,col="blue")
#' par(mar = c(4,4,0.5,0.5))
#' plot(df$dist,df$z,type="l",xlab="Distance along river (m)",ylab="Elevation (m)")
#' lines(df$dist,df$z2,col="firebrick")
#' lines(df$dist,df$z3,col="dodgerblue")
#' legend("topleft",c("Raw","Smoothed","Filtered"),lty=1,col=c("black","firebrick","dodgerblue"))
#' plot(df$dist,df$gradient,type="l",xlab="Distance along river (m)",ylab="Gradient (m/m)")
#' plot(df$chi,df$z,type="l",xlab="Chi (m)",ylab="Elevation (m)")
profile_gradient<-function(distance,elevation,window){
  tmp=as.data.frame(cbind(distance,elevation))
  colnames(tmp)<-c("distance","elevation")
  rk=rank(tmp$distance) # retain memory of the initial order
  tmp=tmp[with(tmp, order(distance)), ] #order data
  tmp$grad=NA
  for(i in 1:nrow(tmp)) {
    abs0=tmp$distance[i]
    tmp2=subset(tmp,abs(distance-abs0)<=window/2)
    fit=lm(tmp2$elevation~tmp2$distance)
    tmp$grad[i]=abs(fit$coefficients[2])
  }
  return(tmp$grad[rk])
}


#' Compute chi parameter along a single river
#'
#' @param area drainage area (m2)
#' @param distance distance along stream (m)
#' @param aref reference drainage area (m2)
#' @param mn m/n ratio
#'
#' @return a vector (same length as the inputs) with chi values along the river
#' @export
#'
#' @examples
#' # get river data
#' # see DEM processing vignette for example of the extraction of such data from a DEM
#' data("rivers_cevennes",package = "gtbox")
#'
#' # get DEM (for plotting)
#' data("dem_cevennes",package = "gtbox")
#' dem = rast(dem_cevennes)
#' rm(dem_cevennes)
#'
#' df = rivers_cevennes[rivers_cevennes$river==1,] # select a river
#' df = df[order(df$dist),] # order data according to distance along river
#' df$area = df$acc*30^2 # compute drainage area
#' df$chi = profile_chi(df$area,df$dist,aref=1000,mn=0.5) # compute chi
#' df$z2 = profile_smooth(df$dist,df$z,span=0.1,radius=10e3) # smooth profile
#' df$z3 = profile_gaussian(df$dist,df$z2,span=0.1,sigma=500) # filter profile
#' df$gradient = profile_gradient(df$dist,df$z3,window=200) # compute gradient
#'
#' par(mfrow=c(2,2))
#' plot(dem)
#' lines(df$x,df$y,lwd=3,col="blue")
#' par(mar = c(4,4,0.5,0.5))
#' plot(df$dist,df$z,type="l",xlab="Distance along river (m)",ylab="Elevation (m)")
#' lines(df$dist,df$z2,col="firebrick")
#' lines(df$dist,df$z3,col="dodgerblue")
#' legend("topleft",c("Raw","Smoothed","Filtered"),lty=1,col=c("black","firebrick","dodgerblue"))
#' plot(df$dist,df$gradient,type="l",xlab="Distance along river (m)",ylab="Gradient (m/m)")
#' plot(df$chi,df$z,type="l",xlab="Chi (m)",ylab="Elevation (m)")
profile_chi<-function(area,distance,aref,mn){
  if(sum(area<0)){warning("chi calculation : some areas are negative, we use the absolute value")}
  tmp=as.data.frame(cbind(area,distance))
  colnames(tmp)<-c("area","distance")
  tmp$param=(aref/(abs(tmp$area)))^mn
  rk=rank(tmp$distance) # retain memory of the initial order
  tmp=tmp[with(tmp, order(distance)), ] #order data
  tmp$chi=pracma::cumtrapz(tmp$distance,tmp$param)
  tmp$chi=tmp$chi[rk]
  return(tmp$chi)
}





