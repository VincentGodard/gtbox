
#' Smooth channel profile by removing artifacts
#'
#' @param distance distance along channel network (m)
#' @param elevation elevation along channel network (m)
#' @param span span parameter for the loess fit used to detrend the profile (fraction of the total points used for the local fits)
#' @param radius radius of the rolling circle used to smooth the profile (m)
#'
#' @return
#' @export
#'
#' @examples
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
#' @return
#' @export
#'
#' @examples
profile_gaussian<-function(distance,elevation,span,sigma){
  tmp=as.data.frame(cbind(distance,elevation))
  colnames(tmp)<-c("distance","elevation")
  rk=rank(tmp$distance) # retain memory of the intial order
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
#' @return
#' @export
#'
#' @examples
profile_gradient<-function(distance,elevation,window){
  tmp=as.data.frame(cbind(distance,elevation))
  colnames(tmp)<-c("distance","elevation")
  rk=rank(tmp$distance) # retain memory of the intial order
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







