


#' Project points along a section.
#'
#' Return a dataframe with two columns (same row order as the input)
#' `distance_along` : distance along  the profile of the projections
#' `distance_to` : perpendicular distance to the profile of the points
#' All input should be in cartesian coordinates with the same units
#'
#' @param x X coordinates of the points to be projected (vector)
#' @param y Y coordinates of the points to be projected (vector)
#' @param X1 X coordinate of starting point defining the projection line
#' @param Y1 Y coordinate of starting point defining the projection line
#' @param X2 X coordinate of end point defining the projection line
#' @param Y2 Y coordinate of end point defining the projection line
#' @param width half width of the projection swath (optional)
#' @param inc number of bins along the profile for the stats (optional)
#' @param value parameter for which the statistics should be returned (optional, vector of same length and order as `x` and `y`)
#'
#' @return A dataframe with the distance along and to the sections of the input points. If `width`, `inc` and `value` are supplied a dataframe containing the binned (`inc` bins) statistics of `value` along the section is returned instead
#' @export
#'
#' @examples
project_points<-function(x,y,X1,Y1,X2,Y2,width=NULL,inc=NULL,value=NULL){
  if (length(x) != length(y)) {stop("x and y should have same length")}
  if ((length(X1) != 1)) {stop("X1 should be a scalar")}
  if ((length(Y1) != 1)) {stop("Y1 should be a scalar")}
  if ((length(X2) != 1)) {stop("X2 should be a scalar")}
  if ((length(Y2) != 1)) {stop("Y2 should be a scalar")}
  #if (!is.null(width) | !is.null(inc) | !is.null(value)) {stop("all three of width, inc and value should be supplied for binned statistics along the section")}
  # notations http://www.sunshine2k.de/coding/java/PointOnLine/PointOnLine.html#step5
  e1_x = X2-X1
  e1_y = Y2-Y1
  e2_x = x-X1
  e2_y = y-Y1
  e1_n = sqrt((e1_x)^2+(e1_y)^2) #norme
  e2_n = sqrt((e2_x)^2+(e2_y)^2) #norme
  ps_e1e2 = e1_x*e2_x+e1_y*e2_y #produit scalaire
  v1pp_n = ps_e1e2/e1_n
  xp = X1+v1pp_n/e1_n*e1_x
  yp = Y1+v1pp_n/e1_n*e1_y
  l = ((xp-X1)*e1_x+(yp-Y1)*e1_y)/e1_n
  w = ((X2-X1)*(y-Y1)-(x-X1)*(Y2-Y1))/sqrt((X2-X1)^2+(Y2-Y1)^2)
  if (is.null(width)){
    res = cbind(l,w)
    colnames(res) <- c("distance_along","distance_to")
  return(as.data.frame(res))
  } else {
    index = abs(w)<width
    len1 = sqrt((X1-X2)^2+(Y1-Y2)^2)
    bk = seq(0,len1,length.out = inc+1)
    bn = fields::stats.bin(l[index],value[index],breaks=bk)
    res = cbind(bn$centers,as.data.frame(t(bn$stats)))
    colnames(res)<-c("distance","n","mean","sd","min","q1","median","q3","max","missing_values")
    return(res)
  }

}



#' Project raster data along a swath profile
#'
#' The profile is defined by its end points and its width.
#' The results are returned as summary statistics for equally spaced bins along the profile.
#'
#' @param rast input raster
#' @param X1 X coordinate of starting point defining the projection line
#' @param Y1 Y coordinate of starting point defining the projection line
#' @param X2 X coordinate of end point defining the projection line
#' @param Y2 Y coordinate of end point defining the projection line
#' @param width one-side width of the swath profile
#' @param inc number of bins along the profile for the stats
#'
#' @return a data frame with statistics computed for each increment along the profile
#' @export
#'
#' @examples
swath_profile <- function(rast,X1,Y1,X2,Y2,width,inc){
  # notations http://www.sunshine2k.de/coding/java/PointOnLine/PointOnLine.html#step5
  # create line
  ln1 = terra::vect(matrix(c(x1,y1,x2,y2),ncol =2,byrow=T),
                    type="lines",crs=crs(rast))
  # create buffer
  bf1 = terra::vect(sf::st_buffer(sf::st_as_sf(l1),dist=width,
                                 endCapStyle = "FLAT",joinStyle = "ROUND"))

  tmp = na.omit(terra::as.data.frame(terra::mask(rast,bf1),xy=TRUE))
  colnames(tmp) <- c("x","y","val")
  # project
  e1_x = X2-X1
  e1_y = Y2-Y1
  e2_x = tmp$x-X1
  e2_y = tmp$y-Y1
  e1_n = sqrt((e1_x)^2+(e1_y)^2) #norme
  e2_n = sqrt((e2_x)^2+(e2_y)^2) #norme
  ps_e1e2 = e1_x*e2_x+e1_y*e2_y #produit scalaire
  v1pp_n = ps_e1e2/e1_n
  xp = X1+v1pp_n/e1_n*e1_x
  yp = Y1+v1pp_n/e1_n*e1_y
  tmp$l = ((xp-X1)*e1_x+(yp-Y1)*e1_y)/e1_n
  # stats
  len1 = sqrt((X1-X2)^2+(Y1-Y2)^2)
  bk = seq(0,len1,length.out = inc+1)
  bn = fields::stats.bin(tmp$l,tmp$val,breaks=bk)
  res = cbind(bn$centers,as.data.frame(t(bn$stats)))
  colnames(res)<-c("distance","n","mean","sd","min","q1","median","q3","max","missing_values")
  return(res)
}


