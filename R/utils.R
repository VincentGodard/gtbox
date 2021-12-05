
#' Initiate a grass throwaway session with proper region definition
#'
#'
#' @param rast  Raster to import (RasterLayer)
#' @param name for the raster in the location
#' @param gisBase The directory path to GRASS binaries and libraries, containing bin and lib subdirectories among others
#'
#' @return
start_grass<-function(rast,name,gisBase){

  rgrass7::use_sp()

  Sys.setenv("GRASS_VERBOSE"=0)

  # create Grass session
  #grid =  SpatialGrid(GridTopology(c(extent(rast)[1]+0.5*res(rast)[1],extent(rast)[3]+0.5*res(rast)[2]),res(rast),c(rast@ncols,rast@nrows)))
  #rgrass7::initGRASS(gisBase,SG=grid,override = TRUE)

  # init
  #rgrass7::initGRASS(gisBase,SG=as(rast, 'SpatialGrid'),override = TRUE)

  # start grass
  rgrass7::initGRASS(gisBase,override = TRUE)
  # shift to PERMANENT
  rgrass7::execGRASS("g.mapset",parameters=list(mapset="PERMANENT"))
  # set projection
  crs = as.character(crs(rast))
  rgrass7::execGRASS("g.proj",flags=c("c"),
                     parameters=list(proj4=crs))
  # set location properties
  ext=extent(rast)
  rgrass7::execGRASS("g.region",parameters=list(rows=dim(rast)[1],
                                                cols=dim(rast)[2],
                                                ewres=as.character(res(rast)[1]),
                                                nsres=as.character(res(rast)[2]),
                                                w=as.character(ext[1]),
                                                e=as.character(ext[2]),
                                                s=as.character(ext[3]),
                                                n=as.character(ext[4])))

  # write raster into location
  write_raster_to_grass(rast,name)

}


write_raster_to_grass <- function(rast,name){
  if(!(is.character(name))){stop("second argument must be a character string")}
  if(class(rast)[1]!="SpatRaster"){stop("first argument must be a SpatRaster")}
  rgrass7::writeRAST(as(raster::raster(rast), 'SpatialGridDataFrame'),vname=c(name))
}

read_raster_from_grass <- function(name){
  if(!(is.character(name))){stop("argument must be character string")}
  rast = terra::rast(raster::raster(rgrass7::readRAST(c(name))))
  return(rast)
}

write_vector_to_grass <- function(vect,name){
  if(!(is.character(name))){stop("second argument must be a character string")}
  if(class(vect)[1]!="SpatVector"){stop("first argument must be a SpatVector")}
  rgrass7::writeVECT(as(vect,"Spatial"),vname=c(name),v.in.ogr_flags=c("overwrite","o"))
}

read_vector_from_grass <- function(name){
  if(!(is.character(name))){stop("argument must be character string")}
  vector = terra::vect(rgrass7::readVECT(c(name)))
  return(vector)
}






#' Dowload a STRM raster from srtm.csi.cgiar.org
#'
#' Compute for each pixel of a raster the value of the following pixel according to flow direction
#'
#' @param sp  an object for which a spatial extent can be extracted
#' @param buf a buffering distance (in degrees) to expand the extent
#'
#' @return
#' @export
#' @examples
get_srtm <- function(sp,buf=0.1) {
  ext = extend(extent(sp),buf)
  lon1 = ext[1]
  lon2 = ext[2]
  lat1 = ext[3]
  lat2 = ext[4]
  stopifnot(lon1 >= -180 & lon1 <= 180 & lon2 >= -180 & lon2 <= 180)
  stopifnot(lat1 >= -60 & lat1 <= 60 & lat2 >= -60 & lat2 <= 60)
  rs <- raster(nrows=24, ncols=72, xmn=-180, xmx=180, ymn=-60, ymx=60 )
  rowTile1 <- rowFromY(rs, lat2)
  colTile1 <- colFromX(rs, lon1)
  rowTile2 <- rowFromY(rs, lat2)
  colTile2 <- colFromX(rs, lon2)
  rowTile3 <- rowFromY(rs, lat1)
  colTile3 <- colFromX(rs, lon1)
  rowTile4 <- rowFromY(rs, lat1)
  colTile4 <- colFromX(rs, lon2)
  if (rowTile1 < 10) { rowTile1 <- paste('0', rowTile1, sep='') }
  if (colTile1 < 10) { colTile1 <- paste('0', colTile1, sep='') }
  if (rowTile2 < 10) { rowTile2 <- paste('0', rowTile2, sep='') }
  if (colTile2 < 10) { colTile2 <- paste('0', colTile2, sep='') }
  if (rowTile3 < 10) { rowTile3 <- paste('0', rowTile3, sep='') }
  if (colTile3 < 10) { colTile3 <- paste('0', colTile3, sep='') }
  if (rowTile4 < 10) { rowTile4 <- paste('0', rowTile4, sep='') }
  if (colTile4 < 10) { colTile4 <- paste('0', colTile4, sep='') }
  f1 <- paste('srtm_', colTile1, '_', rowTile1, sep="")
  f2 <- paste('srtm_', colTile2, '_', rowTile2, sep="")
  f3 <- paste('srtm_', colTile3, '_', rowTile3, sep="")
  f4 <- paste('srtm_', colTile4, '_', rowTile4, sep="")
  lf = unique(c(f1,f2,f3,f4))
  lr = list()
  for (i in 1:length(lf)){
    f = lf[[i]]
    theurl <- paste("https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/", f, ".zip", sep="")
    utils::download.file(url=theurl, destfile=paste("/tmp/",f,".zip",sep=""), method="auto", quiet = FALSE, mode = "wb", cacheOK = TRUE)
    system(paste("unzip -o -d /tmp /tmp/",f,".zip",sep=""))
    lr[[i]] = raster(paste("/tmp/",f,".tif",sep=""))
  }
  if (length(lr)==1){
    return(raster::crop(lr[[1]],ext))
  } else {
    rs = lr[[1]]
    for (i in 2:length(lr)){rs=raster::merge(lr[[i]],rs)}
    return(raster::crop(rs,ext))
  }
}



#' Trim a raster to remove NA values on the periphery
#'
#' The trim function from the raster package ony remove cols and rows full of NA
#' Useful to get a nicely formatted raster grid after reprojection
#'
#' @param rast input raster
#'
#' @return a trimmed raster
#' @export
#'
#' @examples
trim_na <- function(rst0){
  # rst0 = terra::trim(rst0)
  # #
  # s1 = rst0[terra::nrow(rst0),]
  # a1 = (1:terra::ncol(rst0))[!is.na(s1)]
  # s3 = rst0[1,]
  # a3 = (1:terra::ncol(rst0))[!is.na(s3)]
  # s2 = rst0[,1]
  # a2 = (1:terra::nrow(rst0))[!is.na(s2)]
  # s4 = rst0[,terra::ncol(rst0)]
  # a4 = (1:terra::nrow(rst0))[!is.na(s4)]
  # #
  # lc = range(a1,a3)
  # lr = range(a2,a4)
  # #
  # e = terra::ext(rst0)
  # xmin = e[1] + lc[1]*terra::res(rst0)[1]
  # xmax = e[2] - (terra::ncol(rst0)-lc[2])*terra::res(rst0)[1]
  # ymin = e[3] + (terra::nrow(rst0)-lr[2])*terra::res(rst0)[2]
  # ymax = e[4] - lr[1]*terra::res(rst0)[2]
  #
  # rst2 = terra::crop(rst0,terra::ext(c(xmin,xmax,ymin,ymax)))
  # return(rst2)
}




#' Get value of next pixel following flow direction.
#'
#' Compute for each pixel of a raster the value of the following pixel according to flow direction
#'
#' @param rast  Raster (RasterLayer)
#' @param dir Flow direction raster (RasterLayer)
#'
#' @return
#'
#' @examples
get_next<-function(rast,dir){
  rast0 = rast
  dir = raster::as.matrix(dir)
  rast = raster::as.matrix(rast)
  nc = ncol(dir)
  nr = nrow(dir)
  res = matrix(NA,nrow = nrow(dir), ncol = ncol(dir))
  # flow toward E, offset input matrix toward W
  res <- ifelse(abs(dir)==8,cbind(rast[1:nr,2:nc],rep(NA,nr)),res)
  # flow toward W, offset input matrix toward E
  res <- ifelse(abs(dir)==4,cbind(rep(NA,nr),rast[1:nr,1:(nc-1)]),res)
  # flow toward N, offset input matrix toward S
  res <- ifelse(abs(dir)==2,rbind(rep(NA,nc),rast[1:(nr-1),1:nc]),res)
  # flow toward S, offset input matrix toward N
  res <- ifelse(abs(dir)==6,rbind(rast[2:nr,1:nc],rep(NA,nc)),res)
  # flow toward the NE, offset input matrix toward the SW
  res <- ifelse(abs(dir)==1,cbind(rbind(rep(NA,nc-1),rast[1:(nr-1),2:nc]),rep(NA,nr)),res)
  # flow toward the NW, offset input matrix toward the SE
  res <- ifelse(abs(dir)==3,cbind(rep(NA,nr),rbind(rep(NA,nc-1),rast[1:(nr-1),1:(nc-1)])),res)
  # flow toward the SE, offset input matrix toward the NW
  res <- ifelse(abs(dir)==5,cbind(rep(NA,nr),rbind(rast[2:nr,1:(nc-1)],rep(NA,nc-1))),res)
  # flow toward the SW, offset input matrix toward the NE
  res <- ifelse(abs(dir)==7,cbind(rbind(rast[2:nr,2:nc],rep(NA,nc-1)),rep(NA,nr)),res)
  return(raster(res,template=rast0))
}





