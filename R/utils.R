
#' Initiate a grass throwaway session with proper region definition
#'
#'
#' @param rast  Raster to import (SpatRaster)
#' @param name for the raster in the location
#' @param gisbase The directory path to GRASS binaries and libraries, containing bin and lib subdirectories among others
#'
#' @return
start_grass<-function(rast,name,gisbase){
  if(class(rast)[1]!="SpatRaster"){stop("first argument must be a SpatRaster")}
  if(!(is.character(name))){stop("second argument must be a character string")}

  rgrass::use_sp()

  Sys.setenv("GRASS_VERBOSE"=0)

  # start grass
#  rgrass::initGRASS(gisbase,home=tempdir(),mapset="PERMANENT",override = TRUE)

  rgrass::initGRASS(gisbase,home=tempdir(),mapset="PERMANENT",
                     addon_base=paste(Sys.getenv("HOME"), "/.grass8/addons", sep=""),override = TRUE)


  # # start grass
  # rgrass::initGRASS(gisbase,mapset="PERMANENT",override = TRUE)

  # set projection
  crs0 = as.character(terra::crs(rast,proj=TRUE))
  rgrass::execGRASS("g.proj",flags=c("c"),
                     parameters=list(proj4=crs0))
  # set location properties
  ext0 = terra::ext(rast)
  rgrass::execGRASS("g.region",parameters=list(rows=terra::nrow(rast),
                                                cols=terra::ncol(rast),
                                                ewres=as.character(terra::res(rast)[1]),
                                                nsres=as.character(terra::res(rast)[2]),
                                                w=as.character(ext0[1]),
                                                e=as.character(ext0[2]),
                                                s=as.character(ext0[3]),
                                                n=as.character(ext0[4])))

  # write raster into location
  write_raster_to_grass(rast,name)

}


write_raster_to_grass <- function(rast,name,warnings=F){
  opt0 = getOption("warn")
  if(warnings){options(warn=0)}else{options(warn=-1)}
  if(!(is.character(name))){stop("second argument must be a character string")}
  if(class(rast)[1]!="SpatRaster"){stop("first argument must be a SpatRaster")}
  rgrass::writeRAST(as(raster::raster(rast), 'SpatialGridDataFrame'),vname=c(name))
  options(warn=opt0)
}

read_raster_from_grass <- function(name,warnings=F){
  opt0 = getOption("warn")
  if(warnings){options(warn=0)}else{options(warn=-1)}
  if(!(is.character(name))){stop("argument must be character string")}
  rast = terra::rast(raster::raster(rgrass::readRAST(c(name))))
  return(rast)
  options(warn=opt0)
}

write_vector_to_grass <- function(vect,name,warnings=F){
  opt0 = getOption("warn")
  if(warnings){options(warn=0)}else{options(warn=-1)}
  if(!(is.character(name))){stop("second argument must be a character string")}
  if(class(vect)[1]!="SpatVector"){stop("first argument must be a SpatVector")}
  rgrass::writeVECT(as(vect,"Spatial"),vname=c(name),v.in.ogr_flags=c("overwrite","o"))
  options(warn=opt0)
}

read_vector_from_grass <- function(name,layer=1,warnings=F){
  opt0 = getOption("warn")
  if(warnings){options(warn=0)}else{options(warn=-1)}
  if(!(is.character(name))){stop("argument must be character string")}
  vector = terra::vect(rgrass::readVECT(c(name),layer=layer))
  return(vector)
  options(warn=opt0)
}









#' Trim a raster to remove NA values on the periphery
#'
#' The trim function from the terra package only remove cols and rows full of NA.
#' This function remove all cols and rows with NA on the periphery.
#' Useful to get a nicely formatted raster grid after reprojection.
#'
#' @param rast input raster (SpatRaster)
#' @param verbose flag (default FALSE) to switch
#'
#' @return a trimmed raster (SpatRaster)
#' @export
#'
#' @examples
trim_na <- function(rast,verbose=F){
  if(class(rast)[1]!="SpatRaster"){stop("argument must be a SpatRaster")}
  tmp = terra::trim(rast)
  val = sum(is.na(tmp[]))
  if(val==0){print("Raster already clean, nothing removed")}
  while(val!=0){
    n = max(1,round(0.5*val/(nrow(tmp)+ncol(tmp))))
    tmp = crop(tmp,ext(ext(tmp)[1]+n*res(tmp)[1],
                       ext(tmp)[2]-n*res(tmp)[1],
                       ext(tmp)[3]+n*res(tmp)[2],
                       ext(tmp)[4]-n*res(tmp)[2]))
    val = sum(is.na(tmp[]))
    if(verbose){
    print(paste("width in pixels of the group of rows and lines removed :",n))
    print(paste(val,"remaining NA pixels"))
    }
  }
  return(tmp)
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
  if(class(rast)[1]!="SpatRaster"){stop("1st argument must be a SpatRaster")}
  if(class(dir)[1]!="SpatRaster"){stop("2nd argument must be a SpatRaster")}
  rast0 = rast
  dir = terra::as.matrix(dir,wide=TRUE)
  rast = terra::as.matrix(rast,wide=TRUE)
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
  return(terra::rast(res,extent=terra::ext(rast0),crs=terra::crs(rast0,proj=TRUE)))
}

#' Raster projection utility
#'
#' Project a raster according to a specified CRS and resolution.
#' It will adjust the output extent according to the resolution.
#'
#' @param rast  input raster (`SpatRaster`)
#' @param epsg the EPSG code to use for projecting the raster
#' @param resolution the resolution of the projected raster (pixel size in m)
#' @param method the projection method (default `bilinear`, see `project` from `terra` package for other available options)
#' @param trim flag (default `FALSE`) to trim the periphery of the rast after projection and remove all NA pixels (use with caution as it can shrink substantially the output raster)
#'
#' @return a projected raster
#'
#' @export
#'
#' @examples
project_raster <- function(rast,epsg,resolution,method="bilinear",trim=FALSE){
  if(class(rast)[1]!="SpatRaster"){stop("1st argument must be a SpatRaster")}
  tmp <- terra::project(rast,paste("epsg:",epsg,sep=""))
  ext1 <- terra::ext(tmp)
  x1 <- as.numeric(ext1[1] -ext1[1]%%resolution)
  x2 <- as.numeric(ext1[2] -ext1[2]%%resolution + resolution)
  y1 <- as.numeric(ext1[3] -ext1[3]%%resolution)
  y2 <- as.numeric(ext1[4] -ext1[4]%%resolution + resolution)
  ext2 <- terra::ext(c(x1,x2,y1,y2))
  temp2 = rast(ext2,crs=terra::crs(tmp),resolution=c(resolution,resolution))
  rst1 = terra::resample(tmp,temp2,method="bilinear")
  rst2 = terra::trim(rst1)
  if(trim){rst2 = gtbox::trim_na(rst2)}
  return(rst2)
}


#' Easily obtain UTM EPSG code
#'
#' Determine the UTM (WGS84) EPSG code of the center of the extent of a `SpatRaster` or `SpatVector`
#'
#' @param obj  input raster or vector (`SpatRaster` or `SpatVector`)
#'
#' @return EPSG code
#'
#' @export
#'
#' @examples
get_utm<-function(obj){
  if(!class(obj)[1]%in%c("SpatRaster","SpatVector")){stop("Argument must be a SpatRaster or SpatVector")}
  if(is.lonlat(obj)){
    ext0 = terra::ext(obj)
  }else{
    ext0 <- terra::ext(terra::project(dem0,"epsg:4326"))
  }
  epsg=as.numeric(32700-round((45+(ext0[3]+ext0[4])/2)/90,0)*100+round((183+(ext0[1]+ext0[2])/2)/6,0))
  epsg_e = as.numeric(32700-round((45+(ext0[3]+ext0[4])/2)/90,0)*100+round((183+ext0[1])/6,0))
  epsg_w = as.numeric(32700-round((45+(ext0[3]+ext0[4])/2)/90,0)*100+round((183+ext0[2])/6,0))
  if ((epsg!=epsg_e) | (epsg!=epsg_w)){warning(print("Object probably extends over several UTM zones"))}
  return(epsg)
}


#' Compute continuous maps for variable defined along the stream network
#'
#' Transfer statistic values from stream segments to the corresponding contributing basins.
#' It allows to provide a continuous representation of variables defined along the stream network such as chi, for visualization purpose.
#'
#' @param st A SpatRaster object from `dem_process`
#' @param rst A SpatRaster with the variable of interest (defined along streams), which will be transferred to the corresponding contributing areas.
#' @param FUN The function used to summarize the values for each stream segment (character string : mean, min, max, median)
#'
#' @return
#' @export
#'
#' @examples
streams2basins<-function(st,rst,FUN="mean"){
  if(!is.character((FUN))){stop("FUN must be character string")}
  if(!(FUN%in%c("mean","max","min","median"))){stop("FUN must be one of : mean, max, min, median")}
  FUN <- match.fun(FUN)
  names(rst)<-"rst"
  tmp = terra::as.data.frame(c(st$st_id,rst))
  b = as.matrix(aggregate(tmp[,"rst"],by=list(tmp$st_id),FUN=FUN))
  res = terra::classify(st$bs_id,b,othersNA=TRUE)
  return(res)
}



#' #' Dowload a STRM raster from srtm.csi.cgiar.org
#' #'
#' #' Compute for each pixel of a raster the value of the following pixel according to flow direction
#' #'
#' #' @param sp  an object for which a spatial extent can be extracted
#' #' @param buf a buffering distance (in degrees) to expand the extent
#' #'
#' #' @return
#' #' @export
#' #' @examples
#' get_srtm <- function(sp,buf=0.1) {
#'   ext = extend(extent(sp),buf)
#'   lon1 = ext[1]
#'   lon2 = ext[2]
#'   lat1 = ext[3]
#'   lat2 = ext[4]
#'   stopifnot(lon1 >= -180 & lon1 <= 180 & lon2 >= -180 & lon2 <= 180)
#'   stopifnot(lat1 >= -60 & lat1 <= 60 & lat2 >= -60 & lat2 <= 60)
#'   rs <- raster(nrows=24, ncols=72, xmn=-180, xmx=180, ymn=-60, ymx=60 )
#'   rowTile1 <- rowFromY(rs, lat2)
#'   colTile1 <- colFromX(rs, lon1)
#'   rowTile2 <- rowFromY(rs, lat2)
#'   colTile2 <- colFromX(rs, lon2)
#'   rowTile3 <- rowFromY(rs, lat1)
#'   colTile3 <- colFromX(rs, lon1)
#'   rowTile4 <- rowFromY(rs, lat1)
#'   colTile4 <- colFromX(rs, lon2)
#'   if (rowTile1 < 10) { rowTile1 <- paste('0', rowTile1, sep='') }
#'   if (colTile1 < 10) { colTile1 <- paste('0', colTile1, sep='') }
#'   if (rowTile2 < 10) { rowTile2 <- paste('0', rowTile2, sep='') }
#'   if (colTile2 < 10) { colTile2 <- paste('0', colTile2, sep='') }
#'   if (rowTile3 < 10) { rowTile3 <- paste('0', rowTile3, sep='') }
#'   if (colTile3 < 10) { colTile3 <- paste('0', colTile3, sep='') }
#'   if (rowTile4 < 10) { rowTile4 <- paste('0', rowTile4, sep='') }
#'   if (colTile4 < 10) { colTile4 <- paste('0', colTile4, sep='') }
#'   f1 <- paste('srtm_', colTile1, '_', rowTile1, sep="")
#'   f2 <- paste('srtm_', colTile2, '_', rowTile2, sep="")
#'   f3 <- paste('srtm_', colTile3, '_', rowTile3, sep="")
#'   f4 <- paste('srtm_', colTile4, '_', rowTile4, sep="")
#'   lf = unique(c(f1,f2,f3,f4))
#'   lr = list()
#'   for (i in 1:length(lf)){
#'     f = lf[[i]]
#'     theurl <- paste("https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/", f, ".zip", sep="")
#'     utils::download.file(url=theurl, destfile=paste("/tmp/",f,".zip",sep=""), method="auto", quiet = FALSE, mode = "wb", cacheOK = TRUE)
#'     system(paste("unzip -o -d /tmp /tmp/",f,".zip",sep=""))
#'     lr[[i]] = raster(paste("/tmp/",f,".tif",sep=""))
#'   }
#'   if (length(lr)==1){
#'     return(raster::crop(lr[[1]],ext))
#'   } else {
#'     rs = lr[[1]]
#'     for (i in 2:length(lr)){rs=raster::merge(lr[[i]],rs)}
#'     return(raster::crop(rs,ext))
#'   }
#' }


