
#' Compute zonal statistics from a raster layer
#'
#' This function is a frontend for GRASS module `v.rast.stats`
#'
#' @param vect the vector layer used for computation
#' @param rast a `RasterLayer` representing the raster
#' @param name a prefix for the statistics columns
#' @param gisBase The directory path to GRASS binaries and libraries, containing bin and lib sub-directories among others
#'
#'
#' @return a vector with attribute table filled with statistics from the raster layer
#' @export
#'
#' @examples
compute_zonal_stats <- function(vect,rast,name,gisbase){
  start_grass(rast,"rast",gisbase)
  rgrass7::writeVECT(vect,vname=c("vect"),v.in.ogr_flags=c("overwrite","o"))
  rgrass7::execGRASS("v.rast.stats", flags=c("c"),
                     parameters=list(map="vect",
                                     raster="rast",
                                     column_prefix=name))
  res = rgrass7::readVECT(c("vect"))
  crs(res)<-crs(vect)
  return(res)
}



#' Compute curvature over a raster
#'
#' Fits quadratic surfaces at each pixel
#'
#' The surface is parametrized as :
#' z = ax^2 + by^2 +cxy +dx + ey + f
#'
#' curvature is C =  2a + 2b
#'
#' the computation relies on GRASS GIS `r.param.scale` function
#'
#' @param rast a `RasterLayer` representing the raster (usually a Digital Elevation Model)
#' @param win an odd integer (>=3) specifying the size of the computation window (in pixels)
#' @param gisBase The directory path to GRASS binaries and libraries, containing bin and lib sub-directories among others
#' @param central Constrain fitted surface through central window cell (default TRUE)
#' @param planar Compute planar curvature instead (default FALSE)
#'
#' @return a curvature raster
#' @export
#'
#' @examples
compute_curvature <- function(rast,win,gisBase,central=TRUE,planar=FALSE){

  if ((win<3) & (win%%2==0)){stop("win must be an odd integer >=3")}

  if (central) {
    fl = c("overwrite","c")
  }else{
    fl = c("overwrite")
  }

  # start grass session
  start_grass(rast,"rast",gisBase)
  Sys.setenv(OMP_NUM_THREADS=1)

  if (!planar){

    rgrass7::execGRASS("r.param.scale", flags=fl, parameters=list(input="rast",
                                                                  output="r_maxic",
                                                                  method="maxic",
                                                                  size=win))

    rgrass7::execGRASS("r.param.scale", flags=fl, parameters=list(input="rast",
                                                                  output="r_minic",
                                                                  method="minic",
                                                                  size=win))

    rgrass7::execGRASS("r.mapcalc", expression="curv = r_maxic + r_minic")

  }else{
    rgrass7::execGRASS("r.param.scale", flags=fl, parameters=list(input="rast",
                                                                  output="curv",
                                                                  method="planc",
                                                                  size=win))
  }

  curv = raster::raster(rgrass7::readRAST(c("curv")))

  crs(curv) <- crs(rast)
  extent(curv) <- extent(rast)

  return(curv)
}



#' Compute slope over a raster
#'
#' Fits quadratic surfaces at each pixel
#'
#' The surfaces are parametrized as :
#' z = ax^2 + by^2 +cxy +dx + ey + f
#'
#'
#'
#' the computation relies on GRASS GIS `r.param.scale` function
#'
#' @param rast a `RasterLayer` representing the raster (usually a Digital Elevation Model)
#' @param win an odd integer (>=3) specifying the size of the computation window (in pixels)
#' @param central Constrain fitted surface through central window cell (default TRUE)
#' @param gisBase The directory path to GRASS binaries and libraries, containing bin and lib sub-directories among others
#'
#'
#' @return a slope raster
#' @export
#'
#' @examples
compute_slope <- function(rast,win,gisBase,central=TRUE){

  if ((win<3) & (win%%2==0)){stop("win must be an odd integer >=3")}

  if (central) {
    fl = c("overwrite","c")
  }else{
    fl = c("overwrite")
  }

  # start grass session
  start_grass(rast,"rast",gisBase)
  Sys.setenv(OMP_NUM_THREADS=1)

  rgrass7::execGRASS("r.param.scale", flags=c("overwrite"), parameters=list(input="rast",
                                                                            output="slope",
                                                                            method="slope",
                                                                            size=win))

  slope = raster::raster(rgrass7::readRAST(c("slope")))

  crs(slope) <- crs(rast)
  extent(slope) <- extent(rast)

  return(slope)

}
