
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



# TODO check win is odd and >1
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
#' @param planar flag to compute planar curvature instead (default FALSE)
#'
#' @return a curvature raster
#' @export
#'
#' @examples
compute_curvature <- function(rast,win,gisBase,planar=FALSE){

  # start grass session
  start_grass(rast,"rast",gisBase)
  Sys.setenv(OMP_NUM_THREADS=1)

  if (!planar){

    pars <- list(input="rast",output="r_maxic",method="maxic",size=win)
    rgrass7::execGRASS("r.param.scale", flags=c("overwrite"), parameters=pars)
    #rgrass7::execGRASS("r.univar", parameters=list(map="r_maxic"))

    pars <- list(input="rast",output="r_minic",method="minic",size=win)
    rgrass7::execGRASS("r.param.scale", flags=c("overwrite"), parameters=pars)
    #rgrass7::execGRASS("r.univar", parameters=list(map="r_minic"))

    rgrass7::execGRASS("r.mapcalc", expression="curv = r_maxic + r_minic")
    #rgrass7::execGRASS("r.univar", parameters=list(map="curv"))

  }else{
    pars <- list(input="rast",output="curv",method="planc",size=win)
    rgrass7::execGRASS("r.param.scale", flags=c("overwrite"), parameters=pars)
  }

  #  curv = raster(as.matrix(raster(rgrass7::readRAST(c("curv")))),template=rast)
  curv = raster(rgrass7::readRAST(c("curv")))

  crs(curv) <- crs(rast)
  extent(curv) <- extent(rast)


  return(curv)

}



# TODO check win is odd and >1
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
#' @param gisBase The directory path to GRASS binaries and libraries, containing bin and lib sub-directories among others
#'
#'
#' @return a slope raster
#' @export
#'
#' @examples
compute_slope <- function(rast,win,gisBase){

  # start grass session
  start_grass(rast,"rast",gisBase)
  Sys.setenv(OMP_NUM_THREADS=1)

  pars <- list(input="rast",output="slope",method="slope",size=win)
  rgrass7::execGRASS("r.param.scale", flags=c("overwrite"), parameters=pars)

  #slope = raster(as.matrix(rgrass7::readRAST(c("slope"))),template=rast)
  slope = raster::raster(rgrass7::readRAST(c("slope")))

  crs(slope) <- crs(rast)
  extent(slope) <- extent(rast)

  return(slope)

}
