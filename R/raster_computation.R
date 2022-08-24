
#' Compute zonal statistics from a raster layer
#'
#' This function is a frontend for GRASS module `v.rast.stats`
#'
#' @param vect a `SpatVector` defining the locations/areas where the calculation will take place
#' @param rast a `SpatRaster` representing the raster
#' @param name a prefix for the statistics columns
#' @param gisbase The directory path to GRASS binaries and libraries, containing bin and lib sub-directories among others
#'
#'
#' @return a `SpatVector` with attribute table filled with statistics from the raster layer
#' @export
#'
#' @examples
compute_zonal_stats <- function(vect,rast,name,gisbase){
  start_grass(rast,"rast",gisbase)
  write_vector_to_grass(vect,"vect")
  rgrass::execGRASS("v.rast.stats", flags=c("c"),
                     parameters=list(map="vect",
                                     raster="rast",
                                     column_prefix=name))
  res = read_vector_from_grass("vect")
  return(res)
}



#' Compute curvature over a raster
#'
#' Fits quadratic surfaces at each pixel
#'
#' The surface is parametrized as :
#' z = ax^2 + by^2 +cxy +dx + ey + f
#'
#' Curvature is calculated as C =  2a + 2b
#'
#' the computation relies on GRASS GIS `r.param.scale` function
#'
#' @param rast a `SpatRaster` representing the raster (usually a Digital Elevation Model)
#' @param win an odd integer (>=3) specifying the size of the computation window (in pixels)
#' @param gisbase The directory path to GRASS binaries and libraries, containing bin and lib sub-directories among others
#' @param central Constrain fitted surface through central window cell (default TRUE)
#' @param planar Compute planar curvature instead (default FALSE)
#'
#' @return a curvature `SpatRaster` (in 1/m)
#' @export
#'
#' @examples
compute_curvature <- function(rast,win,gisbase,central=TRUE,planar=FALSE){

  if ((win<3) & (win%%2==0)){stop("win must be an odd integer >=3")}

  if (central) {
    fl = c("overwrite","c")
  }else{
    fl = c("overwrite")
  }

  # start grass session
  start_grass(rast,"rast",gisbase)
  Sys.setenv(OMP_NUM_THREADS=1)

  if (!planar){

    rgrass::execGRASS("r.param.scale", flags=fl,
                       parameters=list(input="rast",
                       output="r_maxic",
                       method="maxic",
                       size=win))

    rgrass::execGRASS("r.param.scale", flags=fl,
                       parameters=list(input="rast",
                       output="r_minic",
                       method="minic",
                       size=win))

    rgrass::execGRASS("r.mapcalc", expression="curv = r_maxic + r_minic")

  }else{
    rgrass::execGRASS("r.param.scale", flags=fl,
                       parameters=list(input="rast",
                       output="curv",
                       method="planc",
                       size=win))
  }

  curv = read_raster_from_grass("curv")
  return(curv)
}



#' Compute slope over a raster
#'
#' Fits quadratic surfaces at each pixel
#'
#' The surfaces are parametrized as :
#' z = ax^2 + by^2 +cxy +dx + ey + f
#'
#' Slope is calculated as S = (d^2 + e^2)^0.5
#'
#' the computation relies on GRASS GIS `r.param.scale` function
#'
#' @param rast a `SpatRaster` representing the raster (usually a Digital Elevation Model)
#' @param win an odd integer (>=3) specifying the size of the computation window (in pixels)
#' @param central Constrain fitted surface through central window cell (default TRUE)
#' @param gisbase The directory path to GRASS binaries and libraries, containing bin and lib sub-directories among others
#'
#'
#' @return a gradient `SpatRaster` (in m/m)
#' @export
#'
#' @examples
compute_slope <- function(rast,win,gisbase,central=TRUE){

  if ((win<3) & (win%%2==0)){stop("win must be an odd integer >=3")}

  if (central) {
    fl = c("overwrite","c")
  }else{
    fl = c("overwrite")
  }

  # start grass session
  start_grass(rast,"rast",gisbase)
  Sys.setenv(OMP_NUM_THREADS=1)

  rgrass::execGRASS("r.param.scale", flags=c("overwrite"),
                     parameters=list(input="rast",
                     output="slope",
                     method="slope",
                     size=win))

  slope = read_raster_from_grass("slope")
  slope = tan(slope/180*pi)
  return(slope)
}
