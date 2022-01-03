#' Compute chi
#'
#'
#' @param st  A SpatRaster object from `dem_process`
#' @param starting_points starting points (SpatVector) for the integration (outlets of basins)
#' @param a0 reference accumulation in pixels
#' @param mn value of the m/n ratio
#' @param gisbase The directory path to GRASS binaries and libraries, containing bin and lib subdirectories among others
#'
#' @return a SpatRaster with chi values along streams
#' @export
#'
#' @examples
compute_chi <- function(st,starting_points,a0,mn,gisbase){
  # get basins
  bas = get_basins(st,starting_points,gisbase=gisbase)
  ext0 = terra::ext(st)
  st = terra::mask(st,bas)
  st = trim(st)

  # compute chi (raster based)
  prm = (a0/st$acc)^mn*terra::res(st)[1]
  prm[is.na(st$st_id)] <- NA
  start_grass(prm,"prm",gisbase)
  write_vector_to_grass(starting_points,"st_pts")
  rgrass7::execGRASS("r.cost", flags=c("overwrite","n"), parameters=list(input="prm",
                                                                         output="chi",
                                                                         start_points="st_pts"))
  chi = read_raster_from_grass("chi")
  chi = terra::extend(chi,ext0)
  return(chi)
}
