library("terra")
library("pryr")
library("gtbox")
library("tictoc")

# test2 -------------
read_raster_from_grass2 <- function(name,warnings=F){
  opt0 = getOption("warn")
  if(warnings){options(warn=0)}else{options(warn=-1)}
  if(!(is.character(name))){stop("argument must be character string")}
  rast = terra::rast(raster::raster(rgrass7::readRAST(c(name))))*1
  #rast = as(raster::raster(rgrass7::readRAST(c(name))),"SpatRaster")
  #tmp = rgrass7::readRAST(c(name))
  #rast = terra::rast(matrix(tmp@data[[1]],nrow=tmp@grid@cells.dim[2],byrow=T),
  #                   extent=ext(tmp@bbox[1], tmp@bbox[3], tmp@bbox[2], tmp@bbox[4]),
  #                   crs=crs(tmp@proj4string))
  return(rast)
  options(warn=opt0)
}


# test 1
gisbase = "/usr/lib/grass78/"
th_px = 2000 # stream initiation threshold in pixels
dem1 = rast("/home/vincent/jobs/gabilan_mesa3/LIDAR_GROUND.tif")
start_grass(dem1,"dem",gisbase)

n=40
mem = rep(NA,n)
for(i in 2:n){
  print(i)
  tic()
  assign(paste("dem",i, sep=""),read_raster_from_grass2("dem"))
  toc()
  mem[i] = free_RAM()
}

