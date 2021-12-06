


#' Perform initial computations on a Digital Elevation Model
#'
#' Perform initial computations on a Digital Elevation Model
#' This function relies on the following GRASS GIS modules :
#' r.watershed, r.stream.distance (add-on to be installed)
#'
#' @param dem  Digital Elevation Model (SpatRaster)
#' @param th_px Accumulation threshold in pixels to initiate a stream (integer)
#' @param gisBase The directory path to GRASS binaries and libraries, containing bin and lib subdirectories among others
#' @param to_net Flag to compute distance to network (`dtn`) and elevation above network (`ean`) (optional)
#' @param precip precipitation raster (m/a) to be used for annual discharge computation
#'
#' @return A SpatRaster with the following layers: `z` (DEM), `acc` (flow accumulation), `dir` (flow direction), `dist` (distance along network), `st_id` (stream id), `nxt_id` (id of the next downstream pixel), `bs_id` (sub basin id)
#' Optionally
#'
#' @export
#'
#' @examples
process_dem<-function(dem,th_px,gisBase,precip=NULL,to_net=FALSE){
  if(class(dem)[1]!="SpatRaster"){stop("dem must be a SpatRaster")}

  # start grass session and import dem
  start_grass(dem,"dem",gisBase)

  dem2 = read_raster_from_grass("dem")
  names(dem2) <- "z"

  # compute accumulation, flow direction and streams/basins ids
  rgrass7::execGRASS("r.watershed", flags=c("overwrite","s"),
                     parameters=list(elevation="dem",
                     threshold=th_px,
                     accumulation="accumulation",
                     drainage="direction",
                     stream="streams",
                     basin="basins"))
  acc = read_raster_from_grass("accumulation")
  names(acc) <- "acc"
  dir = read_raster_from_grass("direction")
  names(dir) <- "dir"
  st_id = read_raster_from_grass("streams")
  names(st_id) <- "st_id"
  bs_id = read_raster_from_grass("basins")
  names(bs_id) <- "bs_id"

  # compute strahler orders
  rgrass7::execGRASS("r.stream.order", flags=c("overwrite"),
                     parameters=list(stream_rast="streams",
                                     direction="direction",
                                     strahler="strahler"))
  sto = read_raster_from_grass("strahler")
  names(sto) <- "sto"

  # get next pixel
  #nxt_id = get_next(st_id,dir)
  #nxt_id@data@names <- "nxt_id"

  # compute distance along network
  rgrass7::execGRASS("r.stream.distance", flags=c("overwrite","o"),
                     parameters=list(stream_rast="streams",
                     direction="direction",
                     distance="distance"))
  dist = read_raster_from_grass("distance")
  names(dist) <- "dist"

  rast_list = c(dem2,acc,dir,dist,st_id,bs_id,sto)

  # distance  to network and  elevation above network
  if (to_net){
    rgrass7::execGRASS("r.stream.distance", flags=c("overwrite"),
                       parameters=list(stream_rast="streams",
                                       direction="direction",
                                       elevation="dem",
                                       distance="dtn",
                                       difference="ean"))
    dtn = read_raster_from_grass("dtn")
    names(dtn) <- "dtn"
    ean = read_raster_from_grass("ean")
    names(ean) <- "ean"
    rast_list =  c(rast_list,dtn,ean)
  }

  # compute discharge if precipitation raster is provided
  if (!(is.null(precip))){
    write_raster_to_grass(precip,"precip")

    rgrass7::execGRASS("r.watershed", flags=c("overwrite","s"),
                       parameters=list(elevation="dem",
                                       threshold=th_px,
                                       accumulation="accumulation",
                                       drainage="tmp1",
                                       stream="tmp2",
                                       flow="precip"))
    discharge = read_raster_from_grass("accumulation")
    discharge = discharge * terra::res(dem)[1]*terra::res(dem)[2]
    names(discharge) <- "discharge"
    rast_list =  c(rast_list,discharge)
  }

  return(rast_list)
}








# TODO NA to NULL?
#' Extract basins outlets locations with various approaches. These various options are mutually exclusives.
#'
#'
#' option `strahler` extract all the basins of the given strahler order
#'
#' option `large` extract the largest possible basins from the DEM.
#' The indicated minimum size should be larger than the accumulation threshold used to delineate the network.
#'
#' option `elevation` extract all outlets at a given elevation
#'
#' @param st  A SpatRaster multi-layers object from `dem_process`
#' @param strahler The Strahler order of the outlets to be extracted
#' @param large The minimum size (in pixels) for the extraction of the largest possible set of basins
#' @param elevation The elevation of the outlets
#'
#' @return a SpatVector with identified outlets
#' @export
#'
#' @examples
get_outlets <- function(st,strahler=NA,large=NA,elevation=NA,gisBase=NA){
  if(class(st)[1]!="SpatRaster"){stop("1st argument must be a SpatRaster")}
  if (is.na(strahler) & is.na(large) & is.na(elevation)){stop("Must specify one outlet selection option : strahler, large or elevation")}
  if (sum(c(!is.na(strahler),!is.na(large),!is.na(elevation)))>1){stop("Must specify only one outlet selection option")}
  #
  if (!is.na(strahler)) {
  #
  next_sto = get_next(st$sto,st$dir)
  next_id = get_next(st$st_id,st$dir)
  pos = st$st_id # initiate the raster
  pos[st$st_id!=next_id] <- 1 # select the end of stream segment
  pos[st$st_id==next_id] <- NA # remove the rest of the segment
  pos[st$sto==next_sto] <- NA # remove the segments (last pixel) that do not see a change in strahler order
  pos[st$acc<0] <- NA # remove negative accumulation
  pos = pos*st$sto
  names(pos) <- "sto"
  pts = terra::as.points(pos)
  pts = pts[pts$sto==strahler,]
  return(pts)
  }
  #
  if (!is.na(large)) {
    th = min(abs(st$acc[])*(st$st_id[]/st$st_id[]),na.rm=T)
   if(large<th){
    stop(paste("For option large : number of pixels must be larger than accumulation threshold (",th," pixels)",sep=""))}
    tmp1 = st$st_id
    tmp1[1,] <- 0
    tmp1[nrow(st),] <- 0
    tmp1[,1] <- 0
    tmp1[,ncol(st)] <- 0
    tmp1[st$acc<=0] <- 0
    tmp2 = get_next(tmp1,st$dir)
    tmp2[tmp2!=0] <- NA
    tmp2[tmp2==0] <- 1
    tmp2[st$acc<large] <- NA
    pts = na.omit(terra::as.data.frame(tmp2*st$acc,xy=TRUE))
    colnames(pts) <- c("x","y","acc")
    pts = pts[order(pts$acc,decreasing = TRUE),]
    row.names(pts) <- 1:nrow(pts)
    tab = as.data.frame(pts$acc)
    colnames(tab)<-"acc"
    pts = terra::vect(as.matrix(pts[,c("x","y")]),type="points",
                      atts=tab,crs=terra::crs(st,proj=TRUE))
    return(pts)
  }
  #
  if (!is.na(elevation)) {
    nxt_z =  get_next(st$z,st$dir)
    st_pts = st$z>=elevation & nxt_z<elevation & st$st_id>0
    st_pts[is.na(st_pts)] <- 0
    #
    start_grass(st$z,"dem",gisbase)
    write_raster_to_grass(st_pts,"flow")
    rgrass7::execGRASS("r.watershed", flags=c("overwrite","s"),
                       parameters=list(elevation="dem",
                                       accumulation="accumulation",
                                       drainage="direction",
                                       flow="flow"))
    acc = read_raster_from_grass("accumulation")
    #
    nxt_acc =  get_next(acc,st$dir)
    nxt_start =  get_next(st_pts,st$dir)
    st_pts = abs(nxt_acc)==1 & abs(nxt_start)==1 & st$st_id>0 & st$acc>0
    st_pts[st_pts==0]<-NA
    pts = terra::as.points(st_pts)
    #
    return(pts)
  }

}



# TODO : add a switch for overlaying basins
#' Compute basins corresponding to outlets
#'
#'
#' Follow the largest accumulation tributary at confluences when moving upstream
#' This function relies on the following GRASS GIS modules :
#' r.stream.basins (add-on to be installed)
#'
#' @param st  A RasterStack object from dem_process
#' @param outlets  Outlets points (SpatialPointsDataFrame)
#' @param gisBase The directory path to GRASS binaries and libraries, containing bin and lib sub-directories among others
#'
#' @return A SpatialPolygonDataFrame with the corresponding basins. The attribute table of the outlets vector is transferred.
#' @export
#'
#' @examples
get_basins<-function(st,outlets,gisbase){
  if(class(st)[1]!="SpatRaster"){stop("1st argument must be a SpatRaster")}
  if(class(outlets)[1]!="SpatVector"){stop("2nd argument must be a SpatVector")}

  # start grass session
  start_grass(st$dir,"direction",gisbase)

  # write outlets into location
  write_vector_to_grass(outlets,"outlets")

  # compute basins
  rgrass7::execGRASS("r.mapcalc", expression="direction0 = round(direction)")
  rgrass7::execGRASS("r.stream.basins", flags=c("overwrite"),
                     parameters=list(direction="direction0",
                                     points="outlets",
                                     basins="basins"))

  # vectorize basins
  rgrass7::execGRASS("r.to.vect", flags=c("overwrite","v"),
                     parameters=list(input="basins",
                                     output="basins",
                                     type="area"))

  # attach attribute table
  rgrass7::execGRASS("v.db.join", parameters=list(map="basins",
                                                  column="cat",
                                                  other_table="outlets",
                                                  other_column="cat"))

  # # clean basins
  # px_area = ((direction@extent[2]-direction@extent[1])/direction@ncols * (direction@extent[4]-direction@extent[3])/direction@nrows ) # pixel area
  # pars <- list(input="temp",output="basins",type="area",tool="rmarea",threshold=px_area*10)
  # execGRASS("v.clean", flags=c("overwrite"), parameters=pars)

  # export vector from Grass
  basins = read_vector_from_grass("basins")

  return(basins)
}

#' Compute river network vector based on a SpatRast multi-layers raster
#'
#'
#' This function relies on the following GRASS GIS modules :
#' r.stream.order (add-on to be installed)
#'
#' @param st  A SpatRast object from `process_dem`
#' @param gisBase The directory path to GRASS binaries and libraries, containing bin and lib sub-directories among others
#' @param clip Remove parts with incomplete contributing areas (logical, default FALSE)
#'
#' @return A SpatVector with the corresponding network. The attribute table contains information about network topology.
#' @export
#'
#' @examples
get_network<-function(st,gisbase,clip=FALSE){
  if(class(st)[1]!="SpatRaster"){stop("1st argument must be a SpatRaster")}
  #
  if(clip){st$st_id[st$acc<0]<-NA}
  #
  start_grass(st$z,"dem",gisbase)
  write_raster_to_grass(st$st_id,"streams")
  rgrass7::execGRASS("r.mapcalc", expression="streams0 = round(streams)")
  write_raster_to_grass(st$acc,"acc")
  write_raster_to_grass(st$dir,"dir")
  rgrass7::execGRASS("r.mapcalc", expression="dir0 = round(dir)")
  rgrass7::execGRASS("r.stream.order", flags=c("overwrite"),
                     parameters=list(stream_rast="streams0",
                                     accumulation="acc",
                                     elevation="dem",
                                     direction="dir0",
                                     stream_vect="temp"))
  network = read_vector_from_grass("temp")
  return(network)
  }



#' Extract portions of the stream network for river long-profile analysis
#'
#' Get ids for all streams segments moving up (default) or down from a starting stream segment (identified with `id`).
#' Follow the largest accumulation tributary at confluences when moving upstream.
#'
#' @param net  stream network from `get_network` (SpatVector)
#' @param id id of the starting stream segment (integer)
#' @param mode type of extraction : upstream from id (`up`, default), upstream including all tributaries (`up_all`) or downstream from id (`down`).
#'
#' @return
#' If `mode` is `up` (default) returns a vector with the upstream ids.
#' If `mode` is `up_all` returns a dataframe with the first line corresponding to the main stream, and following lines to tributaries.
#' The column `ids` contains the ids of the stream segments, `length` the total length of the stream and `n` the number of segments.
#' If `mode` is `down` returns a vector with the downstream ids.
#'
#' @export
#'
#' @examples
get_streams<-function(net,id,mode="up"){
  if(class(net)[1]!="SpatVector"){stop("1st argument must be a SpatVector")}
  if (!(mode %in% c("up","up_all","down"))){stop("mode must be one of up, up_all or down")}
  net = terra::as.data.frame(net)
  if (mode=="up"){
  list_streams = id
  while(net$prev_str01[net$stream==id]!=0 |
        net$prev_str02[net$stream==id]!=0 |
        net$prev_str03[net$stream==id]!=0){
        prev_str = c(net$prev_str01[net$stream==id],
                     net$prev_str02[net$stream==id],
                     net$prev_str03[net$stream==id])
       prev_acc = c(0,0,0)
       if(prev_str[1]!=0){prev_acc[1] = net$flow_accum[net$stream==prev_str[1]]}
       if(prev_str[2]!=0){prev_acc[2] = net$flow_accum[net$stream==prev_str[2]]}
       if(prev_str[3]!=0){prev_acc[3] = net$flow_accum[net$stream==prev_str[3]]}
       imax = which.max(prev_acc)
       id = prev_str[imax]
    list_streams = c(list_streams,id)
    }
  return(list_streams)
  }
  #
  if (mode=="up_all"){
    tribs = c()
    streams = data.frame(stream=NA,ids=NA)

    k = 1
    while(k<=length(tribs) | k==1 ){

    ll = get_up(net,id)
    streams[k,"ids"][[1]]<-list(ll[[1]])
    streams[k,"stream"] = k
    tribs = c(tribs,ll[[2]])
    id = tribs[k]
    k = k+1
    }

   streams$n = NA
   streams$length = NA
   for  (i in 1:nrow(streams)){
     ids = unlist(streams[i,"ids"])
     streams$n[i] = length(ids)
     streams$length[i] = sum(net$length[net$stream%in%ids])
   }
  return(streams)
  }
 #
  if (mode=="down"){
    list_streams = id
    while(net$next_stream[net$stream==id]!=-1){
      id = net$next_stream[net$stream==id]
      list_streams = c(list_streams,id)
    }
    return(list_streams)
  }

}

get_up<-function(net,id){
  list_streams = id
  list_tribs = c()
  while(net$prev_str01[net$stream==id]!=0 |
        net$prev_str02[net$stream==id]!=0 |
        net$prev_str03[net$stream==id]!=0){
    prev_str = c(net$prev_str01[net$stream==id],
                 net$prev_str02[net$stream==id],
                 net$prev_str03[net$stream==id])
    prev_acc = c(0,0,0)
    if(prev_str[1]!=0){prev_acc[1] = net$flow_accum[net$stream==prev_str[1]]}
    if(prev_str[2]!=0){prev_acc[2] = net$flow_accum[net$stream==prev_str[2]]}
    if(prev_str[3]!=0){prev_acc[3] = net$flow_accum[net$stream==prev_str[3]]}
    imax = which.max(prev_acc)
    id = prev_str[imax]
    tribs = prev_str[-imax]
    tribs = tribs[tribs>0]
    list_streams = c(list_streams,id)
    list_tribs =  c(list_tribs,tribs)
  }
  return(list(list_streams,list_tribs))
}
