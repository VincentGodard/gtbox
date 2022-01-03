library("gtbox")
library("terra")
library("tictoc")
gisbase = "/usr/lib/grass78/"

# get DEM
data("dem_cevennes",package = "gtbox")
dem = rast(dem_cevennes)
rm(dem_cevennes)

st = process_dem(dem,th_px=2000,gisbase=gisbase)
net = get_network(st,gisbase)
out = get_outlets(st,large = 10e4)
bas = get_basins(st,out,gisbase)


# compute chi
tic()
chi = compute_chi(st,starting_points=out,a0=1000,mn=0.5,gisbase)
toc()



# continuous chi map
tmp1 = terra::mask(st$bs_id,bas)
tmp2 = terra::mask(chi,bas)
tmp = as.data.frame(c(tmp1,tmp2))
b = as.matrix(aggregate(tmp[,"chi"],by=list(tmp$bs_id),FUN=mean))
chi_map = terra::classify(tmp1,b)

# network topology
nxt_id = get_next(st$st_id,st$dir)
seg = as.data.frame(c(st$st_id,nxt_id,st$acc))
colnames(seg)<-c("id","down","acc")
seg = seg[seg$id!=seg$down,]
# those that drain outside
ids = unlist(unique(st$st_id))
tmp = data.frame(id=ids[!ids%in%seg$id],down=-1,acc=NA)
seg = rbind(seg,tmp)

seg$up1 = NA
seg$up2 = NA
seg$up3 = NA
for(i in 1:nrow(seg)){
  up = seg[seg$down==seg$id[i],]$id
  acc = seg[seg$down==seg$id[i],]$acc
  up= up[order(acc)]
  if(length(up)==1){
    seg$up1[i]=up[1]
    }
  if(length(up)==2){
    seg$up1[i]=up[1]
    seg$up2[i]=up[2]
  }
    if(length(up)==3){
      seg$up1[i]=up[1]
      seg$up2[i]=up[2]
      seg$up3[i]=up[3]

  }
  }


