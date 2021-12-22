library("terra")
library("gdalUtils")
library("gtbox")
dem_vrt="/home/vincent/data/SRTM_GL1_srtm.vrt"

v_lat = c(45,17,-22,-40,  38,  58,48)
v_lon = c( 5, 8,130,-71,-116,-125,97)
ext=0.5
resolution = 100


for (i in 1:length(v_lat)){
#for (i in 1:1){
  print(i)
  lat=v_lat[i]
  lon=v_lon[i]
  gdal_translate(src_dataset = dem_vrt,dst_dataset = "/tmp/mnt_geo.tif",
                 projwin = c(lon-ext,lat+ext,lon+ext,lat-ext))
  dem0 = terra::rast("/tmp/mnt_geo.tif")
  epsg = get_utm(dem0)
  dem = project_raster(dem0,epsg,resolution,trim=F) # get a projected, clean and nicely aligned DEM
  dem1 = trim_na(dem,verbose = T)
  plot(dem,main=epsg)
  lines(ext(dem1))
}


# dem0=rast("/home/vincent/data/chantiers/france/rasters/rge_1m/valensole/1_fedes_dem_denoised.tif")
# epsg=32631
# resolution = 2







