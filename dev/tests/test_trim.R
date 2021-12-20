library("terra")
library("gdalUtils")
library("terra")
dem_vrt="/home/vincent/data/SRTM_GL1_srtm.vrt"


lat=-10
lon=-60
ext=0.2
resolution = 100


gdal_translate(src_dataset = dem_vrt,dst_dataset = "/tmp/mnt_geo.tif",
               projwin = c(lon-ext,lat+ext,lon+ext,lat-ext))
dem0 = terra::rast("/tmp/mnt_geo.tif")

# dem0=rast("/home/vincent/data/chantiers/france/rasters/rge_1m/valensole/1_fedes_dem_denoised.tif")
# epsg=32631
# resolution = 2



dem = project_raster(dem0,epsg,resolution,trim=F) # get a projected, clean and nicely aligned DEM
dem1 = trim_na(dem,verbose = T)




