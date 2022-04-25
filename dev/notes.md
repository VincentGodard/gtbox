- compute_chi : raster based -> option to pass basins to speed up calculation


# rgrass7 problem
I am using rgrass7 to define throwaway GRASS locations. 
GRASS version : 7.8.6
R version : 3.6.3
Ubuntu 20.04.4 

With rgrass7 0.2.6 I am routinely using :
gisbase = "/usr/lib/grass78/"
rgrass7::initGRASS(gisbase,home=tempdir(),mapset="PERMANENT",override = TRUE)
and then calls to g.proj and g.region to set the location properties (through execGRASS). 

I know I could use the SG option to set everything at once, but  I try to move my R workflow to terra  and would rather avoid dependencies on sp and raster, which would be necessary to convert a terra SpatRaster to the required SpatialGrid if I am not mistaken. 

Now with rgrass7 0.2.8 :

rgrass7::initGRASS(gisbase,home=tempdir(),mapset="PERMANENT",override = TRUE)

returns the following error :

Error in execGRASS("g.proj", flags = c("w"), intern = TRUE, ignore.stderr = ignore.stderr) : 
  The command:
g.proj -w
produced an error (1) during execution:
ATTENTION: fichier <PROJ_INFO> introuvable pour le secteur
           <file3f47d24ad6933>
ATTENTION: fichier <PROJ_UNITS> introuvable pour le secteur
           <file3f47d24ad6933>
ERREUR : Fichiers de projection manquant

Due to the non-existing projection information at this stage in the location (which apparently was not an issue with 0.2.6).

I tried to use the SG option but I also get an error in this case (rgrass7 0.2.8), where it seems to expect a latlon system whereas dem is UTM31N.

dem0 <- as(raster(dem), "SpatialGrid")
rgrass7::initGRASS(gisbase,SG=dem0,home=tempdir(),mapset="PERMANENT",override = TRUE)
Error in execGRASS("g.region", save = "input", flags = "overwrite") : 
  The command:
g.region --overwrite save=input
produced an error (1) during execution:
ERREUR : Lattitude Nord illégale : 4.91676e+06


> dem
class       : SpatRaster 
dimensions  : 1596, 1654, 1  (nrow, ncol, nlyr)
resolution  : 30, 30  (x, y)
extent      : 539757.5, 589377.5, 4868879, 4916759  (xmin, xmax, ymin, ymax)
coord. ref. : WGS 84 / UTM zone 31N (EPSG:32631) 
source      : memory 
name        : srtm_cevennes_utm 
min value   :                85 
max value   :              1672 

