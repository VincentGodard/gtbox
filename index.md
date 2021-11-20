

## Installation
The simplest way to install the package is using `devtools` :

`devtools::install_github("VincentGodard/gtbox")`


### Dependance on Grass

Some of the functions will use the `rgrass7` package to connect to [GRASS GIS](https://grass.osgeo.org), which must be installed on your system.
These functions will need the directory path to GRASS binaries and libraries, containing bin and lib subdirectories (`gisBase` argument). 
On a Linux system this can determined by running `echo $GISBASE` inside a GRASS session.

You must also install a few GRASS addons, which can be done easily in a GRASS session through the module `g.extension` :

- `r.stream.basins`
- `r.stream.distance`
- `r.stream.order`

