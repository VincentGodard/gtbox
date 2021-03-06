% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dem_processing.R
\name{get_network}
\alias{get_network}
\title{Compute river network vector based on a SpatRast multi-layers raster}
\usage{
get_network(st, gisbase, clip = FALSE)
}
\arguments{
\item{st}{A SpatRast object from `process_dem`}

\item{gisbase}{The directory path to GRASS binaries and libraries, containing bin and lib sub-directories among others}

\item{clip}{Remove parts with incomplete contributing areas (logical, default FALSE)}
}
\value{
A SpatVector with the corresponding network. The attribute table contains information about network topology.
}
\description{
This function relies on the following GRASS GIS modules :
r.stream.order (add-on to be installed)
}
