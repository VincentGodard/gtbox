% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_utm}
\alias{get_utm}
\title{Easily obtain UTM EPSG code}
\usage{
get_utm(obj)
}
\arguments{
\item{obj}{input raster or vector (`SpatRaster` or `SpatVector`)}
}
\value{
EPSG code
}
\description{
Determine the UTM (WGS84) EPSG code of the center of the extent of a `SpatRaster` or `SpatVector`
}
