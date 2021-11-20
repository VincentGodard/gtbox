#' Digital Elevation Model of the Cevennes area SE France
#'
#' The Shuttle Radar Topography Mission (SRTM) obtained elevation data on a near-global scale to generate the most complete high-resolution digital topographic database of Earth. SRTM consisted of a specially modified radar system that flew onboard the Space Shuttle Endeavour during an 11-day mission in February of 2000. SRTM is an international project spearheaded by the National Geospatial-Intelligence Agency (NGA) and the National Aeronautics and Space Administration (NASA).
#'
#' NASA Shuttle Radar Topography Mission (SRTM)(2013). Shuttle Radar Topography Mission (SRTM) Global.
#'
#' Distributed by OpenTopography.
#'
#'  https://doi.org/10.5069/G9445JDF
#'
#'  Accessed: 2020-12-04
#'
#' @format A Raster Layer
#' \describe{
#'   \item{dimensions}{1568 rows and 2258 columns}
#'   \item{resolution}{0.0002777778}
#'   \item{extent}{E3.495139, E4.122361, N43.96764, N44.40319}
#'   \item{CRS}{+proj=longlat +datum=WGS84 +no_defs}
#' }
#' @source \url{https://opentopography.org/}
"dem_cevennes"

#' Basins outlets positions in the Cevennes area SE France
#'
#' @format A SpatialPointsDataFrame
#' \describe{
#'   \item{features}{5}
#'   \item{variables}{2}
#'   \item{extent}{565391.4, 584650.7, 4881077, 4892994}
#'   \item{CRS}{+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs}
#' }
#'
"outlets_cevennes"

#' Profile line across the Cevennes area SE France
#'
#' @format A SpatialLinesDataFrame
#' \describe{
#'   \item{features}{1}
#'   \item{variables}{0}
#'   \item{extent}{552118.7, 587431.2, 4880279, 4914009}
#'   \item{CRS}{+proj=utm +zone=31 +datum=WGS84 +units=m +no_defs}
#' }
#'
"profile_cevennes"


