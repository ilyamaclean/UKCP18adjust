#' 5 km raster of Cornwall
#'
#' A 5 km resolution raster of Cornwall with land coded as 1 and sea as NA.
#'
#' @format A raster with 23 rows and 34 columns:
"msk"
#' General Additive Model cofficients for correcting UKCP18 radiation data
#'
#' An object of class `gam` containing coefficients for correcting UKCP18 radiation
#' data, as returned by [gamcorrect()].
"rad_gam"
#' Hourly surface incoming shortwave radiation (observed).
#'
#' An object of class `spatialarray` containing hourly radiation values in Cornwall for
#' the whole of 2005
#'
#' @format A list with the following components:
#' \describe{
#'   \item{arraydata}{a three-dimensional array of hourly radiation values}
#'   \item{times}{`POSIXlt` object of times associated with `arraydata`}
#'   \item{crs}{`crs` object of coordinate reference system associated with `arraydata`}
#'   \item{extent}{`extent` object giving extent covered by `arraydata`}
#'   \item{units}{units of `arraydata`}
#'   \item{description}{character description of `arraydata`}
#' }
#' @source \url{https://www.cmsaf.eu/}
"sis_2005"
#' Daily surface incoming shortwave radiation (UKCP18).
#'
#' An object of class spatialarray containing daily radiation values in Cornwall for
#' the period 2000-12-01 to 2010-11-30
#'
#' @format A list with the following components:
#' \describe{
#'   \item{arraydata}{a three-dimensional array of daily radiation values}
#'   \item{times}{`POSIXlt` object of times associated with `arraydata`}
#'   \item{crs}{`crs` object of coordinate reference system associated with `arraydata`}
#'   \item{extent}{`extent` object giving extent covered by `arraydata`}
#'   \item{units}{units of `arraydata`}
#'   \item{description}{character description of `arraydata`}
#' }
#' @source \url{http://data.ceda.ac.uk/badc/ukcp18/}
"sis_ukcpdaily"