#' Digital elevation dataset for Cornwall
#'
#' A 5 km resolution raster of elevations in Cornwall
#'
#' @format A raster with 23 rows and 34 columns:
"dem5k"
#' Digital elevation dataset for UK
#'
#' A 60 km resolution raster of elevations in the UK
#'
#' @format A raster with 21 rows and 12 columns:
"dem60k"
#' Elevation differences between 60 km and 5km resolution data
#'
#' A 5 km resolution raster of elevation differences between 60 km and 5km resolution data
#'
#' @format A raster with 23 rows and 34 columns:
"demdif"
#' General Additive Model cofficients for correcting UKCP18 diurnal temperature range data
#'
#' An object of class `gam` containing coefficients for correcting UKCP18 diurnal temperature range
#' data, as returned by [gamcorrect()].
"dtr_gam"
#' General Additive Model cofficients for correcting UKCP18 specific humidity data
#'
#' An object of class `gam` containing coefficients for correcting UKCP18 specific
#' humidity data, where specific humdidity in Kg / Kg has been multiplied by 100,000,
#' as returned by [gamcorrect()].
"hus_gam"
#' General Additive Model cofficients for correcting UKCP18 sea-level pressure data
#'
#' An object of class `gam` containing coefficients for correcting UKCP18 diurnal temperature range
#' data, as returned by [gamcorrect()].
"psl_gam"
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
#' Daily maximum temperature (UKCP18).
#'
#' An object of class spatialarray containing daily maximum temperature values in Cornwall for
#' the period 2000-12-01 to 2010-11-30
#'
#' @format A list with the following components:
#' \describe{
#'   \item{arraydata}{a three-dimensional array of daily maximum temperature values}
#'   \item{times}{`POSIXlt` object of times associated with `arraydata`}
#'   \item{crs}{`crs` object of coordinate reference system associated with `arraydata`}
#'   \item{extent}{`extent` object giving extent covered by `arraydata`}
#'   \item{units}{units of `arraydata`}
#'   \item{description}{character description of `arraydata`}
#' }
#' @source \url{http://data.ceda.ac.uk/badc/ukcp18/}
"tmx_ukcpdaily"
#' General Additive Model cofficients for correcting UKCP18 temperature data
#'
#' An object of class `gam` containing coefficients for correcting UKCP18 temperatures
#' data, as returned by [gamcorrect()].
"tc_gam"
#' Daily minimum temperature (UKCP18).
#'
#' An object of class spatialarray containing daily minimum temperature values in Cornwall for
#' the period 2000-12-01 to 2010-11-30
#'
#' @format A list with the following components:
#' \describe{
#'   \item{arraydata}{a three-dimensional array of daily minimum temperature values}
#'   \item{times}{`POSIXlt` object of times associated with `arraydata`}
#'   \item{crs}{`crs` object of coordinate reference system associated with `arraydata`}
#'   \item{extent}{`extent` object giving extent covered by `arraydata`}
#'   \item{units}{units of `arraydata`}
#'   \item{description}{character description of `arraydata`}
#' }
#' @source \url{http://data.ceda.ac.uk/badc/ukcp18/}
"tmn_ukcpdaily"
#' General Additive Model cofficients for correcting UKCP18 wind speed data
#'
#' An object of class `gam` containing coefficients for correcting UKCP18 wind speed
#' data, as returned by [gamcorrect()].
"ws_gam"
