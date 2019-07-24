#' Calculate hourly corrected sea-level pressure from daily data
#'
#' @param ukcp_psl a `spatialarray` of daily sea-level pressure (Pa)
#' @param psl_gam a `gam` object of correction coefficients to apply to
#' specific humdidity data as returned by [gamcorrect()]
#' @param dem an optional digital elevation dataset covering the same extent and
#' of the same resolution as data in ukcp_psl. Can be supplied as a raster or matrix
#' (see details)
#' @details This function applies a correction to `ukcp_psl` using
#' parameters in `psl_gam` and then interpolates daily data to hourly.
#' If digital elevation data are supplied, sea-level pressure is elevation corrected.
#'
#' @return a spatial array of hourly pressure in Pa
#'
#' @seealso [gamcorrect()]
#' @import raster microclima mgcv
#' @export
hourlypsl <- function(ukcp_psl, psl_gam, dem = NA) {
  cat("Applying GAM correction to pressure \n")
  ukcp_psl <- .applygam(ukcp_psl, psl_gam)
  tme <- ukcp_psl$times + 12 * 3600
  tme <- as.POSIXlt(tme)
  sel <- .timesel(tme$year[1] + 1900)
  hrs <- sel[length(sel)] * 24 - 24
  tme2 <- as.POSIXlt(c(0:hrs) * 3600, origin = tme[1], tz = "GMT")
  dst <- "Sea-level pressure"
  a <- ukcp_psl$arraydata
  if(class(dem)[1] != "logical") {
    dd <- is_raster(dem)
    aj <- ((293 - 0.0065 * dd) / 293)^5.26
    dst <- "Elevation adjusted pressure"
  } else aj <- array(1, dim = dim(a)[1:2])
  aout <- array(NA, dim = c(dim(a)[1:2], length(tme2)))
  cat("Spline interpolating to hourly\n")
  for (i in 1:dim(a)[1]) {
    for (j in 1:dim(a)[2]) {
      tst <- mean(a[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        psl <- a[i,j,]
        psl <- psl * aj[i,j]
        xy <- spline(as.numeric(tme), psl, n = length(tme2))
        aout[i,j,] <- xy$y
      }
    }
  }
  xx <- round(aout, 0)
  xx <- array(as.integer(xx), dim = dim(aout))
  ao <- list(arraydata = xx, times = tme2[sel], crs = ukcp_psl$crs,
             extent = ukcp_psl$extent, units = ukcp_psl$units,
             description = dst)
  class(ao) <- "spatialarray"
  return(ao)
}
#' Calculate hourly corrected specific humidity from daily data
#'
#' @param ukcp_huss a `spatialarray` of daily specific humdidity (Kg / Kg)
#' @param hus_gam a `gam` object of correction coefficients to apply to
#' pecific humdidity data as returned by [gamcorrect()]
#' @param htemps optional vector of hourly temperatures (deg C - see details)
#' @param hpre optional vector of hourly pressure (Pa - see details)
#' @details This function applies a correction to `ukcp_huss` using
#' parameters in `hus_gam` and then interpolates daily data to hourly.
#' If temperature data are supplied, hourly specific humdidities are conveted
#' to relative humidity, which is capped at 100 prior to converting back to
#' specific humdidity. Hourly pressure data are used in converting to relative
#' humidity. If not supplied, then pressure is assumed to be 101300 Pa. Not that
#' returned values are multiplied by 100,000 to renable storing as integers
#'
#' @return a spatial array of hourly specific humidities in Kg / Kg x 100,000 stored
#' as an integer.
#'
#' @seealso [gamcorrect()]
#' @import microclima zoo mgcv
#' @export
hourlyhuss <- function(ukcp_huss, hus_gam, htemps = NA, hpre = NA) {
  cat("Applying GAM correction to humidity \n")
  ukcp_huss <- .applygam(ukcp_huss, hus_gam)
  tme <- ukcp_huss$times + 12 * 3600
  tme <- as.POSIXlt(tme)
  sel <- .timesel(tme$year[1] + 1900)
  hrs <- sel[length(sel)] * 24 - 24
  tme2 <- as.POSIXlt(c(0:hrs) * 3600, origin = tme[1], tz = "GMT")
  a <- ukcp_huss$arraydata
  aout <- array(NA, dim = c(dim(a)[1:2], length(tme2)))
  cat("Spline interpolating to hourly\n")
  for (i in 1:dim(a)[1]) {
    for (j in 1:dim(a)[2]) {
      tst <- mean(a[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        hus <- a[i,j,]
        xy <- spline(as.numeric(tme), hus, n = length(tme2))
        hus3 <- xy$y
        if (class(htemps) != "logical") {
          tc <- htemps[i,j,]
          if (class(pre) != "logical") {
            pr <- hpre[i,j,]
          } else pr <- 101300
          rh <- suppressWarnings(humidityconvert(hus3, intype = "specific", tc, pr)$relative)
          rh[rh > 100] <- 100
          hus3 <- humidityconvert(rh, intype = "relative", tc, pr)$specific
        }
        aout[i,j,] <- hus3
      }
    }
  }
  dts <- "Hourly specific humidity"
  if (class(htemps) != "logical") dts <- "Hourly temperature adjusted specific humidity"
  xx <- round(aout, 0)
  xx <- array(as.integer(xx), dim = dim(aout))
  ao <- list(arraydata = xx, times = tme2[sel], crs = ukcp_psl$crs,
             extent = ukcp_huss$extent, units = ukcp_huss$units,
             description = dts)
  class(ao) <- "spatialarray"
  return(ao)
}


