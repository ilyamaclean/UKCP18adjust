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
  a <- ukcp_psl$arraydata
  tme <- ukcp_psl$times + 12 * 3600
  tme <- as.POSIXlt(tme)
  sel <- .timesel(tme$year[1] + 1900)
  hrs <- sel[length(sel)] * 24 - 24
  # Correct times
  tme2 <- as.POSIXlt(c(0:hrs) * 3600, origin = tme[1], tz = "GMT")
  tme2 <- as.numeric(tme2)
  tme2 <- c(tme2[1] - c(24:1) * 3600, tme2, tme2[length(tme2)] + c(1:24) * 3600)
  tme2 <- as.POSIXlt(tme2, origin = "1970-01-01", tz = "GMT")
  tme <- as.numeric(tme)
  tme <- c(tme[1] - 24 * 3600, tme, tme[length(tme)] + 24 * 3600)
  tme <- as.POSIXlt(tme, origin = "1970-01-01", tz = "GMT")
  selt <- c(13:(length(tme2)-13))
  if(class(dem)[1] != "logical") {
    dd <- is_raster(dem)
    aj <- ((293 - 0.0065 * dd) / 293)^5.26
    dst <- "Elevation adjusted pressure"
  } else aj <- array(1, dim = dim(a)[1:2])
  cat("Spline interpolating to hourly\n")
  aout <- array(NA, dim = c(dim(a)[1:2], length(selt)))
  for (i in 1:dim(a)[1]) {
    for (j in 1:dim(a)[2]) {
      tst <- mean(a[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        psl <- a[i,j,]
        psl <- c(psl[1], psl, psl[length(psl)])
        psl <- psl * aj[i,j]
        xy <- spline(as.numeric(tme), psl, n = length(tme2))
        aout[i,j,] <- xy$y[selt]
      }
    }
  }
  xx <- round(aout, 0)
  xx <- array(as.integer(xx), dim = dim(aout))
  ao <- list(arraydata = xx, times = tme2[selt], crs = ukcp_psl$crs,
             extent = ukcp_psl$extent, units = ukcp_psl$units,
             description = dst)
  class(ao) <- "spatialarray"
  return(ao)
}
#' Calculate hourly corrected specific humidity from daily data
#'
#' @param ukcp_huss a `spatialarray` of daily specific humidity (Kg / Kg)
#' @param hus_gam a `gam` object of correction coefficients to apply to
#' specific humdidity data as returned by [gamcorrect()]
#' @param htemps optional spatial array of hourly temperatures (deg C - see details)
#' @param hpre optional spatial array  of hourly pressure (Pa - see details)
#' @details This function applies a correction to `ukcp_huss` using
#' parameters in `hus_gam` and then interpolates daily data to hourly.
#' If temperature data are supplied, hourly specific humdidities are conveted
#' to relative humidity, which is capped at 100 prior to converting back to
#' specific humdidity. Hourly pressure data are used in converting to relative
#' humidity. If not supplied, then pressure is assumed to be 101300 Pa.
#'
#' @return a spatial array of hourly specific humidities in Kg / Kg
#'
#' @seealso [gamcorrect()]
#' @import microclima zoo mgcv
#' @export
hourlyhuss <- function(ukcp_huss, hus_gam, htemps = NA, hpre = NA) {
  cat("Applying GAM correction to humidity \n")
  ukcp_huss <- .applygam(ukcp_huss, hus_gam)
  a <- ukcp_huss$arraydata
  tme <- ukcp_huss$times + 12 * 3600
  tme <- as.POSIXlt(tme)
  sel <- .timesel(tme$year[1] + 1900)
  hrs <- sel[length(sel)] * 24 - 24
  # Correct times
  tme2 <- as.POSIXlt(c(0:hrs) * 3600, origin = tme[1], tz = "GMT")
  tme2 <- as.numeric(tme2)
  tme2 <- c(tme2[1] - c(24:1) * 3600, tme2, tme2[length(tme2)] + c(1:24) * 3600)
  tme2 <- as.POSIXlt(tme2, origin = "1970-01-01", tz = "GMT")
  tme <- as.numeric(tme)
  tme <- c(tme[1] - 24 * 3600, tme, tme[length(tme)] + 24 * 3600)
  tme <- as.POSIXlt(tme, origin = "1970-01-01", tz = "GMT")
  selt <- c(13:(length(tme2)-13))
  aout <- array(NA, dim = c(dim(a)[1:2], length(selt)))
  cat("Spline interpolating to hourly\n")
  for (i in 1:dim(a)[1]) {
    for (j in 1:dim(a)[2]) {
      tst <- mean(a[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        hus <- a[i,j,]
        hus <- c(hus[1], hus, hus[length(hus)])
        xy <- spline(as.numeric(tme), hus, n = length(tme2))
        hus3 <- xy$y
        hus3 <- hus3[selt]
        if (class(htemps) != "logical") {
          tc <- htemps$arraydata[i,j,]
          if (class(hpre) != "logical") {
            pr <- hpre$arraydata[i,j,]
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
  ao <- list(arraydata = aout, times = tme2[selt], crs = ukcp_huss$crs,
             extent = ukcp_huss$extent, units = ukcp_huss$units,
             description = dts)
  class(ao) <- "spatialarray"
  return(ao)
}
#' Calculate hourly temperature from daily max and min
#'
#' @param huss a `spatialarray` of hourly specific humidities (Kg / Kg) as returned by [hourlyhuss()]
#' @param pre a `spatialarray` of hourly pressures (Pa) as returned by [hourlypsl()]
#' @param cfc a `spatialarray` of hourly percentage cloud cover as returned by [radsplit()]
#' @param dni a `spatialarray` of hourly direct normal irradiances (W / m^2) as returned by [radsplit()]
#' @param dif a `spatialarray` of hourly horizontal diffuse radiation (W / m^2) as returned by [radsplit()]
#' @param tmin a `spatialarray` of daily minimum temperature (deg C)
#' @param tmax a `spatialarray` of daily maximum temperature (deg C)
#' @param dem a raster of eleveations of the same resolution and extent as `huss`
#' @param dtr_gam a `gam` object of correction coefficients to apply to
#' diurnal temperature ranges as returned by [gamcorrect()]
#' @param tc_gam a `gam` object of correction coefficients to apply to
#' mean temperatures as returned by [gamcorrect()]
#' @details This function uses the [hourlytemp()] in the `microclima` package to compute
#' hourly temperatures. Percentage cloud cover is automatically converted to fractional
#' cloud cover, and units of radiation are automatically converted to MJ / m^2 / hr
#'
#' @return a `spatialarray` of hourly temperatures (deg C)
#'
#' @seealso [hourlyhuss()] [hourlypsl()] [radsplit()] [gamcorrect()]  [microclima::hourlytemp()]
#' @import microclima zoo mgcv raster
#' @export
hourlytc <- function(huss, pre, cfc, dni, dif, tmin, tmax, dem, dtr_gam, tc_gam) {
  dem60k <- aggregate(dem5k, 12)
  dem60k <- resample(dem60k, dem5k)
  demd <- dem5k - dem60k
  demd <- is_raster(demd)
  dys <- length(huss$times) / 24  - 1
  tme <- as.POSIXlt(c(0:dys) * 3600 * 24, origin = tmin$times[1], tz = "GMT")
  jd <- julday(tme$year + 1900, tme$mon +1, tme$mday)
  r <- raster(tmin$arraydata[,,1])
  extent(r) <- tmin$extent
  crs(r) <- tmin$crs
  lats <- .latsfromr(r)
  lons <- .lonsfromr(r)
  cat("Calculating hourly temperatures\n")
  sel <- .timesel(tme$year[1] + 1900)
  tm <- rep(NA, sel[length(sel)])
  tout <- array(NA, dim = dim(huss$arraydata))
  for (i in 1:dim(lats)[1]) {
    for (j in 1:dim(lats)[2]) {
      tst <- mean(tmin$arraydata[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        xy <- data.frame(x = lons[i,j], y = lats[i,j])
        coordinates(xy) = ~x + y
        proj4string(xy) = crs(r)
        ll <- as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
        tmn <- tm
        tmx <- tm
        tmn[sel] <- tmin$arraydata[i,j,]
        tmx[sel] <- tmax$arraydata[i,j,]
        tmn <- na.approx(tmn)
        tmx <- na.approx(tmx)
        dtr <- tmx - tmn
        tmean <- (tmn + tmx) / 2
        tmean <- predict.gam(tc_gam, newdata = data.frame(v2 = tmean))
        dtr <- predict.gam(dtr_gam, newdata = data.frame(v2 = dtr))
        tmn <- tmean - 0.5 * dtr
        tmx <- tmean + 0.5 * dtr
        ht <- hourlytemp(jd, em = NA, h = huss$arraydata[i,j,],
                         n = cfc$arraydata[i,j,] / 100,
                         p = pre$arraydata[i,j,],
                         dni = dni$arraydata[i,j,] * 0.0036,
                         dif = dif$arraydata[i,j,] * 0.0036,
                         mintemp = tmn,
                         maxtemp =tmx,
                         lat = ll$y, long = ll$x, merid = 0)
        lr <- lapserate(ht, huss$arraydata[i,j,], pre$arraydata[i,j,])
        tout[i,j,] <- ht + lr * demd[i,j]
      }
    }
  }
  ao <- list(arraydata = tout, times = huss$times, crs = huss$crs,
             extent = huss$extent, units = "Deg C",
             description = "Hourly temperature")
  class(ao) <- "spatialarray"
  return(ao)
}
#' Calculate hourly corrected wind speed and direction from daily data
#'
#' @param ukcp_uw a `spatialarray` of u vector of daily wind speed (m / s)
#' @param ukcp_uw a `spatialarray` of v vector of daily wind speed (m / s)
#' @param ws_gam a `gam` object of correction coefficients to apply to
#' wind speed data as returned by [gamcorrect()]
#' @details This function interpolates u and v wind vectors to hourly and then
#' calculates wind speed and direction. It then applies a correction
#' using parameters in `ws_gam`.
#'
#' @return a list of two spatialarrays
#' \describe{
#'   \item{ws}{A `spatialarray` of hourly wind speed (m /s)}
#'   \item{wd}{A `spatialarray` of hourly wind direction (degrees)}
#'}
#' @seealso [gamcorrect()]
#' @import microclima raster mgcv
#' @export
hourlywind <- function(ukcp_uw, ukcp_vw, ws_gam) {
  u <- ukcp_uw$arraydata
  v <- ukcp_vw$arraydata
  tme <- ukcp_uw$times + 12 * 3600
  tme <- as.POSIXlt(tme)
  sel <- .timesel(tme$year[1] + 1900)
  hrs <- sel[length(sel)] * 24 - 24
  # Correct times
  tme2 <- as.POSIXlt(c(0:hrs) * 3600, origin = tme[1], tz = "GMT")
  tme2 <- as.numeric(tme2)
  tme2 <- c(tme2[1] - c(24:1) * 3600, tme2, tme2[length(tme2)] + c(1:24) * 3600)
  tme2 <- as.POSIXlt(tme2, origin = "1970-01-01", tz = "GMT")
  tme <- as.numeric(tme)
  tme <- c(tme[1] - 24 * 3600, tme, tme[length(tme)] + 24 * 3600)
  tme <- as.POSIXlt(tme, origin = "1970-01-01", tz = "GMT")
  selt <- c(13:(length(tme2)-13))
  aout1 <- array(NA, dim = c(dim(u)[1:2], length(selt)))
  aout2 <- aout1
  cat("Spline interpolating to hourly\n")
  for (i in 1:dim(u)[1]) {
    for (j in 1:dim(u)[2]) {
      tst <- mean(u[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        uw <- u[i,j,]
        vw <- v[i,j,]
        uw <- c(uw[1], uw, uw[length(uw)])
        vw <- c(vw[1], vw, vw[length(vw)])
        xy1 <- spline(as.numeric(tme), uw, n = length(tme2))
        xy2 <- spline(as.numeric(tme), vw, n = length(tme2))
        uw2 <- xy1$y[selt]
        vw2 <- xy2$y[selt]
        ws <- sqrt(uw2^2 + vw2^2)
        ws <- predict.gam(ws_gam, newdata = data.frame(v2 = ws))
        wd <- atan2(uw2, vw2) * 180/pi + 180
        wd <- round(wd, 0)%%360
        aout1[i,j,] <- ws
        aout2[i,j,] <- wd
      }
    }
  }
  ao1 <- list(arraydata = aout1, times = tme2[selt], crs = ukcp_uw$crs,
              extent = ukcp_uw$extent, units = "m / s",
              description = "Wind speed")
  ao2 <- list(arraydata = aout2, times = tme2[selt], crs = ukcp_uw$crs,
              extent = ukcp_uw$extent, units = "Degrees",
              description = "Wind direction")
  class(ao1) <- "spatialarray"
  class(ao2) <- "spatialarray"
  return(list(ws = ao1, wd = ao2))
}
#' Saves data by year in specified directory
#' @param sa a `spatialarray` of data spanning the period in `yrs`
#' @param yrs a vector of years
#' @param dirout directory in which to save data. If it does not contain an absolute
#' path, the directory  name is relative to the current working directory, [getwd()].
#' Tilde-expansion is performed where supported.
#' @param filename file name of data to be saved. Each file is saved a sa `spatialarray`
#' with the following name: `filenameyear.R`
savebyyear <- function(sa, yrs, dirout, filename) {
  dir.create(dirout, showWarnings = F)
  tme <- sa$times
  a <- sa$arraydata
  if (length(tme) != dim(a)[3]) {
    stop("time vector length not equal to number of time entries in array")
  }
  for (yr in yrs) {
    sel <- which(tme$year + 1900 == yr)
    sao <- sa
    sao$times <- tme[sel]
    sao$arraydata <- a[,,sel]
    fo <- paste0(dirout, filename, yr, ".R")
    save(sao, file = fo)
  }
}


