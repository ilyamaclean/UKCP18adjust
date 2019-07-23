#' Calculate hourly radiation patchiness coefficients
#'
#' @param sis an object of class spatial array with hourly radiation values
#' @return a list containing the following components:
#' \describe{
#'   \item{m1}{A list of semivariagram model coefficients for the six-hourly deviations from daily values}
#'   \item{m2}{A list of semivariagram model coefficients for the hourly deviations from six-hourly values}
#' }
#' @import raster sp gstat microclima
#' @export
#' @seealso [hourlysis()]
#' @details Using [hourlysis()], hourly radiation values are computed by deriving
#' daily optical depths from daily radiation values, interpolating these to hourly
#' and back-calculating radiation for each hour. In so doing, cloud patchiness is
#' under-estimated, leading to potential errors in the calculation of direct and
#' diffuse radiation. The function [radfit()] derives model coefficients from real
#' hourly radiation data so that spatially and temporally autocorrelated cloud
#' patchiness can be simulated.
#'
#' @examples
#' modfit <- radfit(sis_2005)
#' modfit
#'
radfit <- function(sis) {
  a <- sis$arraydata
  r <- raster(a[,,1])
  extent(r) <- sis$extent
  crs(r) <- sis$crs
  lats <- .latsfromr(r)
  lons <- .lonsfromr(r)
  tme <- sis$times
  jd <- microclima::julday(tme$year + 1900, tme$mon + 1, tme$mday)
  aodh <- array(NA, dim = dim(a))
  aods <- aodh
  aodd <- aodh
  Nm <- lats * NA
  cat("Computing optical depths\n")
  for (i in 1:dim(lats)[1]) {
    for (j in 1:dim(lats)[2]) {
      xy <- data.frame(x = lons[i,j], y = lats[i,j])
      coordinates(xy) = ~x + y
      proj4string(xy) = crs(r)
      ll <- as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
      # Compute optical depth for global radiation
      ex <- a[i,j,] / 1353
      am <- airmasscoef(tme$hour, ll$y, ll$x, jd, merid = 0)
      od <- suppressWarnings(log(ex) / -am)
      od[od > 3] <- NA
      od[od <= 0] <- NA
      # Recalculate optical depth
      odd <- matrix(od, ncol = 24, byrow = T)
      odd <- apply(odd, 1, mean, na.rm = T)
      odd <- rep(odd, each = 24)
      lN <- -0.78287 -0.43582 * log(odd)
      N <- exp(lN)
      od2 <- suppressWarnings(log(ex) / -am ^ N)
      od2[od2 > 10] <- NA
      od2[od2 <= 0] <- NA
      aodh[i,j,] <- suppressWarnings(log(od2))
      odm <- matrix(aodh[i,j,], ncol = 24, byrow = T)
      od24 <- apply(odm, 1, mean, na.rm = T)
      aodd[i,j,] <- rep(od24, each = 24)
      odm <- matrix(aodh[i,j,], ncol = 6, byrow = T)
      od6 <- apply(odm, 1, mean, na.rm = T)
      aods[i,j,] <- rep(od6, each = 6)
    }
  }
  # Compute 6-hourly deviation from daily
  sod <- aodd - aods
  # Compute hourly deviation from six hourly
  hod <- aods - aodh
  cat("Computing semivariogram models\n")
  # Compute model coefficients sod:
  m1 <- .semivar(sod, r)
  m2 <- .semivar(hod, r)
  return(list(m1 = m1, m2 = m2))
}
#' calculate optical depth from radiation
#' @param rad `spatialarray` of total surface incoming radiation
#' @return a `spatialarray` of optical depths
#' @import raster sp microclima
#' @export
#' @examples
#' od <- opticaldepth(sis_2005)
opticaldepth <- function(rad) {
  a <- rad$arraydata
  tme <- rad$times
  r <- raster(a[,,1])
  extent(r) <- rad$extent
  crs(r) <- rad$crs
  lats <- .latsfromr(r)
  lons <- .lonsfromr(r)
  jd <- julday(tme$year + 1900, tme$mon + 1, tme$mday)
  mx <- max(a, na.rm = T)
  if (mx > 10) a <- a * 0.0036
  oda <- array(NA, dim = dim(a))
  for (i in 1:dim(a)[1]) {
    for (j in 1:dim(a)[2]) {
      tst <- mean(a[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        xy <- data.frame(x = lons[i,j], y = lats[i,j])
        coordinates(xy) = ~x + y
        proj4string(xy) = crs(r)
        ll <- as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
        am <- airmasscoef(tme$hour, ll$y, ll$x, jd, merid = 0)
        radmj <- a[i,j,]
        ext <- radmj / 4.87
        # Compute optical depth for global radiation
        od <- suppressWarnings(log(ext) / (- am))
        od[od > 4.5] <- NA
        od[od < 0] <- NA
        # Adjusting optical depths
        odd <- matrix(od, ncol = 24, byrow = T)
        odd <- apply(odd, 1, mean, na.rm = T)
        odd <- rep(odd, each = 24)
        lN <- -0.78287 -0.43582 * log(odd)
        N <- exp(lN)
        od <- suppressWarnings(log(ext) / (- am^N))
        od[od > 4.5] <- NA
        od[od < 1e-6] <- NA
        oda[i,j,] <- od
      }
    }
  }
  lst <- list(arraydata = oda, times = rad$times, crs = rad$crs,
              extent = rad$extent, units = "unitless", description = "Cloud optical depth")
  class(lst) <- "spatialarray"
  return(lst)
}
#' Calculate hourly radiation from daily radiation
#'
#' @param dailysis an object of class spatialarray with hourly radiation values
#' @param modfit a list of semivariagram model coefficients as returned by [radfit()]
#' @param Trace logical value indicating whether to plot progress
#' @return an object of class spatialarray containing the following components:
#' \describe{
#'   \item{arraydata}{a three-dimensional array of hourly radiation values}
#'   \item{times}{`POSIXlt` object of times associated with `arraydata`}
#'   \item{crs}{`crs` object of coordinate reference system associated with `arraydata`}
#'   \item{extent}{`extent` object giving extent covered by `arraydata`}
#'   \item{units}{units of `arraydata`}
#'   \item{description}{character description of `arraydata`}
#' }
#' @import raster microclima zoo
#' @export
#' @seealso [radfit()]
#' @details Hourly radiation values are computed by deriving daily optical depths
#' from daily radiation values, interpolating these to hourly and back-calculating
#' radiation for each hour. To accunt for sub-daily variation in cloud patchiness
#' spatially and temporally autocorrelated cloud patchiness is simulated using
#' paramaters from `modfit`.
#'
#' @examples
#' # Takes a few seconds to run
#' modfit <- radfit(sis_2005)
#' # Takes a few minutes to run
#' ukcpsishourly <- hourlysis(sis_ukcpdaily, modfit, 2000)
hourlysis <- function(dailysis, modfit, startyear, Trace = T) {
  a <- dailysis$arraydata
  tme <- dailysis$times
  r <- raster(a[,,1])
  extent(r) <- dailysis$extent
  crs(r) <- dailysis$crs
  lats <- .latsfromr(r)
  lons <- .lonsfromr(r)
  jd <- rep(julday(tme$year + 1900, tme$mon + 1, tme$mday), each = 24)
  lt <- rep(c(0:23), dim(a)[3])
  sel <- .timesel(startyear)
  nhrs <- sel[length(sel)] * 24
  odm <- array(NA, dim = c(dim(a)[1:2], nhrs))
  cat("Computing optical depths\n")
  for (i in 1:dim(lats)[1]) {
    for (j in 1:dim(lats)[2]) {
      tst <- mean(a[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        # Calculate lat and long
        xy <- data.frame(x = lons[i,j], y = lats[i,j])
        coordinates(xy) = ~x + y
        proj4string(xy) = crs(r)
        ll <- as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
        # Calculate optical depth
        am <- airmasscoef(lt, ll$y, ll$x, jd, merid = 0)
        radmj <- a[i,j,] * 0.0036
        radmj <- rep(radmj, each = 24)
        ext <- radmj / 4.87
        # Compute optical depth for global radiation
        od <- log(ext) / (- am)
        od[od > 10] <- NA
        od[od < 0] <- NA
        # Adjusting optical depths
        odd <- matrix(od, ncol = 24, byrow = T)
        odd <- apply(odd, 1, mean, na.rm = T)
        odd <- rep(odd, each = 24)
        lN <- -0.78287 -0.43582 * log(odd)
        N <- exp(lN)
        od <- log(ext) / (- am^N)
        od[od > 10] <- 10
        od[od < 1e-6] <- 1e-6
        # Calculate for every day
        odd <- matrix(od, ncol = 24, byrow = T)
        odd <- apply(odd, 1, mean, na.rm = T)
        # Set missing days to NA
        odd2 <- rep(NA, sel[length(sel)])
        odd2[sel] <- odd
        odd2 <- rep(odd2, each = 24)
        mult <- rep(c(rep(NA, 12), 1, rep(NA, 11)), sel[length(sel)])
        odd2 <- odd2 * mult
        odd2[1] <- odd2[13]
        odd2[length(odd2)] <- odd2[length(odd2) - 11]
        odd2 <- na.approx(odd2)
        odm[i,j,] <- log(odd2)
      }
    }
  }
  ####
  cat("Computing spatially autocorrelated cloud patchiness grids\n")
  xy <- expand.grid(1:dim(r)[2], 1:dim(r)[1])
  names(xy) <- c('x','y')
  g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                    model=vgm(psill = modfit$m1$psill, range = modfit$m1$rge, model='Sph'), nmax = 10)
  g.dummy2 <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                    model=vgm(psill = modfit$m2$psill, range = modfit$m2$rge, model='Sph'), nmax = 10)
  yy1 <- predict(g.dummy1, newdata=xy, nsim=dim(odm)[3]/6)
  yy2 <- predict(g.dummy2, newdata=xy, nsim=dim(odm)[3])
  cat("Applying spatially autocorrelated cloud patchiness grids\n")
  odm2 <- odm
  mt <- is_raster(r)
  mt <- mt * 0 + 1
  for (i in 1:dim(odm)[3]) {
    j <- floor(((i-1) / 6) + 1)
    xx1 <- data.frame(x = yy1$x, y = yy1$y, z = yy1[,j+2])
    xx2 <- data.frame(x = yy2$x, y = yy2$y, z = yy2[,i+2])
    r1 <- rasterFromXYZ(xx1)
    r2 <- rasterFromXYZ(xx2)
    m1 <- is_raster(r1)
    m1 <- m1 - mean(m1)
    m1 <- m1 * (modfit$m1$sdev / sd(m1))
    m2 <- is_raster(r2)
    m2 <- m2 - mean(m2)
    m2 <- m2 * (modfit$m2$sdev / sd(m2))
    m3 <- m2 + m1
    od <- suppressWarnings(log(odm[,,i]) + m3)
    od[is.na(od)] <- 0
    od <- od * mt
    odm2[,,i] <- exp(od)
    if(i%%2000 == 0 | i == dim(odm)[3] & Trace) {
      p <- round((i / dim(odm)[3]) * 100, 1)
      p <- paste0(p, "%")
      ro <- raster(odm2[,,i], template = r)
      plot(ro, main = p)
    }
  }
  cat("Computing hourly shortwave radiation\n")
  radm <- odm2
  scs <- seq(as.numeric(tme[1]), as.numeric(tme[3600]) + 23 * 3600, by = 3600)
  tme2 <- as.POSIXlt(scs, origin = "1970-01-01 00:00", tz = "GMT")
  jd <- julday(tme2$year + 1900, tme2$mon + 1, tme2$mday)
  for (i in 1:dim(lats)[1]) {
    for (j in 1:dim(lats)[2]) {
      tst <- mean(a[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        # Calculate lat and long
        xy <- data.frame(x = lons[i,j], y = lats[i,j])
        coordinates(xy) = ~x + y
        proj4string(xy) = crs(r)
        ll <- as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
        # Calculate optical depth
        am <- airmasscoef(tme2$hour, ll$y, ll$x, jd, merid = 0)
        odh <- exp(odm2[i,j,])
        # Adjusting optical depths
        odd <- matrix(odh, ncol = 24, byrow = T)
        odd <- apply(odd, 1, mean, na.rm = T)
        odd <- rep(odd, each = 24)
        lN <- -0.78287 -0.43582 * log(odd)
        N <- exp(lN)
        radh <- exp(-am^N * odh) * 4.87
        radh[is.na(radh)] <- 0
        radm[i,j,] <- radh / 0.0036
      }
    }
  }
  xx <- radm
  xx <- round(xx, 0)
  xx <- array(as.integer(xx), dim = dim(radm))
  lst <- list(arraydata = radm, times = tme2, crs = dailysis$crs,
       extent = dailysis$extent, units = "Watts / m^2",
       description = "Total incoming shortwave radiation")
  class(lst) <- "spatialarray"
  return(lst)
}
#' Calculate and save corrected radiation by year
#'
#' @param ukcpsishourly a `spatialarray` of hourly uncorrected
#' surface incoming solar irradiance (W / m^2) as returned by [hourlysis()]
#' @param rad_gam a `gam` object of correction coefficients to apply to
#' radiation data
#' @param years years for which data are required. Must overlap with years
#' in `ukcpsishourly`
#' @param destfolder folder in which to save data (see details).If it does
#' not contain an absolute path, the file name is relative to the
#' current working directory, [getwd()]. Tilde-expansion is performed where
#' supported.
#'
#' @details This function applies a correction to `ukcpsishourly` using
#' parameters in `rad_gam`. Total surface incoming solar irradiance is then
#' aparitioned between direct and diffuse, and direct irradiance normal to
#' solar beam computed. A cloud cover coefficient is also computed from the
#' ratio of received to clearsky radiation. The resulting data are saved seperately
#' by year as `spatialarrays` in three subfolders of `destfolder`.
#'
#' @seealso [hourlysis()] [gamcorrect()] [microclima::difprop()]
#' @import microclima zoo mgcv
#' @export
#'
#' @examples
#' # Takes a few seconds to run
#' modfit <- radfit(sis_2005)
#' # Takes a few minutes to run
#' ukcpsishourly <- hourlysis(sis_ukcpdaily, modfit, 2000)
#' yrs <- c(2001:2009)
#' radbyyear(ukcpsishourly, rad_gam, yrs, destfolder = "C:/radiation/")
radbyyear <- function(ukcpsishourly, rad_gam, years, destfolder) {
  dir.create(destfolder, showWarnings = FALSE)
  a <- ukcpsishourly$arraydata
  tme <- ukcpsishourly$times
  r <- raster(a[,,1])
  extent(r) <- ukcpsishourly$extent
  crs(r) <- ukcpsishourly$crs
  lats <- .latsfromr(r)
  lons <- .lonsfromr(r)
  jd <- julday(tme$year + 1900, tme$mon + 1, tme$mday)
  d.dni <- paste0(destfolder, "dni/")
  d.dif <- paste0(destfolder, "dif/")
  d.cfc <- paste0(destfolder, "cfc/")
  dir.create(d.dni, showWarnings = FALSE)
  dir.create(d.dif, showWarnings = FALSE)
  dir.create(d.cfc, showWarnings = FALSE)
  for (yr in years) {
    sel <- which(tme$year + 1900 == yr)
    ay <- a[,,sel]
    dni <- ay
    dif <- ay
    cfc <- ay
    for (i in 1:dim(a)[1]) {
      for (j in 1:dim(a)[2]) {
        tst <- mean(ay[i,j,], na.rm = T)
        if (is.na(tst) == F) {
          # Calculate lat and long
          xy <- data.frame(x = lons[i,j], y = lats[i,j])
          coordinates(xy) = ~x + y
          proj4string(xy) = crs(r)
          ll <- as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
          rad <- ay[i,j,]
          sis <- predict.gam(rad_gam, newdata = data.frame(v2 = rad))
          sis[sis > 1352.778] <- 1352.778
          dp <- difprop(sis, jd[sel], tme$hour[sel], ll$y, ll$x, hourly = TRUE,
                        watts = TRUE, merid = 0)
          si <- siflat(tme$hour[sel], ll$y, ll$x, jd[sel], merid = 0)
          pdif <- dp * sis
          pdni <- ((1 - dp) * sis) / si
          pdni[is.na(pdni)] <- 0
          pdni[pdni > 1352.778] <- 1352.778
          pdif[pdif > 1352.778] <- 1352.778
          # Compute maximum potential radiation
          # Calculate optical depth
          ext <- sis / 1352.778
          am <- airmasscoef(tme$hour[sel], ll$y, ll$x, jd[sel], merid = 0)
          od <- log(ext) / (- am)
          # Adjusting optical depths
          odd <- matrix(od, ncol = 24, byrow = T)
          odd <- apply(odd, 1, mean, na.rm = T)
          odd <- rep(odd, each = 24)
          lN <- -0.78287 -0.43582 * log(odd)
          N <- exp(lN)
          od <- suppressWarnings(log(ext) / (- am^N))
          od[od > 10] <- NA
          od[od < 0] <- NA
          # Maximum radiation
          mdir <- si * 1352.778
          n <- sis / mdir
          n[n > 1] <- NA
          n[n < 0] <- 0
          le <- length(n)
          n[1] <- mean(n[1:24], na.rm = T)
          n[le] <- mean(n[(le-24):le], na.rm = T)
          n <- na.approx(n)
          dni[i,j,] <- pdni
          dif[i,j,] <- pdif
          cfc[i,j,] <- 1 - n
        }
      }
    }
    fo1 <- paste0(d.dni, "dni_", yr, ".R")
    fo2 <- paste0(d.dif, "dif_", yr, ".R")
    fo3 <- paste0(d.cfc, "cfc_", yr, ".R")
    xx1 <- round(dni, 0)
    xx2 <- round(dif, 0)
    xx3 <- round(cfc * 100, 0)
    xx1 <- array(as.integer(xx1), dim = dim(dni))
    xx2 <- array(as.integer(xx2), dim = dim(dif))
    xx3 <- array(as.integer(xx3), dim = dim(cfc))
    ao1 <- list(arraydata = xx1, times = tme[sel], crs = ukcpsishourly$crs,
                extent = ukcpsishourly$extent, units = "Watts / m^2",
                description = "Direct Normal irradiance")
    ao2 <- list(arraydata = xx2, tme = tme[sel], crs = ukcpsishourly$crs,
                extent = aaukcpsishourly$extent, units = "Watts / m^2",
                description = "Diffuse horizontal irradiance")
    ao3 <- list(arraydata = xx3, tme = tme[sel], crs = ukcpsishourly$crs,
                extent = ukcpsishourly$extent, units = "Percentage",
                description = "Cloud cover")
    class(ao1) <- "spatialarray"
    class(ao2) <- "spatialarray"
    class(ao3) <- "spatialarray"
    save(ao1, file = fo1)
    save(ao2, file = fo2)
    save(ao3, file = fo3)
  }
}
