#' Perform thin-plate spline with rasters
#' @import raster fields
#' @export
.rasterTps <- function(y, rc1, rf1, rc2 = NA, rf2 = NA, rc3 = NA, rf3 = NA) {
  xy <- data.frame(xyFromCell(rc1, 1:ncell(rc1)))
  z1 <- extract(rc1, xy)
  xyz <- cbind(xy, z1)
  if (class(rc2) != "logical")  xyz$z2 <- extract(rc2, xy)
  if (class(rc3) != "logical")  xyz$z3 <- extract(rc3, xy)
  v <- extract(y, xy)
  sel <- which(is.na(v) == F)
  xyz <- xyz[sel,]
  v <- v[sel]
  sel <- which(is.na(xyz$z1) == F)
  xyz <- xyz[sel,]
  v <- v[sel]
  if (class(rc2) != "logical") {
    sel <- which(is.na(xyz$z2) == F)
    xyz <- xyz[sel,]
    v <- v[sel]
  }
  if (class(rc3) != "logical") {
    sel <- which(is.na(xyz$z3) == F)
    xyz <- xyz[sel,]
    v <- v[sel]
  }
  tps <- suppressWarnings(fields::Tps(xyz, v, m = 3))
  xy <- data.frame(xyFromCell(rf1, 1:ncell(rf1)))
  z1 <- extract(rf1, xy)
  xyz <- cbind(xy, z1)
  if (class(rf2) != "logical")  xyz$z2 <- extract(rf2, xy)
  if (class(rf3) != "logical")  xyz$z3 <- extract(rf3, xy)
  sel <- which(is.na(xyz$z1) == F)
  xyz <- xyz[sel,]
  if (class(rf2) != "logical") {
    sel <- which(is.na(xyz$z2) == F)
    xyz <- xyz[sel,]
  }
  if (class(rf3) != "logical") {
    sel <- which(is.na(xyz$z3) == F)
    xyz <- xyz[sel,]
  }
  xy <- data.frame(x = xyz$x, y = xyz$y)
  xy$z <- NA
  xy$z <- fields::predict.Krig(tps, xyz)
  r <- rasterFromXYZ(xy)
  r
}
#' Order rain from lowest to highest
#' @export
.orderain <- function(rd, v1, v2) {
  o0 <- order(rd)
  out <- o0
  sel <- which(rd[o0] == 0)
  if (length(sel) > 1) {
    p1 <- v1[o0[sel]]
    o1 <- order(p1)
    for (k in 1:length(sel)) out[k] <- o0[o1[k]]
  }
  sel <- which(rd[out] == 0 & v1[out] == 0)
  out2 <- out
  if (length(sel) > 1) {
    p2 <- v2[out[sel]]
    o2 <- order(p2)
    for (k in 1:length(sel)) out2[k] <- out[o2[k]]
  }
  out2
}
#' Crop rainfall, convert to daily and resample to desired extent and resolution
#'
#' @param file the name of the .nc file which the data are to be read from.
#' If it does not contain an absolute path, the file name is relative to the
#' current working directory, [getwd()]. Tilde-expansion is performed where
#' supported.
#' @param r an object of type raster of elevations. The returned `spatialarray` has the same
#' resolution and extent as `r`.
#' @param Trace optional logical indicating whether to create and plot
#' a raster of every 100th entry to enable progress to be tracked.
#'
#' @return An object of class `spatialarray` containing the following components:
#' \describe{
#'   \item{arraydata}{A three-dimensional array of daily rainfall values}
#'   \item{times}{An object of class `POSIXlt` of times `arraydata`}
#'   \item{crs}{An object of class `crs` indicating the coordinate reference system associated with `arraydata`}
#'   \item{extent}{An object of class `extent` indicating the geographic extent covered by `arraydata`}
#'   \item{units}{Units of `arraydata`}
#'   \item{description}{Description of `arraydata`}
#' }
#' @details This function crops and resamples the rainfall data in `file` to the extent of
#' `r`. It also performs elevation adjustments by applying a thin-plate spline model with
#' elevation as a covariate to calculate the expected total rainfall and number of zero
#' rainfall days. The resampled data is then adjusted accordingly. The 3600 values for each
#' grid cell location are then spline interpolated to give daily values. Currently, the
#' UKCP18 nc files contain errors in both the times and the defined coordinate reference
#' system. These are corrected in the function, but users are advised to check their data.
#' @import raster ncdf4 zoo microclima
#' @export
cropresamplerain <- function(file, r, Trace = TRUE) {
  rain <- .nctoarray(file)
  nc <- nc_open(file)
  tme <- ncvar_get(nc, varid = 'time')
  tme <- as.POSIXlt(tme * 3600, origin = "1970-01-01", tz = "GMT")
  nc_close(nc)
  tme <- .tmecreate(tme$year[1] + 1900)
  # Calculate total and proportion rain days
  train <- apply(rain, c(1,2), sum)
  rainyn <- ifelse(rain > 0, 1, 0)
  prain <- apply(rainyn, c(1,2), sum) / 3600
  prain[prain == 0] <- NA
  train[train == 0] <- NA
  prain <- na.approx(prain)
  train <- na.approx(train)
  rte <- suppressWarnings(raster(file))
  crs(rte) <- crs(dtm100m)
  train <- if_raster(train, rte)
  prain <- if_raster(prain, rte)
  train <- resample(train, dem60k)
  prain <- resample(prain, dem60k)
  train <- mask(train, dem60k)
  prain <- mask(prain, dem60k)
  train5k <- .rasterTps(train, dem60k, r)
  prain5k <- .rasterTps(prain, dem60k, r)
  aout <- array(NA, dim = c(dim(r)[1:2], dim(rain)[3]))
  # Calculate resampled array for desired area
  for (i in 1:dim(rain)[3]) {
    ro <- if_raster(rain[,,i], rte)
    ro <- aggregate(ro, 5)
    ro <- resample(ro, r)
    ro <- mask(ro, r)
    aout[,,i] <- is_raster(ro)
    if(Trace & i%%100 == 0) plot(ro, main = tme[i])
  }
  # Extend from 3600 to daily
  sq <- .timesel(tme$year[1] + 1900)
  v1 <- apply(aout, 3, mean, na.rm = T)
  v2 <- apply(rain, 3, mean, na.rm = T)
  aout2 <- array(NA, dim = c(dim(aout)[1:2], sq[length(sq)]))
  for (i in 1:dim(aout)[1]) {
    for (j in 1:dim(aout)[2]) {
      tst <- mean(aout[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        rd <- aout[i,j,]
        # Sort out zeros
        pzero <- round((1- (is_raster(prain5k)[i,j])) * 3600, 0)
        azero <- length(which(rd == 0))
        o <- .orderain(rd, v1, v2)
        if (pzero > azero) {
          rd[o[1:pzero]] <- 0
        } else {
          # How many zeros need replacing
          rplc <- azero - pzero
          # Generate random sequence
          rnz <- rd[rd > 0]
          mn <- mean(log(rnz))
          sde <- sd(log(rnz))
          rrd <- exp(rnorm(length(rnz), mn, sde))
          rrd <- rrd[order(rrd)][1:rplc]
          # replace requisite number of zeros with lowest values from random sequence
          rd[o[(pzero + 1):azero]] <- rrd
        }
        # Adjust by total
        adjt <- is_raster(train5k)[i,j] / sum(rd)
        rd <- rd * adjt
        # Lengthen to daily
        rd2 <- rep(NA, sq[length(sq)])
        rd2[sq] <- rd
        aout2[i,j,] <- na.approx(rd2)
      }
    }
  }
  xx <- c(1:sq[length(sq)]) - 1
  tme2 <- as.POSIXlt(xx * 3600 * 24, origin = tme[1], tz = "GMT")
  lst <- list(arraydata = aout2, times = tme2, crs = crs(r),
              extent = extent(r), units = "mm / day", description = "Daily precipitation")
  return(lst)
}
#' Calculates ratio of observed to ukcp18 days with zero precipitation
#' @param observed_pre a `spatialarray` of observed daily precipitation values
#' @param ukcp_pre a `spatialarray` of ukcp daily precipitation values for the period
#' and extent corresponding to `observed`
#' @return a single numeric value of the proportion of days with zero precipitation in
#' the observed dataset divided by the corresponding proportion of days in the UKCP18 data
#' @export
rainzero_adj <- function(observed_pre, ukcp_pre) {
  a <- observed_pre$arraydata
  v <- as.vector(a)
  v <- v[is.na(v) == F]
  sel <- which(v == 0)
  pz1 <- length(sel) / length(v)
  a <- ukcp_pre$arraydata
  v <- as.vector(a)
  v <- v[is.na(v) == F]
  sel <- which(v == 0)
  pz2 <- length(sel) / length(v)
  return(pz1 / pz2)
}
#' Adjusts UKCP18 precipitation data to give expected proportion of days with zero precipitation
#'
#' @param ukcp_pre a `spatialarray` of UKCP18 daily rainfall data as returned by [cropresamplerain()]
#' @param rzero an optional single numeric value corresponding to the  ratio of observed to ukcp18
#' days with zero precipitation as returned by [rainzero_adj()].
#'
#' @details if no value for `rzero` is given, a default value derived for Cornwall, UK for
#' the period 2000-12-01 to 2010-11-30 is used. If `rzero` > 1 (likely to be almost always true),
#' then the requisite number of zeros is assigned to days with the lowest rainfall. If `rzero` < 1
#' then a dataset of precipitation values with the same mean and standard deviation of non-zero
#' precipitations values is generated, and the requisite number of days with zero rainfall are
#' assigned the lowest randomly generated precipitations values.
#'
#' @export
rainapplyzero <- function(ukcp_pre, rzero = 18.59483) {
  a <- ukcp_pre$arraydata
  rd <- as.vector(a)
  s1 <- which(is.na(rd) == F)
  # Calculate order
  v <- apply(a, 3, mean, na.rm = T)
  m <- a[,,1]
  s2 <- which(is.na(m) == F)
  v <- rep(v, each = length(s2))
  # Adjust order
  o0 <- order(rd[s1])
  out <- o0
  sel <- which(rd[s1[o0]] == 0)
  if (length(sel) > 1) {
    p1 <- v[o0[sel]]
    o1 <- order(p1)
    for (k in 1:length(sel)) out[k] <- o0[o1[k]]
  }
  azero <- length(which(rd[s1] == 0))
  przero <- round(azero * rzero, 0)
  if (rzero >= 1) {
    rd[s1[o0[1:przero]]] <- 0
  } else {
    # How many zeros need replacing
    rplc <- azero - przero
    # Generate random sequence
    rnz <- rd[s1]
    rnz <- rnz[rnz > 0]
    mn <- mean(log(rnz))
    sde <- sd(log(rnz))
    rrd <- exp(rnorm(length(rnz), mn, sde))
    rrd <- rrd[order(rrd)][1:rplc]
    # replace requisite number of zeros with lowest values from random sequence
    rx <- rd[s1]
    rx[o0[(przero + 1):azero]] <- rrd
    rd[s1] <- rx
  }
  ao <- array(rd, dim = dim(a))
  lst <- ukcp_pre
  lst$arraydata <- ao
  return(lst)
}
#' Adjusts ukcp daily precipitation data to make it more consistent with observed data
#'
#' @param ukcp_pre a `spatialarray` of UKCP18 daily rainfall data as returned by [cropresamplerain()]
#' @param rzero an optional single numeric value corresponding to the  ratio of observed to ukcp18
#' days with zero precipitation as returned by [rainzero_adj()].
#' @param pre_gam a `gam` object of correction coefficients to apply to
#' precipitation data as returned by [gamcorrect()
#'
#' @details This function applies a correction to `ukcp_pre` using
#' parameters in `pre_gam` and then interpolates daily data to hourly.
#'
#' @return a `spatialarray` of daily rainfall in mm
#'
#' @seealso [gamcorrect()]
#' @export
rain_adj <- function(ukcp_pre, rzero = 18.59483, pre_gam) {
  ukcp <- rainapplyzero(ukcp_pre, rzero)
  lukcp <- ukcp
  lukcp$arraydata <- log(lukcp$arraydata)
  lukcp$arraydata[lukcp$arraydata < -999] <- NA
  lukcp <- .applygam(lukcp, pre_gam)
  lukcp$arraydata <- exp(lukcp$arraydata)
  lukcp$arraydata[ukcp$arraydata == 0] <- 0
  return(lukcp)
}
