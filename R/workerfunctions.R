#' Resample nc file to spatial array
#'
#' @param file the name of the .nc file which the data are to be read from.
#' If it does not contain an absolute path, the file name is relative to the
#' current working directory, [getwd()]. Tilde-expansion is performed where
#' supported.
#' @param r an object of type raster. The returned `spatialarray` has the same
#' resolution and extent as `r`.
#' @param Trace optional logical indicating whether to create and plot
#' a raster of every 100th entry to enable progress to be tracked.
#'
#' @return an object of class `spatialarray`. An object of class `spatialarray` is
#' a list containing the following components:
#' \describe{
#'   \item{arraydata}{A three-dimensional array of values corresponding to values of the resampled .nc file}
#'   \item{times}{An object of class `POSIXlt` of times `arraydata`}
#'   \item{crs}{An object of class `crs` indicating the coordinate reference system associated with `arraydata`}
#'   \item{extent}{An object of class `extent` indicating the geographic extent covered by `arraydata`}
#'   \item{units}{Units of `arraydata`}
#'   \item{description}{Description of `arraydata`}
#' }
#' @import ncdf4 raster
#' @export
cropresamplenc <- function(file, r, Trace = TRUE) {
  nc <- nc_open(file)
  tme <- ncvar_get(nc, varid = 'time')
  tme <- as.POSIXlt(tme * 3600, origin = "1970-01-01", tz = "GMT")
  a <- .nctoarray(file)
  rx <- suppressWarnings(raster(file, band = 1))
  ecrop <- extent(r)
  e <- extent(rx)
  ae <- .croparray(a, e, ecrop)
  ao <- .resamplearray(ae, r, tme, Trace)
  ncv <- nc$var
  varname <- names(ncv)[1]
  nca <- ncatt_get(nc, varid = varname)
  nc_close(nc)
  lst <- list(arraydata = ao, times = tme, crs = crs(r),
              extent = ecrop, units = nca[[3]], description = nca[[2]])
  class(lst) <- "spatialarray"
  return(lst)
}
#' Fit General Additive Model to correct UKCP18 data
#'
#' @param observed object of class `spatialarray` of observed climate data
#' @param ukcp object of class `spatialarray` of UKCP18 climate data for the
#' corresponding period and location
#' @param Trace optional logical indicating whether to plot model fit
#' @param positive optional logical indicating whether all values are positive (as
#' with e.g. radiation)
#' @return an object of class `gam`.
#' @details This function ranks the values of `observed` and `ukcp` and fits
#' a General Additive Model (GAM) to these data. The area and time period covered by
#' `observed` and `ukcp`, and the spatial resolution of the data should be
#' identical. If the datasets contain more than 10,000 values then 10000 values
#' spanning the full range of data are selected. If `positive` is TRUE, then only
#' values greater than 0 are selected when fitting the GAM.
#' @import mgcv
#' @export
#' @examples
#' # Takes a few seconds to run
#' modfit <- radfit(sis_2005)
#' # Takes a few minutes to run
#' ukcpsishourly <- hourlysis(sis_ukcpdaily, modfit, 2000)
#' # Select data for 2005
#' tme <- ukcpsishourly$times
#' sel <- which(tme$year + 1900 == 2005)
#' xx <-  ukcpsishourly$arraydata[,,sel]
#' ukcp_2005 <- list(arraydata = xx, times = tme[sel], crs = ukcpsishourly$crs,
#'                  extent = ukcpsishourly$extent, units = ukcpsishourly$array.units,
#'                  description = ukcpsishourly$array.description)
#' class(ukcp_2005) <- "spatialarray"
#' od_obs <- opticaldepth(sis_2005)
#' od_pre <- opticaldepth(ukcp_2005)
#' gam.mod <- gamcorrect(od_obs, od_pre)
gamcorrect <- function(observed, ukcp, Trace = TRUE, positive = TRUE) {
  v1 <- as.vector(observed$arraydata)
  v2 <- as.vector(ukcp$arraydata)
  v1 <- v1[is.na(v1) == F]
  v2 <- v2[is.na(v2) == F]
  if (positive) {
    v1 <- v1[v1 > 0]
    v2 <- v2[v2 > 0]
  }
  v1 <- v1[order(v1)]
  v2 <- v2[order(v2)]
  le1 <- length(v1)
  le2 <- length(v2)
  lgt <- min(le1, le2)
  if (lgt > 10000) {
    s1 <- floor(seq(101, le1 - 100, length.out = 9800))
    s2 <- floor(seq(101, le2 - 100, length.out = 9800))
    v1 <- c(v1[1:100], v1[s1], v1[(le1-100):le1])
    v2 <- c(v2[1:100], v2[s2], v2[(le2-100):le2])
  }
  m1 <- gam(v1~s(v2))
  pred <- predict.gam(m1, newdata = data.frame(v2 = v2))
  if (Trace) {
    par(mar=c(7,7,5,2))
    plot(v1 ~ v2, xlim = c(0,max(v2)), ylim = c(0,max(v1)),
         xlab = "UKCP18 data",
         ylab = "Observed data", cex.axis = 2, cex.lab = 2, pch = 15,
         cex = 1, col = "gray")
    par(new=T)
    plot(pred ~ v2, xlim = c(0,max(v2)), ylim = c(0,max(v1)),
         type = "l", lwd = 2, col = "red",
         xlab = "", ylab = "", cex.axis = 2, cex.lab = 2)
  }
  return(m1)
}
