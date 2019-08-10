#' Derive parameters relating land to sea temperatures
#'
#' @param temps a data.frame of hourly temperatures with the following columns:
#' \describe{
#'   \item{obs_time}{an object of class `POSIXlt` of observation times}
#'   \item{land}{land temperatures (Deg C)}
#'   \item{sst}{Sea-surface temperatures (Deg C)}
#' }
#' @param Trace an optional logical indicating whether to plot observed versus predicted
#' sea-surface temperatures
#'
#' @return a data frame with the following columns
#' \describe{
#'   \item{x1}{parameter relating land-sea temperature difference to sea-surface temperature change}
#'   \item{x2}{parameter relating mean sea-surfacetemperatures to sea-surface temperature change}
#'   \item{x3}{a constant}
#'   \item{k}{mean sea-surface temperature - mean land temperature}
#'}
#' @seealso [sst_predict()]
#'
#' @examples
#' sstparams <- sst_fit(temperature2017)
#' sstparams
sst_fit <- function(temps, Trace = TRUE) {
  fitone <- function(land, sst, x1, x2, x3, out = "x") {
    pred <- sst[1]
    sstdeep <- mean(sst)
    for (i in 2:length(sst)) {
      heatdown <- (land[i - 1] - pred[i - 1])
      heatup <- (sstdeep - pred[i - 1])
      pred[i] <- pred[i - 1] + x1 * heatdown + x2 * heatup + x3
    }
    ssqr <- sum((sst - pred)^2)
    if (out == "x") {
      return(sum(ssqr))
    } else return(pred)
  }
  fitloop <- function(xs1, xs2, xs3, x1add, x2add, x3add, land, sst) {
    dfa <- data.frame(x1 = 0, x2 = 0, x3 = 0, ssqr = 0)
    for (ii in 1:length(xs1)) {
      for (jj in 1:length(xs2)) {
        for (kk in 1:length(xs3)) {
          x1 <- xs1[ii] + x1add
          x2 <- xs2[jj] + x2add
          x3 <- xs3[kk] + x3add
          ssqr <- fitone(land, sst, x1, x2, x3)
          df1 <- data.frame(x1 = x1, x2 = x2, x3 = x3, ssqr = ssqr)
          dfa <- rbind(dfa, df1)
        }
      }
    }
    dfa <- dfa[-1, ]
    sel <- which(dfa$ssqr == min(dfa$ssqr, na.rm = T))
    dfa[sel, ]
  }
  xs1 <- c(0:10) / 100
  xs2 <- c(0:10) / 100
  xs3 <- c(0:10) / 100
  dfo <- fitloop(xs1, xs2, xs3, 0, 0, 0, temps$land, temps$sst)
  if (dfo$x1 > 0) {
    xs1 <- c(-5:5)/1000
  } else xs1 <- c(0:10)/1000
  if (dfo$x2 > 0) {
    xs2 <- c(-5:5)/1000
  } else xs2 <- c(0:10)/1000
  if (dfo$x3 > 0) {
    xs3 <- c(-5:5)/1000
  } else xs3 <- c(0:10)/1000
  dfo <- fitloop(xs1, xs2, xs3, dfo$x1, dfo$x2, dfo$x3, temps$land, temps$sst)
  if (dfo$x1 > 0) {
    xs1 <- c(-5:5)/10000
  } else xs1 <- c(0:10)/10000
  if (dfo$x2 > 0) {
    xs2 <- c(-5:5)/10000
  } else xs2 <- c(0:10)/10000
  if (dfo$x3 > 0) {
    xs3 <- c(-5:5)/10000
  } else xs3 <- c(0:10)/10000
  dfo <- fitloop(xs1, xs2, xs3, dfo$x1, dfo$x2, dfo$x3, temps$land, temps$sst)
  if (dfo$x1 > 0) {
    xs1 <- c(-5:5)/1e+05
  } else xs1 <- c(0:10)/1e+05
  if (dfo$x2 > 0) {
    xs2 <- c(-5:5)/1e+05
  }  else xs2 <- c(0:10)/1e+05
  if (dfo$x3 > 0) {
    xs3 <- c(-5:5)/1e+05
  }  else xs3 <- c(0:10)/1e+05
  dfo <- fitloop(xs1, xs2, xs3, dfo$x1, dfo$x2, dfo$x3, temps$land, temps$sst)
  if (dfo$x1 > 0) {
    xs1 <- c(-5:5)/1e+06
  } else xs1 <- c(0:10)/1e+06
  if (dfo$x2 > 0) {
    xs2 <- c(-5:5)/1e+06
  }  else xs2 <- c(0:10)/1e+06
  if (dfo$x3 > 0) {
    xs3 <- c(-5:5)/1e+06
  }  else xs3 <- c(0:10)/1e+06
  dfo <- fitloop(xs1, xs2, xs3, dfo$x1, dfo$x2, dfo$x3, temps$land, temps$sst)
  dfo$k <- mean(temps$sst) - mean(temps$land)
  if (Trace) {
    pred <- fitone(temps$land, temps$sst, dfo$x1, dfo$x2, dfo$x3, out = "p")
    tme <- as.POSIXct(tme)
    mn <- min(min(pred), min(temps$sst))
    mx <- max(max(pred), max(temps$sst))
    tme <- as.POSIXct(temps$obs_time)
    par(mar=c(7,7,2,2))
    plot(temps$sst~tme, ylim = c(mn, mx), type = "l", col = "red",
         ylab = "Temperature", xlab = "Month", cex.lab = 2, cex.axis = 2,
         main = "Red: observed, Grey: predicted", lwd = 2)
    par(new=T)
    plot(pred~tme, ylim = c(mn, mx), type = "l", col = "darkgray", lwd = 3,
         ylab = "", xlab = "",  cex.lab = 2, cex.axis = 2)

  }
  return(dfo)
}
#' Derive sea-surface temperatures from land temperatures
#'
#' @param sstparams of parameter values as returned by [sst_fit()]
#' @param htemps a vector of hourly land temperatues, normally for one year (deg C)
#'
#' @return a vector of hourly sea-surface temperatures
#' @seealso [sst_fit()]
#'
#' @examples
#' sstparams <- sst_fit(temperature2017)
#' htemp <- rep(temperature2017$land, each = 24)
#' sstp <- sst_predict(sstparams, htemp)
#' ssto <- rep(temperature2017$sst, each = 24)
#' plot(ssto, type = "l", ylim = c(5, 20), col = "red", ylab = "SST")
#' par(new = T)
#' plot(sstp, type = "l", ylim = c(5, 20), col = "darkgray", ylab = "")
sst_predict <- function(sstparams, htemp) {
  land <- matrix(htemp, ncol = 24, byrow = T)
  land <- apply(land, 1, mean)
  pred <- land[1] + sstparams$k * 1.85
  sstdeep <- mean(land) + sstparams$k
  for (i in 2:length(land)) {
    heatdown <- (land[i - 1] - pred[i - 1])
    heatup <- (sstdeep - pred[i - 1])
    pred[i] <- pred[i - 1] + sstparams$x1 * heatdown +
      sstparams$x2 * heatup + sstparams$x3
  }
  pred <- c(pred[1], pred, pred[length(pred)])
  n <- (length(pred) - 1) * 24 + 1
  x <- c(1:length(pred))
  ssth <- spline(x, pred, n = n)$y
  ssth <- ssth[13:(n - 13)]
  ssth
}
#' Generates spatial array of sea-surface temperatures
#'
#' @param sstv vector of hourly sea surface temperatures over one year
#' @param monthlysst an object of class `spatialarray` of spatial anomolies in
#' sea-surface tempertaures in each month with the following components:
#' \describe{
#'   \item{arraydata}{a three-dimensional array of monthlysea-surface temperatures }
#'   \item{times}{`POSIXlt` object of times associated with `arraydata`}
#'   \item{crs}{`crs` object of coordinate reference system associated with `arraydata`}
#'   \item{extent}{`extent` object giving extent covered by `arraydata`}
#'   \item{units}{units of `arraydata`}
#'   \item{description}{character description of `arraydata`}
#' }
#' @return an object of class `spatialarray` of hourly sea-surface temperatures
#' @example
#' library(raster)
#' # Get vector of hourly sea-surface temperatures
#' sstparams <- sst_fit(temperature2017)
#' htemp <- rep(temperature2017$land, each = 24)
#' sstv <- sst_predict(sstparams, htemp)
#' # Get spatial array of hourly sea-surface temperatures using inbuilt dataset
#' ssta <- sst_spatial(sstv, monthlysst)
#' # Plot mean January sea-surface temperatures
#' sel <- which(ssta$times$mon + 1 == 1)
#' sstjan <- apply(ssta$arraydata[,,sel], c(1, 2), mean, na.rm = T)
#' r <- microclima::if_raster(sstjan, dem5k)
#' plot(r)
sst_spatial <- function(sstv, monthlysst) {
  tme.mon <- monthlysst$times
  hiy <- 365 * 24
  yr <- as.numeric(tme.mon$year[1] + 1900)
  if (yr%%4 == 0) hiy <- 366 * 24
  x0 <- as.numeric(tme.mon[length(tme.mon)]) - hiy * 3600
  x0 <- as.POSIXlt(x0, origin = "1970-01-01 00:00", tz = "GMT")
  x1 <- as.numeric(tme.mon[1]) + hiy * 3600
  x1 <- as.POSIXlt(x1, origin = "1970-01-01 00:00", tz = "GMT")
  tme.mon <- c(x0, tme.mon, x1)
  tmeh <- seq(as.numeric(tme.mon[1]), as.numeric(tme.mon[length(tme.mon)]),
              by = 3600)
  tmeh <- as.POSIXlt(tmeh, origin = "1970-01-01", tz = "GMT")
  ain <- monthlysst$arraydata
  aout <- array(NA, dim = c(dim(ain)[1:2], hiy))
  sel <- which(tmeh$year + 1900 == yr)
  for (i in 1:dim(ain)[1]) {
    for (j in 1:dim(ain)[2]) {
      tst <- mean(ain[i,j,], na.rm = T)
      if (is.na(tst) == F) {
        y <- c(ain[i,j,12], ain[i,j,], ain[i,j,1])
        yhr <- spline(as.numeric(tme.mon), y, n = length(tmeh))$y
        aout[i,j,] <- yhr[sel] + sstv
      }
    }
  }
  lst <- monthlysst
  lst$arraydata <- aout
  lst$times <- tmeh[sel]
  return(lst)
}





