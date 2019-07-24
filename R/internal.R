#' Create POSIXlt object of 3600 daily times for a  10 year UKCP18 dataset
#' @export
.tmecreate <- function(startyear) {
  yrs <- 0
  mts <- 0
  dys <- 0
  i <- 1
  year <- startyear
  mth <- 12
  for (day in 1:30) {
    yrs[i] <- year
    mts[i] <- mth
    dys[i] <- day
    i <- i + 1
  }
  yrx <- c(1:9) + startyear
  for (year in yrx) {
    for (mth in 1:12) {
      for (day in 1:30) {
        yrs[i] <- year
        mts[i] <- mth
        dys[i] <- day
        i <- i + 1
      }
    }
  }
  year <- startyear + 10
  for (mth in 1:11) {
    for (day in 1:30) {
      yrs[i] <- year
      mts[i] <- mth
      dys[i] <- day
      i <- i + 1
    }
  }
  # Shift March by plus one
  sel <- which(mts == 3)
  dys[sel] <- dys[sel] + 1
  # Shift Feb by minus 1
  sel <- which(mts == 2)
  dys[sel] <- dys[sel] - 1
  # 0 of Feb as 31 Jan
  sel <- which(dys == 0)
  dys[sel] <- 31
  mts[sel] <- 1
  # 29 of Feb as 1 Mar
  sel <- which(mts == 2 & dys == 29)
  dys[sel] <- 1
  mts[sel] <- 3
  dys <- ifelse(dys < 10, paste0("0", dys), paste0("", dys))
  mts <- ifelse(mts < 10, paste0("0", mts), paste0("", mts))
  tme <- paste0(yrs, "-", mts, "-", dys)
  tme <- as.POSIXlt(tme, tz = "GMT")
  tme
}
#' Convert nc file to array
#' @import ncdf4
#' @export
.nctoarray <- function(filein, varid = NA) {
  nc <- nc_open(filein)
  a <- aperm(ncvar_get(nc, varid = varid), c(2,1,3))
  a <- apply(a, c(2,3), rev)
  nc_close(nc)
  a
}
#' Crop array of extent e to extent ecrop
#' @import raster
#' @export
.croparray <- function(a, e, ecrop) {
  r <- raster(a[,,1])
  extent(r) <- e
  xs <- seq(e@xmin + 0.5 * res(r)[1], e@xmax - 0.5 * res(r)[1], by = res(r)[1])
  ys <- seq(e@ymin + 0.5 * res(r)[2], e@ymax - 0.5 * res(r)[2], by = res(r)[2])
  ys <- rev(ys)
  xmins <- abs(ecrop@xmin - xs)
  xmaxs <- abs(ecrop@xmax - xs)
  ymins <- abs(ecrop@ymin - ys)
  ymaxs <- abs(ecrop@ymax - ys)
  xmn <- which(xmins == min(xmins))
  xmx <- which(xmaxs == min(xmaxs))
  ymn <- which(ymins == min(ymins))
  ymx <- which(ymaxs == min(ymaxs))
  aout <- a[ymx[1]:ymn[1], xmn[1]:xmx[1], ]
  eout <-extent(xs[xmn[1]] - 0.5 * res(r)[1],
                xs[xmx[1]] + 0.5 * res(r)[1],
                ys[ymn[1]] - 0.5 * res(r)[2],
                ys[ymx[1]] + 0.5 * res(r)[2])
  return(list(aout = aout, eout = eout))
}
#' resample spatial array object to resolution of r
#' @export
.resamplearray <- function(ae, r, tme, Trace) {
  a <- ae$aout
  ao <- array(NA, dim = c(dim(r)[1:2], dim(a)[3]))
  for (i in 1:dim(a)[3]) {
    ro <- raster(a[,,i])
    extent(ro) <- ae$eout
    ro <- resample(ro, r)
    ro <- mask(ro, r)
    ao[,,i] <- is_raster(ro)
    if(Trace & i%%100 == 0) plot(ro, main = tme[i])
  }
  ao
}
# Calculate latitudes of raster cells
#' @export
.latsfromr <- function(r) {
  e <- extent(r)
  lts <- rep(seq(e@ymax - res(r)[2] / 2, e@ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}
# Calculate longitudes of raster cells
#' @export
.lonsfromr <- function(r) {
  e <- extent(r)
  lns <- rep(seq(e@xmin + res(r)[1] / 2, e@xmax - res(r)[1] / 2, length.out = dim(r)[2]), dim(r)[1])
  lns <- lns[order(lns)]
  lns <- array(lns, dim = dim(r)[1:2])
  lns
}
# Calculate semivariogram
#' @export
.semivar <- function(a, r) {
  xx <- apply(a, 3, mean)
  sel <- which(is.na(xx) == F)
  ns <- floor(seq(1, length(sel), length.out = 20))
  psill <- 0
  rge <- 0
  sdev <- 0
  for (i in 1:20) {
    m <- a[,,sel[ns[i]]]
    sdev[i] <- sd(m)
    rs <- raster(m, template = r)
    xy <- data.frame(xyFromCell(rs, 1:ncell(rs)))
    xy$z <- getValues(rs)
    coordinates(xy) = ~x+y
    v <- variogram(z~1, xy)
    vf <- tryCatch(fit.variogram(v, vgm("Sph")),error=function(e) e,
                   warning=function(w) w)
    if (class(vf)[1] == "variogramModel") {
      psill[i] <- vf[2,2]
      rge[i] <- vf[2,3]
    } else {
      psill[i] <- NA
      rge[i] <- NA
    }
  }
  sel <- which(is.na(psill) == F)
  return(list(psill = mean(psill, na.rm = T), rge = mean(rge, na.rm = T),
              sdev = mean(sdev)))
}
#' creates sequence matching 10 years of UKCP18 data, with missing days ommited
#' @export
.timesel <- function(startyear) {
  tme2 <- .tmecreate(startyear)
  scs <- seq(as.numeric(tme2[1]), as.numeric(tme2[3600]), by = 24 * 3600)
  tme <- as.POSIXlt(scs, origin = "1970-01-01 00:00", tz = "GMT")
  tme2 <- round(as.numeric(tme2) / (3600 * 24), 0)
  tme <- round(as.numeric(tme) / (3600 * 24), 0)
  sel <- 0
  for (i in 1:length(tme)) {
    xx <- which(tme[i] == tme2)
    if (length(xx) > 0) {
      sel[i] <- i
    } else sel[i] <- NA
  }
  sel <- sel[is.na(sel) == F]
  sel
}
#' Apply gam model to spatial array to adjust data
#' @export
.applygam <- function(a, mod_gam) {
  x1 <- a$arraydata
  x1 <- as.vector(x1)
  p1 <- predict.gam(mod_gam, newdata = data.frame(v2 = x2))
  p1 <- array(p1, dim = dim(a$arraydata))
  a$arraydata <- p1
  return(a)
}
