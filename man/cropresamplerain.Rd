% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/precipitation.R
\name{cropresamplerain}
\alias{cropresamplerain}
\title{Crop rainfall, convert to daily and resample to desired extent and resolution}
\usage{
cropresamplerain(file, r, Trace = TRUE)
}
\arguments{
\item{file}{the name of the .nc file which the data are to be read from.
If it does not contain an absolute path, the file name is relative to the
current working directory, \code{\link[=getwd]{getwd()}}. Tilde-expansion is performed where
supported.}

\item{r}{an object of type raster of elevations. The returned \code{spatialarray} has the same
resolution and extent as \code{r}.}

\item{Trace}{optional logical indicating whether to create and plot
a raster of every 100th entry to enable progress to be tracked.}
}
\value{
An object of class \code{spatialarray} containing the following components:
\describe{
\item{arraydata}{A three-dimensional array of daily rainfall values}
\item{times}{An object of class \code{POSIXlt} of times \code{arraydata}}
\item{crs}{An object of class \code{crs} indicating the coordinate reference system associated with \code{arraydata}}
\item{extent}{An object of class \code{extent} indicating the geographic extent covered by \code{arraydata}}
\item{units}{Units of \code{arraydata}}
\item{description}{Description of \code{arraydata}}
}
}
\description{
Crop rainfall, convert to daily and resample to desired extent and resolution
}
\details{
This function crops and resamples the rainfall data in \code{file} to the extent of
\code{r}. It also performs elevation adjustments by applying a thin-plate spline model with
elevation as a covariate to calculate the expected total rainfall and number of zero
rainfall days. The resampled data is then adjusted accordingly. The 3600 values for each
grid cell location are then spline interpolated to give daily values. Currently, the
UKCP18 nc files contain errors in both the times and the defined coordinate reference
system. These are corrected in the function, but users are advised to check their data.
}
