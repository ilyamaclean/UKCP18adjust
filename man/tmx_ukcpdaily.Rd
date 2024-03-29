% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{tmx_ukcpdaily}
\alias{tmx_ukcpdaily}
\title{Daily maximum temperature (UKCP18).}
\format{A list with the following components:
\describe{
\item{arraydata}{a three-dimensional array of daily maximum temperature values}
\item{times}{\code{POSIXlt} object of times associated with \code{arraydata}}
\item{crs}{\code{crs} object of coordinate reference system associated with \code{arraydata}}
\item{extent}{\code{extent} object giving extent covered by \code{arraydata}}
\item{units}{units of \code{arraydata}}
\item{description}{character description of \code{arraydata}}
}}
\source{
\url{http://data.ceda.ac.uk/badc/ukcp18/}
}
\usage{
tmx_ukcpdaily
}
\description{
An object of class spatialarray containing daily maximum temperature values in Cornwall for
the period 2000-12-01 to 2010-11-30
}
\keyword{datasets}
