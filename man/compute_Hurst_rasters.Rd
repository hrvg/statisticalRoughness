% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_Hurst_rasters.R
\name{compute_Hurst_rasters}
\alias{compute_Hurst_rasters}
\title{Wrapper function around the computation of Hurst coefficients from box-counting and log-log regression}
\usage{
compute_Hurst_rasters(dem, L, R)
}
\arguments{
\item{dem}{initial DEM}

\item{L}{upper scale region size}

\item{R}{list of scales (aggregation factors)}
}
\value{
a RasterStack with four layers q, H, p-value and R squared values
}
\description{
Wrapper function around the computation of Hurst coefficients from box-counting and log-log regression
}
\keyword{Hurst}
