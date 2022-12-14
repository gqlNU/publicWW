% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processingData.R
\name{approximate.catchment.using.circle}
\alias{approximate.catchment.using.circle}
\title{Approximating catchments using circles with varying radius to achieve a required number of LSOA centroids to fall within}
\usage{
approximate.catchment.using.circle(
  site.coords,
  pc,
  initial.radius = 1,
  nlsoas.req = 10,
  max.nsearches = 40,
  r.inc = 0.2,
  plot = FALSE
)
}
\description{
Approximating catchments using circles with varying radius to achieve a required number of LSOA centroids to fall within
}
