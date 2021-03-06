% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infinitesimalJackknife.R
\name{randomForestInfJack}
\alias{randomForestInfJack}
\title{The infinitesimal jackknife for random forests}
\usage{
randomForestInfJack(rf, newdata, calibrate = TRUE)
}
\arguments{
\item{rf}{A random forest trained with replace = TRUE and keep.inbag = TRUE}

\item{newdata}{A set of test points at which to evaluate standard errors}

\item{calibrate}{whether to apply calibration to mitigate Monte Carlo noise
warning: if calibrate = FALSE, some variance estimates may be negative
         due to Monte Carlo effects if the number of trees in rf is too small}
}
\description{
The infinitesimal jackknife for random forests
}
