% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova.R
\name{arrayTtestAnova}
\alias{arrayTtestAnova}
\title{United function for \code{\link{arrayTtest}} and \code{\link{arrayAnova}}}
\usage{
arrayTtestAnova(
  test,
  .arraydat,
  .arraydat2 = NULL,
  factordef = NULL,
  bwdat = NULL,
  paired = FALSE,
  groups = NULL,
  mu = 0,
  var_equal = FALSE,
  id_dim = "id",
  verbose = TRUE,
  perm = permParams(n = 0L),
  tfce = NULL,
  parallel = NULL,
  seed = NULL
)
}
\description{
United function for \code{\link{arrayTtest}} and \code{\link{arrayAnova}}
}
\keyword{internal}
