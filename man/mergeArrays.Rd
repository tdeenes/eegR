% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/array_transform.R
\name{mergeArrays}
\alias{mergeArrays}
\title{Merge arrays having common dimension identifiers}
\usage{
mergeArrays(..., base_value = NA, sort_dims = FALSE, sort_dimlevels = FALSE)
}
\arguments{
\item{...}{numeric arrays with named dimension names or a list of such
arrays. All arrays must have identically named dimensions, but the order of
dimensions does not need to be identical. Duplicated dimension levels are not
allowed.}

\item{base_value}{while setting up the resulting array, what value should be
given as default (e.g. NA, 0, "", etc.)}

\item{sort_dims}{logical. If FALSE (default), the order of dimensions follows
the first array's dimension order; if TRUE, lexical sorting is applied.}

\item{sort_dimlevels}{logical; should dimension levels be sorted for each
dimension (default: FALSE)}
}
\value{
The resulting array has identical dimension identifiers as the input
arrays, and for each dimension, as many dimension levels as the union of
the dimension levels of the input arrays.
}
\description{
\code{mergeArrays} merges multiple arrays or a list of arrays into one large
array. It can be regarded as the inverse of \code{\link{splitArray}}.
}
