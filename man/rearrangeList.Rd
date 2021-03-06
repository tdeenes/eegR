% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/array_transform.R
\name{rearrangeList}
\alias{rearrangeList}
\title{Rearrange two-level list}
\usage{
rearrangeList(dat, name_listdim, name_datadim = NULL)
}
\arguments{
\item{dat}{the list to be rearranged}

\item{name_listdim}{character string; the name of the dimension which the
list represents}

\item{name_datadim}{character vector or a list of character vectors (for
two-level list input), providing the name of the dimensions for each list
element. The default is NULL, meaning that the original dimension
names will be used.}
}
\value{
A matrix or array, if the input is a one-level list, and a one-level
list, if the input is a two-level list
}
\description{
\code{rearrangeList} reshapes a special type of one- or two-level lists.
}
\details{
One or two level lists which contain identically shaped elements at
the base level, and elements at the base level are vectors, matrices or
arrays, can be rearranged to a matrix/array (from a one-level list) or to a
one-level list (from a two-level list) by binding identical elements.
This way the second level of the list will be represented as an additional
dimension in the ground level elements (vectors will be matrices, matrices
will become arrays, arrays recieve an extra dimension)
}
