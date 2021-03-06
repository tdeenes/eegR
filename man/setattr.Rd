% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{setattr}
\alias{setattr}
\title{Set attributes of objects by reference}
\usage{
setattr(x, name, value, check_for = NULL)
}
\arguments{
\item{x}{an object}

\item{name}{the name of the attribute (a character string)}

\item{value}{the value to assign to the attribute. If value is NULL, 
\code{setattr} removes the attribute (if present).}

\item{check_for}{an object which might have the same address in RAM as
'x' but should not be affected by \code{setattr}}
}
\description{
\code{setattr} changes the attribute of an object by reference, that is, 
without making any copy. This is an imported and re-exported function from
the \bold{data.table} package with two minor enhancements (see Details).
}
\details{
The original \code{\link[data.table]{setattr}} function has an
undesirable consequence on \code{!FALSE} or \code{!TRUE} if the 'x' object 
is a logical scalar (see References). This function corrects this bug. 
Additionally, the standard \code{attr} function is applied for such an 'x' 
object, which is referenced to the 'check_for' object (has the same address 
in RAM, see \code{\link{address}}).
}
\references{
\url{https://github.com/Rdatatable/data.table/issues/1281}
}
\seealso{
\code{\link[data.table]{setattr}} for the original version in the data.table
package
}
