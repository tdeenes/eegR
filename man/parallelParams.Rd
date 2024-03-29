% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param_setters.R
\name{parallelParams}
\alias{parallelParams}
\title{Setting the parameters for parallel computation}
\usage{
parallelParams(
  cl = NULL,
  method = c("auto", "snow", "multicore"),
  ncores = parallel::detectCores() - 1L,
  ...
)
}
\arguments{
\item{cl}{an object of class "cluster" (default: NULL)}

\item{method}{character string of the chosen parallelization method; if 
'auto' (default), 'snow' is chosen on Windows and 'multicore' otherwise}

\item{ncores}{integer; the number of cores}

\item{...}{options to be passed to the function spawning the workers. See 
'Details' in \code{\link[parallel]{makeCluster}} and 
\code{\link[parallel]{mclapply}} for snow- and multicore-parallelization,
respectively (for multicore, one can use the "preschedule", "set.seed", 
"silent" options, e.g. \code{options = list(preschedule = FALSE)}.}
}
\value{
\code{parallelParams} returns a list with five elements:
\itemize{
\item{cl: }{the cluster if a snow-type cluster was requested}
\item{cl_new: }{a logical whether the cluster was created by 
\code{parallelParams}}
\item{snow_options: }{a list of options for snow-clusters}
\item{mc_options: }{a list of options for multicore-clusters}
}
}
\description{
\code{parallelParams} sets the parameters for parallel computation and 
registers the cluster. Newly registered clusters should be stopped after
the computations are performed.
}
\note{
If \code{cl_new = TRUE}, you should stop the cluster after the 
computations (see Examples).
}
\examples{
# create a function  which computes the range of values in each column of
# a matrix and can run in single-core or parallel mode
rangeColumns <- function(x, parallel = FALSE) {
    # 
    stopifnot(require(doParallel))
    ob <- getDoBackend() # store active backend
    #
    # parallel argument can be a logical or a direct call to parallelParams
    # or a .(key = value)-type call
    parallel <- argumentDeparser(substitute(parallel), "parallelParams", 
                                 null_params = list(ncores = 0L))
    #
    # stop cluster on exit if it was created by parallelParams
    if (parallel$cl_new) {
        on.exit(stopCluster(parallel$cl))
    }
    on.exit(setDoBackend(ob), add = TRUE)
    #
    # call foreach and compute range
    out <- foreach(xi = iter(x, by = "col"), .combine = "cbind",
                   .options.snow = parallel$snow_options,
                   .options.multicore = parallel$mc_options) \%dopar\% 
                   range(xi)
    #
    # return with sensible dimension names
    dimnames(out) <- list(range = c("min", "max"), colnames(x))
    out
}
    
# create a toy data matrix
mat <- matrix(rnorm(100), 25, 4)
colnames(mat) <- paste0("column", 1:4)

# compute the range of values in each column and print to the console
ranges_parallel <- rangeColumns(mat, 
                                parallel = .(method = "snow", ncores = 2L))
ranges_parallel

# compare to single-core calculation
ranges_single <- apply(mat, 2, range)
stopifnot(identical(unname(ranges_parallel), 
                    unname(ranges_single)))

}
