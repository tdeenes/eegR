#
# <<< statistics on matrices >>> --------
#

#' Find local peaks (maxima) and valleys (minima) of a vector
#' 
#' \code{findExtremes} identifies local peaks and valleys of a vector or 
#' the columns of a matrix. Endpoints are considered peaks/valleys if they are 
#' part of a plateau.
#' @param x an integer or numeric vector or matrix (or any object which can
#' be coerced to a matrix)
#' @param n the number of neighbouring points (default: 1L) to the left and to 
#' the right of each data point; a data point is a local minimum/maximum if it 
#' is below/above all data points in its neighbourhood
#' @param has_NA if FALSE, x is not checked for missing values, thereby speeding
#' up the computations; if has_NA is NULL (default), a fast check is performed 
#' and if x has missing values, special corrections are applied (see Details).
#' @details There are two special aspects of \code{findExtremes}. First, 
#' endpoints (the first and last elements of the vector or column) are 
#' considered as local minima/maxima if they are part of a plateau (i.e. if 
#' the endpoint is equal to its nearest neighbour). Second, neighbours of a 
#' missing value are returned as missing values because if a given data window
#' has at least one missing value, no minimum or maximum can be computed.
#' @export
#' @return An integer vector or matrix of the same length as x, with values 1 
#' (local maximum), -1 (local minimum) and 0 (neither minimum nor maximum), or 
#' NA (not available). 
#' @examples
#' # create a vector with two local minima and one local maximum
#' x <- c(10, 3, 1, 2, -1, 0, 4, 5)
#' 
#' # find local minima/maxima
#' (x_extr <- findExtremes(x))
#' 
#' # the same with a more stringent criterion
#' (x_extr2 <- findExtremes(x, 2L))
#' 
#' # note that findExtremes always returns an integer vector or matrix
#' stopifnot(is.integer(x_extr))
#' 
#' # check results
#' stopifnot(identical(x[x_extr < 0L], c(1, -1)))
#' stopifnot(identical(x[x_extr > 0L], 2))
#' stopifnot(identical(x[x_extr2 != 0L], -1))
#' 
#' # modify the vector to have a plateau at the start, and a missing value at
#' # the position 7; consider only the nearest neighbours
#' x <- c(10, x)
#' x[7] <- NA
#' x
#' 
#' # now the first two elements should also be identified as local maxima,
#' # but the value -1 is not a local minimum any more because there is a 
#' # missing value in its neighbourhood
#' (x_extr <- findExtremes(x))
#' stopifnot(all(x_extr[1:2] == 1L))
#' stopifnot(all(is.na(x_extr[6:8])))
#' 
#' # visualize the results (blue = local min., red = local max.)
#' plot(x, type = "l", lty = 3)
#' points(x, pch = 16, col = c("blue", "grey", "red")[findExtremes(x) + 2L])
#' 
#' # x can be a matrix
#' x <- cbind(sin(seq(0, 3*pi, pi/4)), cos(seq(0, 3*pi, pi/4)))
#' (x_extr <- findExtremes(x))
#' matplot(x, type = "l", lty = 3, col = 1)
#' points(x[, 1L], pch = 16, 
#'        col = c("blue", "grey", "red")[findExtremes(x[, 1L]) + 2L])
#' points(x[, 2L], pch = 16, 
#'        col = c("blue", "grey", "red")[findExtremes(x[, 2L]) + 2L])
findExtremes <- function(x, n = 1L, has_NA = NULL) {
    # argument checks
    if (!is.atomic(x) | length(dim(x)) > 2) 
        stop("Provide a vector or matrix as input")
    if (is.null(has_NA)) has_NA <- anyNA(x)
    n <- as.integer(n)
    # store original dimensions before transforming to matrix
    origdim <- dim(x)
    origdimnames <- dimnames(x)
    orignames <- names(x)
    x <- as.matrix(x)
    # main computations
    out <- matrixIP(0L, nrow(x), ncol(x))
    out[x == rollFun(x, 2*n + 1L, max)] <- 1L
    out[x == rollFun(x, 2*n + 1L, min)] <- -1L
    # external elements
    out[1L, which(x[1L, ] != x[2L, ])] <- 0L
    ind <- which(abs(x[1L, ] - x[2L, ]) < .Machine$double.eps ^ 0.5)
    out[1L, ind] <- out[2L, ind]
    nr <- nrow(out)
    out[nr, which(out[nr, ] != out[nr - 1L, ])] <- 0L
    ind <- which(abs(x[nr, ] - x[nr - 1L, ]) < .Machine$double.eps ^ 0.5)
    out[nr, ind] <- out[nr - 1L, ind]
    # neighbourhood of NA values
    if (has_NA) {
        naind <- which(is.na(x), arr.ind = TRUE)
        out[naind] <- NA
        for (i in 1:n) {
            ind <- cbind(pmax(1L, naind[, 1L] - i), naind[, 2L])
            out[ind] <- NA
            ind <- cbind(pmin(nrow(x), naind[, 1L] + i), naind[, 2L])
            out[ind] <- NA
        }
    }
    # original attributes
    setattr(out, "dim", origdim)
    setattr(out, "dimnames", origdimnames)
    setattr(out, "names", orignames)
    # return
    out
}

#' Compute rolling (a.k.a. moving) window statistics
#' 
#' \code{rollFun} computes rolling window statistics on vectors or matrices.
#' @param dat a numeric vector, matrix or data.frame. In the latter cases
#' rolling statistics are computed column-wise.
#' @param width width of moving window; can be an integer value or vector.
#' @param FUN the function to be applied to compute moving window statistics. 
#' See details.
#' @param force_rollapply logical variable; if yes, \code{zoo::rollapply} is 
#' called (default = FALSE).
#' @param ... optional arguments to the corresponding function in \pkg{caTools}
#' or \code{zoo::rollapply}
#' @details If FUN is one of \code{min}, \code{max}, \code{mean}, \code{sd}, 
#' \code{mad}, \code{quantile} (OR "min", "max", "mean", etc.) \code{rollFun} 
#' calls the corresponding function from the \pkg{caTools} package (e.g. 
#' \code{caTools::runmin}). Otherwise, or if \code{force_rollapply} is TRUE,
#' \code{zoo::rollapply} is called.
#' @export
#' @return An object having the same attributes as dat.
rollFun <- function(dat, width, FUN, force_rollapply = FALSE, ...) {
    dims <- dim(dat)
    attribs <- attributes(dat)
    if (is.data.frame(dat)) dat <- as.matrix(dat)
    if (is.list(dat)) 
        stop("Provide a vector, a matrix or a data.frame as input!")
    if (is.vector(dat)) {
        dat <- matrix(dat, ncol = 1)
        dims <- dim(dat)
    }
    if ((min(width) < 1) || (max(width) > nrow(dat)))
        stop("width must be an integer value or vector between 1 and 
             length(dat) or nrow(dat)")
    width <- as.integer(width)
    FUN <- match.fun(FUN)
    if (!force_rollapply) {
        funlist <- list(min = min, max = max, mean = mean, 
                        sd = sd, mad = mad, quantile = quantile)
        funind <- sapply(funlist, identical, FUN)
        if (any(funind)) {
            if (requireNamespace("caTools", quietly = TRUE)) {
                FUN <- names(funlist)[funind]
                FUN <- getExportedValue("caTools", paste0("run", FUN))
                if (identical(FUN, caTools::runquantile)) {
                    probs <- list(...)$probs
                    dims <- c(dim(dat), length(probs))
                    attribs$dim <- dims
                    if (!is.null(attribs$dimnames)) {
                        attribs$dimnames <- c(attribs$dimnames, 
                                              list(quantiles = probs))
                    }
                }
                if (length(width) == 1) {
                    out <- FUN(dat, width, ...)
                } else {
                    out <- arrayIP(0, dims)
                    if (identical(FUN, caTools::runquantile)) {
                        for (i in unique(width)) {
                            ind <- width == i
                            out[ind, , ] <- FUN(dat, i, ...)[ind, , , drop = F]
                        }
                    } else {
                        for (i in unique(width)) {
                            ind <- width == i
                            out[ind, ] <- FUN(dat, i, ...)[ind, , drop = F]
                        }
                    }
                }
            } else {
                force_rollapply <- TRUE
                warning("Install package:caTools to speed up computations!")
            }
        } else {
            force_rollapply <- TRUE
        }
    } 
    if (force_rollapply) {
        if (reqFn("zoo")) {
            out <- zoo::rollapply(dat, width, FUN, partial = TRUE, 
                                  by.column = TRUE, ...)
        } else {
            stop("Package zoo is not installed but called by rollFun().")
        }
    }
    if (!is.null(attribs$class) && attribs$class == "data.frame") {
        out <- as.data.frame(out)
    }
    attributes(out) <- attribs
    # return
    out
}

#' Low-level TFCE function which calls C++ code
#' 
#' \code{tfceFn} performs TFCE correction. This function is not 
#' intended for direct use.
#' @param x numeric matrix or array
#' @param chn channel neighbourhood matrix
#' @param eh numeric vector of E and H parameters
#' @param nr_steps number of threshold steps (default: 50L)
#' @param channel_dim the dimension of x which represents channels (default: 1L)
#' @export
#' @keywords internal
#' @return numeric matrix or array of the same dimensions as x
#' @author The original C code was written by Christian Gaser, further modified
#' by Pau Coma and adopted to EEG analysis by Armand Mensen. Denes Toth ported
#' the C code to C++ (Rcpp) with minor improvements and added the R wrapper.
tfceFn <- function(x, chn, eh, nr_steps = 50L, channel_dim = 1L) {
    chan_dim = which(names(dimnames(x)) == "chan")
    if (length(chan_dim)==0) {
        stopifnot(dim(x)[channel_dim] == nrow(chn))
        chan_dim <- 1L
    }
    out <- arrayIP(0, dim(x))
    ind.neg <- x < 0
    ind.pos <- x > 0
    if (any(ind.pos)) {
        sdat <- x
        sdat[ind.neg] <- 0
        out <- tfce(inData = sdat, chan_dim = chan_dim, 
                    ChN = chn, EH = eh, numSteps = nr_steps)
    }
    if (any(ind.neg)) {
        sdat <- -x
        sdat[ind.pos] <- 0
        out <- out - tfce(inData = sdat, chan_dim = chan_dim, 
                          ChN = chn, EH = eh, numSteps = nr_steps)
    }
    # return
    out
}