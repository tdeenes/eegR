#
# <<< statistics on matrices >>> --------
#

#' Find local/global peaks (maxima) and valleys (minima) of a vector
#' 
#' \code{findExtrema} identifies local (or global) peaks and valleys of a 
#' vector (or subviews of a matrix or array).
#' @param x an integer or numeric vector, matrix, or array
#' @param n the number of neighbouring points (default: 1L) to the left and to 
#' the right of each data point; a data point is a local minimum/maximum if it 
#' is below/above all data points in its neighbourhood. See also Details.
#' @param global logical value whether global extrema should be identified
#' instead of local extrema (default: FALSE). If TRUE, \code{n} is ignored.
#' See Note.
#' @param along_dim the dimension of \code{x} which defines the vectors that
#' should be tested for local extrema. If \code{along_dim} is of type character,
#' \code{x} must have named dimnames.
#' @param tail character string defining how tails should be handled (might be
#' abbreviated): "if_plateau" (the default) means extrema at the tails of the 
#' vector are only valid if they are part of a plateau; "never" means no 
#' extrema at the tails; "do_not_care" means no special treatment of extrema 
#' at the tails
#' @param topN a positive integer scalar; if not NULL (the default), the local 
#' extrema must be among the top N highest/lowest extrema to be marked as such. 
#' If there are ties, all of them are returned.
#' @param has_NA if FALSE, x is not checked for missing values, thereby speeding
#' up the computations; if has_NA is NULL (default), a fast check is performed 
#' and if x has missing values, special corrections are applied (see Details).
#' @details There are three special features of \code{findExtrema}. First, 
#' repeated neighbouring values ('plateaus') are treated as if they were a 
#' single data point. This has two consequences: 1) If the plateau is a local 
#' extrema, all points of the plateau are considered as extrema. 2) The 
#' argument 'n' of \code{findExtrema} is applied on the 'unitized' (de-repeated)
#' vector.\cr
#' Second, extrema at the tails (the first and last elements of the vector) 
#' are often problematic, because for those data points, only one-sided 
#' comparisons are available. By default, \code{findExtrema} considers such
#' endpoints as local minima/maxima if they are part of a plateau (i.e. if 
#' the endpoint is equal to its nearest neighbour). However, this behaviour
#' can be changed by setting the argument 'tail' to "never" or "do_not_care".\cr
#' Third, neighbours of a missing value are returned as missing values because 
#' if a given data window has at least one missing value, no minimum or 
#' maximum can be computed.
#' @note Instead of setting 'global' to TRUE, you can also set 'n' large 
#' enough (e.g. length(x)) to achieve the same effect.
#' @export
#' @return \code{findExtrema} returns an object of the same shape and length 
#' as \code{x}, recoding the original values in \code{x} to integer values 1 
#' (local maximum), -1 (local minimum), and 0 (neither minimum nor maximum), or 
#' NA (not available).
#' @examples
#' # create a vector with two local minima (which are equal) and two 
#' # local maxima (which are different)
#' x <- c(10, 7, -1, 6, 2, -1, 5, 4, 3)
#' 
#' # find local minima/maxima
#' (x_extr <- findExtrema(x))
#' 
#' # the same with a more stringent criterion
#' (x_extr2 <- findExtrema(x, 2L))
#' 
#' # return only the top 1 extrema; note that the local minima are equal,
#' # so both of them are returned, but only the higher local maximum is kept
#' (x_extr3 <- findExtrema(x, topN = 1L))
#' 
#' # findExtrema() always returns an integer vector or matrix
#' stopifnot(is.integer(x_extr))
#' 
#' # check results
#' stopifnot(identical(x[x_extr < 0L], c(-1, -1)))
#' stopifnot(identical(x[x_extr > 0L], c(6, 5)))
#' stopifnot(identical(x[x_extr2 != 0L], c(-1, -1, 5)))
#' stopifnot(identical(x[x_extr3 != 0L], c(-1, 6, -1)))
#' 
#' # look for global extrema
#' (x_global <- findExtrema(x, global = TRUE))
#' stopifnot(identical(x[x_global != 0L], c(-1, -1)))
#' 
#' # the same with large 'n'
#' (x_global2 <- findExtrema(x, length(x)))
#' stopifnot(identical(x_global, x_global2))
#' 
#' # modify the vector to have a plateau at the start, and a missing value at
#' # the position 8; consider only the nearest neighbours
#' x <- c(10, x)
#' x[8] <- NA
#' x
#' 
#' # now the first two elements should also be identified as local maxima,
#' # but the second local minimum is not a local minimum any more because 
#' # there is a missing value in its neighbourhood
#' (x_extr <- findExtrema(x))
#' stopifnot(all(x_extr[1:2] == 1L))
#' stopifnot(all(is.na(x_extr[7:9])))
#' 
#' # visualize the results (blue: local minimum, red: local maximum)
#' plot(x, type = "l", lty = 3)
#' points(x, pch = 16, col = c("blue", "grey", "red")[findExtrema(x) + 2L])
#' 
#' # however, if 'tail' is set to "never", the first two elements are not
#' # extrema
#' (x_extr <- findExtrema(x, tail = "n"))
#' stopifnot(all(x_extr[1:2] == 0L))
#' 
#' # x can be a matrix (or even an array)
#' x <- cbind(sin(seq(0, 3*pi, pi/4)), cos(seq(0, 3*pi, pi/4)))
#' (x_extr <- findExtrema(x))
#' matplot(x, type = "l", lty = 3, col = 1)
#' points(x[, 1L], pch = 16, 
#'        col = c("blue", "grey", "red")[x_extr[, 1L] + 2L])
#' points(x[, 2L], pch = 16, 
#'        col = c("blue", "grey", "red")[x_extr[, 2L] + 2L])
#'        
findExtrema <- function(x, n = 1L, global = FALSE, along_dim = 1L, 
                        tail = c("if_plateau", "never", "do_not_care"), 
                        topN = NULL, has_NA = NULL) {
    #
    # argument checks
    assertAtomic(x, .var.name = "x")
    assertInt(n, lower = 1, .var.name = "n")
    assertFlag(global, .var.name = "global")
    if (global) n <- length(x)
    tail <- match.arg(tail)
    if (!is.null(topN) && !identical(topN, Inf)) {
        assertCount(topN, .var.name = "topN")
        topN <- as.integer(topN)
    }
    if (is.null(has_NA)) {
        has_NA <- anyNA(x) 
    } else {
        assertFlag(has_NA, .var.name = "has_NA")
    }
    n <- as.integer(n)
    # return
    if (is.vector(x)) {
        out <- as.vector(
            findExtremaMatrix(as.matrix(x), n, tail, topN, has_NA))
        setattr(out, "names", names(x))
        out
    } else {
        fnDims(x, along_dim, findExtremaMatrix, 
               arg_list = list(n, tail, topN, has_NA),
               vectorized = TRUE, keep_dimorder = TRUE)
    }
}

#' Find local extrema in column vectors
#' 
#' \code{findExtremaMatrix} looks for local minima and maxima in the
#' columns of a matrix. This is an internal function which supports the
#' exported function \code{\link{findExtrema}}.
#' @param x a numeric matrix
#' @inheritParams findExtrema
#' @seealso \code{\link{findExtrema}}
findExtremaMatrix <- function(x, n, tail, topN, has_NA) {
    # helper function to test near equality
    test_equal <- function(dat, ref) {
        if (!identical(dim(dat), dim(ref))) {
            abs(sweepMatrix(dat, 2, ref, "-")) < .Machine$double.eps^0.5
        } else {
            abs(dat - ref) < .Machine$double.eps^0.5
        }
    }
    # find local extrema
    if (n >= (nrow(x) - 1L)) {
        # if n >= (nrow(x) - 1), rolling statistics are not needed
        out <- matrix_(0L, nrow(x), ncol(x))
        stats <- colRanges(x)
        out[test_equal(x, stats[, 1L])] <- -1L
        out[test_equal(x, stats[, 2L])] <- 1L
        # this is needed later to handle external elements
        if (tail != "do_not_care") {
            rle_x <- matrixRle(x)
            ind <- cumsum(c(1L, rle_x$lengths[-length(rle_x$length)]))
            rle_x$values <- out[ind]
            # set out to NULL to signal that we want to use rle_x
            out <- NULL; ind <- NULL
        }
    } else {
        out <- NULL
        # main computation (need run-length encoding to handle repeated values
        # in x)
        rle_x <- matrixRle(x)
        # replicate external elements if x is a multi-column matrix
        if (ncol(x) > 1L) {
            ins <- diff(rle_x$matrixcolumn)
            ins <- 1L + n * c(0L, ins) + c(ins, 0L)
            rle_x_ins <- lapply(rle_x[2:3], rep.int, ins)
            ins_ind <- which(diff(rle_x_ins$matrixcolumn) > 0L)
            ins_ind <- as.vector( outer((1-n):n, ins_ind, "+") )
        } else {
            rle_x_ins <- rle_x
            ins_ind <- NULL
        }
        # store values and remove temporary variables
        values <- rle_x_ins$values
        rle_x_ins <- NULL; ins <- NULL
        # find extrema
        temp_values <- integer(length(values))
        temp_values[test_equal(values, rollFun(values, 2 * n + 1L, max))] <- 1L
        temp_values[test_equal(values, rollFun(values, 2 * n + 1L, min))] <- -1L
        # drop inserted elements
        if (!is.null(ins_ind)) temp_values <- temp_values[-ins_ind]
        # assign the new values
        rle_x$values <- temp_values
        temp_values <- NULL; ins_ind < NULL
    }
    # handle external elements
    if (tail != "do_not_care") {
        border_ind <- c(1L, length(rle_x$values))
        if (ncol(x) > 1L) {
            temp_ind <- which(diff(rle_x$matrixcolumn) > 0L)
            border_ind <- c(border_ind,
                            outer(0:1, temp_ind, "+"))
        }
        if (tail == "if_plateau") {
            border_ind <- border_ind[rle_x$lengths[border_ind] == 1L]
        }
        rle_x$values[border_ind] <- 0L
    }
    # reshape to a matrix (unless n >= nrow(x), and tail == "do_not_care")
    if (is.null(out)) out <- inverse.matrixRle(rle_x)
    #
    # handle the neighbourhood of NA values
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
    # return only the topN extrema
    if (is.integer(topN)) {
        temp_x <- x
        temp_x[out != -1L] <- NA
        out[which(colRanks(temp_x, 
                           ties.method = "min",
                           preserveShape = TRUE) > topN)] <- 0L
        temp_x <- -x
        temp_x[out != 1L] <- NA
        out[which(colRanks(temp_x, 
                           ties.method = "min",
                           preserveShape = TRUE) > topN)] <- 0L
    }
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
#' @examples
#' # either caTools or zoo must be installed before using this function;
#' # here follows a timing comparison for caTools and zoo, so we need both
#' if (require(caTools) && require(zoo)) {
#'     # create a matrix
#'     x <- matrix_(rnorm(2e4), 1e2, 2e2)
#' 
#'     # compute rolling mean for each columns, set the width of the 
#'     # sliding window to 5
#'     system.time(roll_mean_catools <- rollFun(x, 5, mean))
#'     system.time(roll_mean_zoo <- rollFun(x, 5, mean, force_rollapply = TRUE))
#'     
#'     # caTools is much faster for the standard statistics, and the results
#'     # are the same
#'     stopifnot(all.equal(roll_mean_catools, roll_mean_zoo))
#' }
#'
rollFun <- function(dat, width, FUN, force_rollapply = FALSE, ...) {
    # check arguments and store attributes
    dims <- dim(dat)
    attribs <- attributes(dat)
    if (is.data.frame(dat)) dat <- as.matrix(dat)
    if (is.list(dat)) 
        stop("Provide a vector, matrix or data.frame as input!")
    if (is.vector(dat)) {
        dat <- matrix(dat, ncol = 1)
        dims <- dim(dat)
    }
    if (length(dims) > 2)
        stop("dat can not be a multidimensional array. Consider the combination of fnDims and rollFun with vectorized = TRUE")
    assertIntegerish(width, lower = 1, upper = dims[1], min.len = 1L,
                     any.missing = FALSE, .var.name = "width")
    width <- as.integer(width)
    if (length(width) > 1L) width <- rep_len(width, dims[1])
    FUN <- match.fun(FUN)
    assertFlag(force_rollapply, .var.name = "force_rollaply")
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
                if (length(width) == 1L) {
                    out <- FUN(dat, width, ...)
                } else {
                    out <- array_(0, dims)
                    if (identical(FUN, caTools::runquantile)) {
                        for (i in unique(width)) {
                            ind <- width == i
                            out[ind, , ] <- FUN(dat, i, ...)[ind, , , drop = FALSE]
                        }
                    } else {
                        for (i in unique(width)) {
                            ind <- width == i
                            out[ind, ] <- FUN(dat, i, ...)[ind, , drop = FALSE]
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
#' @param has_neg logical value if there are (or can be) negative values in x 
#' (default: TRUE)
#' @param has_pos logical value if there are (or can be) positive values in x
#' (default: TRUE)
#' @export
#' @keywords internal
#' @return numeric matrix or array of the same dimensions as x
#' @author The original C code was written by Christian Gaser, further modified
#' by Pau Coma and adopted to EEG analysis by Armand Mensen. Denes Toth ported
#' the C code to C++ (Rcpp) with minor improvements and added the R wrapper.
tfceFn <- function(x, chn, eh, nr_steps = 50L, channel_dim = 1L, 
                   has_neg = TRUE, has_pos = TRUE) {
    chan_dim <- which(names(dimnames(x)) == "chan")
    if (length(chan_dim) == 0) {
        stopifnot(dim(x)[channel_dim] == nrow(chn))
        chan_dim <- channel_dim
    }
    # return
    .Call('eegR_tfce', PACKAGE = 'eegR', 
          x, chan_dim, chn, eh, nr_steps, has_neg, has_pos)
}
