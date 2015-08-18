
#' Compute time-warping weights
#' 
#' \code{warp} computes time-warping weights which can be used for 
#' time-alignment of ERP curves (see \code{\link{align}}). It is a wrapper 
#' around \code{\link[fdasrvf]{time_warping}}. See its documentation for 
#' references to the time-warping algorithm.
#' @param erp a numeric matrix or array of ERPs
#' @param align_dim numeric or character index of the dimension which groups
#' the individual time series. If missing, and 'erp' is a matrix, it is the
#' dimension which is not 'time_dim'. See Details.
#' @param time_dim numeric or character index of the time dimension (default: 
#' "time")
#' @param lambda numeric scalar which controls the elasticity (default: 0; 
#' higher lambda means less elastic curves)
#' @param method warp and calculate to Karcher Mean ("mean", the default) or 
#' Karcher Median ("median"). The latter option seems buggy; do not use it
#' in serious analyses.
#' @param smooth_data logical value whether to smooth data using box filter 
#' (default: FALSE)
#' @param smooth_times numeric scalar; number of times to apply box filter if 
#' 'smooth_data' is TRUE (default: 25L)
#' @param verbose logical value whether iteration messages should be displayed
#' and plots should be produced
#' @param time_points numeric vector of time points at which the ERPs were 
#' sampled. If NULL (default), it is derived from the dimnames attribute of 
#' time dimension of the ERP array.
#' @param ... further arguments passed to \code{\link{fnDims}} (e.g. parallel)
#' @note The defaults of \code{warp} are not necessarily identical to the 
#' defaults of the underlying \code{\link[fdasrvf]{time_warping}} function and
#' are subject to change in the future.
#' @details In the standard case, \code{warp} takes a matrix as input which 
#' represent multiple time series (where the individual time series are either 
#' the columns or the rows of the matrix). However, \code{warp} also accepts 
#' multidimensional array as input which is treated as a collection of 
#' independent time X align_dim matrices (the order of the dimensions is not
#' relevant).\cr
#' \code{warp} returns an array of the same size as \code{erp} consisting
#' of the warping weights. The weights can be used for time-alignment of the 
#' same or an other ERP array in \code{\link{align}}.
#' @export
#' @seealso \code{\link{align}} for time-aligning transformation and examples, 
#' and \code{\link{fnDims}} for additional parameters
#' 
warp <- function(erp, align_dim, time_dim = "time", lambda = 0, 
                 method = c("mean", "median"), smooth_data = FALSE, 
                 smooth_times = 25L, verbose = FALSE, time_points = NULL, ...) {
    # helper function
    checkdim <- function(d, dim., dimid., arg_name) {
        assertScalar(d, .var.name = arg_name)
        if (is.character(d) && (!d %in% dimid.)) {
            stop(sprintf("%s: 'erp' does not have '%s' dimension identifier",
                         arg_name, d))
        } else if (is.numeric(d) && !d %in% seq_along(dim.)) {
            stop(sprintf("%s: wrong dimension identifier; must be between 1 and %d",
                         arg_name, length(dim.)))
        }
        if (is.character(d))
            d <- match(d, dimid.)
        # return
        d
    }
    # argument checks
    assertArray(erp, mode = "numeric", min.d = 2L, 
                .var.name = "erp")
    assertNumber(lambda, .var.name = "lambda")
    method <- match.arg(method)
    assertLogical(smooth_data, any.missing = FALSE, len = 1L,
                  .var.name = "smooth_data")
    assertIntegerish(smooth_times, any.missing = FALSE, len = 1L, 
                     .var.name = "smooth_times")
    assertLogical(verbose, any.missing = FALSE, len = 1L,
                  .var.name = "verbose")
    dims <- dim(erp)
    dimn <- dimnames(erp)
    dimid <- names(dimnames(erp))
    names(dims) <- dimid
    time_dim <- checkdim(time_dim, dims, dimid, "time_dim")
    if (missing(align_dim)) {
        if (length(dims) == 2L) {
            align_dim <- setdiff(seq_along(dims), time_dim)
        } else {
            stop(paste0("the argument 'align_dim' must be provided ",
                        "if 'erps' is a multidimensional array"))
        }
    }
    align_dim <- checkdim(align_dim, dims, dimid, "align_dim")
    if (is.null(time_points)) {
        time_points <- as.integer(dimn[[time_dim]])
        if (is.null(time_points) || anyNA(time_points)) {
            stop(paste0("The time dimension does not have proper dimnames; ",
                        "provide time_points explicitly"))
        }
    }
    tw_params <- list(
        time = time_points, lambda = lambda, method = method, 
        showplot = verbose, smooth_data = smooth_data, sparam = smooth_times,
        parallel = FALSE
    )
    # sink time_warping output if not verbose
    if (!verbose) {
        tmp <- tempfile()
        sink(tmp)
        on.exit({sink(file = NULL);unlink(tmp)})
    }
    # main part
    ncols <- dims[align_dim]
    out <- fnDims(erp, c(time_dim, align_dim),
                  function(x, ncols, params) {
                      do("time_warping", matrix_(x, ncol = ncols),
                         arg_list = params)$gam
                  },
                  arg_list = list(ncols = ncols, params = tw_params),
                  ...)
    # return
    setattr(out, "time_points", time_points)
    out
}

#' Adjust ERP curves by time warping
#' 
#' \code{align} transforms the ERP curves by aligning them on the time 
#' dimension.
#' @param erp the numeric vector or array of ERPs
#' @param w the numeric vector or array of the warping weights (see 
#' \code{\link{warp}})
#' @param time_dim either a numeric or a character index indicating the time
#' dimension of the ERP array (default: "time")
#' @details This function is useful if the warping weights are computed on a 
#' specific subset of the ERP data, but the whole array should be aligned using 
#' these pre-computed weights. So the standard analysis pipeline is to call 
#' \code{\link{warp}} on the subset of the ERPs, store the computed weights in 
#' an array, and call \code{align} on the original (non-subsetted) array of ERPs
#' and the stored warping weights.\cr
#' In the simplest case, \code{erp} and \code{w} are vectors. They must be of 
#' equal length. In the second case, \code{erp} and \code{w} are arrays of 
#' equal dimensions. In the third case, \code{erp} is an array with more
#' dimensions than \code{w}. In this scenario \code{erp} and \code{w} must 
#' have named dimensions, and all dimension identifiers in \code{w} must be
#' present in \code{erp}.
#' @export
#' @seealso \code{\link{warp}}
#' @examples
#' # load example data
#' data(erps)
#' 
#' # collapse stimclass and pairtype dimensions + downsample
#' x <- avgDims(erps, c("stimclass", "pairtype"))
#' x <- avgBin(x, "time", bin_length = 3L)
#' str(x) 
#' 
#' # let's assume we want to align the individual ERPs so that the individual
#' # variation in the peak latencies can be controlled for;
#' # we will use the Global Field Power instead of warping the ERPs for each
#' # channel separately
#' w <- warp(compGfp(x), align_dim = "id", verbose = TRUE) # takes some time
#' 
#' # now perform the alignment on the original data
#' x_aligned <- align(x, w)
#' 
#' # bind the original and the aligned arrays
#' x_all <- bindArrays(raw = x, aligned = x_aligned, along_name = "warp")
#' 
#' # plot the GFP of the raw and aligned curves
#' plotERParray(compGfp(x_all), sepdim = "id", 
#'              minus_up = FALSE, grid = c(2, 1))
#' 
#' # repeat the above analysis with lambda = 1 to demonstrate how 
#' # the lambda parameter affects the results
#' w <- warp(compGfp(x), align_dim = "id", lambda = 1) # takes some time
#' 
#' # now perform the alignment on the original data
#' x_aligned <- align(x, w)
#' 
#' # bind the new results to the previous array
#' x_all <- bindArrays(x_all, aligned_lambda1 = x_aligned, along_name = "warp")
#' 
#' # plot the GFP of the raw, the default aligned, and the new aligned curves
#' plotERParray(compGfp(x_all), sepdim = "id", 
#'              minus_up = FALSE, grid = c(3, 1))
#'              
align <- function(erp, w, time_dim = "time") {
    #
    # helper function
    checkdim <- function(d, dim., dimid., arg_name) {
        assertScalar(d, .var.name = arg_name)
        if (is.character(d) && (!d %in% dimid.)) {
            stop(sprintf("%s: 'erp' does not have '%s' dimension identifier",
                         arg_name, d))
        } else if (is.numeric(d) && !d %in% seq_along(dim.)) {
            stop(sprintf("%s: wrong dimension identifier; must be between 1 and %d",
                         arg_name, length(dim.)))
        }
        if (is.character(d))
            d <- match(d, dimid.)
        # return
        d
    }
    # function for interpolation
    warpVector <- function(erp_vec, w_vec, time) {
        approx(
            time, erp_vec, 
            xout = (time[length(time)] - time[1]) * w_vec + time[1])$y
    }
    # function to prepare interpolation
    warpArray <- function(y, w, time, return_array = TRUE) {
        dims <- dim(w)
        dim(y) <- c(dims[1L], prod(dims[-1L]))
        dim(w) <- c(dims[1L], prod(dims[-1L]))
        for (i in 1:ncol(y)) y[, i] <- warpVector(y[, i], w[, i], time)
        # return
        if (return_array) {
            array_(y, dims, arg_check = FALSE)
        } else {
            as.vector(y)
        }
    }
    #
    # check arguments (if vectors, return early)
    time_points <- attr(w, "time_points")
    if (is.null(time_points))
        stop("w must have 'time_points' attribute")
    #
    assertNumeric(erp, .var.name = "erp")
    assertNumeric(w, .var.name = "w")
    #
    if (is.vector(erp) && is.vector(w)) {
        if (length(erp) != length(w)) {
            stop("If erp and w are vectors, they must have equal length")
        }
        return(warpVector(erp, w, time_points))
    }
    #
    dims_erp <- dim(erp); dims_w <- dim(w);
    dimn_erp <- dimnames(erp); dimn_w <- dimnames(w);
    dimid_erp <- names(dimn_erp); dimid_w <- names(dimn_w);
    names(dims_erp) <- dimid_erp; names(dims_w) <- dimid_w;
    #
    # if identical dimensions, return early
    if (identical(dim(erp), dim(w)) && identical(length(erp), length(w))) {
        time_dim <- checkdim(time_dim, dims_erp, dimid_erp, "time_dim")
        reord <- order(c(time_dim, seq_along(dims_erp)[-time_dim]))
        erp <- apermArray(erp, first = time_dim)
        w <- apermArray(w, first = time_dim)
        out <- warpArray(erp, w, time_points, TRUE)
        out <- aperm(out, reord)
        return(array_(out, dims_erp, dimn_erp))
    }
    #
    # if erp has more dimensions than w, they must have named dimnames
    if (is.null(dimid_erp) || is.null(dimid_w)) {
        stop("If not vectors, both 'erp' and 'w' must have named dimension names")
    }
    if (!all(dimid_w %in% dimid_erp)) {
        stop("w may not have any dimension identifier which is not present in erp")
    }
    if (!time_dim %in% dimid_w) {
        stop("If not vectors, both erp and w must have time dimension")
    }
    if (!identical(dimnames(w), dimn_erp[dimid_w])) {
        stop("Common dimensions in erp and w must have identical dimension names")
    }
    # change dimension order for efficiency: 
    # first time, then commons, then others
    dimord <- c(time_dim, setdiff(dimid_w, time_dim))
    w <- apermArray(w, first = dimord)
    erp <- apermArray(erp, first = dimord)
    # computations
    out <- fnDims(
        erp, dimid_w, warpArray,
        arg_list = list(w = w, time = time_points, 
                        return_array = FALSE))
    # return
    apermArray(out, dimid_erp)
}
