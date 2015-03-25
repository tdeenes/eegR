
#
# <<< general functions on arrays >>> --------
#

#' Call a function on array chunks
#' 
#' \code{chunkify} creates subsets of an array according to various chunking
#' schemes and performs the given function on these data chunks. The resulting
#' array is identical to the result of a function call on the whole array.
#' @param dat matrix or array with named dimnames
#' @param fun name of the function (can be a character string). See details.
#' @param arg_list a list of parameters passed to \code{fun}
#' @param chunks a named list which defines the chunking scheme. See details.
#' @details The main purpose of \code{chunkify} is twofold: 1) it can decrease
#' memory load if a memory-intensive function is called on a large array; 2) it
#' enables within-dimension computations (e.g. perform baseline correction 
#' separately for various levels of the trial dimension). The parameter 
#' \code{chunks} is a named list; each element must correspond to a named 
#' dimension of \code{dat}. If the list element contains one integer, it gives
#' the number of chunks for the given dimension. If it contains as many elements
#' as the length of the given dimension, it defines chunk code; i.e., 
#' which chunk the given level of the dimension belongs to. Otherwise, the
#' list elements are expanded (recycled to the dimensions' length) and 
#' interpreted as chunk codes. \code{fun} must be a function which operates on
#' arrays (or matrices), and returns a matrix or array with named dimensions.
#' Chunking can be based only on those dimensions which are not affected by
#' \code{fun}. 
#' @export
#' @return An array (or vector/matrix); its shape and dimension names are 
#' identical to the return value of \code{fun} called on \code{dat} without 
#' chunking
chunkify <- function(dat, fun, arg_list = NULL, chunks = NULL) {
    chunkCheckFn <- function(x, name) {
        dimlen <- dimlens[name]
        if (is.character(x)) {
            x <- match(x, unique(x))
        } else if (is.logical(x)) {
            x <- as.integer(x) + 1L            
        } else if (!is.integer(x)) {
            x <- as.integer(x)
        } 
        out <- 
            if (length(x) == 1) {
                rep(seq_len(x), 
                    vapply(parallel::splitIndices(dimlen, x), length, 0L))
            } else if (length(x) == dimlen) {
                as.integer(x)
            } else {
                rep(as.integer(x), length.out = dimlen)
            }
        out
    }
    # 
    if (is.null(chunks)) {
        return( do.call(fun, append(list(dat), arg_list)) )
    }
    if (!is.list(chunks) && is.null(names(chunks))) {
        stop("chunks must be a named list")
    }
    dimlens <- dim(dat)
    names(dimlens) <- names(dimnames(dat))
    #
    chunks <- setNames(
        mapply(chunkCheckFn, chunks, names(chunks),
               SIMPLIFY = FALSE, USE.NAMES = FALSE),
        names(chunks))
    chunkgrid <- expand.grid(lapply(chunks, unique), 
                             KEEP.OUT.ATTRS = FALSE,
                             stringsAsFactors = FALSE)
    message(sprintf("\n==== Analyze Chunk 1/%i ====", nrow(chunkgrid)))
    subs <- mapply("==", chunks, chunkgrid[1, ], SIMPLIFY = FALSE)
    x <- subsetArray(dat,
                     subsets = subs,
                     drop = FALSE)
    message()
    outpart <- do.call(fun, append(list(x), arg_list))
    #
    if (nrow(chunkgrid) == 1) return( outpart )
    #
    out.dims <- dim(outpart)
    out.dims[match(names(chunks), names(dimnames(outpart)))] <- 
        dimlens[names(chunks)]
    out.dimnames <- setNames(
        lapply(names(dimnames(outpart)), function(n) {
            if (n %in% names(chunks)) {
                dimnames(dat)[[n]]
            } else {
                dimnames(outpart)[[n]]
            }}),
        names(dimnames(outpart))
    )
    na <- NA
    storage.mode(na) <- storage.mode(outpart)
    out <- arrayIP(na, out.dims, out.dimnames)
    for (i in 1:nrow(chunkgrid)) {
        if (i > 1) {
            message(sprintf("\n==== Analyze Chunk %i/%i ====", 
                            i, nrow(chunkgrid)))
            subs <- mapply("==", chunks, chunkgrid[i, ], SIMPLIFY = FALSE)
            x <- subsetArray(dat,
                             subsets = subs,
                             drop = FALSE)
            outpart <- do.call(fun, append(list(x), arg_list))
        }
        subsetArray(
            out,
            subsets = subs) <- outpart
        gc()
    }
    message("\n==== Analysis completed ====")
    # return
    out
}


#' Collapse data in array format across specific dimensions
#' 
#' \code{avgDims} collapses data in array format across specific dimensions by 
#' simple averaging
#' @param dat matrix or array to average
#' @param dims character vector of dimension names or numeric vector of 
#' dimension indices across which averaging should ocuur. If character vector, 
#' dat must have named dimnames. 
#' @param na_rm logical variable; if TRUE (default), missing values are removed 
#' before averaging
#' @export
#' @return An array (or matrix) with averaged data (the number of dimensions of 
#' dat is decreased by the length of dims)
avgDims <- function(dat, dims, na_rm = TRUE) {
    if (is.character(dims)) dims <- match(dims, names(dimnames(dat)))
    out <- aperm(dat, c(dims, seq_along(dim(dat))[-dims]))
    out <- colMeans(out, na.rm = na_rm, dims = length(dims))
    if (is.vector(out)) {
        dim(out) <- c(length(out), 1)
        dimnames(out) <- c(dimnames(dat)[-dims], list(NULL))
    }
    # return
    out
}

#' Call a user-specified function on whatever dimension of an array
#'
#' \code{fnDims} calls a user-specified function on whatever dimension(s) of a 
#' matrix or array. Function can be a vectorized function in which case fnDims 
#' is very fast.
#' @param dat matrix or array on which function should be called
#' @param target_dim The dimension(s) of dat on which the function should be 
#' performed. Can be an integer or numeric vector of dimension index(es) or a 
#' character vector if dat has named dimension names.
#' @param target_fn function definition
#' @param arg_list a list of arguments in addition to dat for target_fn 
#' (default: NULL)
#' @param newdims a named list of character vector(s) defining the dimension 
#' names of the output of target_fn (see Details)
#' @param vectorized logical value; should be set to TRUE if target_fn expects 
#' two-dimensional input (e.g. colSums)
#' @param columnwise logical value; if vectorized is TRUE, columnwise indicates 
#' if target_fn operates column-wise (TRUE, default) or row-wise (FALSE)
#' @param keep_dimorder logical value; if TRUE, and the returned matrix or array 
#' is of the same size as dat, the order of dimensions in the returned matrix 
#' or array will be the same as in dat (default = FALSE)
#' @param useparallel logical value; if TRUE, a snow-type parallel backend is 
#' used (default: FALSE)
#' @param cl cluster definition if useparallel is TRUE
#' @param ncores number of cpus if useparallel is TRUE and clust is NULL
#' @details This function works by transforming dat to a matrix with a row 
#' dimension corresponding to target_dim and calling target_fn on each column 
#' of the matrix repeatedly, or if vectorized is TRUE, by feeding the transformed
#' matrix directly to target_fn. After calling target_fn, the results are
#' back-transformed to the original array format. If target_fn produces one 
#' value per input vector, target_dim is dropped from the resulting array. 
#' If target_fn produces an equal-length vector as the input, and newdims 
#' parameter is not provided, fnDims assumes that dimensions correspond to the 
#' original input. If target_fn outputs a list, no back-transformation is done.
#' @export
#' @return see Details
#' @seealso \code{\link{avgDims}} if target_fn is one of mean/colMeans/rowMeans,
#' avgDims might be a better choice
fnDims <- function(dat, target_dim, target_fn, arg_list = NULL, 
                   newdims = list(), vectorized = FALSE, columnwise = TRUE, 
                   keep_dimorder = FALSE, 
                   useparallel = FALSE, cl = NULL, ncores = NULL) {
    if (is.data.frame(dat)) dat <- as.matrix(dat)
    if (is.character(target_dim)) 
        target_dim <- match(target_dim, names(dimnames(dat)))
    dimord <- c(target_dim, seq_along(dim(dat))[-target_dim])
    dims <- dim(dat)[dimord]
    dims.n <- dimnames(dat)[dimord]
    out <- 
        if (length(target_dim) == 1L) {
            array2mat(dat, target_dim, return_attributes=FALSE, 
                      keep_dimnames = FALSE)
        } else {
            mergeDims(dat, list(target_dim, seq_along(dim(dat))[-target_dim]),
                      return_attributes=FALSE, keep_dimnames = FALSE)
        }
    if (vectorized) {
        out <- do.call(target_fn, append(list(out), arg_list))
    } else {
        if (!columnwise) out <- t(out)
        if (!useparallel) {
            out <- sapply(1:ncol(out), function(i) 
                do.call(target_fn, append(list(out[, i]), arg_list)))
        } else {
            if (is.null(ncores)) ncores <- detectCores()
            if (is.null(cl)) cl <- makePSOCKcluster(ncores)
            tempfn <- function(i) 
                do.call(target_fn, append(list(out[, i]), arg_list))
            clusterExport(cl, 
                          varlist = c("tempfn","out","arg_list","target_fn"),
                          envir = environment())
            out <- parLapply(cl, 1:ncol(out), tempfn)
            stopCluster(cl)
            rm(cl)
            if (is.vector(out[[1]])) {
                tempdim <- c(length(out[[1]]), length(out))
            } else {
                tempdim <- c(dim(out[[1]]), length(out))
            }
            out <- drop(arrayIP(unlist(out, use.names = FALSE), tempdim))
        }
    }
    if (length(out) == 1 | is.list(out)) {
        return(out)
    } else if (length(dim(out)) <= 1) {
        if (columnwise) {
            dims <- dims[-seq_along(target_dim)]
            dims.n[seq_along(target_dim)] <- NULL
        } else {
            dims <- dims[seq_along(target_dim)]
            dims.n <- dims.n[seq_along(target_dim)]
        }
    } else if (length(out) != length(dat) || length(newdims) > 0) {
        if (length(newdims) == 0)
            newdims <- list(function.values = 1:nrow(out))
        if (columnwise) {
            dims <- c(vapply(newdims, length, 0L), 
                      dims[-seq_along(target_dim)])
            dims.n <- append(newdims, 
                             dims.n[-seq_along(target_dim)])
        } else {
            dims <- c(vapply(newdims, length, 0L), 
                      dims[seq_along(target_dim)])
            dims.n <- append(newdims, 
                             dims.n[seq_along(target_dim)])
        }
    } 
    arrayIP(out, dims, dims.n)
    if (length(out) == length(dat) && keep_dimorder) {
        out <- aperm(out, names(dimnames(dat)))
    }
    # return
    out
}

#' Compute averages of bins on a dimension of an array
#'
#' \code{avgBin} computes averages of bins on a dimension of an array. The
#' size of bins or directly the bins themselves can be specified by the user.
#' @param dat matrix or array
#' @param target_dim name of binned dimension of dat
#' @param bin_length number of data points in one bin. If rolling parameter is 
#' FALSE (default), bin_length must match the length of the target dimension
#' @param bin_ind an abbreviation for bin indices; a numeric or character vector
#' or a factor which provides the bin membership for each data point of 
#' target_dim. If bin_ind is provided, bin_length is not used. 
#' @param rolling logical variable; if yes, avgBin computes rolling (or moving) 
#' window averages (the size of the moving window is determined by bin_length)
#' @param newnames one of "avgs" (default), NULL, or a character vector defining
#' the new dimension names of the binned dimension. If "avgs", the original
#' names of target_dim are averaged per each bin. Not relevant if rolling is 
#' TRUE (which reuses the original labels). 
#' @export
#' @return A matrix or array, depending on the input
avgBin <- function(dat, target_dim, bin_length = NULL, bin_ind = NULL, 
                   rolling = FALSE, newnames = "avg") {
    if (rolling) {
        if (is.null(bin_length)) {
            stop("Provide bin_length if rolling_avg is set to TRUE")
        }
        out <- fnDims(dat, target_dim, rollFun, 
                      arg_list = list(width = bin_length, FUN = mean), 
                      vectorized = TRUE)
    } else {
        avgfn <- function(x, y) y %*% x
        if (is.null(bin_length) & is.null(bin_ind)) {
            stop("Either bin_length or bin_ind has to be provided!")
        } else if (!is.null(bin_ind)) {
            bin_ind <- as.numeric(factor(bin_ind))
            multmat <- diag(max(bin_ind))[, bin_ind]
            multmat <- multmat/colSums(multmat)
            out <- fnDims(dat, target_dim, avgfn, arg_list = list(y = multmat), 
                          vectorized = TRUE)
        } else {
            dimlen <- length(dimnames(dat)[[target_dim]])
            rest <- dimlen %% bin_length
            bins <- dimlen %/% bin_length
            if (rest > 0) {
                keep <- rep(TRUE, dimlen)
                keep[(dimlen - rest + 1L):dimlen] <- FALSE
                dat <- subsetArray(dat, listS(.target.dim = keep))
                warning("Target dimension length is not multiple of bin_length")
            }
            multmat <- diag(bins)[, rep(1:bins, each = bin_length)]/bin_length
            out <- fnDims(dat, target_dim, avgfn, arg_list = list(y = multmat), 
                          vectorized = TRUE)
        }
        if (newnames == "avg") {
            newdim <- list(as.character(
                multmat %*% as.numeric(dimnames(dat)[[target_dim]])))
        } else if (is.null(newnames)) {
            newdim <- dimnames(out)[1]
        } else {
            newdim <- list(newnames)
        }
        names(newdim) <- target_dim
        dimnames(out) <- c(newdim, dimnames(out)[-1])
    }
    out <- aperm(out, names(dimnames(dat)))
    # return
    out
}

#' Scale array across arbitrary dimension(s)
#' 
#' \code{scaleArray} centers and possibly scales the dimensions of an array 
#' across the dimension(s) chosen by the user
#' @param dat numeric matrix or array
#' @param by_dims,along_dims integer or character vector referring to the 
#' dimension(s) of dat across/along which centering/scaling should occur. See
#' Details.
#' @param center logical value if means should be subtracted (default: TRUE)
#' @param scale logical value if scaling should occur (default: TRUE)
#' @param base_subset argument passed to \code{\link{subsetArray}} before
#' computing the mean and/or the scale
#' @param center_subset the same as base_subset but only applies for the 
#' centering step
#' @param scale_subset the same as base_subset but only applies for the 
#' scaling step
#' @param keep_dimorder logical value; if TRUE (default), the order of the 
#' dimensions is kept intact, otherwise the channel dimension will be 
#' the first dimension in the resulting matrix or array
#' @details \code{scaleArray} is theoretically similar to \code{\link{scale}},
#' and can be conceived of as if one merges the \code{by_dims} dimensions and
#' also the \code{along_dims} dimensions of an array, and calls \code{scale}
#' on the resulting matrix. Thus, the arguments \code{by_dims} and 
#' \code{along_dims} are redundant: use whatever is more convenient in the 
#' given use case (\code{by_dims} defines the rows, \code{along_dims} the 
#' columns).
#' @export
#' @return A numeric matrix or array of the same size as the input object
#' @examples
#' # example dataset
#' data(erps)
#' 
#' # center and scale the data to normalize the subjects (coded with id in the 
#' # array), i.e. all subjects should have mean = 0 and SD = 1
#' scaled_erps <- scaleArray(erps, along_dims = "id")
#' means <- apply(scaled_erps, "id", mean)
#' sds <- apply(scaled_erps, "id", sd)
#' stopifnot(all.equal(unname(means), rep(0, length(means))))
#' stopifnot(all.equal(unname(sds), rep(1, length(means))))
#' 
#' # you might want to normalize each subject in each experimental condition 
#' # based only on the pre-stimulus activity
#' prestim <- as.numeric(dimnames(erps)$time) < 0
#' base <- list(time = prestim)
#' scaled_erps <- scaleArray(erps, by_dims = c("chan", "time"),
#'                           base_subset = base)
#' means <- apply(subsetArray(scaled_erps, base), 
#'                c("stimclass", "pairtype", "id"), 
#'                mean)
#' sds <- apply(subsetArray(scaled_erps, base), 
#'              c("stimclass", "pairtype", "id"), 
#'              sd)
#' stopifnot(all.equal(as.vector(means), rep(0, length(means))))
#' stopifnot(all.equal(as.vector(sds), rep(1, length(sds))))
scaleArray <- function(dat, by_dims = NULL, along_dims = NULL,
                       center = TRUE, scale = TRUE,
                       base_subset = NULL, 
                       center_subset = NULL, scale_subset = NULL,
                       keep_dimorder = TRUE) {
    sweepFn <- function(x, stats, FUN) {
        setattr(stats, "dim", NULL)
        FUN <- match.fun(FUN)
        FUN(x, stats)
    }
    #
    if (!center & !scale) return(dat)
    if (!is.null(base_subset) & is.null(center_subset)) {
        center_subset <- base_subset
    }
    if (!is.null(base_subset) & is.null(scale_subset)) {
        scale_subset <- base_subset
    }
    #
    na_rm <- anyNA(dat)
    if (is.character(by_dims))
        by_dims <- match(by_dims, names(dimnames(dat)))
    if (is.character(along_dims)) 
        along_dims <- match(along_dims, names(dimnames(dat)))
    if (is.null(by_dims) & is.null(along_dims)) 
        stop("Either by_dims or along_dims must be provided")
    if (!is.null(along_dims)) {
        if (is.null(by_dims)) {
            by_dims <- setdiff(seq_along(dim(dat)), along_dims)
        } else if (!setequal(c(by_dims, along_dims), seq_along(dim(dat)))) {
            stop("Both by_dims and along_dims are provided, but do not match with the dimensions of dat")
        }
    }
    if (anyNA(by_dims) || !all(by_dims %in% seq_along(dim(dat)))) 
        stop("Wrong by_dims or along_dims argument")
    along_dims <- setdiff(seq_along(dim(dat)), by_dims)
    out <- aperm(dat, c(along_dims, by_dims))
    if (center) {
        stat <- 
            if (!is.null(center_subset)) {
                rowMeans(subsetArray(out, center_subset, drop = FALSE),
                         na.rm = na_rm, dims = length(along_dims))
            } else {
                rowMeans(out, na.rm = na_rm, dims = length(along_dims))
            }
        out <- sweepFn(out, stat, "-")
    }
    if (scale) {
        if (!is.null(scale_subset)) {
            subs <- subsetArray(out, scale_subset, drop = FALSE)
            Sum <- rowSums(subs^2, na.rm = na_rm, dims = length(along_dims))
            N <- 
                if (na_rm) {
                    rowSums(!is.na(subs), na.rm = FALSE, 
                            dims = length(along_dims))
                } else {
                    prod(dim(subs)[-seq_along(along_dims)])
                }
        } else {
            Sum <- rowSums(out^2, na.rm = na_rm, dims = length(along_dims))
            N <- 
                if (na_rm) {
                    rowSums(!is.na(out), na.rm = FALSE, 
                            dims = length(along_dims))
                } else {
                    prod(dim(out)[-seq_along(along_dims)])
                }
        }
        stat <- sqrt(Sum / (N - 1))
        out <- sweepFn(out, stat, "/")
    }
    # return
    if (keep_dimorder) {
        aperm(out, order(c(along_dims, by_dims)))
    } else {
        out
    }
}

