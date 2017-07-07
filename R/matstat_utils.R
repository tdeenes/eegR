#
# <<< utilities for matrix computations >>> --------
#

#' Fast version of unique for matrices
#'
#' \code{fastUnique} finds the unique rows or columns of a matrix. It is 
#' much faster and also more memory-efficient than 
#' \code{\link{unique.matrix}}, especially if the input matrix is not a 
#' character matrix.
#' @param x matrix
#' @param margin 1 (default) or 2, referring to rows and columns of the 
#' matrix, respectively. Can be a character string which is interpreted
#' as the name of the dimension if \code{x} has named dimnames. 
#' @param freq logical value if the frequency of the given rows or columns 
#' should also be returned (default: FALSE)
#' @export
#' @return \code{fastUnique} returns a matrix without duplicated units. If 
#' requested, the matrix has a \code{freq} attribute containing the frequency
#' of the given rows or columns.  
#' @examples
#' # create a matrix of random integer values ranging from 1 to 3
#' x <- matrix_(sample(1:3, 4e5, TRUE), 1e5, 4)
#' 
#' # compare computation times
#' system.time(un_x <- unique(x))
#' system.time(un_x_fast <- fastUnique(x))
#' 
#' # the same if x has dimension names
#' decorateDims_(x)
#' system.time(un_x <- unique(x))
#' system.time(un_x_fast <- fastUnique(x))
#' 
#' # the results are identical
#' stopifnot(identical(un_x, un_x_fast))
#' 
#' # fastUnique can also return frequency counts
#' system.time(un_x_withfreq <- fastUnique(x, freq = TRUE))
#' 
#' # check that frequencies add up to the number of rows of x
#' stopifnot(sum(attr(un_x_withfreq, "freq")) == nrow(x))
fastUnique <- function(x, margin = 1L, freq = FALSE) {
    x_class <- names(
        which(vapply(c("data.table", "data.frame", "matrix"), 
                     function(cl) inherits(x, cl), 
                     FALSE))[1]
        )
    if (is.na(x_class)) 
        stop("Provide a matrix, data.frame, or data.table as input!")
    if (is.character(margin) && !is.null(names(dimnames(x)))) {
        margin <- match(margin, names(dimnames(x)))
    }
    margin <- as.integer(margin)
    if (length(na.omit(margin)) != 1 || margin < 0L || margin > 2L) {
        stop("Wrong margin argument")
    }
    #
    if ((margin == 1L & nrow(x) == 1L) | (margin == 2L & ncol(x) == 1L))
        return(x)
    #
    dimn <- dimnames(x)
    if (margin > 1) x <- t(x)
    x <- as.data.table(x)
    setkey(x, NULL)
    #
    if (freq) {
        if (is.null(dimn[[margin]])) {
            out <- x[, list(X_counts_X = .N), by = names(x)]
            counts <- out[, X_counts_X]
            out[, X_counts_X := NULL]
        } else {
            out <- x[, list(X_counts_X = .N, X_idx_X = min(.I)), by = names(x)]
            counts <- out[, X_counts_X]
            idx <- out[, X_idx_X]
            out[, X_counts_X := NULL]
            out[, X_idx_X := NULL]
        }
    } else {
        if (is.null(dimn[[margin]])) {
            out <- unique(x)
        } else {
            out <- x[, list(X_idx_X = min(.I)), by = names(x)]
            idx <- out[, X_idx_X]
            out[, X_idx_X := NULL]
        }
    }
    if (margin > 1) out <- t(out)
    out <- do.call(paste0("as.", x_class), list(out))
    if (!is.null(dimn[[margin]])) dimn[[margin]] <- dimn[[margin]][idx]
    setattr(out, "dimnames", dimn)
    if (freq) setattr(out, "freq", counts)
    # return
    out
}

#' Maximum number of consecutive TRUE values in each column/row of a matrix
#' 
#' \code{consectrue} is helper function of \code{\link{tanova}}. It 
#' computes the longest streak of TRUE values in each column (or row) of a 
#' two-dimensional input.  
#' @param x logical vector, matrix or data.frame
#' @param col logical value if the search shall be performed column-wise (TRUE,
#' the default). Set to FALSE for row-wise counting.
#' @export
#' @keywords internal
consectrue <- function(x, col = TRUE, any_NA = NULL) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.list(x) || length(dim(x)) > 2L) {
        stop("x must be a vector, matrix or data.frame")    
    }
    if (is.null(any_NA)) any_NA <- anyNA(x)
    if (any_NA) stop("Missing values in x are not allowed")
    if (!is.logical(x)) stop("x must contain logical values")
    if (length(dim(x)) < 2L) {
        x <- matrix(x)
        col <- TRUE
    }
    if (!col) x <- t(x)
    out <- consectrueRcpp(x)
    names(out) <- colnames(x)
    # return
    out
}

#' Collapse row or column elements of a matrix
#' 
#' \code{collapse} is a fast version of 
#' \code{apply(x, along_dim, paste, collapse = "")}. It also works on vectors.
#' @param x vector, matrix or data.frame
#' @param col logical value; if TRUE (default), elements are collapsed 
#' column-wise.
#' @export
#' @return \code{collapse} always returns a character vector.
collapse <- function(x, col = TRUE) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.list(x)) stop("Provide a vector, matrix or data.frame")
    if (is.null(dim(x))) x <- matrix(x)
    if (length(dim(x)) > 2L) stop("Provide a vector, matrix or data.frame")
    if (!is.character(x)) storage.mode(x) <- "character"
    charmatCollapse(x, as.integer(col))
}

#' Fast replacement for \code{sweep}
#' 
#' \code{sweepMatrix} is a faster and less memory-hungry version of 
#' \code{sweep(x, MARGIN, STATS, FUN = c("-", "+", "*", "/", "<", ...))}.
#' @param x numeric matrix
#' @param MARGIN an integer index of the dimension to sweep over (1 for rows,
#' 2 for columns). Can be a character value if the input matrix has named 
#' dimension names.
#' @param STATS a numeric vector; the summary statistic which is to be swept 
#' out. It is recycled to match the size of \code{MARGIN}.
#' @param FUN character value, one of "+", "-", "*", "/", "^", "==", "!=", "<", 
#' "<=", ">=", ">"
#' @param has_NA logical value indicating if there is any missing value in 
#' \code{x}. If NULL (default), it is explicitly checked by calling 
#' \code{\link{anyNA}} on \code{x}.
#' @export
#' @return \code{sweepMatrix} returns a matrix of the same type as
#' \code{FUN(x[1, 1], y[1])}
#' @examples
#' # create a matrix and choose its first row as a statistic
#' mat <- matrix(1:20, 4, 5)
#' stat <- mat[1,]
#' 
#' # suppose we have missing values in the matrix and the statistic
#' mat[2, 3] <- NA
#' stat[1] <- NA
#' 
#' # compute and print the result
#' ( result <- sweepMatrix(mat, 2, stat, "-") )
#' 
#' # check
#' stopifnot(all.equal(result, 
#'                     sweep(mat, 2, stat, "-")))
#' stopifnot(identical(result[1, ], c(NA, rep(0L, ncol(mat) - 1))))
sweepMatrix <- function(x, MARGIN, STATS, 
                        FUN = c("-", "+", "*", "/", "^", 
                                "==", "!=", "<", "<=", ">", ">="),
                        has_NA = NULL) {
    assertMatrix(x, mode = "numeric")
    if (is.character(MARGIN)) MARGIN <- match(MARGIN, names(dimnames(x)))
    MARGIN <- as.integer(MARGIN[[1]])
    if (length(MARGIN) < 1L || !MARGIN %in% 1:2) stop("Wrong MARGIN provided")
    assertNumeric(STATS)
    STATS <- rep_len(STATS, dim(x)[MARGIN])
    FUN <- FUN[[1]]
    Rfun <- match.fun(FUN)
    if (MARGIN == 1L) {
        out <- Rfun(x, STATS)
    } else {
        if (is.double(x) && !is.double(STATS)) {
            storage.mode(STATS) <- "double"
        } else if (is.double(STATS) && !is.double(x)) {
            storage.mode(x) <- "double"
        }
        if (is.null(has_NA) || has_NA) {
            has_NA_x <- anyNA(x)
            has_NA_stats <- anyNA(STATS)
            has_NA <- if (has_NA_x || has_NA_stats) TRUE else FALSE
        } else {
            has_NA_x <- has_NA_stats <- FALSE 
        }
        if (has_NA_stats) {
            indNA <- rep(is.na(STATS), each = nrow(x))
            if (has_NA_x) {
                indNA <- indNA | is.na(x)
            }
        } else if (has_NA_x) {
            indNA <- is.na(x)
        } 
        cppFn <- 
            if (FUN %in% c("+", "-", "*")) {
                "eegR_sweepcol_multitype_cpp"
            } else if (FUN %in% c("/", "^")) {
                "eegR_sweepcol_double_cpp"
            } else if (FUN %in% c("==", "!=", "<", "<=", ">", ">=")) {
                "eegR_sweepcol_logical_cpp"
            } else {
                warning("FUN not implemented - fallback to sweep()")
                return(sweep(x, MARGIN, STATS, Rfun))
            }
        out <- .Call(eval(cppFn), PACKAGE = 'eegR', x, STATS, FUN)
        out_mode <- storage.mode(Rfun(x[1,1], STATS[1]))
        if (storage.mode(out) != out_mode) {
            warning("conversion was needed - check results")
            storage.mode(out) <- out_mode
        }
        if (has_NA) out[indNA] <- NA
    }
    # return
    out
}

# rle on the columns or rows of a matrix
matrixRle <- function(x, col = TRUE) {
    if (length(dim(x)) > 2 || (!is.atomic(x) && !is.data.frame(x)) )
        stop("Provide vector, matrix, or data.frame as input!")
    if (length(dim(x)) <= 1) {
        col <- TRUE
        x <- matrix(x)
    }
    out <- if (col) rleRcpp(x) else rleRcpp(t(x))
    if (!col) names(out)[3] <- "matrixrow"
    # return
    out
}

# inverse of matrixRle
inverse.matrixRle <- function(x) {
    out <- matrix_(rep(x$values, x$lengths), ncol = max(x$matrix))
    if ("matrixrow" %in% names(x)) {
        out <- t(out)
    }
    # return
    out
}

