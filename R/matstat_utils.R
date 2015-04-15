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
#' x <- matrixIP(sample(1:3, 4e5, TRUE), 1e5, 4)
#' 
#' # compare computation times
#' system.time(un_x <- unique(x))
#' system.time(un_x_fast <- fastUnique(x))
#' 
#' # the same if x has dimension names
#' decorateDims(x, in_place = TRUE)
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
    out <- matrixIP(rep(x$values, x$lengths), ncol = max(x$matrix))
    if ("matrixrow" %in% names(x)) {
        out <- t(out)
    }
    # return
    out
}

