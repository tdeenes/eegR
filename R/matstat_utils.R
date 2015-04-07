#
# <<< utilities for matrix computations >>> --------
#

#' Fast version of unique for matrices
#'
#' \code{fastUnique} finds the unique rows or columns of a matrix.
#' @param x matrix
#' @param units_in_rows units are in rows (TRUE, default) or columns of x
#' @export
#' @return matrix without duplicated units
fastUnique <- function(x, units_in_rows = TRUE) {
    if (!is.matrix(x)) stop("Provide a matrix as input!")
    if (units_in_rows) x <- t(x)
    dupl <- duplicated(str_collapse(lapply(1:ncol(x), function(i) x[,i])))
    x <- x[, !dupl, drop = F]
    if (units_in_rows) x <- t(x)
    # return
    x
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

