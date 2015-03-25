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

# TODO: add documentation
consectrue <- function(x, col = TRUE) {
    if (is.vector(x)) {
        x <- matrix(x, ncol = 1)
        col <- T
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

