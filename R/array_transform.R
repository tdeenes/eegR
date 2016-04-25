#
# <<< array reformatting functions >>> --------
#

#' Array transposition
#'
#' Transpose an array by permuting its dimensions and optionally resizing it.
#' This function is a simple wrapper around \code{\link{aperm.default}} with two
#' minimal enhancements: 1) It checks if the permutation order is the same as
#' the original order of the dimensions. In this case it simply returns a copy
#' of the input array without any unnecessary time-consuming manipulation. See
#' also NOTE. 2) It allows for partial definition of the permutation order.
#' See argument 'first'.
#' @param a the array to be transposed
#' @param perm the subscript permutation vector, usually a permutation of the
#' integers 1:n, where n is the number of dimensions of a. When a has named
#' dimnames, it can be a character vector of length n giving a permutation of
#' those names. The default (used whenever perm has zero length) is to reverse
#' the order of the dimensions.
#' @param resize a flag indicating whether the vector should be resized as well
#' as having its elements reordered (default: TRUE)
#' @param first dimension(s) which should come first; either numeric index(es),
#' or if the array has named dimensions, 'first' can be a character vector.
#' Ignored if 'perm' is not NULL.
#' @param keep_attributes. see NOTE (default: FALSE)
#' @export
#' @return A transposed version of array \code{a}, with subscripts permuted as
#' indicated by the array \code{perm}. If \code{resize} is TRUE, the array is
#' reshaped as well as having its elements permuted, the dimnames are also
#' permuted; if \code{resize = FALSE} then the returned object has the same
#' dimensions as \code{a}, and the dimnames are dropped. In each case other
#' attributes are dropped unless \code{keep_attributes. = TRUE}. See NOTE!
#' @note The documentation of \code{aperm} wrongly states that "other
#' attributes are copied from \code{a}". Instead, \code{aperm} drops all
#' attributes other than \code{dim} and \code{dimnames}. To make
#' \code{apermArray} a replacement function of \code{aperm}, this unexpected
#' behaviour is also reproduced unless explicitly required by the user not to
#' do so.
#' @examples
#' # example data
#' data(erps)
#' str(erps) # five dimensions
#'
#' # a use case where argument 'first' is helpful: suppose we need 'time'
#' # as the first dimension
#' time_first <- apermArray(erps, first = "time")
#'
#' # compare to the more cumbersome perm = setdiff(...) solution
#' time_first_aperm <- aperm(erps,
#'                           perm = c("time",
#'                                    setdiff(names(dimnames(erps)), "time")))
#' stopifnot(identical(time_first, time_first_aperm))
#'
#' # create permutation order by simply replicating the original order;
#' # note that in this special case it is more efficient to request
#' # keep_attributes. = TRUE
#' perm <- seq_along(dim(erps))
#'
#' # timing
#' if (require(microbenchmark)) {
#'     microbenchmark(
#'         aperm(erps, perm),
#'         apermArray(erps, perm),
#'         apermArray(erps, names(dimnames(erps))),
#'         apermArray(erps, perm, keep_attributes. = TRUE),
#'         times = 200L
#'     )
#' }
#'
#' # check if identical
#' aperm_orig <- aperm(erps, perm)
#' stopifnot(identical(aperm_orig,
#'                     apermArray(erps, perm)))
#' stopifnot(identical(aperm_orig,
#'                     apermArray(erps, names(dimnames(erps)))))
#' stopifnot(all.equal(aperm_orig,
#'                     apermArray(erps, perm, keep_attributes. = TRUE),
#'                     check.attributes = FALSE))
#' stopifnot(identical(erps,
#'                     apermArray(erps, perm, keep_attributes. = TRUE)))
#'
#' # if resize = FALSE, dimension names are dropped
#' noresize <- apermArray(erps, perm, resize = FALSE)
#' stopifnot(is.null(dimnames(noresize)))
#' stopifnot(identical(unname(aperm_orig),
#'                     noresize))
#'
apermArray <- function(a, perm = NULL, resize = TRUE, first = perm,
                       keep_attributes. = FALSE) {
    assertArray(a, min.d = 1L, .var.name = "input array ('a')")
    dims <- dim(a)
    num_dims <- length(dims)
    orig_attribs <- attributes(a)
    orig_attribs <- orig_attribs[setdiff(names(orig_attribs),
                                         c("dim", "dimnames"))]
    if (is.null(perm) && !is.null(first)) {
        if (is.numeric(first)) {
            assertNumeric(first, lower = 1, upper = num_dims,
                          any.missing = FALSE, max.len = num_dims,
                          .var.name = "first")
        } else if (is.character(first)) {
            dimn <- names(dimnames(a))
            if (is.null(dimn))
                stop("Input array ('a') does not have named dimnames")
            assertCharacter(first,
                            any.missing = FALSE, max.len = num_dims,
                            .var.name = "first")
            first <- match(first, dimn)
            if (anyNA(first))
                stop(sprintf("Not all elements in 'first' match the dimenion names: %s",
                             paste(dimn, collapse = ", ")))
        }
        perm <- c(first, seq_len(num_dims)[-first])
    }
    if (!is.null(perm) &&
        (
         (is.character(perm) && identical(perm, names(dimnames(a)))) ||
         (is.numeric(perm) && identical(as.integer(perm), seq_along(dim(a))))
        )
       ) {
        x <- copy(a)
        if (keep_attributes.) {
            if (!resize) setattr(x, "dimnames", NULL)
            x
        } else {
            dims <- dim(a)
            if (!resize) {
                attributes(x) <- NULL
                setattr(x, "dim", dims)
            } else {
                dimn <- dimnames(a)
                attributes(x) <- NULL
                setattr(x, "dim", dims)
                setattr(x, "dimnames", dimn)
            }
            x
        }
    } else {
        x <- aperm.default(a, perm, resize)
        if (keep_attributes.) {
            for (i in names(orig_attribs))
                setattr(x, i, orig_attribs[[i]])
        }
        x
    }
}


#' Add dimension names to a matrix or array with (partially) missing dimension
#' names
#'
#' \code{decorateDims} is not unlike \code{\link{provideDimnames}}: it provides
#' dimension names for dimensions which have 'missing' dimension names
#' ('missing' in a broad sense: either NA, NULL or ""). The main difference is
#' that it has a more convenient interface for providing names of dimensions,
#' and via \code{decorateDims_} the user can modify the input data in place
#' (without creating even a shallow copy).
#' @name decorateDims
#' @param dat a matrix or array
#' @param .names either 1) a logical value whether dimension identifiers (names
#' of dimnames) should be checked and corrected (default: TRUE), or 2) a
#' character vector of new dimension identifiers for missing, NA or "" dimension
#' identifiers, or 3) a single character value after which 1, 2, ... is
#' appended to create unique identifiers (see Examples)
#' @param .dimnames either 1) a logical value whether dimnames should be checked
#' and corrected (default: TRUE), or 2) a list of character vectors of new
#' dimnames for dimensions with missing, NA or "" dimnames, or 3) a single
#' character value after which 1, 2, ... is appended to create unique dimension
#' levels (see Examples). If '.names' is not explicitly provided but '.dimnames'
#' is a named list, its names are considered as '.names'.
#' @note Use \code{decorateDims_} with extra care because it modifies in place
#' all objects which 'dat' refers to.
#' @return \code{decorateDims} returns a \emph{copy} of the original matrix or
#' array with non-null dimension names. \code{decorateDims} invisibly returns
#' the \emph{original} matrix or array with modified dimnames attribute.
#' @export
#' @seealso \code{\link{provideDimnames}} for a slightly different solution
#' @examples
#' # create a matrix without dimension names
#' ( mat <- matrix(letters[1:8], 2, 4) )
#'
#' # add default dimension names
#' ( mat2 <- decorateDims(mat) )
#'
#' # remove the second dimension identifier (name of the 'dimnames' attribute)
#' names(dimnames(mat2))[2] <- ""
#' mat2
#'
#' # create a new variable which is referenced to 'mat2'
#' mat3 <- mat2
#'
#' # modify the names of dimnames in place
#' decorateDims_(mat3, .names = "Dimension")
#' mat3
#' stopifnot(identical(names(dimnames(mat3)),
#'                     c("_Dim1", "Dimension2")))
#'
#' # NOTE that the attributes of 'mat2' has been also changed!
#' mat2
#' stopifnot(identical(dimnames(mat2), dimnames(mat3)))
#'
#' # decorate only with dimension identifiers, but no dimension names
#' ( mat4 <- decorateDims(mat,
#'                        .names = c("rows", "columns"),
#'                        .dimnames = FALSE) )
#' stopifnot(identical(dimnames(mat4), list(rows = NULL, columns = NULL)))
#'
decorateDims <- function(dat, .names = TRUE, .dimnames = TRUE) {
    assertArray(dat, min.d = 1L, .var.name = "dat")
    new_dimn <- fillMissingDimnames(dimnames(dat), dim(dat),
                                    .names, .dimnames)
    # return
    if (identical(new_dimn, dimnames(dat))) {
        dat
    } else {
        dimnames(dat) <- new_dimn
        dat
    }
}

#' @rdname decorateDims
#' @export
decorateDims_ <- function(dat, .names = TRUE, .dimnames = TRUE) {
    assertArray(dat, min.d = 1L, .var.name = "dat")
    new_dimn <- fillMissingDimnames(dimnames(dat), dim(dat),
                                    .names, .dimnames)
    setattr(dat, "dimnames", new_dimn)
    # return
    invisible(dat)
}

#' Find and fill missing dimension names and identifiers
#'
#' \code{fillMissingDimnames} is the workhorse function for
#' \code{\link{decorateDims}} and \code{\link{decorateDims_}}. It replaces
#' missing dimension names and identifiers and returns the modified 'dimnames'
#' attribute.
#' @param dimn a list of the original dimension names
#' @param .dim an integer vector of dimensions
#' @inheritParams decorateDims
#' @keywords internal
fillMissingDimnames <- function(dimn, .dim, .names = TRUE, .dimnames = TRUE) {
    # helper function
    checkForMissing <- function(x) {
        if (is.null(x)) {
            TRUE
        } else {
            is.na(x) | x == ""
        }
    }
    # argument checks and fast return if no missings are found
    check_names <- if (!identical(.names, FALSE)) TRUE else FALSE
    check_dimnames <- if (!identical(.dimnames, FALSE)) TRUE else FALSE
    check <- FALSE
    if (check_names && any(checkForMissing(names(dimn)))) {
        check <- TRUE
    }
    if (check_dimnames) {
        missings <- 
            if (is.null(dimn)) {
                as.list(rep_len(TRUE, length(.dim)))
            } else {
                lapply(dimn, checkForMissing)    
            }
        decor_dim <- vapply(missings, any, FALSE)
        if (any(decor_dim)) check <- TRUE
    }
    if (!check) return(dimn)
    #
    new_dimnames <-
        if (check_names && !check_dimnames) {
            vector("list", length(.dim))
        } else if (identical(.dimnames, TRUE) || is.null(.dimnames)) {
            rep(list(""), sum(decor_dim))
        } else if (!is.list(.dimnames)) {
            rep(list(.dimnames), sum(decor_dim))
        } else {
            rep_len(.dimnames, sum(decor_dim))
        }
    new_names <- if (is.logical(.names)) "_Dim" else .names
    if (is.null(new_names)) new_names <- names(new_dimnames)
    #
    if (!is.null(new_names))
        assertCharacter(new_names, any.missing = FALSE, .var.name = "new_names")
    if (!is.null(new_dimnames))
        assertList(new_dimnames, types = c("null", "character"),
                   any.missing = FALSE, .var.name = "new_dimnames")
    # main part
    if (is.null(dimn)) dimn <- vector("list", length(.dim))
    if (check_dimnames) {
        tempfn <- function(i, counter) {
            newdimn <- new_dimnames[[counter]]
            if (is.null(newdimn)) return(NULL)
            miss <- missings[[i]]
            # return
            if (length(miss) == 1L && length(newdimn) == 1L) {
                paste0(newdimn, seq_len(.dim[i]))
            } else if (length(newdimn) == 1L) {
                out <- dimn[[counter]]
                out[miss] <- paste0(newdimn, which(miss))
                make.unique(out, "__")
            } else {
                out <- rep_len(newdimn, .dim[i])
                out[!miss] <- dimn[[counter]][!miss]
                make.unique(out, "__")
            }
        }
        counter <- 1L
        for (i in which(decor_dim)) {
            dimn[i] <- list(tempfn(i, counter))
            counter <- counter + 1L
        }
    }
    dimn.n <- names(dimn)
    if (check_names) {
        if (is.null(dimn.n)) dimn.n <- rep("", length(dimn))
        ind <- dimn.n == "" | is.na(dimn.n)
        dimn.n[ind] <-
            if (length(new_names) == 1L) {
                paste0(new_names, which(ind))
            } else {
                out <- rep_len(new_names, sum(ind))
                make.unique(as.character(out), "__")
            }
        dimn.n <- make.unique(dimn.n, "__")
    }
    names(dimn) <- dimn.n
    # return
    dimn
}


#' Recode dimension names
#'
#' \code{recodeDimnames} renames dimension identifiers and dimension levels
#' (useful before plotting). \code{recodeDimnames_} does the same in place, 
#' thus without making any copy. 
#' @param dat an array with dimension names
#' @param dim_levels a named list of named character vectors. The name of the 
#' list element identifies the dimension, the names of the character vector 
#' refer to the original dimension levels, and the values of the character
#' vector correspond to the new dimension levels. See also Details.
#' @param dim_ids a named character vector which is used to rename the dimension
#' identifiers
#' @details \code{recodeDimnames} goes through the 'dim_levels' argument first.
#' If a character vector in 'dim_levels' is not named, it must be of the same
#' length as the length of the given dimension. 
#' @note Use \code{recodeDimnames_} with extra care because it modifies in place
#' all objects which 'dat' refers to.
#' @export
#' @examples
#' ## load example data
#' data(erps)
#' 
#' ## recode the 'subst', 'ident', and 'transp' levels of the 'stimclass' 
#' ## dimension to 'Substitution', 'Identical', and 'Transposition', 
#' ## respectively; also rename the 'stimclass' dimension name to 
#' ## 'Stimulus class'
#' erps2 <- recodeDimnames(erps, 
#'                         list(stimclass = c(subst = "Substitution",
#'                                            ident = "Identical",
#'                                            transp = "Transposition")),
#'                         c(stimclass = "Stimulus class"))
#'  
recodeDimnames <- function(dat, dim_levels = NULL, dim_ids = NULL) {
    dat <- copy(dat)
    recodeDimnames_(dat, dim_levels, dim_ids)
    # return
    dat
}

#' @describeIn recodeDimnames Modify by reference
#' @export
recodeDimnames_ <- function(dat, dim_levels = NULL, dim_ids = NULL) {
    dimn <- dimnames(dat)
    dimid <- names(dimn)
    dim_levels <- dim_levels[names(dim_levels) %in% dimid]
    dim_ids <- dim_ids[names(dim_ids) %in% dimid]
    if (!is.null(dim_levels)) {
        dimn[names(dim_levels)] <- lapply(
            names(dim_levels), 
            function(n) {
                old <- names(dim_levels[[n]])
                if (is.null(old)) old <- dimn[[n]]
                Replace(dimn[[n]], old, as.vector(dim_levels[[n]])) 
            })
    }
    if (!is.null(dim_ids)) {
        dimid <- Replace(dimid, names(dim_ids), as.vector(dim_ids))    
    }
    setattr(dimn, "names", dimid)
    setattr(dat, "dimnames", dimn)
    invisible(dat)
}


#' Drop singleton dimensions of an array
#'
#' \code{dropDims} drops singleton dimensions (whose lengths is 1) of a
#' multidimensional array. Compared to \code{\link[base]{drop}} and
#' \code{\link[abind]{adrop}} this function gives more control over which
#' singleton dimensions to drop.
#' @param x an array (or a matrix, or a list with \code{dim} attribute)
#' @param drop 1) either a single logical whether all singleton dimensions
#' should be dropped (default: TRUE); or 2) a logical, integer or character
#' vector indicating which dimensions to drop (if character, \code{x} must have
#' named dimensions). If \code{drop} is a vector (case 2), the referred
#' dimensions must be present in x.
#' @param keep an integer or character vector indicating those dimensions which
#' must remain in the returned value even if they are singletons (if character,
#' \code{x} must have named dimensions). Note that \code{keep} has a higher
#' priority than \code{drop}. Also note that to-be-kept dimensions which are
#' actually not present in \code{x} are simply ignored (instead of resulting in
#' error).
#' @param return_array logical value whether \code{dropDims} should return a
#' one-dimensional array instead of a vector even if all or all but one
#' dimension of \code{x} is dropped (default: TRUE)
#' @param named_vector logical value whether a vector result should be named
#' (TRUE, the default). Ignored if \code{return_array} is TRUE.
#' @param stop_if_missing logical value whether dropping or keeping
#' non-existent dimensions should result in error (TRUE, default), or should
#' be ignored (FALSE)
#' @export
#' @examples
#' # create example data
#' x <- array(1:4, c(2, 1, 2, 1),
#'            dimnames = list(dimA = letters[1:2],
#'                            dimB = "a",
#'                            dimC = LETTERS[1:2],
#'                            dimD = "z"))
#'
#' # drop all singleton dimensions
#' x0 <- dropDims(x)
#' stopifnot(identical(dim(x0), c(2L, 2L)))
#'
#' # drop all singleton dimensions but always keep dimD
#' x1 <- dropDims(x, keep = "dimD")
#' stopifnot(identical(dim(x1), c(2L, 2L, 1L)))
#'
#' # create a new example; a list with dim attribute
#' ( x <- array(list(1:2, letters[1:2]), c(2, 1, 1),
#'              dimnames = list(type = c("numeric", "character"),
#'                              single1 = "a",
#'                              single2 = "b")) )
#'
#' # drop the single1 dimension
#' x0 <- dropDims(x, "single1")
#' stopifnot(identical(dim(x0), c(2L, 1L)))
#'
#' # wrong dimension in drop: by default, it results in an error
#' # (note that x has only 3 dimensions, not 4)
#' x1 <- try(dropDims(x, drop = 1:4), silent = TRUE)
#' stopifnot(inherits(x1, "try-error"))
#'
#' # however, you can also ask for skipping those missing dimensions in 'drop'
#' # and 'keep'
#' x2 <- dropDims(x, drop = 1:4, stop_if_missing = FALSE)
#' stopifnot(identical(dim(x2), 2L))
#'
#' # by default, dropDims returns an array, even if it has only one dimension
#' ( x <- array(1:3, c(3, 1), list(dimA = letters[1:3], dimB = "A")) )
#' ( x0 <- dropDims(x) )
#' stopifnot(is.array(x0))
#'
#' # you can change this behaviour
#' ( x1 <- dropDims(x, return_array = FALSE) )
#' stopifnot(!is.array(x1))
#'
#' # vector results are named by default...
#' stopifnot(identical(names(x1), letters[1:3]))
#'
#' # ...but not necessarily
#' ( x2 <- dropDims(x, return_array = FALSE, named_vector = FALSE) )
#' stopifnot(is.null(names(x2)))
#' stopifnot(identical(unname(x1), x2))
#'
dropDims <- function(x, drop = TRUE, keep = NULL,
                     return_array = TRUE, named_vector = TRUE,
                     stop_if_missing = TRUE) {
    # helper function
    checkStop <- function(subset, fullset, dropkeep = c("drop", "keep"),
                          st = stop_if_missing) {
        if (is.logical(subset)) {
            if (length(subset) != length(fullset)) {
                stop(sprintf("if '%s' is a logical vector, its length must be equal to the number of dimensions in 'x'",
                             dropkeep))
            } else {
                return(rep_len(subset, length(fullset)))
            }
        } else {
            ind <- subset %in% fullset
            if (st && !all(ind)) {
                stop(sprintf("the dimensions of 'x' and the ones listed in '%s' do not match",
                             dropkeep))
            }
            fullset %in% subset
        }
    }
    # check if array
    assertArray(x, .var.name = "x")
    # return if drop is FALSE
    if (identical(drop, "FALSE")) return(x)
    # dimensions
    dims <- dim(x)
    dimn <- dimnames(x)
    dimid <- names(dimn)
    # check drop
    singleton <- dims == 1L
    drop <-
        if (isTRUE(drop)) {
            singleton
        } else if (is.character(drop)) {
            if (is.null(dimid)) {
                stop("x must have named dimensions if 'drop' is a character vector")
            }
            singleton & checkStop(drop, dimid, "drop")
        } else if (is.numeric(drop)) {
            singleton & checkStop(as.integer(drop), seq_along(dims), "drop")
        } else if (is.logical(drop)) {
            singleton & checkStop(drop, logical(length(dims)), "drop")
        } else {
            stop("x is of wrong type. Must be a single logical, or a logical/integer/character vector")
        }
    keep <-
        if (is.null(keep)) {
            !drop
        } else if (is.character(keep)) {
            if (is.null(dimid)) {
                stop("x must have named dimensions if 'keep' is a character vector")
            }
            !drop | checkStop(keep, dimid, "keep")
        } else if (is.numeric(keep)) {
            !drop | checkStop(keep, seq_along(dims), "keep")
        } else if (is.logical(keep)) {
            !drop | checkStop(keep, logical(length(dims)), "keep")
        }
    # return
    nrkeep <- sum(keep)
    asvec <- if (nrkeep <= 1L & !return_array) TRUE else FALSE
    if (asvec) {
        out <- copy(x)
        setattr(out, "dim", NULL)
        if (named_vector && nrkeep > 0L) {
            setattr(out, "names", dimn[[which(keep)]])
        } else {
            out
        }
    } else if (nrkeep == length(dims)) {
        x
    } else {
        out <- copy(x)
        dims <- dims[keep]
        dimn <- dimn[keep]
        setattr(out, "dim", dims)
        setattr(out, "dimnames", dimn)
    }
}

#' Fast in-place transformation to a matrix (without copy)
#'
#' \code{matrix_} transforms its data argument to a matrix by reference. If the
#' length of data remains the same, no copy is made at all, and it invisibly
#' returns the matrix. This is mostly useful for manipulating interim objects
#' in functions, and should not be used for interactive analyses.
#' @param x a data vector, matrix or array
#' @param nrow the desired number of rows
#' @param ncol the desired number of columns
#' @param byrow logical. If FALSE (the default), the matrix is filled by
#' columns, otherwise the matrix is filled by rows. If TRUE, it results in an
#' error!
#' @param dimnames A dimnames attribute for the matrix: NULL or a list of
#' length 2 giving the row and column names respectively. The list can be
#' named, and the list names will be used as names for the dimensions. An empty
#' list is treated as NULL, and a list of length one as row names.
#' @param force_length logical. If TRUE (the default), \code{matrix_} checks if
#' length(x)==nrow*ncol. If not, x is recycled or subsetted to the desired
#' length.
#' @param arg_check logical indicating if argument checks should be performed
#' (TRUE, the default). Do not set to FALSE unless you really know what you
#' are doing!
#' @note Use \code{matrix_} with extra care because it modifies in place all
#' objects which x refers to. If you want to avoid this, call
#' \code{x <- copy(x)} before calling \code{matrix_} or use the standard way as
#' described in the Note section of \code{\link{matrix}}. However, for input
#' objects created on-the-fly (e.g. a temporary vector), \code{matrix_} is safe
#' and more compact than the latter solution, and can be many times faster than
#' \code{matrix}.
#' @return \code{matrix_} invisibly returns a matrix without duplicating the
#' input values.
#' @export
#' @seealso \code{\link{array_}} for in-place transformation of x to an array
#' (without copy) and \code{\link{matrix}} for creating a new matrix without
#' modifying the original input
#' @examples
#' # create two vectors
#' x <- y <- 1:10
#'
#' # suppose 'x' should be a 2x5 matrix
#' matrix_(x, 2, 5)
#' str(x)
#'
#' # however, since 'x' was referenced to 'y', 'y' has been changed, too
#' str(y)
#'
#' # compare the timing for matrix creation
#' if (require(microbenchmark)) {
#'     microbenchmark(
#'         matrix = matrix(0L, 1e3, 1e3),
#'         matrix_ = matrix_(0L, 1e3, 1e3),
#'         times = 100L
#'     )
#' }
#'
matrix_ <- function(x, nrow, ncol, byrow = FALSE,
                    dimnames = NULL, force_length = TRUE, arg_check = TRUE) {
    if (arg_check) {
        # argument checks
        assertLogical(byrow, len = 1, .var.name = "byrow")
        if (byrow) stop("'byrow' must be FALSE")
        assertAtomic(x, .var.name = "x")
        if (missing(nrow) && missing(ncol)) {
            nrow <- length(x)
            ncol <- 1L
        } else if (missing(nrow)) {
            assertCount(ncol, .var.name = "ncol")
            nrow <- length(x)/as.integer(ncol)
        } else if (missing(ncol)) {
            assertCount(nrow, .var.name = "nrow")
            ncol <- length(x)/as.integer(nrow)
        }
        nrow <- as.integer(nrow)
        ncol <- as.integer(ncol)
        # dimension names
        if (!is.null(dimnames)) {
            if (grepl(deparse(substitute(x)),
                      deparse(substitute(dimnames)))) {
                dimnames <- copy(dimnames)
            }

        }
        # check length
        if (force_length && length(x) != nrow*ncol) {
            x <- rep_len(x, nrow*ncol)
        }
    }
    # set attributes
    setattr(x, "dim", c(nrow, ncol))
    setattr(x, "dimnames", dimnames)
    # return
    invisible(x)
}

#' Fast in-place transformation to an array (without copy)
#'
#' \code{array_} transforms its data argument to an array by reference. No
#' copy is made at all, and it invisibly returns the array.
#' @param x a data vector, matrix or array
#' @param dim the dim attribute for the array to be created, that is an integer
#' vector of length one or more giving the maximal indices in each dimension
#' @param dimnames either NULL or the names for the dimensions. This must a
#' list (or it will throw an error) with one component for each dimension,
#' either NULL or a character vector of the length given by dim for that
#' dimension. The list can be named, and the list names will be used as names
#' for the dimensions. If the list is shorter than the number of dimensions, it
#' is extended by NULLs to the length required.
#' @param force_length logical. If TRUE (the default), \code{array_} checks if
#' length(x)==nrow*ncol. If not, x is recycled or subsetted to the desired
#' length.
#' @param arg_check logical indicating if argument checks should be performed
#' (TRUE, the default). Do not set to FALSE unless you really know what you
#' are doing!
#' @note Use \code{array_} with extra care because it modifies in place all
#' objects which x refers to. See \code{\link{matrix_}} for further hints.
#' @return This function (invisibly) returns an array (or a matrix if
#' \code{length(dim)==2L}) without duplicating the input values.
#' @export
#' @seealso \code{\link{matrix_}} for in-place transformation of x to a matrix
#' (without copy) and \code{\link{array}} for creating a new array without
#' modifying the original input
array_ <- function(x, dim, dimnames = NULL,
                   force_length = TRUE, arg_check = TRUE) {
    if (arg_check) {
        # check arguments
        assertAtomic(x, .var.name = "x")
        assertNumeric(dim, .var.name = "dim")
        if (missing(dim)) {
            dim <- length(x)
        }
        if (!is.null(dimnames)) {
            if (grepl(deparse(substitute(x)),
                      deparse(substitute(dimnames)))) {
                dimnames <- copy(dimnames)
            }
        }
        dim <- as.integer(dim)
        # force length
        if (force_length && length(x) != prod(dim)) {
            x <- rep_len(x, prod(dim))
        }
    }
    # set attributes
    setattr(x, "dim", dim)
    setattr(x, "dimnames", dimnames)
    # return
    invisible(x)
}

#' Create subset indices for subsetArray
#'
#' \code{subsetIndices} check arguments and creates subset indices for
#' subsetArray and `subsetArray<-`.
#' @param subset. a possibly named list of subset indices (either logical,
#' integer or is* functional indices)
#' @param which_dims. the indices of the dimensions which should be subsetted
#' @param dim. the dimensions of the data to subset on
#' @param dimnames. the dimension names of the data to subset on
#' @keywords internal
subsetIndices <- function(subset., which_dims., dim., dimnames.) {
    dimid. <- names(dimnames.)
    if (is.null(subset.)) subset. <- list()
    assertList(subset.,
               types = c("logical", "integerish", "character", "function"),
               any.missing = FALSE, max.len = length(dim.),
               .var.name = "subset.")
    if (any(duplicated(names(subset.))))
        stop("Duplicated elements in the joint list of 'subset.' and '...'")
    if (length(subset.) == 0L)
        stop(paste0("The joint list of 'subset.' and '...' is empty. ",
                    "Provide either 'subset.' or use '...' or both."))
    if (is.null(which_dims.)) {
        if (anyDuplicated(dimid.)) {
            stop(paste0(
                "'dat' has duplicated dimension identifiers. ",
                "Provide 'which_dims.' and do not rely on the names of ",
                "the subset."))
        }
        which_dims. <- match(names(subset.), dimid.)
        if (anyNA(which_dims.)) {
            stop(paste0(
                "Not all names in the joint list of 'subset.' and '...' ",
                "correspond to existing dimensions in 'dat'. Provide ",
                "either 'which_dim' to disambiguate or set the names ",
                "in 'subset.' and/or '...' properly."))
        }
    } else {
        if (is.null(which_dims.))
            stop(paste0(
                "If 'subset.' is not a named list and/or the arguments ",
                "in '...' are not named, 'which_dims.' must be provided."))
        assertVector(which_dims., strict = TRUE, any.missing = FALSE,
                     len = length(subset.), unique = TRUE,
                     .var.name = "which_dims.")
        if (is.character(which_dims.)) {
            which_dims. <- match(which_dims., dimid.)
        }
    }
    #
    # subsetting indices
    ind <- as.list(rep(TRUE, length(dim.)))
    for (i in seq_along(subset.)) {
        subset_ <- subset.[[i]]
        dimid_ <- which_dims.[i]
        dimn_ <- dimnames.[[dimid_]]
        ind[[dimid_]] <-
            if (is.function(subset_)) {
                subset_(dimn_)
            } else {
                subset_
            }
    }
    #
    # return
    setattr(ind, "which_dims.", which_dims.)
    ind
}



#' Extract or replace a part of an array
#'
#' \code{subsetArray} is a convenience function for extracting or replacing a
#' part of an array which has dimension names
#' @name subsetArray
#' @usage
#' subsetArray(dat, subset.=list(), which_dims.=NULL, drop.=NULL, keep_attributes.=TRUE, ...)
#' @param dat array to be subsetted
#' @param subset. a (named) list of character, numeric, or logical vectors or
#' a subsetting function (see \code{\link[eegR]{is}}) indicating which levels
#' of which dimensions to subset (see Details). If 'subset.' is an unnamed
#' list, the argument 'which_dims.' must be provided.
#' @param which_dims. numeric or character indices of the dimensions which
#' should be subsetted. If 'which_dims.' is not NULL, 'which_dims.' is used and
#' the names of 'subset.' is ignored.
#' @param drop. either 1) NULL (the default), or 2) a logical value
#' (TRUE or FALSE), or 3) numeric or character indices of the dimensions which
#' should be dropped if they become singleton dimensions (i.e. have only one
#' level) after subsetting (see Details)
#' @param keep_attributes. a logical variable which determines if the result
#' inherits the custom attributes of the input (TRUE, default) or not
#' @param ... an alternative specification of the subsetting rule; one can
#' provide the subsetting vectors as named arguments, where the argument name
#' refers to the name of the subsetted dimension. For programmatic use,
#' 'subset.' is preferred because the names of dimensions might interfere with
#' the default argument names.
#' @details Names of 'subset.' or the indices of dimensions as given in
#' 'which_dim.' indicate which dimensions are to be subsetted in
#' the input array, and each list element indicates which levels of the given
#' dimension will be selected. If a list element is an empty vector, all levels
#' of the correspondig dimension will be selected. Further possibilities for
#' subsetting vectors:
#' \itemize{
#' \item{function: }{a function which returns a logical vector if called on
#' the dimension levels of the given dimension (see Examples for a use case
#' of \code{\link{isBetween}} and \code{\link{isPattern}}). Note that 'dat'
#' must have non-NULL dimension names for functional subsetting.}
#' \item{logical: }{logical vector of the same length as the given dimension,
#' denoting whether the given level of the dimension should be included in the
#' subset (TRUE) or not (FALSE)}
#' \item{integer: }{numeric indices of the dimension levels which should be
#' included in the subset}
#' \item{character: }{character vector of the dimension levels which should be
#' included in the subset. Note that 'dat' must have non-NULL dimension names
#' for character subsetting.}
#' }
#' The argument 'drop.' defines the procedure if a dimension becomes a singleton
#' dimension after subsetting. The default behaviour (\code{drop. = NULL}) is to
#' drop all subsetting dimensions but no others. If 'drop.' is FALSE, all
#' dimensions are kept, if TRUE, all singleton dimensions are dropped. If
#' 'drop.' is a numeric or character vector, its elements define which
#' dimensions to drop.
#' @export
#' @seealso See the \code{\link[eegR]{is}} family of functions for subsetting
#' dimensions by functional expressions; see also \code{\link[abind]{asub}}
#' for a less general approach.
#' @return The function returns a subset of the array or the array with replaced
#' values.
#' @export
#' @examples
#' # example data (see ?erps)
#' data(erps)
#' str(erps)
#'
#' # subsetting without knowing the exact order of dimensions and using
#' # various subsetting schemes
#' sub1 <- subsetArray(erps,
#'                     time = isBetween(0, 10),
#'                     chan = isPattern("Fp"),
#'                     stimclass = c(TRUE, FALSE, FALSE),
#'                     keep_attributes. = FALSE)
#'
#' # traditional subsetting
#' sub2 <- erps["A", , c("Fp1", "Fp2"), c("0", "2", "4", "6", "8", "10"), ]
#'
#' # the results are identical
#' stopifnot(identical(sub1, sub2))
#'
#' # the same for replacement
#' subsetArray(sub1, list(id = 1, pairtype = isSame("transp"))) <- NA
#' sub2["transp", , , 1] <- NA
#' stopifnot(identical(sub1, sub2))
#'
subsetArray <- function(dat, subset. = list(), which_dims. = NULL,
                        drop. = NULL, keep_attributes. = TRUE, ...) {
    # if NULL, return
    if (is.null(dat)) return(NULL)
    #
    # check arguments (dat, subset., which_dims.) and find indices
    assertArray(dat, min.d = 1L, .var.name = "dat")
    dat_d <- dim(dat)
    dat_dn <- dimnames(dat)
    dat_dimid <- names(dat_dn)
    subset. <- c(subset., list(...))
    ind <- subsetIndices(subset., which_dims., dat_d, dat_dn)
    which_dims. <- attr(ind, "which_dims.")
    #
    # do subsetting
    out <- do("[", dat, arg_list = c(ind, list(drop = FALSE)))
    #
    # drop singleton dimensions if requested
    if (!identical(drop., FALSE)) {
        out_d <- dim(out)
        out_dn <- dimnames(out)
        keep_dim <- out_d > 1L
        if (is.null(drop.)) {
            keep_dim[-which_dims.] <- TRUE
        } else if (is.logical(drop.)) {
            if (length(drop.) != 1L) {
                stop("If 'drop.' is logical, it must be a single value")
            }
        } else {
            assertVector(drop., strict = TRUE, any.missing = FALSE,
                         max.len = length(out_d), unique = TRUE,
                         .var.name = "drop.")
            if (is.character(drop.)) {
                drop. <- match(drop., dat_dimid)
                if (anyNA(drop.))
                    stop("Not all elements of 'drop.' correspond to existing dimensions in 'dat'")
                keep_dim[-drop.] <- TRUE
            } else if (is.numeric(drop.)) {
                assertIntegerish(drop., lower = 1L, upper = length(dim(dat)),
                                 any.missing = FALSE, .var.name = "drop.")
                keep_dim[-drop.] <- TRUE
            } else {
                stop("If 'drop.' is not NULL, TRUE or FALSE, it must be a character or numeric vector")
            }
        }
        if (sum(keep_dim) > 0L) {
            setattr(out, "dim", out_d[keep_dim])
            setattr(out, "dimnames", out_dn[keep_dim])
        } else {
            setattr(out, "dim", NULL)
            setattr(out, "dimnames", NULL)
        }
    }
    #
    # reattach attributes if requested
    if (keep_attributes.) {
        a <- attributes(dat)
        a2keep <- setdiff(names(a),
                          c("class", "comment", "dim", "dimnames", "names",
                            "row.names", "tsp"))
        a <- a[a2keep]
        for (i in a2keep) setattr(out, i, a[[i]])
        # this is just to handle the temporary "factor_level" attribute of
        # the eeg data
        if (!is.null(attr(out, "factors")) &&
                ("factor_level" %in% names(subset.)) ) {
            tempa <- splitMarker(dimnames(out)$factor_level,
                                 colnames(attr(out, "factors")),
                                 splitchar = "\\|")
            setattr(out, "factors", tempa)
        }
    }
    #
    # return
    out
}

#' @rdname subsetArray
#' @usage
#' subsetArray(dat, subset. = list(), which_dims. = NULL, drop. = NULL, ...) <- value
#' @param value a vector, matrix or array of the new values
#' @export
# Replace a part of an array
`subsetArray<-` <- function(dat, subset. = list(), which_dims. = NULL, 
                            drop. = NULL, ..., value) {
    #
    # check arguments (dat, subset., which_dims.) and find indices
    assertArray(dat, min.d = 1L, .var.name = "dat")
    dat_d <- dim(dat)
    dat_dn <- dimnames(dat)
    dat_dimid <- names(dat_dn)
    subset. <- c(subset., list(...))
    ind <- subsetIndices(subset., which_dims., dat_d, dat_dn)
    #
    # do subsetting
    value_dnn <- names(dimnames(value))
    if (!is.null(value_dnn) && !anyDuplicated(value_dnn) &&
        !anyNA(value_dnn)) {
        value <- apermArray(value, na.omit(match(dat_dimid, value_dnn)),
                            keep_attributes. = FALSE)
    }
    dat <- do("[<-", dat, arg_list = c(ind, list(value = value)))
    #
    # return
    invisible(dat)
}

#' Expand array
#'
#' \code{expandInto} expands an array to a larger array (either to an array
#' with more dimensions or an array with longer dimensions or both).
#' @param dat an array
#' @param new_dat an array to expand 'dat' into
#' @param expand_levels a list of vectors which define the expanding scheme
#' for each dimension of 'dat' or a named list of vectors where the names refer
#' to selected dimensions in 'dat' (in this case 'dat' must have named
#' dimnames). The length of each vector in 'expand_levels' must match the
#' corresponding dimension size in 'new_dat'. The vectors must contain either
#' numeric or character indices of the levels of the given dimension in 'dat'.
#' @param safe_mode a logical value whether the expansion of non-singleton
#' dimensions is not allowed if the corresponding vectors in 'expand_levels' are
#' not provided (default: TRUE). If 'safe_mode' is TRUE, and both 'dat' and
#' 'new_dat' has dimension names, non-expanded dimensions are checked if the
#' order of levels should be adjusted for the given dimension. See Examples.
#' @param fill a logical value if the 'new_dat' should be filled with the 
#' corresponding values in 'dat' (TRUE, the default). In this case the values 
#' of 'dat' are coerced to match the type of 'new_dat' and the returned array 
#' inherits all attributes of 'new_dat'. Otherwise, only the dimensions and
#' dimension names are preserved. 
#' @export
#' @examples
#' # load example data
#' data(erps)
#'
#' # -----
#' # solve the following task: find all data points for which the amplitudes
#' # in the "Fz" channel are negative, and return TRUE for all corresponding
#' # data points in the other channels as well
#' # -----
#'
#' # subset the data and return TRUE if the values are negative
#' x_Fz <- subsetArray(erps, chan = "Fz") < 0
#' str(x_Fz)
#'
#' # expand this array into the original array
#' result <- expandInto(x_Fz, erps)
#' str(result)
#'
#' # check on a random channel that the results are really fine
#' x_Cz <- subsetArray(result, chan = "Cz")
#' # -->
#' # all TRUEs in x_Fz are also TRUE in x_Cz, and vica versa
#' ( tab <- table(x_Fz, x_Cz) )
#' stopifnot(identical(sum(diag(tab)), length(x_Fz)))
#'
#' # -----
#' # the function is clever enough to reorder the levels for those
#' # dimensions as well, which should not be expanded, but the order
#' # of levels is different in 'new_dat'
#' # -----
#' # reorder the 'stimclass' dimension in the original ERP array
#' erps2 <- subsetArray(erps, stimclass = c("C", "B", "A"))
#'
#' # expand x_Fz again
#' result2 <- expandInto(x_Fz, erps2)
#'
#' # turn 'safe_mode' off
#' result2_notsafe <- expandInto(x_Fz, erps2, safe_mode = FALSE)
#'
#' # check the results -> result2 is fine
#' x_Cz_2 <- subsetArray(result2,
#'                       chan = "Cz",
#'                       stimclass = c("A", "B", "C"))
#' ( tab <- table(x_Fz, x_Cz_2) )
#' stopifnot(identical(sum(diag(tab)), length(x_Fz)))
#'
#' # check the results -> result2_notsafe is wrong
#' x_Cz_2w <- subsetArray(result2_notsafe,
#'                        chan = "Cz",
#'                        stimclass = c("A", "B", "C"))
#' ( tab <- table(x_Fz, x_Cz_2w) )
#' stopifnot(!identical(sum(diag(tab)), length(x_Fz)))
#'
#' # -----
#' # the safest way is to provide 'expand_levels' explicitly for all
#' # dimensions where the order or number of levels do not match;
#' # using this argument it is also possible to copy the values of a given
#' # level to an other one
#' # -----
#' # suppose we want stimclass C to be copied from stimclass B while expanding
#' # to all channels (on the original ERP array)
#' result <- expandInto(x_Fz, erps,
#'                      expand_levels = list(stimclass = c("A", "B", "B")))
#'
#' # check the results
#' x_Cz_B <- subsetArray(result,
#'                       chan = "Cz",
#'                       stimclass = "B")
#' x_Cz_C <- subsetArray(result,
#'                       chan = "Cz",
#'                       stimclass = "C")
#' # --> they are identical:
#' stopifnot(identical(x_Cz_B, x_Cz_C))
#' # --> compared to the results on Fz, stimclass B remained the same:
#' x_Fz_B <- subsetArray(x_Fz, stimclass = "B")
#' ( tab <- table(x_Fz_B, x_Cz_B) )
#' stopifnot(identical(sum(diag(tab)), length(x_Fz_B)))
#'
#' # -----
#' # it is possible that we want to expand an array into a larger array,
#' # but the types do not match; consider the 'fill' argument depending on
#' # your needs
#' # -----
#' # create a logical matrix
#' ( from_logical <- matrix(c(TRUE, FALSE), 2, 1,
#'                          dimnames = list(observation = c("a", "b"),
#'                                          measure = "width")) )
#' 
#' # it should be expanded to a larger, integer matrix with a special class
#' ( to_integer <- matrix(1:4, 2, 2,
#'                        dimnames = list(observation = c("a", "b"),
#'                                        measure = c("height", "width"))) )
#' class(to_integer) <- "mySpecialClass"
#' 
#' # perform to expansions
#' ( res_int <- expandInto(from_logical, to_integer) )
#' ( res_log <- expandInto(from_logical, to_integer, fill = FALSE) )
#' 
#' # res_int is integer, and preserves the class, res_log not
#' stopifnot(is.integer(res_int))
#' stopifnot(inherits(res_int, "mySpecialClass"))
#' stopifnot(is.logical(res_log))
#' stopifnot(!inherits(res_log, "mySpecialClass"))
#' 
#' # however, the dimnames are preserved
#' stopifnot(identical(
#'     dimnames(res_log),
#'     dimnames(to_integer)
#' ))
#' 
expandInto <- function(dat, new_dat, expand_levels = NULL, safe_mode = TRUE,
                       fill = TRUE) {
    # check arguments
    assertArray(dat, min.d = 1L, .var.name = "dat")
    orig_dim <- dim(dat)
    orig_dimn <- fillMissingDimnames(dimnames(dat), orig_dim)
    orig_dimid <- names(orig_dimn)
    assertArray(new_dat, min.d = length(orig_dim), .var.name = "new_dat")
    new_dim <- dim(new_dat)
    new_dimn <- fillMissingDimnames(dimnames(new_dat), new_dim)
    new_dimid <- names(new_dimn)
    if (!is.null(new_dimid) && any(!orig_dimid %in% new_dimid)) {
        stop(paste0("if 'new_dat' has dimension identifiers, ",
                    "it must contain all dimension identifiers of 'dat'"))
    }
    # match dimension order
    shared_dimid <- intersect(new_dimid, orig_dimid)
    dat <-
        if (identical(shared_dimid, orig_dimid)) {
            copy(dat)
        } else {
            apermArray(dat, shared_dimid)
        }
    # insert new dimensions
    orig_dimn2 <- setNames(rep(list("1"), length(new_dim)), new_dimid)
    orig_dimn2[orig_dimid] <- orig_dimn
    array_(dat, vapply(orig_dimn2, length, integer(1L)), orig_dimn2)
    # consider expand_levels
    if (!is.null(expand_levels)) {
        assertList(expand_levels, types = c("character", "numeric"),
                   any.missing = FALSE,
                   min.len = 1L, max.len = length(orig_dim),
                   .var.name = "expand_levels")
        if (is.null(names(expand_levels))) {
            if (length(expand_levels) != length(orig_dim)) {
                stop(paste0("if the length of 'expand_levels' does not ",
                            "match the number of dimension of 'dat', ",
                            "'expand_levels' must be named"))
            } else {
                names(expand_levels) <- orig_dimid
            }
        }
        expand_levels <- fillMissingDimnames(
            expand_levels, vapply(expand_levels, length, integer(1L)))
    }
    # create subsetting indices
    sub_indices <- mapply(
        function(old, new, nam) {
            if (safe_mode && is.null(expand_levels) &&
                length(old) > 1L && length(old) != length(new)) {
                stop(sprintf(
                    paste0("the '%s' dimension is not a singleton dimension ",
                           "in 'dat', but it should be expanded to match ",
                           "the corresponding dimension in 'new_dat'. ",
                           "Please provide 'expand_levels' or if you really ",
                           "know what you are doing, set 'safe_mode' to FALSE."
                           ), nam))
            }
            # return
            if (is.null(exp <- expand_levels[[nam]])) {
                if (safe_mode) {
                    if (all(new %in% old)) {
                        new
                    } else if (length(old) == 1L) {
                        rep_len(old, length(new))
                    } else {
                        stop(sprintf(
                            paste0("the '%s' dimension is neither a singleton ",
                                   "dimension in 'dat', nor does it contain ",
                                   "all levels of the corresponding dimension ",
                                   "in 'new_dat'. Please provide ",
                                   "'expand_levels' or if you really know ",
                                   "what you are doing, set 'safe_mode' to ",
                                   "FALSE."
                            ), nam))
                    }
                } else {
                    rep_len(old, length(new))
                }
            } else if (length(exp) == length(new)) {
                exp
            } else {
                stop(sprintf(
                    paste0("the size of the '%s' dimension in 'new_dat' ",
                           "and the length of the corresponding vector in ",
                           "'expand_levels' must be equal"), nam))
            }
        },
        orig_dimn2, new_dimn, new_dimid, SIMPLIFY = FALSE
    )
    # expand array and return
    if (fill) {
        new_dat[] <- subsetArray(dat, sub_indices)
        new_dat
    } else {
        out <- subsetArray(dat, sub_indices)
        setattr(out, "dim", dim(new_dat))
        setattr(out, "dimnames", dimnames(new_dat))
        out
    }
}


#' Reshape array to matrix with specified row dimension
#'
#' \code{array2mat} reshapes an array to a matrix
#' @param dat an array to reshape
#' @param row_dim name or index of dimension which should be the row dimension
#' of the returned matrix
#' @param return_attributes logical value; if TRUE (default), the attributes
#' of the array is saved as an additional attribute ("array_atributes") of
#' the returned matrix including the row_dim parameter
#' @param keep_dimnames logical value; if TRUE (default), the function not
#' only reshapes the array but retains the dimension names separated by "|"
#' @export
#' @return A matrix
#' @seealso \code{\link{mat2array}} which is the inverse of array2mat
array2mat <- function(dat, row_dim, return_attributes = TRUE,
                      keep_dimnames = TRUE) {
    if (is.character(row_dim)) {
        row_dim <- which(names(dimnames(dat)) == row_dim)
    }
    col_dims <- seq_along(dim(dat))[-row_dim]
    out <- apermArray(dat, c(row_dim, col_dims),
                      keep_attributes. = TRUE)
    if (!keep_dimnames || is.null(dimnames(dat))) {
        out_dimnames <- NULL
    } else {
        dat <- decorateDims(dat)
        out_dimnames <- c(dimnames(dat)[row_dim],
                          list(interaction(expand.grid(dimnames(dat)[col_dims]),
                                           sep = "|")))
        setattr(out_dimnames, "names",
                c(names(dimnames(dat))[row_dim],
                  paste(names(dimnames(dat))[col_dims], collapse = "|")))
    }
    matrix_(out, nrow(out), dimnames = out_dimnames)
    if (return_attributes) {
        setattr(out, "array_attributes",
                c(attributes(dat), list(row_dim = row_dim)))
    }
    # return
    out
}

#' (Back-)Transforms a matrix to an array. See array2mat as a related function.
#'
#' \code{mat2array} reshapes a matrix to an array. It is essentially the inverse
#' function of array2mat.
#' @param dat a matrix to reshape
#' @param dims a numeric vector indicating the dimensions of the resulting array
#' or a (potentially named) list indicating the dimnames attribute of the
#' resulting array. If NULL (default), dat is assumed to be the result of a call
#' to array2mat().
#' @param row_dim name or index of dims which corresponds to the row dimension
#' of dat. If NULL (default), dat is assumed to be the result of a call to
#' array2mat().
#' @export
#' @seealso \code{\link{array2mat}}
mat2array <- function(dat, dims = NULL, row_dim = NULL) {
    dimn <- NULL
    if (is.null(dims)) {
        array_attribs <- attr(dat, "array_attributes")
        if (!is.null(array_attribs)) {
            dims <- array_attribs$dim
            row_dim <- array_attribs$row_dim
            dimn <- array_attribs$dimnames
        } else {
            stop("Provide dims and row_dim parameters!")
        }
    } else {
        if (is.character(row_dim)) row_dim <- which(names(dims) == row_dim)
        if (is.list(dims)) {
            dimn <- dims
            dims <- vapply(dims, length, 0L)
        }
    }
    matdims <- c(dims[row_dim], dims[-row_dim])
    dimord <- order( c(row_dim, seq_along(dims)[-row_dim]) )
    out <- apermArray(array(dat, matdims), dimord,
                      keep_attributes. = TRUE)
    dimnames(out) <- dimn
    # return
    out
}

#' Transforms an array to a data.frame
#'
#' \code{array2df} transforms an array to a data.frame
#' @param dat a matrix or array
#' @param value_name the name of the variable in the data.frame which
#' contains the values of the input (default: "values"). If 'value_dim' is 
#' not NULL or \code{character(0L)}, the name of the value-columns will start 
#' with 'value_name'.
#' @param value_dim numeric or character index(es) of dimension(s) whose
#' levels should be handled as separate measures. They appear as separate
#' variables in the returned data.frame. The default is NULL, which means
#' that one variable will contain all data points of the array.
#' @param value_type a character value which specifies the type of the
#' measure (default: \code{typeof(dat)}). Possible types are "logical",
#' "character", "numeric", "double", "integer", "factor", "raw", "complex".
#' @param auto_convert a logical value whether automatic conversion of
#' dimension names (i.e., characters to numeric (if possible) or to factors)
#' should be performed (default: FALSE). 
#' @param na_omit if TRUE, omit all rows from the data.frame which have missing
#' values (default: FALSE)
#' @param ... named arguments passed to \code{\link{autoConvert}} if 
#' 'auto_convert' is TRUE
#' @details By default, this function returns a data.frame with as many 
#' variables as the number of dimensions of the array, coding the levels
#' of each dimension, plus an extra variable containing the values of the array.
#' However, 'dat' might be an array of different measures, e.g. one of the
#' dimensions might have two levels: 'weight' and 'height' which should 
#' appear as separate variables in the returned data.frame. If this is the case,
#' one should provide the name or index of the given dimension in the argument
#' 'value_dim'. If more dimensions are given in 'value_dim', all combinations 
#' of their levels will be returned as separate variables. The names of these
#' variables will contain the character string as given in 'value_name', plus
#' the concatenated names of the corresponding dimensions and dimension levels.
#' @export
#' @return The function returns a data.frame.
#' @seealso \code{\link{autoConvert}} for coercion, and 
#' \code{\link{transformArray}} for a higher-level version of \code{array2df}
#' @examples
#' # example data
#' data(erps)
#'
#' # transform to data.frame, change the default value name to "amplitudes"
#' str(array2df(erps, value_name = "amplitudes"))
#'
#' # treat all dimensions as factors, except for the time dimension, which
#' # should be integer, and return the three levels of the stimulus class
#' # dimension in wide format (as separate variables)
#' str(array2df(erps, value_name = "amplitudes", value_dim = "stimclass",
#'              auto_convert = TRUE))
array2df <- function(dat, value_name = "values", value_dim = NULL,
                     value_type = typeof(dat), 
                     auto_convert = FALSE,
                     na_omit = FALSE, ...) {
    # argument checks
    assertArray(dat, min.d = 1L, .var.name = "dat")
    assertString(value_name, .var.name = "value_name")
    assertChoice(value_type, 
                 c(typeof(dat), "logical", "character", "integer",
                   "numeric", "double", "factor", "complex", "raw"), 
                 .var.name = "value_type")
    assertFlag(auto_convert, .var.name  = "auto_convert")
    assertFlag(na_omit, .var.name = "na_omit")
    orig_address <- address(dat)
    # if value_dim is not NULL, dat must be permuted
    if (length(value_dim)) {
        if (is.character(value_dim)) {
            value_dim <- match(value_dim, names(dimnames(dat)), 
                                  nomatch = 0L)
        }
        if (!all(value_dim %in% seq_along(dim(dat)))) {
            stop(paste0(
                "array2df: not all dimensions in 'value_dim' ",
                "are present in 'dat'"), call. = FALSE)
        }
        value_dim <- sort(value_dim)
        # permute dat
        dimord <- c(seq_along(dim(dat))[-value_dim], value_dim)
        # modify select (a possible argument to autoConvert)
        if (exists("select", inherits = FALSE) && is.numeric(select)) 
            select <- match(select, dimord)
        dat <- apermArray(dat, dimord)
    }
    row_dim <- seq_len(length(dim(dat)) - length(value_dim))
    var_dim <- seq_along(dim(dat))[-row_dim]
    # coerce dat if needed
    if (value_type != typeof(dat)) {
        dat <- do(paste0("as.", value_type), dat)
    }
    # prepare dimension names
    dimn <- fillMissingDimnames(dimnames(dat), dim(dat), .dimnames = FALSE)
    if (auto_convert) dimn <- autoConvert(dimn, ...)
    dimn_null <- vapply(dimn, is.null, logical(1L))
    dimn[dimn_null] <- lapply(dim(dat)[dimn_null], seq_len)
    # prepare the data.frame
    out <- expand.grid(dimn[row_dim], KEEP.OUT.ATTRS = FALSE, 
                       stringsAsFactors = FALSE)
    # cbind the values
    if (is.null(value_dim)) {
        out[[value_name]] <- as.vector(dat)
    } else {
        if (identical(orig_address, address(dat))) {
            dim(dat) <- c(nrow(out), length(dat)/nrow(out))
        } else {
            matrix_(dat, nrow = nrow(out))
        }
        tempn <- 
            if (length(var_dim)) {
                do.call(
                    "paste",
                    c(expand.grid(mapply(function(x, y) paste(x, y, sep = "_"),
                                         names(dimn)[var_dim], dimn[var_dim],
                                         SIMPLIFY = FALSE), 
                                  KEEP.OUT.ATTRS = FALSE, 
                                  stringsAsFactors = FALSE),
                      sep = "."))
            } else {
                NULL
            }
        colnames(dat) <- 
            if (length(value_name) & length(tempn)) {
                paste(value_name, tempn, sep = ".")
            } else if (length(tempn)) {
                tempn
            } else if (length(value_name)) {
                value_name
            } else {
                "values"
            }
        out <- cbind(out, dat)
    }
    # remove rows with NA if requested
    if (na_omit) out <- na.omit(out)
    # return
    out
}

#' Reshape array by merging specific dimensions
#'
#' \code{mergeDims} reshapes an array by merging user-specified dimensions
#' @param dat array to reshape
#' @param dims list of dimension names or dimension indices to merge. Can be a vector if
#' only one group of dimensions should be merged. If any of the list elements is a
#' character vector, dat must have named dimnames.
#' @param return_attributes logical value (default: TRUE) whether original
#' attributes should be appended to the resulting array
#' @param keep_dimnames logical; if TRUE (default), dimension names are also
#' merged and attached as dimnames attribute
#' @param sep character which separates dimension names after merging (default ".")
#' @export
#' @return Array with merged dimensions
mergeDims <- function(dat, dims, return_attributes = TRUE,
                      keep_dimnames = TRUE, sep = ".") {
    if (return_attributes) {
        attribs <- attributes(dat)
        dimattribs <- list(dim = setNames(attribs$dim, names(attribs$dimnames)),
                           dimnames = attribs$dimnames)
        attribs <- attribs[setdiff(names(attribs), names(dimattribs))]
    }
    if (!is.list(dims)) dims <- list(dims)
    dims <- lapply(dims,
                   function(x) {
                       if (is.character(x)) {
                           match(x, names(dimnames(dat)))
                       } else {
                           x
                       }
                   })
    if (anyDuplicated(unlist(dims)))
        stop("Duplicated dimensions in the dims parameter are not allowed")
    if (length(unlist(dims)) < length(dim(dat))) {
        dims <- c(dims, as.list(setdiff(seq_along(dim(dat)), unlist(dims))))
    }
    if (keep_dimnames) {
        dat <- decorateDims(dat)
        dimn <- lapply(dims, function(i)
            do.call(paste,
                    c(expand.grid(dimnames(dat)[i],
                                  stringsAsFactors = FALSE),
                      list(sep = sep))))
        setattr(dimn, "names",
                vapply(relist(names(dimnames(dat))[unlist(dims)], dims),
                       paste, "", collapse = sep))
    } else {
        dimn <- NULL
    }
    udims <- unlist(dims, use.names = FALSE)
    out <-
        if (identical(udims, seq_along(dim(dat)))) {
            array_(as.vector(dat),
                   vapply(relist(dim(dat)[udims], dims), prod, 0),
                   dimn)
        } else {
            array_(aperm(dat, udims),
                   vapply(relist(dim(dat)[udims], dims), prod, 0),
                   dimn)
        }
    if (return_attributes) {
        setattr(out, "orig_dimattributes",
                c(dimattribs, list(sep = sep, merged_dims = dims)))
        if (length(attribs) > 0) {
            for (i in names(attribs)) setattr(out, i, attribs[[i]])
        }
    }
    # return
    out
}

#' Reverse mergeDims transformation
#'
#' \code{revMergeDims} sets back the original array transformed by
#' \code{\link{mergeDims}}
#' @param dat numeric matrix or array with merged dimensions
#' @export
#' @return An array of the same dimension attributes as the array which
#' \code{\link{mergeDims}} was called on
revMergeDims <- function(dat) {
    attribs <- attributes(dat)
    if (!"orig_dimattributes" %in% names(attribs))
        stop("The input is not a result of mergeDims(..., return_attributes=TRUE)")
    dimattribs <- attr(dat, "orig_dimattributes")
    attribs <- attribs[setdiff(names(attribs),
                               c("orig_dimattributes", names(dimattribs)))]
    orig_dim <- dimattribs$dim
    orig_dimnames <- dimattribs$dimnames
    merged_dims <- unlist(dimattribs$merged_dims, use.names = FALSE)
    if (is.character(merged_dims)) {
        merged_dims <- match(merged_dims, names(orig_dimnames))
    }
    dimord <- c(merged_dims,
                setdiff(seq_along(orig_dim), merged_dims))
    out <- array(dat, orig_dim[dimord])
    out <- aperm(out, order(dimord))
    dimnames(out) <- orig_dimnames
    if (length(attribs) > 0) {
        for (i in names(attribs)) setattr(out, i, attribs[[i]])
    }
    # return
    out
}


#' Transform a specific dimension of an array into a multidimensional array
#'
#' \code{dim2multidim} transforms a user-specified dimension of an array into a
#' multidimensional array while keeping other dimensions intact. It can be
#' conceived of as the inverse of expand.grid on the specified dimension.
#' @param dat is the input matrix or array
#' @param whichdim is the target dimension to be expanded (can be a numerical
#' index, a character string specifying the dimension name or a logical vector).
#' If whichdim is a string, dat must have named dimnames.
#' @param datfr a data.frame or matrix with the combination of factor levels
#' (e.g., the result of a call to expand.grid)
#' @export
#' @return An array
#' @seealso \code{\link{expand.grid}}
dim2multidim <- function(dat, whichdim, datfr) {
    if (!is.data.frame(datfr))
        datfr <- as.data.frame(datfr, stringsAsFactors = FALSE)
    datfr <- lapply(datfr, as.character)
    orig_dimnames <- dimnames(dat)
    orig_dims <- dim(dat)
    add_dimnames <- lapply(datfr, unique)
    if (!identical(datfr,
                   as.list(expand.grid(add_dimnames,
                                       KEEP.OUT.ATTRS = FALSE,
                                       stringsAsFactors = FALSE)))) {
        stop("datfr is not commensurate with the result of expand.grid()")
    }
    add_dims <- vapply(add_dimnames, length, 0L)
    if (is.character(whichdim))
        whichdim <- which(names(orig_dimnames) == whichdim)
    if (is.logical(whichdim))
        whichdim <- which(whichdim)
    dim(dat) <- append(orig_dims[-whichdim], add_dims, whichdim - 1)
    dimnames(dat) <- append(orig_dimnames[-whichdim], add_dimnames,
                            whichdim - 1)
    # return
    dat
}

#' Splits an array along a given dimension
#'
#' \code{splitArray} splits an array along given dimension(s) into a list of
#' sub-arrays
#' @param dat numeric array (preferably with named dimnames)
#' @param whichdim numeric or character vector, the dimension(s) of the array
#' to split along
#' @param f a list of factors in the sense that
#' \code{lapply(f, as.factor)} defines the grouping to split along for
#' each dimensions in \code{whichdim}. If \code{NULL}, \code{splitArray}
#' splits all levels of the \code{whichdim} dimensions. If not \code{NULL},
#' the length of \code{f} must match the length of \code{whichdim}. If \code{f}
#' is a named list, the names are reflected in the \code{dimnames} attribute of
#' the resulting list (see Value section).
#' @param drop logical; should singleton dimensions (dimensions with only
#' one level) be deleted (TRUE) or not (FALSE, the default)
#' @export
#' @return A list of subsets of the original data matrix/array with \code{dim}
#' and \code{dimnames} attributes. The dimensions of the list correspond to the
#' length of each element in \code{f} (after replacing NULL values with correct
#' vectors).
#' @examples
#' # load example data
#' data(erps)
#' 
#' # get the reading group membership of the subjects
#' dat_id <- attr(erps, "id")
#' 
#' # split on the basis of the reading group membership
#' groups <- splitArray(erps, "id", list(readgroup = dat_id$group))
#' 
#' # check
#' str(groups)
#' 
#' \dontshow{
#' stopifnot(identical(groups$control, 
#'                     subsetArray(erps, id = dat_id$group == "control",
#'                                 keep_attributes. = FALSE)))
#' stopifnot(identical(groups$dl, 
#'                     subsetArray(erps, id = dat_id$group == "dl",
#'                                 keep_attributes. = FALSE)))
#' stopifnot(identical(names(groups), c("control", "dl")))
#' stopifnot(identical(dim(groups), 2L))
#' stopifnot(identical(dimnames(groups), list(readgroup = c("control", "dl"))))
#' }
#' 
splitArray <- function(dat, whichdim, f = NULL, drop = FALSE) {
    subFn <- function(ind) {
        abind::asub(dat, lapply(ind, unlist), whichdim_num, drop = drop)
    }
    if (is.character(whichdim)) {
        whichdim_num <- match(whichdim, names(dimnames(dat)))
        if (anyNA(whichdim_num))
            stop("splitArray: wrong dimension name(s) provided", call. = FALSE)
    } else {
        whichdim_num <- whichdim
    }
    assertArray(dat, min.d = 1L, .var.name = "dat")
    dimn <- fillMissingDimnames(dimnames(dat), dim(dat))
    if (is.null(f)) {
        f <- dimn[whichdim_num]
    } else if (!is.list(f)) {
        f <- list(f)
        names(f) <- names(dimn)[whichdim_num[1]]
    }
    if (is.character(whichdim)) {
        ind <- names(f) == "" & vapply(f, is.null, logical(1L))
        names(f)[ind] <- whichdim[ind]
    }
    if (length(f) != length(whichdim_num)) {
        stop("splitArray: length of f must match the length of whichdim", 
             call. = FALSE)
    }
    for (i in seq_along(f)) {
        if (is.null(f[[i]])) {
            #f[[i]] <- seq_len(dim(dat)[whichdim_num[i]])
            f[[i]] <- dimn[[whichdim_num[i]]]
        } else if (is.list(f[[i]])) {
            f[[i]] <- do.call(paste, list(f[[i]], sep="_"))
        }
    }
    f <- lapply(f, function(x) split(seq_along(x), x))
    out.dimnames <- lapply(f, names)
    out.dim <- vapply(out.dimnames, length, integer(1L), USE.NAMES = FALSE)
    f <- expand.grid(f, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    out <- lapply(1:nrow(f), function(i) subFn(f[i, ]))
    setattr(out, "names",
            do.call("paste", c(lapply(f, names), list(sep = "."))))
    setattr(out, "dim", out.dim)
    setattr(out, "dimnames", out.dimnames)
    # return
    out
}

#' Combine several arrays into one large array
#'
#' \code{bindArrays} is just a wrapper around \code{\link[abind]{abind}}.
#' @details This function calls \code{\link[abind]{abind}} and adds
#' the names of the dimension names of the arrays. If the inputs have named
#' dimension names, their dimensions are permuted before feeding to \code{abind}.
#' \code{bindArrays} has the same arguments as \code{abind} except for
#' \code{along_name}, as documented below.
#' @param ... Any number of vectors, matrices, arrays, or data frames. If the
#' objects have named dimension names, their dimensions are permuted before
#' feeding to \code{abind}. Otherwise, the dimensions of all the arrays must
#' match, except on one dimension (specified by along=). If these arguments are
#' named, the name will be used for the name of the dimension along which the
#' arrays are joined. Vectors are treated as having a dim attribute of length
#' one.
#' Alternatively, there can be one (and only one) list argument supplied, whose
#' components are the objects to be bound together. Names of the list components
#' are treated in the same way as argument names.
#' @param along_name a character version of \code{along} in
#' \code{abind}. Only considered if input arrays have named
#' dimension names, and if conflicts with \code{along}, \code{along_name}
#' overrides \code{along}. If \code{along_name} is not among the dimension
#' names, it will be the name of the new dimension name. In this case the
#' position of the new dimension can be controlled by \code{along}.
#' @inheritParams abind::abind
#' @export
#' @seealso \code{\link[abind]{abind}} for the original version in package
#' \pkg{abind}; \code{\link{mergeArrays}} if you have a list of arrays created
#' by \code{\link{splitArray}} or you want to bind on multiple dimensions;
#' \code{\link{rearrangeList}} if you want to bind arrays in two-level lists
bindArrays <- function(..., along = NULL, rev.along = NULL, new.names = NULL,
                       force.array = TRUE, make.names = use.anon.names,
                       use.anon.names = FALSE,
                       use.first.dimnames = FALSE, hier.names = FALSE,
                       along_name = NULL)  {
    dat <- list(...)
    if (length(dat) > 1L && any(vapply(dat, is.list, logical(1L)))) {
        stop("Only one list is allowed as an argument")
    }
    if (is.list(dat[[1L]])) dat <- dat[[1L]]
    if (is.character(along)) stop("Provide along_name if you want to refer to a dimension by its name")
    #
    dimn <- lapply(dat, function(x) names(dimnames(x)))
    if (length(udimn <- unique(dimn)) > 1L) {
        if (is.null(along_name)) {
            stop("The order of dimensions differs. Provide along_name instead of along.")
        }
        alldimn <- unique(unlist(dimn, use.names = FALSE))
        dat <- lapply(dat, function(x) {
            if (length(dim(x)) < length(alldimn)) {
                if (identical(names(dimnames(x)),
                              setdiff(alldimn, along_name))) {
                    x
                } else {
                    aperm(x, setdiff(alldimn, along_name))
                }
            } else {
                if (identical(names(dimnames(x)), alldimn)) {
                    x
                } else {
                    aperm(x, alldimn)
                }
            }
        })
    }
    new_ndimn <- udimn[[1L]]
    N <- max(1L, vapply(dat, function(x) length(dim(x)), integer(1L)))
    if (is.null(along))
        along <- N
    if (!is.null(rev.along))
        along <- N + 1L - rev.along
    if (!is.null(along_name)) {
        ind <- match(along_name, new_ndimn)
        if (!is.na(ind)) {
            along <- ind
        } else {
            if (along >= 1L && along <= length(new_ndimn))
                along <- length(new_ndimn) + 1L
            new_ndimn <- append(new_ndimn, along_name, along)
        }
    } else if (along < 1) {
        new_ndimn <- c("", new_ndimn)
    } else if (along > N) {
        new_ndimn <- c(new_ndimn, "")
    }
    out <- abind(dat, along = along, rev.along = NULL, new.names = new.names,
                 force.array = force.array, make.names = make.names,
                 use.anon.names = use.anon.names,
                 use.first.dimnames = use.first.dimnames,
                 hier.names = hier.names)
    names(dimnames(out)) <- new_ndimn
    # return
    out
}

#' Merge arrays having common dimension identifiers
#'
#' \code{mergeArrays} merges multiple arrays or a list of arrays into one large
#' array. It can be regarded as the inverse of \code{\link{splitArray}}.
#' @param ... numeric arrays with named dimension names or a list of such
#' arrays. All arrays must have identically named dimensions, but the order of
#' dimensions does not need to be identical. Duplicated dimension levels are not
#' allowed.
#' @param base_value while setting up the resulting array, what value should be
#' given as default (e.g. NA, 0, "", etc.)
#' @param sort_dims logical. If FALSE (default), the order of dimensions follows
#' the first array's dimension order; if TRUE, lexical sorting is applied.
#' @param sort_dimlevels logical; should dimension levels be sorted for each
#' dimension (default: FALSE)
#' @export
#' @return The resulting array has identical dimension identifiers as the input
#' arrays, and for each dimension, as many dimension levels as the union of
#' the dimension levels of the input arrays.
mergeArrays <- function(..., base_value = NA,
                        sort_dims = FALSE, sort_dimlevels = FALSE) {
    is_bad <- function(x) {
        anyNA(x) || is.null(x) || any(x == "") || anyDuplicated(x)
    }
    dimnCheck <- function(x) {
        if (is_bad(names(x))) return(FALSE)
        if (any(vapply(x, is_bad, logical(1L)))) return(FALSE)
        TRUE
    }
    #
    dat <- list(...)
    if (length(dat) > 1L && any(vapply(dat, is.list, logical(1L)))) {
        stop("mergeArrays: only one list is allowed as an argument", 
             call. = FALSE)
    }
    if (is.list(dat[[1L]])) dat <- dat[[1L]]
    #
    dimn <- lapply(dat, function(x) dimnames(x))
    if (any(checks <- !vapply(dimn, dimnCheck, logical(1L)))) {
        bad <- which(checks)
        stop(paste0(
            "mergeArrays: dimension names of the arrays ",
            paste(bad, collapse = ", "),
            " are not appropriate (see Arguments in help('mergeArrays')"), 
            call. = FALSE)
    }
    ndimn <- lapply(dimn, names)
    if (length(unique(lapply(ndimn, sort))) > 1L)
        stop(paste0(
            "mergeArrays: at least one list element has a ",
            "unique dimension identifier", 
             call. = FALSE))
    all_dimn <- if (sort_dims) dimn[[1L]][order(ndimn[[1L]])] else dimn[[1L]]
    for (i in names(all_dimn)) {
        for (j in dimn[-1L]) {
            all_dimn[[i]] <- union(all_dimn[[i]], j[[i]])
        }
        if (sort_dimlevels) all_dimn[[i]] <- sort(all_dimn[[i]])
    }
    out <- base_value[1L]
    storage.mode(out) <- typeof(dat[[1L]])
    out <- array_(out, vapply(all_dimn, length, integer(1L)),
                  all_dimn)
    out_touched <- array_(FALSE, dim(out), dimnames(out))
    for (i in seq_along(dat)) {
        x <- dat[[i]]
        if (any(subsetArray(out_touched, dimnames(x)))) {
            stop(paste0(
                "mergeArrays: ",
                "the ", i, ". list element has overlapping dimension ",
                "combination(s) with at least one of the previous list ",
                "elements"), call. = FALSE)
        }
        subsetArray(out_touched, dimnames(x)) <- TRUE
        subsetArray(out, dimnames(x)) <- x
    }
    # return
    out
}


#' Rearrange two-level list
#'
#' \code{rearrangeList} reshapes a special type of one- or two-level lists.
#' @param dat the list to be rearranged
#' @param name_listdim character string; the name of the dimension which the
#' list represents
#' @param name_datadim character vector or a list of character vectors (for
#' two-level list input), providing the name of the dimensions for each list
#' element. The default is NULL, meaning that the original dimension
#' names will be used.
#' @details One or two level lists which contain identically shaped elements at
#' the base level, and elements at the base level are vectors, matrices or
#' arrays, can be rearranged to a matrix/array (from a one-level list) or to a
#' one-level list (from a two-level list) by binding identical elements.
#' This way the second level of the list will be represented as an additional
#' dimension in the ground level elements (vectors will be matrices, matrices
#' will become arrays, arrays recieve an extra dimension)
#' @export
#' @return A matrix or array, if the input is a one-level list, and a one-level
#' list, if the input is a two-level list
rearrangeList <- function(dat, name_listdim, name_datadim = NULL) {
    stopifnot(!missing(name_listdim))
    if (identical(unlist(dat, recursive = FALSE, use.names = FALSE),
                  unlist(dat, recursive = TRUE, use.names = FALSE))) {
        dat <- lapply(dat, as.array)
        newdims <- c(dim(dat[[1]]), length(dat))
        newdimns <- c(dimnames(dat[[1]]), listS(.name_listdim = names(dat)))
        if (!is.null(name_datadim)) {
            names(newdimns) <- c(name_datadim, name_listdim)
        }
        dat <- unlist(dat, recursive = TRUE, use.names = FALSE)
        dim(dat) <- newdims
        dimnames(dat) <- newdimns
    } else {
        names_at_level2 <- names(dat[[1]])
        dat <- lapply(dat, function(x) lapply(x, as.array))
        dat <- lapply(seq_along(names_at_level2), function(i) {
            out <- abind(lapply(dat, function(x) x[[i]]),
                         along = length(dim(dat[[1]][[i]])) + 1)
            dimn <-
                if (!is.null(name_datadim[[i]])) {
                    name_datadim[[i]]
                } else {
                    names(dimnames(dat[[1]][[i]]))
                }
            if (is.null(dimn)) {
                dimn <- character(length(dim(dat[[1]][[i]])))
            }
            names(dimnames(out)) <- c(dimn, name_listdim)
            return( out )
        })
        names(dat) <- names_at_level2
    }
    # return
    dat
}

#' Data preparation mainly aimed at facilitating plotting in lattice or ggplot2
#'
#' \code{prepare2plot} is deprecated. Use \code{\link{transformArray}} instead.
#' @param dat an array of ERP data. Must have named dimnames, one of which must
#' be id (corresponding to participants' identification codes)
#' @param datid data.frame consisting of identification codes ("id") and
#' subject-level factors
#' @param bwFac named list of between-subject factors (for splitting ERP data)
#' @param wiFac named list of within-subject factors (for subsetting ERP data)
#' @param collFac character vector of dimension names: average the ERP data
#' across these dimensions
#' @param diffFac a character vector indicating which levels of which dimensions
#' should be subtracted. The first element is the dimension name, the 2nd and 3rd
#' are the levels to be compared (2nd-3rd), and the 4th element corresponds to
#' the label of the new level.
#' @param compGFP a logical scalar; if TRUE (default), Global Field Power is
#' also computed.
#' @param keep_channels a logical scalar; if Global Field Power is requested,
#' shall individual channels be included in the result (default: FALSE, if GFP
#' is requested, ignored otherwise)
#' @param sc if not NULL, must be a named list indicating how scaling by GFP
#' should be done - either individually (default) or group-based, and either
#' timepoint-by-timepoint or for averaged segments
#' @param datfr a logical scalar (default: TRUE) determining if the resulting
#' array shall be transformed to a data.frame
#' @param iaFac a character vector indicating which dimensions should
#' be combined (ia is the abbreviation of interaction), if datfr is TRUE
#' @param ... additional parameters to be passed to \code{\link{array2df}}
#' @export
#' @return A data.frame if datfr is TRUE, and an array if datfr is FALSE
prepare2plot <- function(dat, datid,
                         bwFac = NULL, wiFac = NULL,
                         collFac = NULL, diffFac = NULL,
                         compGFP = TRUE, keep_channels = !compGFP,
                         sc = NULL,
                         datfr = TRUE, iaFac = NULL, ...) {
    #
    .Deprecated("transformArray")
    #
    bool_collFac <- !is.null(collFac) && collFac %in% names(dimnames(dat))
    if (!is.null(sc)) {
        sc.options.default <- list(method = "phase", indiv = TRUE, div = NULL)
        if (with(sc, exists("sc.options", inherits = FALSE))) {
            sc.options <- lapply(1:length(sc.options.default), function(i) {
                tt <- sc$sc.options[[names(sc.options.default)[i]]]
                if (is.null(tt)) sc.options.default[[i]] else tt
            })
            names(sc.options) <- names(sc.options.default)
            sc$sc.options <- NULL
        } else {
            sc.options <- sc.options.default
        }
    }
    facs <- nofac_facs <- data.frame(id = "all")
    if (!is.null(bwFac)) {
        facs <- do.call(expand.grid, bwFac)
        if (nrow(facs) == 0) facs <- nofac_facs
    }
    out <- lapply(1:nrow(facs), function(i) {
        if (identical(nofac_facs, facs)) {
            ids <- datid$id
        } else {
            ids <- datid$id[rowSums(sapply(1:length(bwFac), function(f)
                datid[, names(bwFac)[f]] == facs[i, f])) == length(bwFac)]
        }
        subFac <- append(wiFac, list(id = ids))
        temp <- subsetArray(dat, subFac)
        if (!is.null(sc) && sc.options$method != "fix" && sc.options$indiv == TRUE) {
            if (bool_collFac) {
                tempcoll <- setdiff(collFac, c("chan", "id", "time"))
                temp <- avgDims(temp, tempcoll)
                collFac <- setdiff(collFac, tempcoll)
            }
            tt0 <- subsetArray(temp, sc)
            tempcoll <- setdiff(names(dimnames(tt0)), c("chan", "id", "time"))
            if (length(tempcoll) > 0) tt0 <- avgDims(tt0, tempcoll)
            tt0 <- compGfp(tt0)
            if (sc.options$method == "phase") {
                tt0 <- apply(tt0, which(names(dimnames(tt0)) == "id"),
                             mean, na.rm = TRUE)
                temp <- sweep(temp, which(names(dimnames(temp)) == "id"),
                              tt0, "/")
            } else {
                temp <- sweep(temp,
                              match(names(dimnames(tt0)), names(dimnames(temp))),
                              tt0, "/")
            }
            if (bool_collFac) temp <- avgDims(temp, collFac)
        } else {
            if (bool_collFac) temp <- avgDims(temp, collFac)
        }
        if (!is.null(sc) && sc.options$method != "fix" && sc.options$indiv == FALSE) {
            tt0 <- subsetArray(temp, sc)
            tempcoll <- setdiff(names(dimnames(tt0)), c("chan", "id", "time"))
            if (length(tempcoll) > 0) tt0 <- avgDims(tt0, tempcoll)
            tt0 <- compGfp(tt0)
            if (sc.options$method == "phase") {
                tt0 <- mean(tt0)
            }
            if (is.vector(tt0)) {
                temp <- temp / tt0
            } else {
                temp <- sweep(temp,
                              match(names(dimnames(tt0)), names(dimnames(temp))),
                              tt0, "/")
            }
        }
        if (!is.null(sc) && sc.options$method == "fix") {
            temp <- temp/sc.options$div[i]
        }
        if (!is.null(diffFac)) {
            if (!is.list(diffFac)) {
                diffFac <- matrix(diffFac, 1, 4)
            } else {
                diffFac <- matrix(unlist(diffFac), length(diffFac), 4, TRUE)
            }
            for (ii in 1:nrow(diffFac)) {
                c1 <- parse(text = paste("list(", diffFac[ii,1], "='",
                                         diffFac[ii,2], "')", sep = ""))
                c2 <- parse(text = paste("list(", diffFac[ii,1], "='",
                                         diffFac[ii,3], "')", sep = ""))
                tempd <- subsetArray(temp, eval(c1)) -
                    subsetArray(temp, eval(c2))
                dimn.orig <- dimnames(temp)
                temp <- aperm(temp, c(setdiff(names(dimn.orig),
                                              diffFac[ii,1]),
                                      diffFac[ii,1]))
                dimn.perm <- dimnames(temp)
                temp <- c(temp,tempd)
                dimn.perm[[diffFac[ii, 1]]] <- dimn.orig[[diffFac[ii, 1]]] <-
                    c(dimn.orig[[diffFac[ii, 1]]], diffFac[ii, 4])
                array_(temp, vapply(dimn.perm, length, 0L), dimn.perm)
                temp <- aperm(temp, names(dimn.orig))
            }
        }
        if (compGFP) {
            temp <- compGfp(temp, keep_channels = keep_channels)
            if (datfr) {
                temp <- array2df(temp, ...)
                temp[colnames(facs)] <- facs[i, ]
            }
        } else if (datfr) {
            temp <- array2df(temp, ...)
            temp[colnames(facs)] <- facs[i, ]
        }
        return(temp)
    })
    if (datfr) {
        out <- do.call(rbind, out)
        if (!is.null(iaFac)) {
            out[, paste(iaFac, collapse = ".")] <- interaction(out[, iaFac])
        }
    } else {
        names.mat <- sapply(1:ncol(facs),
                            function(cc) sapply(1:nrow(facs), function(rr)
                                paste(colnames(facs)[cc],
                                      facs[rr, cc], sep = "-")))
        if (nrow(facs) > 1) {
            names(out) <- apply(names.mat, 1, paste, collapse = "_")
        } else if (ncol(facs) > 1) {
            names(out) <- paste(names.mat, collapse = "_")
        } else {
            names(out) <- names.mat
        }
    }
    # return
    out
}

#' Data preparation mainly aimed at facilitating plotting in lattice or ggplot2
#'
#' \code{transformArray} provides several options to transform an array to a
#' a data.frame which enables direct plotting in lattice or ggplot2
#' afterwards. It can also be used for analyses purposes without data.frame
#' conversion if compact code is desirable.
#' @param formula an object of class "formula" (or one that can be coerced to
#' that class, e.g., a character string): a symbolic description of the 
#' transformation steps before converting the \code{array} to a 
#' \code{data.frame}. See Details.
#' @param data a matrix or an array. Must have named dimnames.
#' @param group a list of grouping factors in the order of appearance in the
#' transformation formula (see 'Details'). If a named list is provided, those
#' names are used as the names of dimnames for the given grouping dimensions.
#' It can be a simple vector if there is only one splitting factor.
#' @param group_fun a function (or symbol or character string 
#' naming a function) which should be performed on the groups (i.e., on the 
#' list of arrays after splitting). \code{group_fun} must be a function which 
#' expects an array and a vector of dimension names as input and returns an 
#' array (or vector). Defaults to \code{avgDims}, which collapses (averages 
#' over) the grouping dimensions.
#' @param subset a list of subsetting vectors on the input array passed to
#' \code{\link{subsetArray}} before any transformation steps
#' @param datfr a logical value (default: TRUE) if the resulting array shall be
#' transformed to a data.frame
#' @param auto_convert a logical value whether automatic conversion of
#' dimension names (i.e., characters to numeric (if possible) or to factors)
#' should be performed (default: TRUE). Set to FALSE and call 
#' \code{\link{autoConvert}} directly on the returned data frame if you need
#' more control.
#' @param ... additional parameters to be passed to \code{\link{array2df}}
#' @details The formula interface of \code{transformArray} shall be given in
#' the form of\cr
#' 
#' \code{fun(y[d1, d2], fun_args) ~ . - d3 | d4 + d5}\cr
#' 
#' where
#' \describe{
#' \item{\code{fun }}{optional; an arbitrary function whose first argument is 
#' the data, and returns an array (optional)}
#' \item{\code{y }}{the name of the variable which holds the values in the 
#' returned data.frame (if 'datfr' is TRUE). See also \code{\link{array2df}}.}
#' \item{\code{[d1,d2] }}{optional; the dimensions of the data array whose 
#' levels should be treated as separate value-variables in the returned 
#' data.frame should be listed between squared brackets after the general 
#' name of the value variable. If you do not want to have a general name, 
#' place a dot (\code{.}) before the brackets. See also 
#' \code{\link{array2df}}.}
#' \item{\code{. }}{a dot on the right-hand side [RHS] of the formula means
#' 'all dimensions of the data array which are not explicitly mentioned in
#' the formula'. The dimension names can be explicitly provided as well, 
#' separated by \code{+}.}
#' \item{\code{d3 }}{optional; any dimension of the data array which is 
#' preceeded by a minus sign or any dimension which is not present in the 
#' formula will be collapsed (averaged over)}
#' \item{\code{d4,d5 }}{dimensions after the \code{|} sign are treated as
#' conditioning (grouping) dimensions, and shall be separated by \code{+} or 
#' \code{*}.}
#' }
#' \code{transformArray} performs the following actions:
#' \enumerate{
#' \item Takes the input array ('data') and subsets it if 'subset' is not NULL
#' or an empty list.
#' \item Calls \code{\link{avgDims}} on the (subsetted) data, and collapses
#' over all dimensions which are preceeded by \code{-} in the formula or are
#' not present in any other part of the formula.
#' \item Calls \code{\link{splitArray}} on the averaged data with the
#' conditioning dimensions in the formula. The 'group' argument is 
#' passed to the \code{\link{splitArray}} as the grouping argument ('f' in 
#' \code{splitArray}). For each data array which is returned after splitting,
#' \code{group_fun} is called with the character vector of the grouping
#' dimension names as its second argument. The resulting arrays are merged 
#' back to form one array.
#' \item If the left-hand side of the formula contains a function (see 
#' \code{fun} above), this function is called on the merged array with its
#' arguments as given in \code{fun_args}.
#' \item If 'datfr' is TRUE (the default), the resulting array is transformed
#' to a data.frame by calling \code{array2df}.
#' }
#' @export
#' @return The function returns a data.frame if 'datfr' is TRUE, and an array 
#' if 'datfr' is FALSE.
#' @examples
#' # example dataset
#' data(erps)
#' dat_id <- attr(erps, "id") # to get reading group memberships
#'
#' # compute simple grand averages (collapse over the 'id' dimension) and
#' # return it as a data.frame
#' DF <- transformArray(~ . - id, erps)
#' head(DF, 10)
#'
#' # compute the grand averages for each level of pairtype in each channel and
#' # time points; return the amplitudes of pairtype as separate variables in 
#' # a data.frame
#' DF <- transformArray(ampl[pairtype] ~ time + chan, erps)
#' head(DF, 10)
#'
#' # compute the grand averages of dyslexic and control subjects, and also 
#' # compute the Global Field Power (and transform to data.frame)
#' res1 <- transformArray(compGfp(ampl, keep_channels = TRUE) ~ . | id,
#'                        erps, list(readgroup = dat_id$group))
#'
#' # the same with much more typing, and it would be even longer to make
#' # it safer (e.g., match the order of dimensions, handle more grouping
#' # dimensions, etc.)
#' res2 <- splitArray(erps, "id", list(readgroup = dat_id$group))
#' res2[] <- lapply(res2, avgDims, "id")
#' res2 <- bindArrays(res2, along_name = "readgroup")
#' res2 <- compGfp(res2, keep_channels = TRUE)
#' res2 <- array2df(res2, value_name = "ampl", auto_convert = TRUE)
#' stopifnot(identical(res1, res2))
transformArray <- function(formula, data, group = NULL, group_fun = "avgDims",
                           subset = NULL, datfr = TRUE, 
                           auto_convert = TRUE, ...) {
    # checks
    assertArray(data, mode = "atomic", min.d = 1L, .var.name = "data")
    group_fun <- match.fun(group_fun)
    opt <- list(...)
    dimn <- dimnames(data)
    dimid <- names(dimn)
    if (is.null(dimid) || any(dimid == ""))
        stop("transformArray: the input array must have named dimnames",
             call. = FALSE)
    #
    if (!is.null(subset)) {
        data <- subsetArray(data, subset)
        dimn <- dimnames(data)
        dimid <- names(dimn)
    }
    #
    formula <- as.formula(formula)
    # LHS
    if (length(formula) == 2L) {
        formula[[3L]] <- formula[[2L]]
        formula[[2L]] <- 
            if (!is.null(opt$value_name)) {
                as.symbol(opt$value_name)
            } else {
                as.symbol(formals(array2df)$value_name)
            }
    }
    LHS <- formula[[2L]]
    value_dims <- 
        regmatches(deparse(LHS), 
                   gregexpr("(?<=\\[).+?(?=\\])", deparse(LHS), perl = TRUE)
                   )[[1L]]
    if (length(value_dims)) value_dims <- strsplit(value_dims, ", *")[[1L]]
    LHS <- parse(text = sub("\\[.*\\]", "", deparse(LHS)))[[1L]]
    # evaluate RHS
    RHS <- deparse(formula[[3L]])
    RHS <- gsub(" ", "", RHS)
    RHS <- gsub("-", "+-", RHS)
    RHS <- strsplit(RHS, "\\|")[[1L]]
    splitdims <- 
        if (length(RHS) > 1L) {
            strsplit(RHS[2L], "\\+")[[1L]] 
        } else {
            NULL
        }
    dims <- strsplit(RHS, "\\+")[[1L]]
    # collapse
    avgdims <- sub("^-", "", dims[grepl("^-", dims)])
    dims <- dims[!grepl("^-", dims)]
    if ("." %in% dims) {
        dims <- setdiff(dimid, c(avgDims, value_dims, splitdims))
    } else if ("." %in% splitdims) {
        splitdims <- setdiff(dimid, c(avgDims, value_dims, dims))
    }
    avgdims <- setdiff(dimid, c(dims, value_dims, splitdims))
    if (length(avgdims)) data <- avgDims(data, avgdims)
    # split
    data <-
        if (!is.null(splitdims)) {
            if (is.null(group)) {
                group <- vector("list", length(splitdims))
            } else if (!is.list(group)) {
                if (length(splitdims) == 1L) {
                    group <- list(group)
                } else {
                    stop(paste0(
                        "The group argument must be a list if there are ",
                        "more than one splitting dimensions"), call. = FALSE)
                }
            }
            ind <-
                if (is.null(names(group))) {
                    !logical(length(splitdims))
                } else {
                    is.na(names(group)) | names(group) == ""
                }
            names(group)[ind] <- splitdims[ind]
            splitArray(data, splitdims, group)
        } else {
            list(data)
        }
    if (!is.null(splitdims)) {
        for (i in seq_along(data)) {
            data[[i]] <- do(group_fun, data[[i]], splitdims)
            if (!is.atomic(data[[i]])) {
                stop(paste0(
                    "transformArray: 'group_fun' must return an atomic object", 
                    "(a vector, matrix, or array)"), call. = FALSE)
            }
        }
    }
    # evaluate LHS
    if (length(LHS) > 1L) {
        value_name <- as.character(LHS[[2]])
        if (identical(value_name, ".")) value_name = ""
        LHS[[2]] <- quote(x)
        data[] <- lapply(data, function(x) eval(LHS))
    } else {
        value_name <- as.character(LHS)
    }
    # back to array
    dimn <- attr(data, "dimnames")
    data <- bindArrays(data, along = 0L)
    data <- dim2multidim(data, 1,
                         expand.grid(dimn,
                                     KEEP.OUT.ATTRS = FALSE,
                                     stringsAsFactors = FALSE))
    # permute to the original dimension order
    dimid[dimid %in% splitdims] <- names(dimn)
    dimid <- intersect(dimid, names(dimnames(data)))
    data <- apermArray(data, dimid)
    # transform to a data.frame
    if (datfr) {
        data <- decorateDims_(data)
        dimn <- dimnames(data)
        singleton <- which(dim(data) == 1L &
                               grepl("_Dim", names(dimn)))
        singleton <- names(dimn)[singleton]
        data <- array2df(data, value_name = value_name, 
                         value_dim = value_dims,
                         auto_convert = auto_convert, 
                         ...)
        if (length(singleton) > 0)
            data <- data[setdiff(colnames(data), singleton)]
    }
    # return
    data
}


#' Convert (coerce) variables according to various schemes
#' 
#' \code{autoConvert} converts the variables which meet specific conditions 
#' on the basis of pre-defined conversion rules.
#' @param dat an object
#' @param select numeric or character indices of the variables or list elements
#' which sould be converted, if 'dat' is a data.frame or a list, respectively
#' @param conversion a character vector referring to the conversion rule
#' which should be applied for the selected variable/list element. If a single 
#' value, the same rule is applied for all variables/list elements. The default
#' is 'ANY', which is handled in a special way (see Details).
#' @param rules a list of conversion rules, see \code{\link{convertParams}}.
#' To save typing, \code{.(key = value)} format is also accepted, and the 
#' arguments inside the parentheses are forwarded to 
#' \code{\link{convertParams}}.
#' @details If 'conversion' is 'ANY' (the default), the 'IF' conditions in 
#' the rule definitions in 'rules' are tested and for the first successful test,
#' the given rule is selected. To avoid these sequential tests, provide the 
#' name of the rule definition explicitly in 'conversion'.
#' The rule definitions can be extended by arbitrary rules; it is suggested to
#' prepare the rules by calling \code{\link{convertParams}} in advance (see
#' Examples).
#' @param keep_dim a logical value whether the dimensions and dimension names
#' should be retained after the conversion (default: TRUE)
#' @export
#' @examples
#' # create an example list with various variable types
#' x <- list(A = c(1L, 0L, 0L),  # integer, but could be simplified to logical 
#'           B = matrix(c(1, 3, 2, 1), 2, 2), # double, but could be integer  
#'           C = c("1.2", "1", "0.92"),   # character, which could be double
#'           D = factor(c("a", "b", "a"), levels = c("b", "a")) 
#'                  # factor, which might be converted to character
#'           )
#' ( x_simplified <- autoConvert(x) )
#' \dontshow{
#' stopifnot(is.logical(x_simplified$A))
#' stopifnot(is.integer(x_simplified$B))
#' stopifnot(identical(dim(x$B), dim(x_simplified$B)))
#' stopifnot(is.double(x_simplified$C))
#' stopifnot(is.character(x_simplified$D))
#' }
#' 
#' # a simple way to convert only those list elements which are factors
#' autoConvert(x, conversion = "factor")
#' 
#' # suppose you do not want to convert characters at all
#' rules <- convertParams()
#' new_rules <- rules[setdiff(names(rules), "character")]
#' ( x_simplified2 <- autoConvert(x, rules = new_rules) )
#' \dontshow{
#' stopifnot(identical(x$C, x_simplified2$C))
#' }
#' 
#' # you can also convert an atomic object
#' autoConvert(c("1.12", "3.4"))
#' \dontshow{
#' stopifnot(identical(autoConvert(c("1.12", "3.4")),
#'                     as.numeric(c("1.12", "3.4"))))
#' }
autoConvert <- function(dat, select = NULL, conversion = "ANY",  
                        rules = convertParams(), keep_dim = TRUE) {
    #
    # workhorse function
    convert <- function(x, conv, rules, x_ind, dat_type) {
        # early return if 'x' does not match 'IF';
        # choose the right element from 'rules' function list
        if (!identical(conv, "ANY") && 
            !isTRUE(rules[[conv]]$IF(x))) {
            return(x)
        } else if (identical(conv, "ANY")) {
            conv_rule <- NULL
            for (r in rules) {
                if (isTRUE(r$IF(x))) {
                    conv_rule <- r
                    break
                }
            }
            if (is.null(conv_rule)) return(x)
        } else {
            conv_rule <- rules[[conv]]
        }
        # informative error message
        err_msg <- "autoConvert: the coercion resulted in an error"
        if (dat_type != "other") {
            tmp <- switch(dat_type,
                          list = " element of the list",
                          data.frame = " variable")
            err_msg <- paste0(err_msg, 
                              sprintf(" for the %d.", x_ind),
                              tmp)
        }
        errFn <- function(e) stop(err_msg, call. = FALSE)
        # 
        for (cr in seq_along(conv_rule$DO)) {
            DO <- conv_rule$DO[[cr]]
            EVAL <- conv_rule$EVAL
            out <- try(suppressWarnings(DO(x)), silent = TRUE)
            # if error, continue with the next cr
            if (inherits(out, "try-error")) {
                out <- x
                next
            }
            # evaluete
            ev <- try(EVAL(x, out), silent = TRUE)
            if (length(eval) > 1L) {
                stop(paste0(
                    "autoConvert: ",
                    "the EVAL function must return a logical value (either ",
                    "TRUE or FALSE)"), call. = FALSE)
            } else if (inherits(ev, "try-error")) {
                stop(paste0(
                    "autoConvert: ",
                    "the EVAL function exited with an error: ",
                    ev), call. = FALSE)
            } else if (isTRUE(ev)) {
                if (keep_dim) {
                    setattr(out, "dim", dim(x))
                    setattr(out, "dimnames", dimnames(x))
                }
                break
            } else {
                out <- x
            }
        }
        # return
        out
    }
    #
    # argument checks
    #
    # deparse 'convert'
    rules <- argumentDeparser(substitute(rules), "convertParams")
    # check 'dat'
    if (is.atomic(dat)) {
        dat_type <- "other"
        dat <- list(dat)
    } else if (is.data.frame(dat)) {
        dat_type <- "data.frame"
    } else if (is.list(dat)) {
        dat_type <- "list"
    } else {
        stop(paste0(
            "autoConvert: ",
            "'dat' must be an atomic or list/data.frame object"), 
            call. = FALSE)
    }
    # check 'select'
    if (!is.null(select)) {
        if (!all(select %in% seq_along(dat)) &&
            !all(select %in% names(dat))) {
            stop("autoConvert: 'select' has invalid element(s)", call. = FALSE)
        }
    } else {
        select <- seq_along(dat)
    }
    select <- select[!vapply(dat[select], is.null, logical(1L))]
    # check 'conversion'
    assertChoice(conversion, c("ANY", names(rules)), .var.name = "conversion")
    conversion <- repLen(conversion, length(select), "conversion")
    # perform coercion
    for (i in seq_along(select)) {
        selind <- select[i]
        dat[[selind]] <- convert(
            dat[[selind]],
            conversion[i], rules,
            selind, dat_type
        )
    }
    # return
    if (dat_type == "other") {
        dat[[1L]]
    } else {
        dat
    }
}
