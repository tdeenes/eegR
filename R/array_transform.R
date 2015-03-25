#
# <<< array reformatting functions >>> --------
#

#' Add dimension names to a matrix or array with (partially) missing dimension 
#' names
#' 
#' \code{decorateDims} is not unlike \code{\link{provideDimnames}}. The main 
#' difference is that \code{decorateDims} can modify the input data in place and
#' it has a more convenient interface for providing names of dimensions
#' @param dat a matrix or array
#' @param names logical; if TRUE (default), named list of dimnames are returned. 
#' Ignored if new_names is provided.
#' @param new_names character vector for named dimnames
#' @param new_dimnames list of new dimension names. If new_names is not 
#' provided but new_dimnames is a named list, its names are considered as
#' new_names
#' @param in_place logical; if TRUE, dat is modified in place (without copy)
#' and it is returned only invisibly (default: FALSE). Use in_place with extra 
#' care because it modifies in place all objects which dat refers to.
#' @return The original matrix or array with non-null dimension names
#' @export
#' @seealso \code{\link{provideDimnames}} for a slightly different solution
decorateDims <- function(dat, names = TRUE, 
                         new_names = NULL, new_dimnames = NULL,
                         in_place = FALSE) {
    dims <- dim(dat)
    dimn <- dimnames(dat)
    if (is.null(dimn)) dimn <- vector("list", length(dims))
    if (!is.null(new_dimnames) && !is.list(new_dimnames)) {
        new_dimnames <- list(new_dimnames)
    }
    if (is.null(new_names) && !is.null(temp <- names(new_dimnames))) {
        new_names <- temp
    }
    counter <- 1L
    tempfn <- 
        if (is.null(new_dimnames)) {
            function() as.character(seq_len(dims[i]))
        } else {
            function() {
                out <- rep_len(new_dimnames[[counter]], dims[i])
                make.unique(as.character(out), "__")
            }
        }
    for (i in which(vapply(dimn, is.null, NA))) {
        dimn[[i]] <- tempfn()
        counter <- counter + 1L
    }
    dimn.n <- names(dimn)
    if (names || !is.null(new_names)) {
        if (is.null(dimn.n)) dimn.n <- rep("", length(dimn))
        ind <- dimn.n == "" | is.na(dimn.n)
        dimn.n[ind] <- 
            if (!is.null(new_names)) {
                out <- rep_len(new_names, sum(ind))
                make.unique(as.character(out), "__")
            } else {
                paste0("_Dim", which(ind))
            }
    }
    if (in_place) {
        setattr(dat, "dimnames", dimn)
        setattr(dimnames(dat), "names", dimn.n)
    } else {
        dimnames(dat) <- dimn
        names(dimnames(dat)) <- dimn.n
        return(dat)
    }
}


#' Fast in-place transformation to a matrix (without copy)
#' 
#' \code{matrixIP} transforms its data argument to a matrix by reference. No
#' copy is made at all, and it invisibly returns the matrix.
#' @param x a data vector, matrix or array
#' @param nrow the desired number of rows
#' @param ncol the desired number of columns
#' @param dimnames A dimnames attribute for the matrix: NULL or a list of 
#' length 2 giving the row and column names respectively. The list can be 
#' named, and the list names will be used as names for the dimensions. An empty 
#' list is treated as NULL, and a list of length one as row names. 
#' @param force_length logical. If TRUE (the default), \code{matrixIP} checks if 
#' length(x)==nrow*ncol. If not, x is recycled or subsetted to the desired 
#' length.
#' @note Use \code{matrixIP} with extra care because it modifies in place all
#' objects which x refers to. If you want to avoid this, call 
#' \code{x <- copy(x)} before calling \code{matrixIP} or use the standard way as
#' described in the Note section of \code{\link{matrix}}. However, for input 
#' objects created on-the-fly (e.g. a temporary vector), \code{matrixIP} is safe 
#' and more compact than the latter solution, and can be many times faster than 
#' \code{matrix}. 
#' @return A matrix (invisibly)
#' @export
#' @seealso \code{\link{arrayIP}} for in-place transformation of x to an array 
#' (without copy) and \code{\link{matrix}} for creating a new matrix without 
#' modifying the original input
matrixIP <- function(x, nrow, ncol, dimnames = NULL, force_length = TRUE) {
    if (missing(nrow) && missing(ncol)) {
        nrow <- length(x)
        ncol <- 1L
    } else if (missing(nrow)) {
        nrow <- length(x)/ncol
    } else if (missing(ncol)) {
        ncol <- length(x)/nrow
    }
    if (!is.null(dimnames) && 
            grepl(deparse(substitute(x)), 
                  deparse(substitute(dimnames)))) {
        dimnames <- copy(dimnames)
    }
    if (force_length && length(x) != nrow*ncol) {
        x <- rep_len(x, nrow*ncol)
    }
    setattr(x, "dim", c(nrow, ncol))
    setattr(x, "dimnames", dimnames)
    invisible(x)
}

#' Fast in-place transformation to an array (without copy)
#' 
#' \code{arrayIP} transforms its data argument to an array by reference. No
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
#' @param force_length logical. If TRUE (the default), \code{arrayIP} checks if 
#' length(x)==nrow*ncol. If not, x is recycled or subsetted to the desired 
#' length.
#' @note Use \code{arrayIP} with extra care because it modifies in place all
#' objects which x refers to. See \code{\link{matrixIP}} for further hints.
#' @return An array (invisibly) or a matrix if length(dim)==2L.
#' @export
#' @seealso \code{\link{matrixIP}} for in-place transformation of x to a matrix 
#' (without copy) and \code{\link{array}} for creating a new array without 
#' modifying the original input
arrayIP <- function(x, dim, dimnames = NULL, force_length = TRUE) {
    if (missing(dim)) {
        dim <- length(x)
    }
    if (!is.null(dimnames) && 
            grepl(deparse(substitute(x)), 
                  deparse(substitute(dimnames)))) {
        dimnames <- copy(dimnames)
    }
    if (force_length && length(x) != prod(dim)) {
        x <- rep_len(x, prod(dim))
    }
    setattr(x, "dim", dim)
    setattr(x, "dimnames", dimnames)
    invisible(x)
}


#' Extract or replace a part of an array
#' 
#' \code{subsetArray} is a convenience function for extracting or replacing a 
#' part of an array which has dimension names
#' @name subsetArray
#' @usage subsetArray(dat, subsets, keep_attr = TRUE, ...)
#' @usage subsetArray(dat, subsets) <- value
#' @param dat array to be subsetted
#' @param subsets a named list of character, numeric, or logical vectors 
#' indicating which levels of which dimensions to subset (see Details)
#' @param keep_attr a logical variable which determines if the result inherits 
#' the custom attributes of the input (TRUE, default) or not  
#' @param value a vector, matrix or array of the new values
#' @param ... further arguments to be passed to \code{\link[abind]{asub}}
#' @details Names of subsets indicate which dimensions are to be subsetted in 
#' the input array, and each list element indicates which levels of the given 
#' dimension will be selected. If a list element is a named empty vector, or 
#' the name of a dimension does not appear in subsets, all levels of the 
#' correspondig dimension will be selected.
#' @export
#' @seealso \code{\link[abind]{asub}}
#' @return A subset of an array or the array with replaced values
#' @rdname subsetArray
#' @export
# Extract a part of an array
subsetArray <- function(dat, subsets, keep_attr = TRUE, ...) {
    dimpos <- match(names(subsets), names(dimnames(dat)))
    out <- asub(dat, subsets, dimpos, ...)
    if (keep_attr) {
        a <- attributes(dat)
        a2keep <- setdiff(names(a), 
                          c("class", "comment", "dim", "dimnames", "names", 
                            "row.names", "tsp"))
        a <- a[a2keep]
        for (i in a2keep) setattr(out, i, a[[i]])
        # this is just to handle the temporary "factor_level" attribute of
        # the eeg data
        if (!is.null(attr(out, "factors")) && 
                ("factor_level" %in% names(subsets)) ) {
            tempa <- splitMarker(dimnames(out)$factor_level,
                                 colnames(attr(out, "factors")),
                                 splitchar = "\\|")
            setattr(out, "factors", tempa) 
        }
    }
    # return
    out
}

#' @rdname subsetArray
#' @export
# Replace a part of an array
`subsetArray<-` <- function(dat, subsets, value) {
    dimn <- dimnames(dat)
    dimn[names(subsets)] <- subsets
    dimvn <- names(dimnames(value))
    if (!is.null(dimvn) && !any(dimvn == "") && !anyNA(dimvn)) {
        value <- aperm(value, na.omit(match(names(dimn), dimvn)))
    }
    eval(parse(
        text = paste("dat[", 
                     paste(paste0("dimn$", names(dimn)), collapse = ","), 
                     "] <- value", sep = ""))
    )
    # return
    dat
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
    out <- aperm(dat, c(row_dim, col_dims))
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
    matrixIP(out, nrow(out), dimnames = out_dimnames)
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
    out <- aperm(array(dat, matdims), dimord)
    dimnames(out) <- dimn
    # return
    out
}

#' Transforms an array to a data.frame
#' 
#' \code{array2df} transforms an array to a data.frame 
#' @param dat a matrix or array
#' @param response_name the name of the variable in the data.frame which 
#' contains the values of the input (default: "values")
#' @param response_type a character value which specifies the type of the 
#' response (default: \code{typeof(dat)}). Possible types are "logical", 
#' "character", "numeric", "double", "integer", "factor", "raw", "complex". 
#' @param dim_types a character value or a named list which specifies 
#' the type of the variable in the data.frame corresponding to the dimension of 
#' the input array. If one character value is given, all variables are 
#' transformed to that type (default: "character"). 
#' @param row_names NULL or a single integer or character string specifying a 
#' column to be used as row names, or a character or integer vector giving the 
#' row names for the data frame
#' @param na_omit if TRUE, omit all rows from the data.frame which have missing 
#' values (default: FALSE)
#' @param ... named arguments passed to \code{\link{decorateDims}}
#' @export
#' @return A data.frame
#' @examples
#' # example data
#' data(erps)
#' 
#' # transform to data.frame, change the default response name to "amplitudes" 
#' str(array2df(erps, response_name = "amplitudes"))
#' 
#' # treat all dimensions as factors, except for the time dimension, which
#' # should be integer
#' str(array2df(erps, dim_types = list(stimclass = "factor", 
#'                                     pairtype = "factor",
#'                                     chan = "factor",
#'                                     time = "integer",
#'                                     id = "factor")))
array2df <- function(dat, response_name = "values", 
                     response_type = typeof(dat), dim_types = "character", 
                     row_names = NULL, na_omit = FALSE, ...) {
    valid_types <- c("logical", "character", "integer", 
                     "numeric", "double", "factor", "complex", "raw")
    if (!response_type %in% valid_types) {
        stop("Invalid response_type argument")
    }
    if (!all(unlist(dim_types, use.names = FALSE) %in% valid_types)) {
        stop("Invalid dim_type argument")
    }
    dimn <- dimnames(decorateDims(dat, ...))
    out <- expand.grid(dimn, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    if (!is.null(row_names)) {
        setattr(out, "row.names", row_names)
    }
    out[[response_name]] <- do.call(paste0("as.", response_type), 
                                    list(dat)) 
    if (!identical(dim_types, "character")) {
        if (!is.list(dim_types)) {
            if (length(dim_types) > 1) {
                stop("dim_types must be a single character value or a named list")
            }
            dim_types <- setNames(as.list(rep(dim_types, ncol(out) - 1)),
                                  setdiff(colnames(out), response_name))
        }
        for (i in names(dim_types)) {
            out[[i]] <- do.call(paste0("as.", dim_types[[i]]), list(out[[i]]))
        }
    }
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
    if (any(duplicated(unlist(dims)))) 
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
    out <- arrayIP(aperm(dat, unlist(dims)),
                   vapply(relist(dim(dat)[unlist(dims)], dims), prod, 0),
                   dimn)
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
#'  \code{revMergeDims} sets back the original array transformed by 
#'  \code{\link{mergeDims}}
#'  @param dat numeric matrix or array with merged dimensions
#'  @export
#'  @return An array of the same dimension attributes as the array which 
#'  \code{\link{mergeDims}} was called on
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
splitArray <- function(dat, whichdim, f = NULL, drop = FALSE) {
    subFn <- function(ind) {
        abind::asub(dat, lapply(ind, unlist), whichdim_num, drop = drop) 
    }
    if (is.character(whichdim)) {
        whichdim_num <- match(whichdim, names(dimnames(dat)))
        if (anyNA(whichdim_num)) 
            stop("Wrong dimension name(s) provided")
    } else {
        whichdim_num <- whichdim
    }
    dat <- decorateDims(dat)
    if (is.null(f)) {
        f <- dimnames(dat)[whichdim_num]
    } else if (!is.list(f)) {
        f <- list(f)
        names(f) <- names(dimnames(dat))[whichdim_num[1]]
    }
    if (length(f) != length(whichdim_num)) {
        stop("Length of f must match the length of whichdim")
    }
    for (i in seq_along(f)) {
        if (is.null(f[[i]])) {
            f[[i]] <- seq_len(dim(dat)[whichdim_num[i]])
        } else if (is.list(f[[i]])) {
            f[[i]] <- do.call(paste, list(f[[i]], sep="_"))
        }
    }
    f <- lapply(f, function(x) split(seq_along(x), x))
    out.dimnames <- lapply(f, names)
    out.dim <- vapply(out.dimnames, length, integer(1))
    setattr(out.dimnames, "names", names(out.dim))
    f <- expand.grid(f, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    out <- lapply(1:nrow(f), function(i) subFn(f[i, ]))
    setattr(out, "dim", out.dim)
    setattr(out, "dimnames", out.dimnames)
    setattr(out, "names", 
            do.call("paste", c(lapply(f, names), list(sep = "."))))
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
#' by \code{\link{splitArrays}} or you want to bind on multiple dimensions; 
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
                if (all.equal(names(dimnames(x)), 
                              setdiff(alldimn, along_name))) {
                    x
                } else {
                    aperm(x, setdiff(alldimn, along_name))    
                }
            } else {
                if (all.equal(names(dimnames(x)), alldimn)) {
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
        anyNA(x) || is.null(x) || any(x == "") || any(duplicated(x))
    }
    dimnCheck <- function(x) {
        if (is_bad(names(x))) return(FALSE)
        if (any(vapply(x, is_bad, logical(1L)))) return(FALSE)
        TRUE
    }
    #
    dat <- list(...)
    if (length(dat) > 1L && any(vapply(dat, is.list, logical(1L)))) {
        stop("Only one list is allowed as an argument")
    }
    if (is.list(dat[[1L]])) dat <- dat[[1L]]
    #
    dimn <- lapply(dat, function(x) dimnames(x))
    if (any(checks <- !vapply(dimn, dimnCheck, logical(1L)))) {
        bad <- which(checks)
        stop(paste0(
            "Dimension names of the arrays ",
            paste(bad, collapse = ", "),
            " are not appropriate (see Arguments in help('mergeArrays')"))
    }
    ndimn <- lapply(dimn, names)
    if (length(unique(lapply(ndimn, sort))) > 1L)
        stop("There are unique dimension identifiers")
    all_dimn <- if (sort_dims) dimn[[1L]][order(ndimn[[1L]])] else dimn[[1L]]
    for (i in names(all_dimn)) {
        for (j in dimn[-1L]) {
            all_dimn[[i]] <- union(all_dimn[[i]], j[[i]])
        }
        if (sort_dimlevels) all_dimn[[i]] <- sort(all_dimn[[i]])
    }
    out <- base_value[1L]
    storage.mode(out) <- typeof(dat[[1L]])
    out <- arrayIP(out, vapply(all_dimn, length, integer(1L)),
                   all_dimn)
    for (i in dat) {
        subsetArray(out, dimnames(i)) <- i
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
#' @param splitFac named list of between-subject factors (for splitting ERP data)
#' @param subFac named list of within-subject factors (for subsetting ERP data)
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
                arrayIP(temp, vapply(dimn.perm, length, 0L), dimn.perm)
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
#' that class): a symbolic description of the transformation steps before 
#' converting the \code{array} to a \code{data.frame}. The details of how to
#' specify the transformations are given under ‘Details’.
#' @param data a matrix or an array. Must have named dimnames.
#' @param group a list of grouping factors in the order of appearance in the 
#' transformation formula (see 'Details'). If a named list is provided, those 
#' names are used as the names of dimnames for the given grouping dimensions.
#' It can be a simple vector if there are only one splitting factor.
#' @param subset a list of subsetting vectors on the input array passed to 
#' \code{\link{subsetArray}} before any transformation steps
#' @param datfr a logical value (default: TRUE) if the resulting array shall be 
#' transformed to a data.frame
#' @param auto_dimtype a logical value (default: TRUE) if dimension types 
#' shall be automatically transformed to "logical", "integer", "numeric" or
#' "factor". Ignored if datfr is FALSE or if dim_types is explicitly provided.
#' @param ... additional parameters to be passed to \code{\link{array2df}}
#' @details The formula interface of \code{transformArray} shall be given in
#' the form of 
#' \code{whateverFn(response_name, args) ~ dimA + .(dimC, dimD)}
#' Dimensions which are not present on the right hand side (RHS) of the formula
#' (e.g., dimB) are collapsed by calling \code{\link{avgDims}}. A simple . 
#' can be used in the RHS and means "all dimensions not otherwise present in 
#' the formula". 
#' 
#' Dimensions in parentheses (dimC and dimD, note the dot before 
#' the parenthesis) refer to the to-be-splitted dimensions (the corresponding 
#' list of splitting vectors must be given in the \code{group} argument). 
#' Splitting is performed by calling \code{\link{splitArray}}, and the splitted 
#' dimensions are immediately collapsed in the resulting list by calling 
#' \code{\link{avgDims}}. 
#' 
#' The left-hand side (LHS) of the formula provides the name of the response 
#' variable ("y") if data.frame conversion is requested (which is the default), 
#' and also allows arbitrary computations (e.g., calling \code{\link{compGfp}} 
#' on the array. 
#' 
#' If there were any splitting factors, the resulting list is
#' back-transformed to an array. After the transformations, the array is
#' converted to a data.frame if \code{datfr} is TRUE.
#' @export
#' @return A data.frame if datfr is TRUE, and an array if datfr is FALSE
#' @examples
#' # example dataset
#' data(erps)
#' dat_id <- attr(erps, "id") # to get reading group memberships
#' 
#' # collapse all dimensions except for stimclass and pairtype
#' str(transformArray(y ~ stimclass + pairtype, erps))
#' 
#' # analyze separately dyslexic and control subjects, also compute Global Field 
#' # Power
#' res1 <- transformArray(compGfp(y, keep_channels = TRUE) ~ . + .(id), 
#'                        erps, list(readgroup = dat_id$group))
#' 
#' # the same with much more typing
#' res2 <- splitArray(erps, "id", list(readgroup = dat_id$group))
#' split_dnn <- dimnames(res2)
#' names(split_dnn) <- names(dim(res2))
#' res2 <- lapply(res2, avgDims, "id")
#' res2 <- lapply(res2, compGfp, keep_channels = TRUE)
#' res2 <- bindArrays(res2, along = 0L)
#' res2 <- dim2multidim(res2, 1, expand.grid(split_dnn))
#' res2 <- array2df(res2, response_name = "y", dim_types = "factor")
#' res2$time <- as.integer(as.character(res2$time))
#' stopifnot(identical(res1, res2)) 
transformArray <- function(formula, data, group = NULL, subset = NULL,
                           datfr = TRUE, auto_dimtype = TRUE, ...) {
    # checks
    if (!is.atomic(data) || !is.array(data)) {
        stop("Provide a matrix or array as input")
    }
    dn <- dimnames(data)
    dnn <- names(dn)
    if (is.null(dnn) || any(dnn == "")) 
        stop("The input array must have named dimnames")
    #
    if (!is.null(subset)) {
        data <- subsetArray(data, subset)
        dn <- dimnames(data)
        dnn <- names(dn)
    }
    #
    formula <- as.formula(formula)
    LHS <- formula[[2]]
    RHS <- deparse(formula[[3]])
    # evaluate RHS
    RHS <- gsub(" ", "", RHS)
    dims <- strsplit(RHS, "\\+")[[1]]
    # collapse
    avgdims <- 
        if ("." %in% dims) {
            NULL
        } else {
            setdiff(dnn, gsub("\\(|\\)", "", dims))
        }
    if (!is.null(avgdims)) data <- avgDims(data, avgdims)
    # split
    pd <- dims[grepl("\\.\\(.*\\)", dims)]
    splitdims <- 
        if (length(pd) > 0L) {
            strsplit(gsub("\\.\\(|\\)", "", pd), ",")[[1]]
        } else {
            NULL
        }
    data <- 
        if (!is.null(splitdims)) {
            if (is.null(group)) {
                stop("Provide grouping factors (see group argument)")
            } else if (!is.list(group)) {
                if (length(splitdims) == 1L) {
                    group <- list(group)
                } else {
                    stop("The group argument must be a list if there are more
                         than one splitting dimensions")
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
    if (!is.null(splitdims)) data[] <- lapply(data, avgDims, splitdims)
    # evaluate LHS
    if (length(LHS) > 1L) {
        response_name <- as.character(LHS[[2]])
        LHS[[2]] <- quote(x)
        data[] <- lapply(data, function(x) eval(LHS))
    } else {
        response_name <- as.character(LHS)
    }
    # back to array
    dnn <- attr(data, "dimnames")
    setattr(dnn, "names", names(attr(data, "dim")))
    data <- bindArrays(data, along = 0L)
    data <- dim2multidim(data, 1, 
                         expand.grid(dnn, 
                                     KEEP.OUT.ATTRS = FALSE, 
                                     stringsAsFactors = FALSE))
    # transform to a data.frame
    if (datfr) {
        data <- decorateDims(data)
        dn <- dimnames(data)
        singleton <- which(dim(data) == 1L & 
                               grepl("_Dim", names(dn)))
        singleton <- names(dn)[singleton]
        datfr_options <- list(...)
        if (is.null(datfr_options$dim_types) && auto_dimtype) {
            datfr_options$dim_types <- suppressWarnings(
                lapply(dn, function(x) {
                    if (identical(x, as.character(as.logical(x)))) "logical"
                    else if (identical(x, as.character(as.integer(x)))) "integer"
                    else if (identical(x, as.character(as.numeric(x)))) "numeric"
                    else "factor"
                }))
        }
        if (is.null(datfr_options$response_name)) {
            datfr_options$response_name <- response_name
        }
        data <- do.call("array2df", c(list(data), datfr_options))
        if (length(singleton) > 0) 
            data <- data[setdiff(colnames(data), singleton)]
    }
    # return
    data
}
