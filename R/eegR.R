
#
# <<< initialize >>> --------
#

#' @useDynLib eegR 
#' @import abind matrixStats permute data.table parallel doParallel foreach iterators
#' @importFrom Rcpp evalCpp
NULL

.onAttach <- function(lib, pkg) {
    packageStartupMessage(paste0("*** eegR ",
                                 packageVersion("eegR"),
                                 " loaded ***"), 
                          appendLF = TRUE)
}

.onUnload <- function (libpath) {
    library.dynam.unload("eegR", libpath)
}

#
# <<< document example data >>> --------
#

#' Averaged ERPs in a visual word recognition experiment
#' 
#' A dataset containing averaged event-related potentials (ERPs) from 20 
#' participants. The data were collected in a same-different matching task, and
#' simplified for present purposes. Only brain responses to the target stimuli 
#' from 10 dyslexic and 10 control participants are included. There were 3 types
#' of stimulus classes (A, B, and C) and 3 types of stimulus pairs (identical,
#' substituted, and transposed). The data were downsampled to 500 Hz, and cover 
#' the following time window: -50 ms to 500 ms. 
#' @docType data
#' @keywords datasets
#' @name erps
#' @usage data(erps)
#' @format An array with 5 dimensions: stimulus class (stimclass, 3 level) x 
#' pair type (pairtype, 3 level) x channels (chan, 33 levels) x time (time, 276 
#' levels) x participant (id, 20 levels). Additionally, the array has an 
#' attribute "id" which is a data.frame of the group memberships (dyslexic or 
#' control), and an attribute "chan" which is a data.frame of the electrode 
#' positions in spherical coordinates.
NULL

#
# <<< simple utility functions >>> --------
#

#' Check for availability of packages
#' @keywords internal
reqFn <- function(packages) {
    for (i in packages) {
        if(!requireNamespace(i, quietly = TRUE)) {
            stop("You have to install package:", i ," before using this function")
        }
    }
    TRUE
}


#' Assign the elements of a named list to the enclosing environment
#' 
#' \code{assignList} assigns the elements of a named list to the enclosing environment. 
#' Be aware that the function does not check if an object with the same name exists in
#' the enclosing environment - if it exists, it will be overwritten.
#' @param listdat a list with named elements
#' @param verbose a logical variable (default: TRUE) which determines if 
#' a warning should be sent to the console
#' @param overwriteGlobal a logical variable (default FALSE) which determines 
#' if the function is allowed to write to the global environment
#' @export
#' @return The function is invoked for its side effect, which is assigning list
#' elements to the enclosing environment
assignList <- function(listdat, verbose = TRUE, overwriteGlobal = FALSE) {
    min_calling_frame <- ifelse(overwriteGlobal, 1, 2)
    if (sys.nframe() >= min_calling_frame) {
        if (is.null(names(listdat)) || 
                any(names(listdat) == "")) {
            stop("All elements of the assigned list should have a name!")
        }
        for (i in names(listdat)) {
            assign(i, listdat[[i]], pos = parent.frame())
        }
        if (verbose) {
            warning(
                paste("The following variables were assigned to the environment:", 
                      paste(names(listdat), collapse = " ")))
        }
    }
}

#' Find local peaks (maxima) and valleys (minima) of a vector
#' 
#' \code{findExtremes} identifies local peaks and valleys of a vector.
#' Endpoints are considered peaks/valleys if they are part of a plateau.
#' @param x a vector
#' @export
#' @return A vector of the same length as x, with values 1 (local maximum), 
#' -1 (local minimum) and 0.
findExtremes <- function(x) {
    if (is.null(names(x))) names(x) <- seq_along(x)
    out <- rep(0, length(x))
    names(out) <- names(x)
    ind <- diff(x)
    ind <- ind[ind != 0]
    startval <- sign(ind[1])*Inf
    endval <- sign(-ind[length(ind)])*Inf
    ind <- diff(c(startval, x, endval))
    if (ind[2] == 0) {
        ind[2] <- -sign(ind[1])
    }
    ind <- ind[ind != 0]
    ind[1] <- 0
    ind[-1] <- diff(sign(ind))
    out[names( ind[which(ind > 0)-1] )] <- -1
    out[names( ind[which(ind < 0)-1] )] <- 1
    for (i in 2:length(out)) {
        if (x[i-1] == x[i]) {
            out[i] <- out[i-1]
        }
    }
    if (out[1] != out[2]) out[1] <- 0
    if (out[length(out)] != out[length(out)-1]) out[length(out)] <- 0
    names(out) <- names(x)
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
#' @param force_rollapply logical variable; if yes, zoo:::rollapply is called 
#' (default = FALSE).
#' @param ... optional arguments to the corresponding function in caTools or
#' zoo:::rollapply
#' @details If FUN is one of min, max, mean, sd, mad, quantile (OR "min", "max",
#' "mean", etc.) rollFun calls the corresponding function from the caTools 
#' package (e.g. caTools:::runmin). Otherwise, or if force_rollapply is TRUE,
#' zoo:::rollapply is called.
#' @export
#' @return An object having the same attributes as dat.
#' @import zoo 
#' @import caTools
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
                    if (identical(FUN, runquantile)) {
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

#' Fast version of unique for matrices
#'
#' \code{fastUnique} finds the unique rows or columns of a matrix.
#' @param x matrix
#' @param units_in_rows units are in rows (TRUE, default) or columns of x
#' @export
#' @return matrix without duplicated units
#' @import Kmisc
fastUnique <- function(x, units_in_rows = TRUE) {
    if (!is.matrix(x)) stop("Provide a matrix as input!")
    if (units_in_rows) x <- t(x)
    dupl <- duplicated(Kmisc::str_collapse(lapply(1:ncol(x), function(i) x[,i])))
    x <- x[, !dupl, drop = F]
    if (units_in_rows) x <- t(x)
    # return
    x
}

#' Create list with substituted names
#' 
#' \code{listS} creates a named list where names are substituted with the 
#' content of the referenced variable.
#' @param ... objects; if not named, listS is equilent to \code{\link{list}}. 
#' Names which should be substituted should start with a dot (.) or INDICES has 
#' to be provided. 
#' @param indices character or numeric vector indicating the position of those 
#' list elements whose name should be substituted. If provided, relevant names 
#' in ... should not be dotted.
#' @export
#' @return A list with substituted names.
listS <- function(..., indices = NULL) {
    list_def <- list(...)
    if (is.null(names(list_def))) {
        return( list_def )
    }
    nam <- names(list_def)
    if (is.null(indices)) {
        indices <- grep("^[.]", nam)
        newnam <- sub("^[.]", "", nam[indices])
    } else {
        newnam <- nam[INDICES]
    }
    if (length(newnam) == 0) return( list_def )
    newnam <- structure(
        lapply(newnam, function(x) {
            for (n in 1:sys.nframe()) {
                tempx <- try(get(x, envir = parent.frame(n)), silent = TRUE)
                if (!inherits(tempx, "try-error")) return( tempx )
            }
        }),
        names = newnam)
    for (i in 1:length(newnam)) {
        if (!is.atomic(newnam[[i]]) || length(newnam[[i]]) > 1) {
            stop(
                "The object denoted by '", names(newnam)[i], "' is not a character string!"
            )
        }
    }
    names(list_def)[indices] <- unlist(newnam, use.names = FALSE)
    # return
    list_def
}

#' Correct p-values in time-series
#' 
#' \code{pvalueConsec} sets a minimal consecutive length criterion on the 
#' p-values in a time series
#' @param dat a numeric matrix or array containing p-values, or a list of such 
#' matrices/arrays. dat must have a named time dimension.
#' @param sig_level numeric value (default: 0.05), the level of significance
#' @param min_length numeric value (default: 10), the minimum number of 
#' consecutive significant time points
#' @export
#' @return An object of the same dimension as the input data
pvalueConsec <- function(dat, sig_level = 0.05, min_length = 10) {
    pCorr <- function(x) {
        temp <- rle(x)
        ind <- (temp$values == 1) & (temp$lengths < min_length)
        temp$values[ind] <- 0
        return( inverse.rle(temp) )
    }
    pCorrMat <- function(datarr) {
        datarr <- array2mat(datarr, "time")
        datarr[apply(datarr <= sig_level, 2, pCorr) == 0] <- 1
        datarr <- mat2array(datarr)
        return(datarr)
    }
    if (is.list(dat)) {
        out <- lapply(dat, pCorrMat)
    } else {
        out <- pCorrMat(dat)
    }
    # return
    out
}

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
#' and more compact than the latter solution. 
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
    if (force_length && length(x) != nrow*ncol) {
        x <- rep_len(x, nrow*ncol)
    }
    setattr(x, "dim", c(nrow, ncol))
    if (!is.null(dimnames)) {
        setattr(x, "dimnames", dimnames)
    }
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
    if (force_length && length(x) != prod(dim)) {
        x <- rep_len(x, prod(dim))
    }
    setattr(x, "dim", dim)
    if (!is.null(dimnames)) {
        setattr(x, "dimnames", dimnames)
    }
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
#' @param ... further arguments to be passed to asub
#' @details Names of subsets indicate which dimensions are to be subsetted in 
#' the input array, and each list element indicates which levels of the given 
#' dimension will be selected. If a list element is a named empty vector, or 
#' the name of a dimension does not appear in subsets, all levels of the 
#' correspondig dimension will be selected.
#' @export
#' @return A subset of an array or the array with replaced values
#' @import abind
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
#' \code{array2df} transforms an array with dimnames to a data.frame 
#' @param dat a matrix or array with named dimnames
#' @param response_name the name of the variable in the data.frame which 
#' contains the values of the input (default: "values")
#' @param dim_types a character value or a named list which specifies 
#' the type of the variable in the data.frame corresponding to the dimension of 
#' the input array. If one character value is given, all variables are 
#' transformed to that type (default: "character"). Possible types are 
#' "logical", "character", "numeric", "double", "integer", "factor".
#' @param na_omit if TRUE, omit all rows from the data.frame which have missing 
#' values (default: FALSE)
#' @param ... additional parameters to be passed to \code{\link{as.data.frame.table}}
#' @export
#' @return A data.frame
array2df <- function(dat, response_name = "values", dim_types = "character", 
                     na_omit = FALSE, ...) {
    if ("stringsAsFactors" %in% names(list(...))) {
        warning("The parameter stringsAsFactors was ignored. It was set to FALSE by default.")
    }
    out <- as.data.frame.table(dat, responseName = response_name, 
                               stringsAsFactors = FALSE, ...)
    if (!is.null(dim_types)) {
        if (!is.list(dim_types)) {
            if (length(dim_types) > 1) {
                stop("dim_types must be a single character value or a named list")
            }
            dim_types <- setNames(as.list(rep(dim_types, ncol(out) - 1)),
                                  setdiff(colnames(out), response_name))
        }
        stopifnot(all(unlist(dim_types, use.names = FALSE) %in% 
                          c("logical", "character", "integer", 
                            "numeric", "double", "factor")))
        for (i in names(dim_types)) {
            out[, i] <- do.call(paste0("as.", dim_types[[i]]), list(out[, i]))
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
    out
}

#' Splits an array along a given dimension
#' 
#' \code{splitArray} splits an array along given dimension(s) into a list of 
#' sub-arrays
#' @param dat numeric array (preferably with named dimnames)
#' @param whichdim numeric or character vector, the dimension(s) of the array 
#' to split along
#' @export
#' @return A list of subsets of the original data matrix/array
splitArray <- function(dat, whichdim) {
    if (length(whichdim) > 1) {
        out.dimnames <- dimnames(dat)[whichdim]
        out.dim <- vapply(out.dimnames, length, 0L)
        dat <- mergeDims(dat, whichdim)
        whichdim <- paste(whichdim, collapse=".")
    } else {
        out.dimnames <- NULL
        out.dim <- NULL
    }
    dimlevels <- dimnames(dat)[[whichdim]]
    out <- lapply(dimlevels, function(i) subsetArray(dat, listS(.whichdim = i)))
    dim(out) <- out.dim
    dimnames(out) <- out.dimnames
    names(out) <- dimlevels
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


#
# <<< general functions on arrays >>> --------
#

# iterator for array chunks
# iarray <- function(a, chunks, chunkgrid) {
#     i <- 0
#     nextElem <- function() {
#         if (i >= nrow(chunkgrid)) stop('StopIteration')
#         i <<- i + 1
#         subsetArray(
#             a,
#             subsets = mapply("==", chunks, chunkgrid[i, ]), 
#             drop = FALSE)
#     }
#     structure(list(nextElem=nextElem), class=c('iarray', 'iter'))
# }
# nextElem.iarray <- function(obj) obj$nextElem()

# chunkify <- function(dat, fun, arg_list = NULL, chunks = NULL,
#                      useparallel = FALSE, ncores = NULL, 
#                      par_method = c("snow", "multicore"), cl = NULL) {
#     require(doParallel)
#     if (useparallel) {
#         par_method = match.arg(par_method)
#         if (is.null(ncores)) ncores <- detectCores()
#         if (par_method == "snow" || .Platform$OS.type=="windows") {
#             if (is.null(cl)) cl <- makePSOCKcluster(ncores)
#             registerDoParallel(cl)
#             on.exit(stopCluster(cl))
#         } else {
#             registerDoParallel(ncores = ncores)
#         }
#     } else {
#         registerDoSEQ()
#     }
#     chunkCheckFn <- function(x, name) {
#         dimlen <- dimlens[name]
#         if (is.character(x)) {
#             x <- match(x, unique(x))
#         } else if (is.logical(x)) {
#             x <- as.integer(x) + 1L            
#         } else if (!is.integer(x)) {
#             x <- as.integer(x)
#         } 
#         out <- 
#             if (length(x) == 1) {
#                 rep(seq_len(x), 
#                     sapply(parallel::splitIndices(dimlen, x), length))
#             } else if (length(x) == dimlen) {
#                 as.integer(x)
#             } else {
#                 rep(as.integer(x), length.out = dimlen)
#             }
#         out
#     }
#     # 
#     if (is.null(chunks)) {
#         return( do.call(fun, append(list(dat), arg_list)) )
#     }
#     if (!is.list(chunks) && is.null(names(chunks))) {
#         stop("chunks must be a named list")
#     }
#     dimlens <- dim(dat)
#     names(dimlens) <- names(dimnames(dat))
#     #
#     chunks <- setNames(
#         mapply(chunkCheckFn, chunks, names(chunks),
#                SIMPLIFY = FALSE, USE.NAMES = FALSE),
#         names(chunks))
#     chunkgrid <- expand.grid(lapply(chunks, unique), 
#                              KEEP.OUT.ATTRS = FALSE,
#                              stringsAsFactors = FALSE)
#     outpart <- foreach(x = iarray(dat, chunks, chunkgrid)) %dopar% 
#         do.call(fun, append(list(x), arg_list))
#     out.dims <- dim(outpart[[1]])
#     out.dims[match(names(chunks), names(dimnames(outpart[[1]])))] <- 
#         dimlens[names(chunks)]
#     out.dimnames <- setNames(
#         lapply(names(dimnames(outpart[[1]])), function(n) {
#             if (n %in% names(chunks)) {
#                 dimnames(dat)[[n]]
#             } else {
#                 dimnames(outpart[[1]])[[n]]
#             }}),
#         names(dimnames(outpart[[1]]))
#     )
#     na <- NA
#     storage.mode(na) <- storage.mode(outpart[[1]])
#     out <- array(na, out.dims, out.dimnames)
#     for (i in 1:nrow(chunkgrid)) {
#         subsetArray(
#             out,
#             subsets = mapply("==", chunks, chunkgrid[i, ])) <- outpart[[i]]
#     }
#     return( out )
# }


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
            multmat <- multmat/colSums(multdat)
            out <- fnDims(dat, target_dim, avgfn, arg_list = list(y = multmat), 
                          vectorized = TRUE)
        } else if (length(dimnames(dat)[[target_dim]]) %% bin_length > 0) {
            stop("Target dimension length is not multiple of bin_length")
        } else {
            bins <- length(dimnames(dat)[[target_dim]])/bin_length
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

#
# <<< t-tests on arrays >>>----
#

#' Workhorse function of arrayTest (for independent samples t-test)
#' 
#' \code{indepTtest} computes t-values and additional p-values & degrees of 
#' freedom if requested.
isTtest <- function(dat, groups, id_dim = "id", mu = 0, var_equal = TRUE, 
                    verbose = FALSE, any_NA = NULL) {
    if (is.null(any_NA)) any_NA <- anyNA(dat)
    res <- lapply(1:2, function(g) {
        datx <- subsetArray(dat, listS(.id_dim = groups==g), drop = FALSE)
        mm <- avgDims(datx, id_dim, na_rm = any_NA)
        vv <- fnDims(datx, id_dim, colVars, arg_list = list(na.rm = any_NA), 
                     vectorized = TRUE, columnwise = TRUE)
        if (any_NA) {
            nn <- fnDims(!is.na(datx), id_dim, colSums, 
                         arg_list = list(na.rm = FALSE), 
                         vectorized = TRUE, columnwise = TRUE)
            nn[nn < 1] <- NA
        } else {
            nn <- sum(groups == g)
        }
        list(mm, vv, nn)
    })
    m1 <- res[[1]][[1]]; v1 <- res[[1]][[2]]; n1 <- res[[1]][[3]]
    m2 <- res[[2]][[1]]; v2 <- res[[2]][[2]]; n2 <- res[[2]][[3]]
    out <- 
        if (var_equal) {
            sd12 <- sqrt(
                ((n1-1)*v1 + (n2-1)*v2)/(n1 + n2 - 2)
            )
            (m1 - m2)/(sd12*sqrt(1/n1 + 1/n2))
        } else {
            sd12 <- sqrt(v1/n1 + v2/n2)
            (m1 - m2)/sd12
        }
    if (verbose) {
        df <- 
            if (var_equal) {
                n1 + n2 - 2 
            } else {
                sd12^4/( (v1/n1)^2/(n1-1) + (v2/n2)^2/(n2-1) )
            }
        setattr(out, "Df", df)
        setattr(out, "pvalues", 2 * pt(-abs(out), df))
    }
    out
}

#' Workhorse function of arrayTest (for one-sample or paired-samples t-test)
#' 
#' \code{indepTtest} computes t-values and additional p-values & degrees of 
#' freedom if requested.
osTtest <- function(dat, id_dim = "id", mu = 0, verbose = FALSE, any_NA = NULL) {
    if (is.null(any_NA)) any_NA <- anyNA(dat)
    m1 <- avgDims(dat, id_dim, na_rm = any_NA)
    sd1 <- fnDims(dat, id_dim, colSds, 
                  arg_list = list(na.rm = any_NA), vectorized = TRUE, 
                  columnwise = TRUE)
    if (any_NA) {
        n1 <- fnDims(!is.na(dat), id_dim, colSums, 
                     arg_list = list(na.rm = FALSE), vectorized = TRUE, 
                     columnwise = TRUE)
        n1[n1 < 1] <- NA
    } else {
        if (is.character(id_dim)) id_dim <- match(id_dim, names(dimnames(dat)))
        n1 <- dim(dat)[id_dim]
    }
    out <- (m1 - mu)/(sd1/sqrt(n1))
    if (verbose) {
        pvalues <- 2 * pt(-abs(out), n1 - 1)
        setattr(out, "Df", n1 - 1)
        setattr(out, "pvalues", pvalues)
    }
    out
}


#' Point-to-point t-tests (potentially with TFCE correction) on arrays
#' 
#' \code{arrayTtest} performs point-to-point t-tests on arrays. 
#' Permutation-based p-values and Threshold-free Cluster Enhancement (TFCE) 
#' correction can be requested.
#' @param arraydat a numeric array with named dimnames containing EEG (or 
#' other) data. Missing values are not allowed. Must have at least three
#' dimensions with names "chan", "time", and "id" (see also id_dim)
#' @param arraydat2 a numeric array with named dimnames containing EEG (or 
#' other) data. If provided, see the parameter \code{paired} for running 
#' dependent or independent-samples t-tests
#' @param paired logical scalar, only used if arraydat2 is provided. If paired 
#' is FALSE (default), the function computes independent samples t-tests, 
#' otherwise paired samples t-tests are performed
#' @param groups provides an alternative (and more efficient) way to perform 
#' independent samples t-tests; a character, factor, or integer vector which 
#' defines group membership. Groups is ignored if arraydat2 is not missing. 
#' NA values code subjects to drop.
#' @param mu a numeric scalar indicating the true value of the mean (or 
#' difference between means, if two-sample tests are performed)
#' @param var_equal a logical scalar whether the variances are equal (only 
#' relevant for independent-samples t-tests). If TRUE (default), the pooled 
#' variance is used to estimate the variance, otherwise the Welch (or 
#' Satterthwaite) approximation to the degrees of freedom is used.
#' @param id_dim name of the dimension which identifies the subjects 
#' (default: "id")
#' @param verbose logical value indicating if p-values should be computed for 
#' the traditional t-test results
#' @param nperm integer value giving the number of permutations (default: 999L)
#' @param useparallel logical value; if TRUE (default), computations are done
#' in parallel
#' @param par_method parallelization method; can be set explicitly to "snow" or
#' "multicore" (ignored on Windows OS), or chosen automatically ("default"). 
#' Ignored if cl is provided.
#' @param ncores integer value corresponding to the number of cores; 
#' if NULL (default), it is set to the maximum number of cores available. 
#' Ignored if cl is provided.
#' @param cl a cluster definition; if NULL (default), it is set up automatically
#' @param usetfce logical value whether TFCE (threshold-free cluster enhancement)
#' correction should also be computed (default: TRUE)
#' @param tfce_options a named list containing the channel neighbourhood matrix
#' (named ChN, no default) and the vector of the E/H parameters (named EH, 
#' defaults to c(0.66, 2))
#' @param seed an integer value which specifies a seed (default: NULL)
#' @details The function assumes that the input array contains at least three 
#' named dimensions: chan (corresponding to the channels [electrodes]) and time 
#' (corresponding to time points), and \code{id_dim} (corresponding to subjects). 
#' All other dimensions are treated in a similar way as chan and time, that is 
#' separate t-tests are computed for each level of those dimensions.
#' @export
#' @return A list object with t-values, TFCE-corrected t-values and
#' permutation-based p-values (if requested)
# TODO: compute also effect sizes + make it general for any arrays without chan
# and time dimensions + clear redundancies in the code for the two types of 
# t-tests
arrayTtest <- function(arraydat, arraydat2, paired = FALSE, groups = NULL,
                       mu = 0, var_equal = TRUE, id_dim = "id", verbose = TRUE, 
                       nperm = 999L, useparallel = TRUE, 
                       par_method = c("default", "snow", "multicore"), 
                       ncores = NULL, cl = NULL, usetfce = TRUE, 
                       tfce_options = NULL, seed = NULL) {
    # helper function for back-transform to original
    backFn <- function(x, origdimnames) {
        if (is.null(dim(x))) return(x)
        if ("TEMPX" %in% (dimn <- names(dimnames(x))))
            x <- subsetArray(x, list(TEMPX=1))
        out <- aperm(x, origdimnames[origdimnames %in% dimn])
        attribs <- attributes(x)
        attribs <- attribs[-match(c("dim", "dimnames"), names(attribs))]
        for (i in names(attribs)) setattr(out, i, attribs[[i]])
        # return
        out
    }
    # 
    if (useparallel) {
        par_method <- match.arg(par_method)
        if (newcl <- is.null(cl)) {
            if (is.null(ncores)) ncores <- detectCores()
            cl <- 
                if (.Platform$OS.type=="windows" || par_method == "snow") {
                    makePSOCKcluster(ncores)
                } else {
                    makeForkCluster(ncores)
                }
        }
        registerDoParallel(cl)
        on.exit(if (newcl) stopCluster(cl))
    } 
    else {
        registerDoSEQ()
    }
    if (usetfce) {
        if (is.null(tfce_options)) stop("Provide at least tfce_options$ChN")
        ChN <- tfce_options$ChN
        EH <- if (is.null(tfce_options$EH)) c(0.66, 2) else tfce_options$EH
    }
    if (missing(arraydat2)) {
        if (!is.null(groups)) {
            if (anyNA(groups)) 
                arraydat <- subsetArray(arraydat, listS(.id_dim=!is.na(groups)))
            groups <- as.integer(as.factor(groups))
            if (max(groups) != 2L) 
                stop("groups must contain two distinct values (e.g. 1 or 2, 'a' or 'b')")
        }
    } 
    else {
        if (paired) {
            if (!identical(dimnames(arraydat), dimnames(arraydat2)) ||
                    !identical(dim(arraydat), dim(arraydat2))) {
                stop("Input array dimensions must be identical")
            }
            arraydat <- arraydat2 - arraydat
            groups <- NULL
        } else {
            groups <- rep.int(1:2, 
                              c(length(dimnames(arraydat)[[id_dim]]),
                                length(dimnames(arraydat2)[[id_dim]])))
            dimnms <- names(dimnames(arraydat))
            arraydat <- abind(arraydat, arraydat2, 
                              along=which(names(dimnames(arraydat))==id_dim))
            setattr(dimnames(arraydat), "names", dimnms)
        }
    }
    any_NA <- anyNA(arraydat)
    origdimnames <- names(dimnames(arraydat))
    if (anyNA(origdimnames) || 
            !all(c(id_dim, "chan", "time") %in% origdimnames) ||
            any(duplicated(origdimnames))) {
        stop("Provide array(s) with named dimensions (at least id_dim, 'chan', 'time')")
    }
    arraydat <- aperm(arraydat, c(id_dim, "chan", "time",
                                  setdiff(names(dimnames(arraydat)),
                                          c(id_dim, "chan", "time"))))
    if (length(dim(arraydat)) == 3L) {
        temp <- dimnames(arraydat)
        arrayIP(arraydat, c(dim(arraydat), 1),
                c(temp, list(TEMPX="1")))
        rm(temp)
    }
    # two independent samples
    if (!is.null(groups)) {
        t_obs <- isTtest(arraydat, groups, id_dim, mu, 
                         var_equal, verbose, any_NA)
        if (nperm > 1) {
            grouplen <- length(groups)
            group2len <- sum(groups == 2L)
            maxperm <- choose(grouplen, group2len)
            if (nperm > maxperm) 
                stop("nperm is too high for this sample size (not enough unique combinations exist)")
            if (!is.null(seed)) set.seed(seed)
            if (nperm > maxperm/3) {
                randind <- combn(grouplen, group2len)
                randgroups <- matrix(1L, grouplen, ncol(randind))
                for (i in 1:ncol(randind)) randgroups[randind[,i], i] <- 2L
                randgroups <- randgroups[, sample.int(nrow(randgroups), nperm)]
            } 
            else {
                randgroups <- matrix(groups, grouplen, 1)
                repeat{
                    nc <- ncol(randgroups) - 1
                    if (nc >= nperm) {
                        randgroups <- randgroups[, 2:(nperm + 1)]
                        break
                    } 
                    else {
                        rg <- matrix(1L, grouplen, (nperm - nc) * 2)
                        for (i in 1:ncol(rg)) {
                            rg[sample.int(grouplen, group2len), i] <- 2L
                        }
                        randgroups <- cbind(randgroups, rg)
                        randgroups <- unique(randgroups, MARGIN = 2L)
                    }
                }
            }
            applydims <- seq(3, length(dim(t_obs)))
            if (usetfce) {
                tfce_obs <- arrayIP(0, dim(t_obs), dimnames(t_obs))
                tfce_obs[] <- apply(t_obs, applydims, tfceFn, 
                                    chn = ChN, eh = EH)
                stat_obs <- abs(tfce_obs)
            } 
            else {
                stat_obs <- abs(t_obs)
            }
            chunks <- min(10L, ceiling(nperm/max(1L, length(cl))))
            permres <- 
                foreach(x = iter(randgroups, by = "column", 
                                 chunksize = chunks), 
                        .combine = "+", .inorder=FALSE) %dopar% {
                            
                            out <- 
                                foreach(i = 1:ncol(x), .combine = "+") %do% {
                                    t_perm <- isTtest(arraydat, x[,i], 
                                                      id_dim, mu, var_equal, 
                                                      verbose = FALSE, any_NA)
                                    if (usetfce) {
                                        t_perm[] <- 
                                            apply(t_perm, applydims, 
                                                  tfceFn, chn = ChN, eh = EH)
                                    }
                                    t_perm <- colMaxs(
                                        mergeDims(abs(t_perm), 
                                                  list(1:2, seq_along(dim(t_perm))[-(1:2)]), 
                                                  keep_dimnames = FALSE, 
                                                  return_attributes = FALSE)
                                    )
                                    sweep(stat_obs, applydims, t_perm, "<")
                                }
                        }
            permres <- (permres + 1) / (nperm + 1)
        }
    } 
    # one-sample (can be difference of two paired samples as well)
    else {
        t_obs <- osTtest(arraydat, id_dim, mu, verbose, any_NA)
        if (nperm > 1) {
            idlen <- length(dimnames(arraydat)[[id_dim]])
            maxperm <- 2^idlen
            if (nperm > maxperm) 
                stop("nperm is too high for this sample size (not enough unique combinations exist)")
            if (!is.null(seed)) set.seed(seed)
            if (nperm > maxperm/3) {
                randi <- matrix(unlist(
                    expand.grid(rep(list(c(-1L, 1L)), idlen), 
                                KEEP.OUT.ATTRS = FALSE),
                    use.names = FALSE), idlen, maxperm, TRUE)
                randi <- randi[, sample(ncol(randi), nperm)]
            } 
            else {
                randi <- matrix(rep.int(1L, idlen))
                repeat{
                    nc <- ncol(randi) - 1
                    if (nc >= nperm) {
                        randi <- randi[, 2:(nperm + 1)]
                        break
                    } 
                    else {
                        ri <- matrixIP(sample(c(-1L, 1L), 
                                              temp <- idlen * (nperm - nc) * 2,
                                              TRUE),
                                       idlen, temp)
                        randi <- unique(cbind(randi, ri), MARGIN = 2L)
                    }
                }
            }
            applydims <- seq(3, length(dim(t_obs)))
            if (usetfce) {
                tfce_obs <- arrayIP(0, dim(t_obs), dimnames(t_obs))
                tfce_obs[] <- apply(t_obs, applydims, tfceFn, 
                                    chn = ChN, eh = EH)
                stat_obs <- abs(tfce_obs)
            } 
            else {
                stat_obs <- abs(t_obs)
            }
            chunks <- min(10L, ceiling(nperm/max(1L, length(cl))))
            permres <- 
                foreach(x = iter(randi, by = "column", 
                                 chunksize = chunks), 
                        .combine = "+", .inorder=FALSE) %dopar% {
                            out <- 
                                foreach(i = 1:ncol(x), .combine = "+") %do% {
                                    t_perm <- osTtest(arraydat * x[,i], 
                                                      id_dim, mu, 
                                                      verbose = FALSE, any_NA)
                                    if (usetfce) {
                                        t_perm[] <- 
                                            apply(t_perm, applydims, 
                                                  tfceFn, chn = ChN, eh = EH)
                                    }
                                    t_perm <- colMaxs(
                                        mergeDims(abs(t_perm), 
                                                  list(1:2, seq_along(dim(t_perm))[-(1:2)]), 
                                                  keep_dimnames=FALSE, 
                                                  return_attributes=FALSE)
                                    )
                                    sweep(stat_obs, applydims, t_perm, "<")
                                }
                        }
            permres <- (permres + 1) / (nperm + 1)
        }
    }
    # back-transform to original shape
    if (verbose) {
        setattr(t_obs, "Df", backFn(attr(t_obs, "Df"), origdimnames))
        setattr(t_obs, "pvalues", backFn(attr(t_obs, "pvalues"), origdimnames))
    }
    t_obs <- backFn(t_obs, origdimnames)
    out <- list(call = match.call(), effect_t_obs = t_obs)
    if (usetfce) {
        tfce_obs <- backFn(tfce_obs, origdimnames)
        out <- c(out, list(effect_tfce_obs = tfce_obs))
    }
    if (nperm > 0) {
        permres <- backFn(permres, origdimnames)
        out <- c(out, list(perm_pvalues = permres))
    }
    # return
    out
}

#
# <<< ANOVA on arrays >>> ----
#

#' Compute marginal means in an ANOVA design
#' 
#' \code{marginalMeans} computes marginal means in ANOVA designs
#' @param form formula of the model
#' @param f_dat data.frame of factors
#' @param a_dat matrix which contains the dependent variables. It must have as
#' many rows as f_dat.  Usually a_dat is the result of calling array2mat on the
#' orginal data array.
#' @param dimn the original dimension names for the array corresponding to the
#' column names of a_dat #### SHOULD BE SIMILAR TO MAT2ARRAY ###################
#' @param keep_term_order logival variable whether to keep the order of the 
#' model terms as provided in form (TRUE) or all interaction terms should follow
#' all main effect terms (FALSE, default)
#' @param residualmean logical variable; if TRUE, each factor term is extracted
#' from the data after computing the given marginal mean (default: FALSE)
#' @param whichterm character, numeric or logical vector; indices of model terms
#' which should be analyzed (default: NULL, meaning all terms are included)
#' @param no_icpt logical value; if TRUE, global mean is subtracted before 
#' computing the marginal means (default: FALSE)
#' @export
#' @return A named list containing the marginal means for each model terms
marginalMeans <- function(form, f_dat, a_dat, dimn, keep_term_order = FALSE, 
                          residualmean = FALSE, whichterm = NULL, no_icpt = FALSE) {
    labels <- attr(terms(as.formula(form), keep.order = keep_term_order), 
                   "term.labels")
    termL <- if (is.null(whichterm) || is.na(whichterm)) {
        strsplit(labels, ":")
    } else if (is.character(whichterm)) {    
        strsplit(whichterm[whichterm %in% labels], ":")
    } else {
        strsplit(labels[whichterm], ":")
    }
    marg_means <- vector("list", length(termL))
    if (residualmean && !no_icpt) a_dat <- sweep(a_dat, 2, colMeans(a_dat), "-")
    for (i in 1:length(termL)) {
        groups <- factor(interaction(f_dat[,termL[[i]]], drop = TRUE))
        groupfreq <- tabulate(groups)
        names(groupfreq) <- levels(groups)
        marg_means[[i]] <- arrayIP(rowsum(a_dat, groups)/groupfreq, 
                                   c(length(groupfreq), vapply(dimn, length, 0L)),
                                   c(list(modelterm = names(groupfreq)), dimn))
        setattr(marg_means[[i]], "freq", groupfreq)
        if (residualmean) {
            a_dat <- a_dat - c(marg_means[[i]][as.numeric(groups),])
        }
    }
    names(marg_means) <- labels
    # return
    marg_means
}



#' Compute Sum of Squares
#' @keywords internal
sumsq <- function(form, f_dat, a_dat, dimn, keep_term_order = FALSE, 
                  whichterm = NULL, no_icpt = FALSE, return_means = FALSE,
                  return_array = TRUE) {
    labels <- attr(terms(as.formula(form), keep.order = keep_term_order), 
                   "term.labels")
    termL <- if (is.null(whichterm) || is.na(whichterm)) {
        strsplit(labels, ":")
    } else if (is.character(whichterm)) {    
        strsplit(whichterm[whichterm %in% labels], ":")
    } else {
        strsplit(labels[whichterm], ":")
    }
    ssq <- matrixIP(0, length(termL), ncol(a_dat))
    df <- rep.int(0L, length(termL))
    names(df) <- labels
    if (return_means) {
        marg_means <- vector("list", length(termL))
        names(marg_means) <- labels
    }
    #
    if (!no_icpt) a_dat <- sweep(a_dat, 2, colMeans(a_dat), "-")
    #
    for (i in 1:length(termL)) {
        groups <- factor(interaction(f_dat[,termL[[i]]], drop = TRUE))
        groupfreq <- tabulate(groups)
        names(groupfreq) <- levels(groups)
        tmeans <- rowsum(a_dat, groups)/groupfreq
        a_dat <- a_dat - tmeans[as.numeric(groups),]
        ssq[i, ] <- colSums(groupfreq*tmeans^2)
        if (length(termL[[i]]) == 1) {
            df[i] <- length(groupfreq)-1
        } else {
            ind <- which(sapply(termL[1:(i-1)], 
                                function(x) all(x%in%termL[[i]])))
            df[i] <- prod(df[termL[[i]]])
        }
        if (return_means) {
            marg_means[[i]] <- array(tmeans, 
                                     c(length(groupfreq), vapply(dimn, length, 0L)),
                                     c(list(modelterm = names(groupfreq)), dimn))
        }
    }
    #
    if (return_array) {
        arrayIP(ssq, c(nrow(ssq), vapply(dimn, length, 0L)),
                c(list(modelterm = labels), dimn))
    } else {
        rownames(ssq) <- labels
        names(dimnames(ssq))[1] <- "modelterm"
    }
    setattr(ssq, "Df", df)
    if (return_means) setattr(ssq, "term_means", marg_means)
    # return
    ssq
}

#' Parameter checks and data preparation before arrayAnovaSub
#' @keywords internal
preAnova <- function(arraydat, factordef, bwdat, 
                     useparallel, ncores, par_method, usetfce = FALSE) {
    params <- list()
    # fast check of arguments
    if (is.null(factordef$w_id)) factordef$w_id <- "id"
    if (is.data.frame(arraydat)) arraydat <- as.matrix(arraydat)
    if (!is.array(arraydat) & !is.matrix(arraydat)) {
        stop("Improper input data class (should be data.frame, matrix or array))")
    }
    if (!is.null(factordef$between)) {
        if (!all(factordef$between %in% colnames(bwdat))) {
            stop("Between-subject factors and 'bwdat' dataset do not match.")
        }
        bwdat[,factordef$between] <- lapply(
            data.frame(bwdat[,factordef$between]), factor)
    }
    if (useparallel) {
        params$ncores <- ncores
        if (par_method == "snow") {
            fns <- ls(envir = parent.env(environment()))
            fns <- fns[sapply(fns, function(x) is.function(get(x)))]
            params$varlist2snow <-  
                if (usetfce) {
                    c(fns, "randind", "ChN","EH") 
                } else {
                    c(fns, "randind")
                }
        } 
    }
    #
    op <- options(contrasts = c("contr.sum", "contr.poly"))
    origdimnames <- dimnames(arraydat)
    modeldims <- c(factordef$w_id, factordef$within)
    keepdims <- setdiff(names(origdimnames), modeldims)
    arraydat <- mergeDims(arraydat, list(modeldims, keepdims))
    dat <- data.frame(matrix(unlist(strsplit(rownames(arraydat), "[.]"), 
                                    use.names = FALSE), 
                             nrow = nrow(arraydat), byrow = TRUE))
    colnames(dat) <- modeldims
    dat[,factordef$between] <- bwdat[
        match(dat[,factordef$w_id],bwdat[,factordef$w_id]),factordef$between]
    # return
    list(dat = dat, arraydat = arraydat, factordef = factordef, 
         origdimnames = origdimnames, par_params = params)
}

#' Workhorse function of arrayAnova
#' 
#' \code{arrayAnovaSub} computes F-values in ANOVA designs, and provides 
#' additional p-values, degrees of freedom, generalized effect sizes if 
#' requested.
arrayAnovaSub <- function(a_dat, f_def, d_names, f_dat, verbose = TRUE) {
    origdims <- vapply(d_names, length, 0L)
    modeldims <- c(f_def$w_id, f_def$within)
    keepdims <- setdiff(names(d_names), modeldims)
    #
    # distinct algorithms for between-subject, within-subject, and mixed models
    if (is.null(f_def$within)) {
        # SHOULD BE IMPROVED - use it with one between-subject factor and a 
        # balanced design
        aov_formula <- as.formula(paste(
            as.character(quote(a_dat)), "~", 
            paste(as.character(f_def$between), collapse = "*")))
        results <- summary(aov(aov_formula, data = f_dat))
        results <- arrayIP(
            unlist(results, use.names = FALSE),
            c(dim(results[[1]]), length(results)),
            list(modelterm = gsub(" ", "", rownames(results[[1]])),
                 measures = colnames(results[[1]]),
                 columns = gsub("[*Response ]", "", names(results))))
        termrows <- tolower(rownames(results)) != "residuals"
        Fvals <- arrayIP(results[termrows,"F value",], 
                         c(sum(termrows), origdims[keepdims]),
                         c(list(modelterm = rownames(results)[termrows]),
                           d_names[keepdims]))
        if (verbose) {
            Df <- as.numeric(results[termrows, "Df", 1])
            rDf <- results[!termrows, "Df", 1]
            pvalues <- results[termrows,"Pr(>F)",]
            # effect size (Olejnik and Algina, 2003; Bakeman, 2005)
            results <- subsetArray(results, list(measures = "Sum Sq"), drop = FALSE)
            matrixIP(results, nrow(results), dim(results)[3], 
                     dimnames = dimnames(results)[-2])
            SSr <- results["Residuals", ]
            if (!is.null(f_def$observed)) {
                ind <- apply(sapply(f_def$observed, 
                                    function(x) 
                                        grepl(x, rownames(results))),
                             1, any)
                SS1 <- colSums(results[ind, , drop = F])
                SS2 <- sweep(results, 1, ind, "*")
                ges <- results / sweep(sweep(results - SS2, 2, SSr, "+"), 2,
                                       SS1, "+")
            } else {
                ges <- results / sweep(results, 2, SSr, "+")
            }
            ges <- ges[-nrow(ges), ]
            factor_means <- marginalMeans(aov_formula, f_dat, a_dat, 
                                          d_names[keepdims], residualmean = FALSE)
        }
    } else if (is.null(f_def$between)) {
        # compute marginal means
        if (verbose) {
            aov_formula <- as.formula(paste(
                as.character(quote(a_dat)), "~", 
                paste(as.character(f_def$within), collapse = "*")
            )
            )
            factor_means <- marginalMeans(aov_formula, f_dat, a_dat, 
                                          d_names[keepdims], residualmean = FALSE)
        }
        # model formula
        aov_formula <- as.formula(paste(
            as.character(quote(a_dat)), "~", 
            paste(
                f_def$w_id, 
                paste(as.character(f_def$within), collapse = "*"),
                sep = "*")
        )
        )
        #
        # run model
        #
        a_datsc <- sweep(a_dat, 2, colMeans(a_dat), "-")
        model <- sumsq(aov_formula, f_dat, a_datsc, d_names[keepdims], 
                       no_icpt = TRUE, return_array = FALSE)
        model_df <- attr(model, "Df")
        # compute indices to handle error terms
        modnames <- dimnames(model)$modelterm
        m_ind <- !grepl(f_def$w_id, modnames)
        err_ind <- match(
            paste(f_def$w_id, ":", modnames[m_ind], sep = ""),
            modnames)
        Df <- model_df[m_ind]
        rDf <- model_df[err_ind]
        Fvals <- (model[m_ind, , drop = F] / Df) / 
            (model[err_ind, , drop = F] / rDf)
        if (verbose) {
            pvalues <- mapply(pf, Fvals, Df, rDf, MoreArgs = list(lower.tail = FALSE))
            dim(pvalues) <- dim(Fvals)
            # effect size (Olejnik and Algina, 2003; Bakeman, 2005)
            SSr <- colSums(model[grep(f_def$w_id, 
                                      rownames(model)), ])
            ges <- model[m_ind, , drop = F] / 
                sweep(model[m_ind, , drop = F], 2, SSr, "+")
        }
        arrayIP(Fvals, c(sum(m_ind), origdims[keepdims]),
                c(list(modelterm = modnames[m_ind]), d_names[keepdims]))
    } else {   
        aov_formula1 <- as.formula(paste(
            as.character(quote(a_dat)), "~", 
            paste(
                paste(as.character(f_def$between), collapse = "*"), 
                paste(as.character(f_def$within), collapse = "*"),
                sep = "*")
        )
        )
        aov_formula2 <- as.formula(paste(
            as.character(quote(a_dat)), "~", 
            paste(
                f_def$w_id, 
                paste(as.character(f_def$within), collapse = "*"),
                sep = "*")
        )
        )
        # compute marginal means
        if (verbose) 
            factor_means <- marginalMeans(aov_formula1, f_dat, a_dat, 
                                          d_names[keepdims], residualmean = FALSE)
        #
        # run model
        #
        a_datsc <- sweep(a_dat, 2, colMeans(a_dat), "-")
        model1 <- sumsq(aov_formula1, f_dat, a_datsc, d_names[keepdims], 
                        no_icpt = TRUE, return_array = FALSE)
        model1_df <- attr(model1, "Df")
        model2 <- sumsq(aov_formula2, f_dat, a_datsc, d_names[keepdims], 
                        no_icpt = TRUE, return_array = FALSE)
        model2_df <- attr(model2, "Df")
        # compute indices to handle error terms
        modnames <- dimnames(model2)$modelterm
        m2indices <- grep(f_def$w_id, modnames)
        modnames <- modnames[m2indices]
        bwnames <- unlist(lapply(seq_along(f_def$between), 
                                 function(i) apply(combn(f_def$between, i), 2, 
                                                   paste, collapse = ":")))
        sumindices <- sapply(bwnames, function(i) 
            match(sub(f_def$w_id, i, modnames), dimnames(model1)$modelterm))
        modnames <- sub(paste(f_def$w_id,"\\:?", sep = ""), "", modnames)
        fillindices <- lapply(1:nrow(sumindices), function(i)
            c(sumindices[i,], 
              na.omit(match(modnames[i], 
                            dimnames(model1)$modelterm))))
        # compute F values
        model1_err <- arrayIP(0, dim(model1))
        model1_err_df <- rep.int(0L, length(model1_df))
        for (i in 1:nrow(sumindices)) {
            errs <- model2[m2indices[i],] - colSums(model1[sumindices[i,], , drop = F])
            model1_err[fillindices[[i]], ] <- 
                rep(c(errs), each = length(fillindices[[i]])) 
            model1_err_df[fillindices[[i]]] <- model2_df[m2indices[i]] - 
                sum(model1_df[sumindices[i, ]])
        }
        Fvals <- model1/model1_df*model1_err_df/model1_err
        if (verbose) {
            Df <- as.numeric(model1_df)
            rDf <- model1_err_df
            pvalues <- mapply(pf, Fvals, Df, rDf, MoreArgs = list(lower.tail = FALSE))
            # effect size (Olejnik and Algina, 2003; Bakeman, 2005)
            SSr <- colSums(model1_err[grep(f_def$w_id, 
                                           rownames(model2)), ])
            if (!is.null(f_def$observed)) {
                ind <- apply(sapply(f_def$observed, 
                                    function(x) grepl(x, rownames(model1))),
                             1, any)
                SS1 <- colSums(model1[ind, ])
                SS2 <- sweep(model1, 1, ind, "*")
                ges <- model1 / sweep(sweep(model1 - SS2, 2, SSr, "+"), 2,
                                      SS1, "+")
            } else {
                ges <- model1 / sweep(model1, 2, SSr, "+")
            }
        }
        dimnms <- c(list(modelterm = rownames(Fvals)), d_names[keepdims])
        arrayIP(Fvals, c(nrow(Fvals), origdims[keepdims]), dimnms)
    }
    # rearrange dimensions to be compatible with arrayTtest output
    if (verbose) {
        tempfn <- function(x) aperm(array(x, dim(Fvals), dimnames(Fvals)), 
                                    c(keepdims, "modelterm"))
        pvalues <- tempfn(pvalues)
        ges <- tempfn(ges)
    }
    Fvals <- aperm(Fvals, c(keepdims, "modelterm"))
    if (verbose) {
        setattr(Fvals, "Df.term", Df)
        setattr(Fvals, "Df.resid", rDf)
        setattr(Fvals, "pvalues", pvalues)
        setattr(Fvals, "ges", ges)
        setattr(Fvals, "factor_means", factor_means)
    }
    # return
    Fvals
}

#' Perform ANOVA (potentially with TFCE correction) on arrays
#' 
#' \code{arrayAnova} performs point-to-point ANOVAs on arrays. Permutation-based
#' p-values and Threshold-free Cluster Enhancement (TFCE) correction can be requested.
#' @param arraydat a numeric array with named dimnames containing the EEG (or 
#' other) data. Missing values are not allowed.
#' @param factordef a named list of factor definitions, containing the following 
#' elements:
#' \itemize{
#' \item{between:}{character vector of between-subject factors (default: NULL)}
#' \item{within:}{character vector of within-subject factors (default: NULL)}
#' \item{w_id:}{name of the dimension which identifies the subjects 
#' (default: "id")}
#' }
#' @param bwdat a data.frame which contains the identification codes 
#' (factordef$w_id) and the group memberships (factordef$between) of the 
#' subjects. Missing values are not allowed.
#' @param verbose logical value indicating if p-values and effect sizes should 
#' be computed for the traditional ANOVA results
#' @param nperm integer value giving the number of permutations (default: 999L)
#' @param useparallel logical value; if TRUE (default), computations are done
#' in parallel
#' @param ncores integer value corresponding to the number of cores; 
#' if NULL (default), it is set to the maximum number of cores available
#' @param par_method parallelization method; can be "snow" (default) or "mc"
#' (multicore)
#' @param cl a cluster definition for snow-type parallelization; if NULL 
#' (default), it is set up automatically
#' @param usetfce logical value whether TFCE (threshold-free cluster enhancement)
#' correction should also be computed (default: TRUE)
#' @param tfce_options a named list containing the channel neighbourhood matrix
#' (named ChN) and the vector of the E/H parameters (named EH)
#' @param seed an integer value which specifies a seed (default: NULL)
#' @details The function assumes that the input array contains at least two 
#' named dimensions: chan (corresponding to the channels [electrodes]) and time 
#' (corresponding to time points). All dimensions which are not listed as 
#' within-subject factors are treated in a similar way as chan and time, that is 
#' separate ANOVA-s are computed for each level of those dimensions.
#' @note The function computes type I p-values - this is correct if the design
#' is fully balanced and orthogonal (if the number of between-subject 
#' factors is one, it may have slightly unequal group sizes).
#' @export
#' @return A list object with F values, TFCE-corrected F-values and
#' permutation-based p-values (if requested)
arrayAnova <- function(arraydat, factordef, bwdat = NULL, verbose = TRUE, 
                       nperm = 999L, useparallel = FALSE, ncores = NULL, 
                       par_method = "snow", cl = NULL,
                       usetfce = TRUE, tfce_options = NULL, seed = NULL) {
    # some checks
    stopifnot(is.array(arraydat))
    stopifnot(is.list(factordef))
    if (!is.null(factordef$between) && is.null(bwdat)) {
        stop("No between-participant data provided")
    } 
    if (usetfce) {
        if (is.null(tfce_options)) stop("Provide at least tfce_options$ChN")
        ChN <- tfce_options$ChN
        EH <- if (is.null(tfce_options$EH)) c(0.66, 2) else tfce_options$EH
    }
    if (useparallel) {
        if (is.null(ncores)) ncores <- detectCores()
    }
    #
    out <- list(call = match.call())
    #
    #
    temp <- preAnova(arraydat, factordef, bwdat, 
                     useparallel, ncores, par_method, usetfce)
    # assign variables to this environment, and potentially overwrite existing ones
    assignList(temp, verbose = FALSE)
    rm(temp)
    #
    Fvals_obs <- arrayAnovaSub(arraydat, factordef, origdimnames, 
                               dat, verbose)
    if (usetfce) {
        mergedimnames <- setdiff(names(dimnames(Fvals_obs)), c("chan", "time"))
        tfce_obs <- mergeDims(Fvals_obs, mergedimnames)
        for (i in 1:nrow(tfce_obs)) 
            tfce_obs[i,,] <- tfceFn(tfce_obs[i,,], ChN, EH)
        tfce_obs <- revMergeDims(tfce_obs)
    }
    if (nperm > 1) {
        #
        permfn_f <- function(i) {
            x <- arrayAnovaSub(arraydat[randind[i,], ], 
                               factordef, origdimnames, dat, FALSE)
            x <- mergeDims(x, mergedimnames)
            maxf <- sapply(1:nrow(x), 
                           function(ii) max(x[ii,,])) 
            return(maxf)
        }
        permfn_tfce <- function(i) {
            verb <- ifelse(i == 1, verbose, FALSE)
            x <- arrayAnovaSub(arraydat[randind[i,], ], 
                               factordef, origdimnames, dat, FALSE)
            x <- mergeDims(x, mergedimnames)
            maxf <- sapply(1:nrow(x), 
                           function(ii) max(abs(tfceFn(x[ii,,], ChN, EH))))
            return(maxf)
        }
        permfn <- if (usetfce) permfn_tfce else permfn_f
        #
        # generate random orders (dim(randind) = nperm X nrow(dat))
        if (!is.null(seed)) set.seed(seed)
        randind <- shuffleSet(nrow(dat), nperm, 
                              how(within = Within(type = "free"), 
                                  plots = Plots(strata = dat[,factordef$w_id], 
                                              type = "free")))
        #
        if (!useparallel) {
            maxperm <- lapply(1:nperm, permfn)
        } else if (par_method == "snow") {
            if (is.null(cl)) cl <- makePSOCKcluster(par_params$ncores, 
                                                    outfile = "")
            clusterExport(cl, 
                          varlist = par_params$varlist2snow,
                          envir = environment())
            maxperm <- parLapply(cl, 1:nperm, permfn)
            stopCluster(cl)
            rm(cl)
        } else if (par_method == "mc") {
            maxperm <- mclapply(1:nperm, permfn, mc.cores = par_params$ncores)
        } 
        maxperm <- matrixIP(unlist(maxperm, use.names = FALSE), ncol = nperm)
        sig <- if (usetfce) tfce_obs else Fvals_obs
        sig <- mergeDims(sig, mergedimnames)
        for (i in 1:nrow(sig)) {
            sig[i,,] <- (rowSums(outer(c(sig[i,,]), maxperm[i,], "<="))+1) / 
                (nperm + 1)
        }
        sig <- revMergeDims(sig)
    } 
    out$effect_F_obs <- Fvals_obs
    if (usetfce) out$effect_tfce_obs <- tfce_obs
    if (nperm > 1L) out$perm_pvalues <- sig
    # return
    out
}

#' Extract interaction
#' 
#' \code{extractInteraction} computes the highest-order (interaction) effect 
#' (i.e. the difference between differences) from means of each cell of the 
#' design matrix
#' @param dat a named list with array elemens, see details below
#' @param sep the character which separates the name of factor levels in the 
#' model terms (default: ".")
#' @param sep_fixed logical indicating as sep should be treated "as is" while 
#' splitting the model term names
#' @details This function is not intended for direct use, and it only works for 
#' specifically formatted data. The input data is supposed to be the return 
#' value of \code{\link{marginalMeans}}; a named list which contains the 
#' marginal means for each effect of an ANOVA. The names are the names of the 
#' factor(s) - for interaction effects the factor names are separated by ":". 
#' The list elements are three-dimensional arrays with the following dimension 
#' names: modelterm, chan and time. For pure main effects, this function 
#' computes all pairwise differences between the levels of the factor. For 
#' interaction effects, the function computes all pairwise differences at the 
#' level of the highest order interaction. That means for the A:B interaction 
#' where factor A has levels Aa & Ab and factor B has levels Ba, Bb and Bc, 
#' the function computes the following differences: (Aa-Ab) - (Ba-Bb), 
#' (Aa-Ab) - (Ba-Bc), (Aa-Ab) - (Bb-Bc).
#' @export
#' @return A named list is returned with the same names as the input, but the 
#' list elements are ia_level x chan x time arrays.  
extractInteraction <- function(dat, sep = ".", sep_fixed = TRUE) {
    # function to compute the weights (+1 or -1) of the factor levels (f)
    iasign <- function(f) {
        combin <- lapply(f, function(x) combn(levels(x), 2))
        combin_ind <- expand.grid(lapply(combin, 
                                         function(x) seq.int(ncol(x))))
        combin_names <- expand.grid(lapply(f, function(x) 
            combn(levels(x), 2, FUN = paste, collapse = "-")))
        combin_names <- apply(combin_names, 1, paste, collapse = ".")
        out <- sapply(1:nrow(combin_ind), function(ii) {
            out <- rep(0, nrow(f))
            ind <- sapply(names(combin), function(x) 
                f[, x] %in% combin[[x]][, combin_ind[ii, x]])
            ind <- apply(ind, 1, all)
            if (sum(ind) > 2) {
                tempf <- sapply(droplevels(f[ind, , drop = F]), as.numeric)
                tempf <- sign(1.5 - tempf)
                out[which(ind)] <- apply(tempf, 1, prod)
            } else {
                out[which(ind)] <- c(1, -1)
            }
            return(out)
        })
        colnames(out) <- combin_names
        return(out)
    }
    # function to compute the interaction effects
    compia <- function(x, modterm_name) {
        modterm <- dimnames(x)$modelterm
        facs <- strsplit(modterm, sep, fixed = sep_fixed)
        facs <- data.frame(matrix(unlist(facs), 
                                  ncol = length(facs[[1]]), 
                                  byrow = TRUE))
        colnames(facs) <- unlist( strsplit(modterm_name, ":") )
        signs <- iasign(facs)
        rownames(signs) <- modterm
        avg_ia <- fnDims(x, "modelterm", 
                         function(xx, y) t( crossprod(xx, y) ), 
                         list(y = signs),
                         newdim = list(ia_level = colnames(signs)),
                         vectorized = TRUE)
        return( avg_ia )
    }
    # ------
    # run computations
    out <- mapply(compia, dat, names(dat), SIMPLIFY = FALSE)
    # return
    out
}



#
# <<< channel positions >>> --------
#

#' Transformation between spherical, cartesian and geographical coordinates
#' 
#' \code{sph2cart}, \code{sph2geo}, \code{cart2sph}, \code{cart2geo}, 
#' \code{geo2sph}, \code{geo2cart} transform spherical (sph), cartesian (cart), 
#' or geographical (geo) coordinates into the respective coordinate system.
#' @name coordinates
#' @param ch_pos a data frame or matrix containing the spherical, cartesian,
#' or geographical coordinates of the electrode positions. It should contain 
#' at least the following (named) columns, respectively: theta and phi; 
#' x, y, and z; or lat and long.
#' @param r radius (default = 1)
#' @param deg logical variable indicating whether spherical or geographical 
#' coordinates are or should be given in degrees (TRUE, default)
#' @param long360 logical variable; if TRUE (default), longitudes range from
#' 0 to 360 (otherwise from -180 to 180)
#' @param orient a character value specifying the orientation in the 
#' geographical coordinate system; "northpole" (default) or "equatorial"
#' @return A data.frame with converted coordinates
NULL

#' @rdname coordinates
#' @export
# Transform cartesian to spherical coordinates
sph2cart <- function(ch_pos, r = 1, deg = TRUE) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("theta", "phi") %in% colnames(ch_pos)))
        stop("Either theta or phi angles are missing!")
    if (!is.null(ch_pos$r)) r <- ch_pos$r
    if (deg) 
        ch_pos[, c("theta", "phi")] <- ch_pos[, c("theta", "phi")]/180*pi
    theta <- ch_pos$theta
    phi <- ch_pos$phi
    out <- data.frame(
        x = r * sin(theta) * cos(phi),
        y = r * sin(theta) * sin(phi),
        z = r * cos(theta)
    )
    rownames(out) <- rownames(ch_pos)
    # return
    out
}

#' @rdname coordinates 
#' @export
# Transform cartesian to spherical coordinates
cart2sph <- function(ch_pos, deg = TRUE) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("x", "y", "z") %in% colnames(ch_pos)))
        stop("Either x, y, or z coordinates are missing!")
    x <- ch_pos$x
    y <- ch_pos$y
    z <- ch_pos$z
    r <- sqrt(x^2 + y^2 + z^2)
    theta <- sign(x) * abs(acos(z/r))
    phi <- abs(atan(y/x))
    phi[y == x & x == 0] <- 0
    phi <- phi * sign(x) * sign(y) 
    out <- data.frame(theta, phi, r)
    if (nrow(out) > 1) rownames(out) <- rownames(ch_pos)
    if (deg) out[, 1:2] <- out[, 1:2] * 180/pi
    # return
    out
}

#' @rdname coordinates
#' @export
# Transform spherical to geographical coordinates
sph2geo <- function(ch_pos, r = 1, deg = TRUE, long360 = TRUE, 
                    orient = c("northpole", "equatorial")) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("theta", "phi") %in% colnames(ch_pos)))
        stop("Either theta or phi angles are missing!")
    if (!is.null(ch_pos$r)) r <- ch_pos$r
    if (!deg) ch_pos <- ch_pos * 180/pi
    orient <- match.arg(orient)
    if (orient == "equatorial") {
        temp <- sph2cart(ch_pos)[, c("z", "x", "y")]
        colnames(temp) <- c("x", "y", "z")
        temp <- cart2sph(temp)
        ch_pos$theta <- temp$theta
        ch_pos$phi <- temp$phi
    }
    theta <- ch_pos$theta
    phi <- ch_pos$phi
    long <- ifelse(theta < 0, 180, 0) + phi
    long[long < 0] <- 360 + long[long < 0]
    if (!long360) long[long > 180] <- long[long > 180] - 360
    lat <- 90 - abs(theta)
    out <- data.frame(long = long, lat = lat, r = r)
    rownames(out) <- rownames(ch_pos)
    # return
    out
}

#' @rdname coordinates
#' @export
# Transform geographical to spherical coordinates
geo2sph <- function(ch_pos, r = 1, deg = TRUE) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("lat", "long") %in% colnames(ch_pos)))
        stop("Either latitude [lat] or longitude [long] coordinates are missing!")
    if (!is.null(ch_pos$r)) r <- ch_pos$r
    neglong <- ch_pos$long < 0
    long360 <- if (any(neglong)) F else T
    theta <- (90 - ch_pos$lat) * sign(90 - ch_pos$long) * sign(270 - ch_pos$long)
    theta <- ifelse(ch_pos$lat == 45 & theta == 0, 45, theta)
    phi <- -(ifelse(theta < 0, 180, 0) - ch_pos$long)
    phi <- ifelse(phi > 180, phi - 360, phi)
    out <- data.frame(r = r, theta = theta, phi = phi)
    backcheck <- sum(abs(as.matrix(
        sph2geo(out, long360 = long360)[, c("long", "lat")] - 
            ch_pos[, c("long", "lat")]))) < 1e-8
    if (backcheck) {
        if (!deg) out <- out/180 * pi
        return(out)
    } else {
        stop("Oops, something went wrong with the transformation and I don't know why.")
    }
}

#' @rdname coordinates
#' @export
# Transform geographical to cartesian coordinates
geo2cart <- function(ch_pos, r = 1, deg = TRUE) {
    sph2cart( geo2sph(ch_pos, r, deg) )
}

#' @rdname coordinates
#' @export
# Transform cartesian to geographical coordinates
cart2geo <- function(ch_pos, deg = TRUE) {
    sph2geo( cart2sph(ch_pos, deg) )
}


#' Find channel neighbours
#' 
#' \code{chanNb} finds neighbouring channels
#' @param ch_pos electrode positions
#' @param check_alpha a two-element numeric vector defining the range which is 
#' supposed to contain the optimal value of alpha
#' @param alpha a numeric value which influences the allowed distance between
#' neighbouring electrodes; if other than NULL (the default), check_alpha is 
#' ignored
#' @param ... parameters to \code{\link{sph2cart}}
#' @export
#' @return An electrode neighbourhood matrix
chanNb <- function(ch_pos, check_alpha = c(0.1, 10), alpha = NULL, ...) {
    options(rgl.useNULL = TRUE)
    reqFn(c("alphashape3d", "geometry"))
    if (!all(c("x", "y", "z") %in% colnames(ch_pos))) {
        if (all(c("theta", "phi") %in% colnames(ch_pos))) {
            ch_pos <- sph2cart(ch_pos, ...)
        } else {
            stop("Channel coordinates should be spherical (polar) or cartesian coordinates!")
        }
    } 
    ch_pos <- ch_pos[, c("x", "y", "z")]
    channames <- paste(1:nrow(ch_pos), rownames(ch_pos), sep = ". ")
    if (is.null(alpha)) {
        require("shiny")
        require("shinyRGL")
        alpha <- runApp(list(
            ui = pageWithSidebar(
                # Application title
                headerPanel("Find channel neighbours"),
                # Sidebar with a slider input for number of points
                sidebarPanel(
                    sliderInput("alpha",
                                "Alpha value: ",
                                min = min(check_alpha), 
                                max = max(check_alpha), 
                                value = 1, step = 0.1),
                    actionButton("submit", "Use selected alpha")
                ),
                # Show the generated 3d scatterplot
                mainPanel(
                    webGLOutput("chanPlot", height = "700px")
                )
            ),
            server = function(input, output) {
                output$chanPlot <- renderWebGL({
                    #bg3d("grey40")
                    a <- suppressWarnings(alphashape3d::ashape3d(as.matrix(ch_pos), 
                                                   input$alpha, pert = TRUE))
                    plot(a, walpha = TRUE, transparency = 0.95, 
                         col = c("red", "red2", "red"), shininess = 100)
                    rgl::text3d(ch_pos$x, ch_pos$y, ch_pos$z, cex = 0.6, 
                                channames, col = "black")
                })
                observe({
                    if (input$submit == 0)
                        return()
                    stopApp(input$alpha)
                })
            }
        ))
    }
    a <- suppressWarnings(alphashape3d::ashape3d(as.matrix(ch_pos), 
                                   alpha, pert = TRUE))$edge
    a <- a[,c(1, 2, ncol(a))]
    out <- matrixIP(0, nrow(ch_pos), nrow(ch_pos))
    rownames(out) <- rownames(ch_pos)
    mirror <- out > 0
    for (i in 1:nrow(ch_pos)) {
        nb <- sort(unique(c(
            which(mirror[i,]), 
            i, 
            a[a[,1]==i & a[,3]>0, 2])))
        out[i, 1:length(nb)] <- nb
        mirror[nb, i] <- T
    }
    out <- out[, apply(out > 0, 2, any)]
    # return
    out
}

#' Cosine of the angles between electrodes
#' 
#' \code{cosAngle} computes the cosine of angles between electrodes or 
#' between any N-dimensional vectors, or between two sets of electrodes or 
#' any N-dimensional vectors. 
#' @param x a matrix or data.frame with at least two columns and rows
#' @param y If not NULL (default), a matrix or data.frame with at least two 
#' columns and rows
#' @param coords a logical variable; if TRUE (default), data in x are 
#' coordinates. Coordinates should be named as "x", "y", "z" or "theta", "phi", 
#' and if names are present, only coordinate columns are used in the 
#' computations. 
#' @param units_in_rows a logical variable indicating whether the units (vectors)
#' are in the rows of x (TRUE, default) or in the columns.
#' @param check_params logical; if TRUE (default), the appropriateness of 
#' input data are checked before computing the cosine of angles. Set it to 
#' FALSE if you really know what you are doing. 
#' @export
#' @return A symmetric matrix  
cosAngle <- function(x, y = NULL, coords = TRUE, units_in_rows = TRUE, 
                     check_params = TRUE) {
    checkFn <- function(xx) {
        if (!units_in_rows) xx <- t(xx)
        if (length(dim(xx)) <= 1 || ncol(xx)<2) {
            stop("Provide a matrix or data.frame with at least 2 columns.")
        }
        if (coords) {
            if (!is.null(colnames(xx))) {
                if (!all(c("x", "y", "z") %in% colnames(xx)) && 
                        all(c("theta", "phi") %in% colnames(xx))) {
                    xx <- sph2cart(xx)
                }
                ind <- na.omit(match(c("x", "y", "z"), colnames(xx)))
                if (length(ind) < 2) {
                    stop("Provide a matrix or data.frame with x and y coordinates.")
                } else {
                    xx <- xx[, ind]
                }
            }
        }
        if (is.data.frame(xx)) xx <- as.matrix(xx) 
        return( xx )
    }
    if (is.null(y)) {
        if (check_params) {
            x <- checkFn(x)
        }
        len <- sqrt(rowSums(x^2))
        out <- (x %*% t(x)) / outer(len, len)
        rownames(out) <- colnames(out) <- rownames(x)
    } else {
        if (check_params) {
            x <- checkFn(x)
            y <- checkFn(y)
        }
        lenx <- sqrt(rowSums(x^2))
        leny <- sqrt(rowSums(y^2))
        out <- (x %*% t(y)) / outer(lenx, leny)
        rownames(out) <- rownames(x)
        colnames(out) <- rownames(y)
    }
    # return
    out
}


#' Spherical spline interpolation
#' 
#' \code{chanInterp} performs spherical spline interpolation. It can impute 
#' bad channels (i.e. channels with missing values) or interpolate to 
#' artbitrary positions on the scalp.
#' @param dat a numeric vector, matrix, data.frame or array with named dimnames
#' @param ch_pos channel (electrode) positions; should be a matrix or data.frame
#' with the following column names (order is not important): "theta", "phi" 
#' (spherical coordinates), or "x", "y", "z" (cartesian coordinates)
#' @param interp_pos interpolation positions, if they differ from ch_pos; 
#' should be in the same format as ch_pos (default: NULL)
#' @param maxNA numeric value given as ratio (0 < maxNA < 1) or integer 
#' (maxNA >= 1 ); it determines the maximum ratio or number of non-missing 
#' channels in a time sample. If exceeded, no interpolation occurs in the given
#' time sample.
#' @param m integer value (default: 4); the flexibility of the spline
#' @param N integer value (default: 7); number of terms in the Legendre 
#' polynomial. You should probably increase it (even up to 100L) for 
#' realistic (i.e. non-spherical) electrode arrangements.
#' @param lambda numeric value (default: 1e-10); smoothing factor
#' @param type character value, do not use yet
#' @param alarm_tolerance numeric value (default: 1e-2); if the maximal absolute 
#' interpolation error at any time sample exceeds this limit, a message is
#' shown or an error is thrown depending on \code{error_on_alarm}. If set to 
#' NULL, no check is performed.
#' @param error_on_alarm defaults to TRUE, see \code{alarm_tolerance}
#' @export
#' @return An object having the same attributes as dat
#' @import orthopolynom
#' @import corpcor
#' @import Kmisc
chanInterp <- function(dat, ch_pos, interp_pos = NULL, maxNA = 0.3, 
                       m = 4L, N = 7L, lambda = 1e-10, 
                       type = c("voltage", "laplacian", "scd"),
                       alarm_tolerance = 1e-2, error_on_alarm = TRUE) {
    require(orthopolynom)
    require(corpcor)
    require(Kmisc)
    message("\n****\nStart interpolation...\n")
    # Function to compute G matrices:
    # G = g(cos(ch_pos)) / interpG = g(cos(ch_pos, interp_pos))
    compGmat <- function(cos_angles, m, N) {
        n <- seq.int(N)
        series <- (2*n + 1) / (n^m * (n+1)^m)
        k <- 1/4 / pi
        P <- polynomial.values(
            legendre.polynomials(N, normalized = FALSE)[-1], cos_angles)
        G <- k * Reduce("+", mapply("*", P, series, SIMPLIFY = FALSE))
        return( G )
    }
    # Function to compute C coefficients
    compCmat <- function(x, G, tol) {
        Gx <- arrayIP(1, dim(G) + 1)
        Gx[-1, -1] <- G
        Gx[1, 1] <- 0
        diag(Gx) <- diag(Gx) + tol
        invG <- pseudoinverse(Gx, tol)[, -1]
        return( invG %*% x )
    }
    # Function to perform spherical spline interpolation
    interpFn <- function(y, m, N, lambda, type) {
        na_ind <- is.na(y)
        message("...Looking for missing value patterns - ", appendLF = FALSE)
        if (is.null(interp_pos)) {
            temp <- colMeans(na_ind)
            keep_columns <- temp <= maxNA & temp > 0
            if (!any(keep_columns)) return(y)
            yorig <- y
            y <- y[, keep_columns, drop = FALSE]
            na_ind <- na_ind[, keep_columns, drop = FALSE]
            missing_patterns <- fastUnique(na_ind, units_in_rows = FALSE)
            mischan <- which(rowAnys(missing_patterns))
            mischan <- 
                if (!is.null(rownames(ch_pos))) {
                    paste(rownames(ch_pos)[mischan], collapse = ", ")
                } else {
                    paste(seq_along(mischan), collapse = ", ")
                }
            message("Done")
            message("...There are missing values in the following channels: ",
                    mischan)
        } else {
            missing_patterns <- fastUnique(na_ind, units_in_rows = FALSE)
            message("Done")
            yy <- arrayIP(0, c(nrow(interp_pos), ncol(y)),
                          list(rownames(interp_pos), colnames(y)))
            names(dimnames(yy)) <- names(dimnames(y))
        }
        if (!is.null(alarm_tolerance)) {
            dev <- rep(FALSE, ncol(y))
            names(dev) <- colnames(y)
            maxdev <- 0
        }
        message("...Perform interpolation - ", appendLF = FALSE) 
        G0 <- compGmat(cosAngle(ch_pos, check_params = FALSE), m, N)
        for (i in 1:ncol(missing_patterns)) {
            na_vec <- missing_patterns[, i]
            y_ind <- colSums(na_ind==na_vec) == length(na_vec)
            ch_good <- ch_pos[!na_vec, , drop = F]
            ch_interp <- 
                if (is.null(interp_pos)) ch_pos[na_vec,,drop = F] else interp_pos
            #
            G <- G0[!na_vec, !na_vec]
            Coef <- compCmat(y[!na_vec, y_ind, drop = F], G, lambda)
            interpG <- compGmat(cosAngle(ch_good, ch_interp, check_params = FALSE), 
                                m, N)
            if (is.null(interp_pos)) {
                res <- crossprod(Coef[-1, , drop = F], interpG) + Coef[1, ]
                y[na_vec, y_ind] <- t(res)
            } else {
                yy[, y_ind] <- t(crossprod(Coef[-1, , drop = F], interpG) + 
                                     Coef[1, ])
            }
            if (!is.null(alarm_tolerance)) {
                ch_interp <- ch_pos
                interpG <- compGmat(cosAngle(ch_good, 
                                             ch_interp, check_params = FALSE), 
                                    m, N)
                res <- crossprod(Coef[-1, , drop = F], interpG) + Coef[1, ]
                mdevs <- colMaxs(abs(y[, y_ind] - t(res)), na.rm = TRUE)
                maxdev <- max(c(maxdev, mdevs))
                dev[y_ind] <- temp <- ( maxdev > alarm_tolerance )
                if (error_on_alarm && temp) 
                    stop(sprintf("Deviance %f exceeds threshold %f", 
                                 maxdev, alarm_tolerance))
            }
        }
        if (!is.null(alarm_tolerance)) {
            if (any(dev)) {
                devind <- which(dev)
                devind <- devind[1:min(c(20, length(devind)))]
                message(
                    "****\nThe maximum absolute interpolation error exceeds the limit 
                    at the following time samples (only the first 20 are shown):\n-----\n", 
                    paste(names(dev)[devind], collapse = " ") 
                )
                message(paste("Maximal deviation = ", maxdev, sep = ""))
            }
        }
        message("Done") 
        if (is.null(interp_pos)) {
            yorig[, keep_columns] <- y
            return(yorig)
        } else {
            return(yy)
        }
    }
    #
    type <- match.arg(type)
    if (is.null(lambda)) {
        lambda <- if (type == "voltage") 1e-7 else 1e-5
    }
    #
    # check if ch_pos has an appropriate format
    if (!all(c("x", "y", "z") %in% colnames(ch_pos))) {
        if (all(c("theta", "phi") %in% colnames(ch_pos))) {
            ch_pos <- sph2cart(ch_pos)
        } else {
            stop("Provide channel locations in spherical or cartesian coordinates!")
        }
    }
    ch_pos <- as.matrix(ch_pos[, c("x", "y", "z")])
    #
    # do the same for the interpolation locations
    if (!is.null(interp_pos)) {
        if (!all(c("x", "y", "z") %in% colnames(interp_pos))) {
            if (all(c("theta", "phi") %in% colnames(interp_pos))) {
                interp_pos <- sph2cart(interp_pos)
            } else {
                stop("Provide interpolation locations in spherical or cartesian coordinates!")
            }
        }
        interp_pos <- as.matrix(interp_pos[, c("x", "y", "z")])        
    } 
    #
    # set up interpolation infos
    argnames <- setdiff(names(as.list(args(chanInterp))), "dat")
    argnames <- argnames[argnames  != ""]
    procstep <- list(
        what = "interpolation", call = match.call(),
        options = mget(argnames), nr_of_missings = sum(is.na(dat)))
    #
    # check input data and run interpolation
    dim_names <- dimnames(dat)
    if (is.list(dat) || !is.numeric(dat)) {
        stop("Provide a numeric vector, matrix, data.frame or array as input!")
    } else if (is.null(interp_pos) && 
                   procstep$nr_of_missings == 0) {
        setattr(dat, "processing_steps",
                c(attr(dat, "processing_steps"), list(procstep)))
        message("...No missing data to interpolate.\nInterpolation finished.")
        return(dat)
    } else if (is.vector(dat)) {
        if (length(dat) != nrow(ch_pos)) {
            stop("Number of channels and number of data points do not match!")
        }
        out <- interpFn(as.matrix(dat), m, N, lambda, type)
    } else {
        if (!any(dim(dat) == nrow(ch_pos))) {
            stop("Number of channels and data size do not match!")
        }
        target_dim <- if (!is.null(dim_names$chan)) "chan" else 1
        arg_list <- list(m = m, N = N, lambda = lambda, type = type)
        if (length(arg_list) == 0) arg_list <- NULL
        out <- fnDims(dat, "chan", interpFn, arg_list = arg_list, 
                      vectorized = TRUE)
        out <- aperm(out, names(dim_names))
    }
    #
    if (is.null(interp_pos)) {
        NAs_left <- sum(is.na(out))
        procstep$nr_of_interp <- 
            procstep$nr_of_missings - NAs_left
        message("Number of NAs after bad channel interpolation: ", 
                NAs_left)
    }
    message("Interpolation finished.")
    setattr(out, "processing_steps", 
            c(attr(dat, "processing_steps"), list(procstep)))
    # return
    out
}

#' Project channel positions onto 2D plane
#' 
#' \code{project3dMap} projects 3D coordinates onto a 2D plane or vice versa
#' @param pos electrode positions (matrix or data.frame)
#' @param r radius
#' @param projection character string (default = "laea"). See 
#' \link{http://www.remotesensing.org/geotiff/proj_list/} for common projections
#' @param projref projection reference (pole [ = default] or equator)
#' @param origo a named character vector of lat ( = latitude) and long (longitude)
#' @param inverse if set to TRUE, back-projection is performed (default = FALSE)
#' @import proj4
#' @export
#' @return A data.frame containing the projected coordinates
project3dMap <- function(pos, r = 1,
                         projection = "laea", projref = c("pole", "equator"),
                         origo = c(lat = ifelse(projref == "pole", 90, 0),
                                 long = ifelse(projref == "pole", 270, 0)), 
                         inverse = FALSE) {
    #
    xyfn <- function(dat) data.frame(x = dat[,1], y = dat[,2])
    projfn <- function(p, projcall, ...) { 
        res <- lapply(unique(projcall), function(x) {
            ind <- projcall == x
            out <- project(as.data.frame(p[ind, , drop = F]), proj = x, ...)
            out <- as.data.frame(out)
            rownames(out) <- which(ind)
            return( out )
        })
        res <- Reduce("rbind", res)
        res <- res[order(as.numeric(rownames(res))), , drop = F]
        return( xyfn(res) )
    }
    #
    if (!is.data.frame(pos) & is.list(pos)) {
        stop("Pos can be a vector, matrix, or data.frame.")
    }
    if (is.atomic(pos) && is.vector(pos)) {
        if (length(pos)  != 2) {
            stop("If pos is a vector, it must contain exactly two elements!")
        }
        pos <- cbind(x = pos[1], y = pos[2])
    }
    pos <- as.data.frame(pos)
    if (is.null(pos$r)) pos$r <- 1
    projref <- match.arg(projref)
    proj4call <- paste0("+proj=", projection, 
                       " +lat_0=", origo[1], 
                       " +lon_0=", origo[2],
                       paste0(" +R=", pos$r))
    #
    if (inverse) {
        if (is.null(colnames(pos))) colnames(pos)[1:2] <- c("x", "y")
        if (!all(c("x", "y") %in% colnames(pos))) {
            stop("Provide a matrix or data.frame with x and y coordinates!")
        }
        mapxy <- projfn(as.matrix(pos[, c("x", "y")]), proj4call, inverse = TRUE) 
        colnames(mapxy) <- c("long", "lat")
    } else {
        if (!all(c("long", "lat") %in% colnames(pos))) {
            if (!all(c("theta", "phi") %in% colnames(pos))) 
                pos <- cart2sph(pos)
            if (projref == "equator") {
                temp <- sph2cart(pos)
                colnames(temp) <- c("y", "z", "x")
                pos <- cart2sph(temp)
            }
            pos <- sph2geo(pos, long360 = FALSE)
        } 
        mapxy <- projfn(pos[, c("long", "lat")], proj4call) 
    }
    rownames(mapxy) <- rownames(pos)
    # return
    mapxy
}

#
# <<< ERP functions >>> --------
#

# import functions =========== 

#' Import options for importBVdat
#' 
#' \code{importOptions} allows to set the raw data importation parameters
#' @param eeg_ext character value, the extension of the segmented EEG dataset 
#' file (default: "dat")
#' @param marker_ext character value, the extension of the raw marker file
#' (default: "vmrk")
#' @param info_ext character value, the extension of the information file 
#' (default: "vhdr")
#' @param marker_skip numeric value; number of rows to be skipped while 
#' importing the marker file (default: 14)
#' @param marker_segment character vector; name(s) of stimuli defining the
#' the segments (default: "Target")
#' @param marker_badcat character string(s) identifying bad segments in the 1st 
#' column of the marker file (default: "Bad Interval")
#' @param marker_badstim character string(s) identifying bad stimuli or 
#' responses in the 2nd column of the marker file (default: "resp_false")
#' @param keepstim character vector identifying which stimuli or responses 
#' should be kept by checking the 2nd column of the marker file (default: "*")
#' @param marker_badchan character vector of length 1 or 2; markers identifying 
#' bad intervals of channels. If marker_badchan has only one element, "_start" 
#' and "_stop" are automatically appended to the end of the string.
#' @param marker_time0 character string identifying the marker name at time 0 
#' in the imported EEG dataset
#' @param segment_dpoints numeric vector of data point indices defining a 
#' segment (default: -100:1023)
#' @param marker_regexp logical value; should marker_segment/badcat/.../time0 
#' strings be handled as regular expressions (default) or treated as they are
#' @param marker_ignorecase logical value; should the case of marker_segment/
#' .../time0 definitions be ignored (default)
#' @param {marker_header,marker_fill,marker_asis,marker_sep} parameters 
#' to be passed to \code{\link{read.table}} while importing the marker file
#' @export
#' @return A list with named parameters
importOptions <- function(eeg_ext = "dat", marker_ext = "vmrk", info_ext = "vhdr",
                          marker_skip = 14, marker_segment = "Target", 
                          marker_badcat = "Bad Interval", 
                          marker_badstim = "resp_false",
                          marker_keepstim = "*", 
                          marker_badchan = c("badchan_start",
                                             "badchan_stop"),
                          marker_time0 = "Time 0",
                          segment_dpoints = (-100):1023,
                          marker_regexp = TRUE,
                          marker_ignorecase = TRUE,  
                          marker_header = FALSE, marker_fill = TRUE, 
                          marker_asis = TRUE, marker_sep = ",") {
    #
    if (length(marker_badchan) == 1) {
        marker_badchan <- paste(marker_badchan, c("start", "stop"), sep = "_")
    }
    # return
    mget(ls())
}

#' Import binary file exported from BrainVision
#' 
#' \code{importBVdat} imports a binary file exported from the BrainVision 
#' software
#' @param file_name character string; the name of the input files without
#' extensions
#' @param file_path character string; the path to the files if they are not in 
#' the working directory (default)
#' @param id character string denoting the identification code of the 
#' participant
#' @param import_options a list, which should be given by calling 
#' \code{\link{importOptions}}
#' @note This is a custom function tailored for the special datasets collected 
#' in our lab. Use it with extra care for general purposes!
#' @export
#' @return A list with three named elements: eeg (array), markers (data.frame), 
#' channels (data.frame)
importBVdat <- function(file_name, file_path = getwd(), id = "",
                        import_options = importOptions()) {
    message(paste("\nImport ", file_name, "...", sep = ""))
    mygrepl <- function(patterns, ...) {
        if (is.null(patterns)) {
            rep(FALSE, length(list(...)$x))
        } else {
            rowSums(sapply(patterns, 
                           grepl, 
                           ignore.case = marker_ignorecase, 
                           fixed = !marker_regexp,
                           ...)) > 0
        }
    }
    extractInfo <- function(x, type = "char") {
        out <- tolower(strsplit(info[grep(x, info)], "=")[[1]][2])
        if (type == "num") out <- as.numeric(out)
        return(out)
    }
    chanInfo <- function() {
        chan <- strsplit(info[grep("Ch.{1,3}=", info)], "=")
        chan_names <- sapply(chan[1:(length(chan)/2)],
                             function(x) strsplit(x[[2]], ",")[[1]][1])
        chan_pos <- t(sapply(chan[-(1:(length(chan)/2))],
                             function(x) as.numeric(strsplit(x[[2]], ",")[[1]])))
        colnames(chan_pos) <- c("r", "theta", "phi")
        return( data.frame(chan_pos, row.names = chan_names) )
    }
    #
    assignList(import_options, verbose = FALSE)
    #
    # info from vhdr file
    con1 <- file(file.path(file_path, paste(file_name, info_ext, sep = ".")))
    info <- readLines(con1)
    close(con1)
    eeg_orientation <- extractInfo("DataOrientation=")
    nr_chan <- extractInfo("NumberOfChannels=", "num")
    nr_tpoints <- extractInfo("SegmentDataPoints=", "num")
    eeg_length <- extractInfo("DataPoints=", "num") * nr_chan
    eeg_Hz <- 1e6 / extractInfo("SamplingInterval=", "num")
    chan <- chanInfo()
    #
    # import eeg data
    if (!is.null(eeg_ext)) {
        eeg <- readBin(file.path(file_path, paste(file_name, eeg_ext, sep = ".")),
                       what = "integer", n = eeg_length)
        eeg <- eeg / 1000
        if (eeg_orientation == "vectorized") {
            eeg <- t(matrixIP(eeg, eeg_length/nr_chan, nr_chan))
        } else {
            matrixIP(eeg, nr_chan, eeg_length/nr_chan)
        }
        setattr(eeg, "subject_id", as.character(id))
        procstep <- list(
            what = "import", call = match.call(),
            file_name = file_name, file_path = file_path,
            options = import_options)
    } else {
        eeg <- NULL
    }
    #
    # import markers
    markers.orig <- read.table(
        file.path(file_path, 
                  paste(file_name, import_options$marker_ext, sep = ".")),
        header = marker_header, fill = marker_fill, 
        as.is = marker_asis, sep = marker_sep, skip = marker_skip)
    markers.orig[,3] <- suppressWarnings(as.numeric(markers.orig[, 3]))
    markers.orig[,5] <- suppressWarnings(as.integer(markers.orig[, 5]))
    markers.orig <- markers.orig[!is.na(markers.orig[, 3]), 1:5]
    markers.orig$segmind <- findInterval(
        markers.orig[,3],
        markers.orig[grepl("New Segment", 
                           markers.orig[, 1]), 3])
    if (eeg_length/nr_chan/nr_tpoints  != max(markers.orig$segmind))
        stop("Marker and data files do not match!")
    #
    # remove data in bad channel intervals
    if (!is.null(marker_badchan)) {
        chan_names <- rownames(chan)
        bch <- markers.orig[, 2]
        ind <- markers.orig[, 5] > 0
        bch[ind] <- paste(bch[ind],
                          chan_names[markers.orig[ind, 5]],
                          sep = "_")
        bch <- paste0(bch, "_")
        badchanFn <- function(i) {
            indstart <- grep(paste(marker_badchan[1], ".*_", 
                                   chan_names[i], "_", sep = ""), 
                             bch)
            indstart <- unique(markers.orig[indstart, 3])
            indstop <- grep(paste(marker_badchan[2], ".*_", 
                                  chan_names[i], "_", sep = ""), 
                            bch)
            indstop <- unique(markers.orig[indstop, 3])
            len1 <- length(indstart)
            len2 <- length(indstop)
            if (len1 == 0 & len2 == 0) {
                return(NULL)
            } else {
                lendiff <- len1 - len2
                if (lendiff == 1) {
                    message(
                        paste0("One badchan_stop marker was missing at channel ", chan_names[i], 
                              " and was automatically set to the last sampling point.")
                    )
                    indstop <- c(indstop, ncol(eeg))
                } else if (lendiff == -1) {
                    message(
                        paste0("One badchan_start marker was missing at channel ", chan_names[i],
                              " and was automatically set to the first sampling point.")
                    )
                    indstart <- c(1, indstart)    
                } else if (abs(lendiff) > 1) {
                    stop(
                        paste0("Badchan_start and badchan_stop markers do not match at channel ", chan_names[i], ".")
                    )
                }
                return(cbind(start = indstart, stop = indstop))
            }
            if (any((indstop - indstart) < 0)) {
                stop(
                    paste0("Badchan_start and badchan_stop markers do not match at channel ", chan_names[i], ".")
                )
            }
        }
        badchans <- lapply(seq_along(chan_names), badchanFn)
        names(badchans) <- chan_names
        badchans <- badchans[!sapply(badchans, is.null)]
        if (!is.null(eeg_ext)) {
            for (i in names(badchans)) {
                eeg[i, unlist(mapply(seq, 
                                     badchans[[i]][, "start"], 
                                     badchans[[i]][, "stop"], 
                                     SIMPLIFY = FALSE), 
                              use.names = FALSE)] <- NA
            }
        }
        procstep$bad_channels <- badchans
    } 
    # 
    # remove bad segments
    markers.orig$badsegm <- rep(
        sapply(
            split(markers.orig, markers.orig$segmind), 
            function(x) 
                any(mygrepl(marker_badcat, x = x[, 1]) | 
                        mygrepl(marker_badstim, x = x[, 2]))
        ),
        table(markers.orig$segmind))
    keepind <- 
        mygrepl(marker_segment, x = markers.orig[, 1]) & 
        mygrepl(marker_keepstim, x = markers.orig[, 2]) &
        rep(markers.orig[mygrepl("Time 0", x = markers.orig[, 1]), 3] %in% 
                markers.orig[mygrepl(marker_time0, x = markers.orig[, 1]), 3],
            table(markers.orig$segmind)) &
        !markers.orig$badsegm
    markers <- markers.orig[keepind, ]
    markers <- data.frame(
        segment = sapply(strsplit(markers[, 1], "="), "[[", 2),
        fullcode = markers[, 2],
        dpoint = markers[, 3],
        segmind = markers$segmind
    )
    #
    # format eeg data
    if (!is.null(eeg_ext)) {
        eeg <- eeg[, c(outer(segment_dpoints, markers$dpoint,"+"))]
        arrayIP(eeg, 
                c(nr_chan, length(segment_dpoints), nrow(markers)),
                list(chan = rownames(chan),
                     time = segment_dpoints * 1000 / eeg_Hz,
                     trial = paste(markers$segment, 
                                   markers$fullcode, 
                                   sep = "_")))
    }
    # decorate
    setattr(eeg, "processing_steps", list(procstep))
    # return
    list(eeg = eeg, markers = markers, channels = chan)
}

#' Split concatenated strings. 
#' 
#' \code{splitMarker} is a conveniance function wrapping strsplit. It returns
#' a data.frame with columns of customizable classes. Useful for post-processing
#' markers.
#' @param marker a vector which contains concatenated elements
#' @param header a character vector which determines the column names of the
#' resulting data.frame
#' @param type a character vector which determines the class of each column in
#' the resulting vector. If its length is less then the number of columns, it
#' will be recycled. If set to NULL (default), character vectors are transformed to
#' factors.
#' @param splitchar a character vector of length one indicating the splitting
#' character (default: _). Splitting characters which are special characters in
#' R (e.g. "|", ".", etc.) should be given as in strsplit (e.g. "\\\|")
#' @export
#' @return A data.frame with header and appropriate classes, containing the 
#' spltted substrings of the original vector elements
#' @seealso \code{\link{strsplit}}
splitMarker <- function(marker, header, type = NULL, splitchar = "_") {
    out <- data.frame(
        matrix(unlist(strsplit(as.character(marker), splitchar), 
                      use.names = FALSE), 
               nrow = length(marker), ncol = length(header), byrow = TRUE))
    colnames(out) <- header
    if (!is.null(type)) {
        if (length(type) < length(header)) {
            type <- rep(type, length.out = length(header))
        }
        for (i in 1:length(type)) {
            if (is.factor(out[, i])) {
                out[, i] <- do.call(paste("as", type[i], sep = "."), 
                                    list(as.character(out[, i])))
            } else {
                out[, i] <- do.call(paste("as", type[i], sep = "."), list(out[, i]))
            }
        }
    }
    # return
    out
}

# base functions =========== 

#' Baseline correction 
#' 
#' \code{baselineCorr} performs baseline correction
#' @param dat numeric array containing the ERPs
#' @param basedim character vector identifying the dimensions of dat along which
#' separate baseline averaging should occur (default: "chan")
#' @param baseind numeric or character vector identifying time points which
#' form the baseline
#' @export
#' @return A numeric array with the same attributes as dat
baselineCorr <- function(dat, basedim = "chan", baseind = NULL) {
    message("\n****\nPerform baseline correction ... ", appendLF = FALSE)
    origattr <- attributes(dat)
    if (is.null(baseind)) baseind <- as.numeric(dimnames(dat)$time)<0
    if (length(basedim) == 1) {
        means <- fnDims(subsetArray(dat, list(time = baseind)),
                        basedim, rowMeans, list(na.rm = TRUE), vectorized = TRUE, columnwise = FALSE)
        out <- sweep(dat, which(names(dimnames(dat)) == basedim), means, "-")
    } else {
        dat <- mergeDims(dat, basedim)
        means <- fnDims(subsetArray(dat, list(time = baseind)),
                        names(dimnames(dat))[1], rowMeans, list(na.rm = TRUE), 
                        vectorized = TRUE, columnwise = FALSE)
        out <- arrayIP(
            sweep(dat, 1, means, "-"),
            c(vapply(origattr$dimnames[basedim], length, 0L), dim(dat)[-1]),
            c(origattr$dimnames[basedim], 
              origattr$dimnames[setdiff(names(origattr$dimnames), basedim)]))
        out <- aperm(out, names(origattr$dimnames))
    }
    attributes(out) <- origattr
    setattr(out, "processing_steps",
            c(attr(out, "processing_steps"),
              list(list(what = "baseline correction", 
                        call = match.call(), base_dimensions = basedim, 
                        base_indices = baseind))))
    message("Done")
    # return
    out
}

#' Options for artifact rejection
#' 
#' \code{artrejOptions} allows to set the parameters of the artifact rejection
#' methods.
#' @param sampling_freq numeric value, the sampling frequency of the EEG-data
#' @param channels character vector containing the name or index of channels
#' which are subject to artifact rejection. If set to "all" (default), all 
#' channels are included.
#' @param apply_maxgrad logical value, if set to TRUE (default), the maximum 
#' gradient criterion is applied.
#' @param maxgrad_limit numeric value, the maximum gradient / millisecond 
#' (default: 50)
#' @param maxgrad_mark numeric vector of length 2; the placement of the Bad
#' Interval mark in milliseconds before and after the occurence of maxgrad_limit
#' violation (default: c(-200, 200))
#' @param apply_diffrange logical value, if set to TRUE (default), the difference 
#' range criterion is applied.
#' @param diffrange_limit numeric vector of length 2, the minimum and maximum
#' voltage difference in a given interval (default: 200)
#' @param diffrange_mark numeric vector of length 2; the placement of the Bad
#' Interval mark in milliseconds before and after the occurence of diffrange_limit
#' violation (default: c(-200, 200))
#' @param diffrange_interval numeric value, the length of interval for the 
#' difference range criterion in milliseconds (default: 200)
#' @param apply_amplrange logical value, if set to TRUE (default), the amplitude 
#' range criterion is applied.
#' @param amplrange_limit numeric vector of length 2, the minimum and maximum 
#' voltage in the whole segment (default: c(-200, 200))
#' @param amplrange_mark numeric vector of length 2; the placement of the Bad
#' Interval mark in milliseconds before and after the occurence of 
#' amplrange_limit violation (default: c(-200, 200))
#' @details The short definitions of the possible artifact rejection criteria 
#' are as follows:
#' \itemize{
#' \item{"Maximum gradient:"}{"The absolute difference between the voltages 
#' measured at successive milliseconds."}
#' \item{"Difference range:"}{"The minimum and maximum difference between the 
#' maximum and minimum voltages in a given sampling interval."}
#' \item{"Amplitude range:"}{"The minimum and maximum voltages in the segments."}
#' }
#' @note The algorithm takes care of the sampling frequency for all parameters
#' which are provided in milliseconds (or /ms) and makes adjustments if needed.
#' However, *_mark parameters are not used since only segmented data can be
#' analyzed in the present version of artifactRejection().
#' @export
#' @return A list object with all parameters.
artrejOptions <- function(
    sampling_freq = 1000, channels = "all", do_reject = TRUE,
    apply_maxgrad = TRUE, maxgrad_limit = 50, maxgrad_mark = c(-200, 200),
    apply_diffrange = TRUE, diffrange_limit = c(0.5, 100), diffrange_mark = c(-200, 200), 
    diffrange_interval = 200,
    apply_amplrange = TRUE, amplrange_limit = c(-200, 200), amplrange_mark = c(-200, 200)) {
    #
    opt <- mget(setdiff(ls(), "sampling_freq"))
    freqmod <- sampling_freq/1000
    ind <- grep("mark|interval|limit", names(opt))
    lapply(opt[ind], function(x) is.numeric(x))
    opt[ind] <- lapply(opt[ind], sort)
    ind <- grep("mark|interval", names(opt))
    opt[ind] <- lapply(opt[ind], "*", freqmod)
    opt$maxgrad_limit <- opt$maxgrad_limit / freqmod
    # return
    opt
}

#' Artifact rejection
#' 
#' \code{artifactRejection} performs artifact rejection on segmented data.
#' @param dat numeric array (EEG-data) with the following named dimensions 
#' (dimension order does not matter): chan, time, trial
#' @param markers if not NULL (default), a matrix or data.frame containing the 
#' characteristics of the trials (markers)
#' @param artrej_options a named list containing the parameters for the 
#' artifact rejection criteria. See \code{\link{artrejOptions}} for details.
#' @param return_data logical value, if TRUE (default), dat and markers without
#' rejected trials are returned
#' @param return_details logical value, if TRUE (default), the full array of 
#' results (e.g., bad trials for each channel and for each criterion) is 
#' returned as an attribute of bad_trials (see Values section)
#' @param print_result logical value, if TRUE (default), a summary of the 
#' results is printed to the console
#' @export
#' @return A named list containing bad_trials (trials identified with artifacts)
#' and the modified input data (dat and markers without contaminated trials)
artifactRejection <- function(dat, markers = NULL, artrej_options = artrejOptions(),
                              return_data = TRUE, return_details = TRUE,
                              print_result = TRUE) {
    # main function
    aRej <- function(x, details = return_details) {
        assignList(artrej_options, verbose = FALSE)
        crits <- sub("apply_", "", 
                     names(artrej_options)[grep("apply", 
                                                names(artrej_options))])
        if (length(crits) == 0) {
            stop("No criterion to apply; check the apply_* parameters in artrej_options!")
        }
        attribs <- attr(x, "array_attributes")
        row_dim <- attribs$row_dim
        dims <- attribs$dim[-row_dim]
        dimn <- attribs$dimnames[-row_dim]
        names(dims) <- names(dimn)
        out <- matrixIP(FALSE, ncol(x), length(crits))
        colnames(out) <- crits
        message("\n****\nStart artifact rejection / Criterion: ...\n")
        if (apply_maxgrad) {
            message("... Maximum gradient - ", appendLF = FALSE)
            out[, "maxgrad"] <- colAnys( abs(x[-1,]-x[-nrow(x),]) > maxgrad_limit )
            message(" done\n")
        }
        if (apply_diffrange) {
            message("... Difference range - ", appendLF = FALSE)
            ind <- rollFun(x, diffrange_interval, max, endrule = "NA") - 
                rollFun(x, diffrange_interval, min, endrule = "NA")
            ind <- ind[!rowAlls(is.na(ind)), ]
            out[, "diffrange"] <- colAnys( abs(ind) < diffrange_limit[1] ) |
                colAnys( abs(ind) > diffrange_limit[2] )
            message(" done\n")
        }
        if (apply_amplrange) {
            message("... Amplitude range - ", appendLF = FALSE)
            out[, "amplrange"] <- ( colMaxs(x) > amplrange_limit[2] |
                                        colMins(x) < amplrange_limit[1] )
            message(" done\n")
        }
        arrayIP(out, c(dims, ncol(out)), c(dimn, list(crit = crits)))
        artrej_summary <- matrixIP(0, dims["chan"]+1, length(crits)+1,
                                   list(chan = c(dimn$chan, "all"),
                                        crit = c(crits, "all")))
        dimres <- dim(artrej_summary)
        artrej_summary[-dimres[1], -dimres[2]] <- avgDims(out, "trial")
        artrej_summary[dimres[1], -dimres[2]] <- 
            colMeans(apply(out, c("trial", "crit"), any))
        artrej_summary[-dimres[1], dimres[2]] <- 
            colMeans(apply(out, c("trial", "chan"), any))
        out.details <- out
        out <- apply(out, "trial", any)
        artrej_summary[dimres[1], dimres[2]] <- mean( out )
        setattr(out, "summary", artrej_summary)
        if (details) setattr(out, "details", out.details)
        return( out )
    }
    # input data check
    if (is.null(dimnames(dat)) || is.null(names(dimnames(dat))) ||
            !identical(sort(names(dimnames(dat))), 
                       c("chan", "time", "trial"))) {
        stop("Provide EEG data as an array with the following named 
             dimensions (dimension order does not matter): 
             chan, time, trial.")
    }
    if (is.null(markers)) {
        markers <- data.frame(fullcode = seq_along(dimnames(dat)$trial))
    }
    if (length(dimnames(dat)$trial) != nrow(markers)) {
        stop("EEG and marker data contain different number of trials!")
    }
    if (!is.list(artrej_options)) {
        stop("The artrej_options parameter must be a list; provide it through 
             artrejOptions() to avoid inconsistent results!")
    }
    keepchan <- 
        if (identical(artrej_options$channels, "all")) {
            dimnames(dat)$chan
        } else {
            artrej_options$channels
        }
    # run artifact rejection
    tempdat <- array2mat(subsetArray(dat, list(chan = keepchan), drop = FALSE),
                         "time", keep_dimnames = FALSE)
    badtrials <- aRej(tempdat)
    if (print_result) {
        cat("\n----- Proportion of bad trials -----\n")
        print(attr(badtrials, "summary"))
        cat("\n------------------------------------\n")
    }
    if (return_data) {
        dat <- subsetArray(dat, list(trial = which(!badtrials)))
        setattr(dat, "processing_steps",
                c(attr(dat, "processing_steps"),
                  list(list(
                      what = "artifact rejection",
                      call = match.call(), results = badtrials, 
                      options = artrej_options)))
        )
        markers <- droplevels( markers[!badtrials, ] )
        # return
        list(bad_trials = badtrials, eeg = dat, markers = markers)
    } else {
        setattr(badtrials, "options", artrej_options)
        # return
        list(bad_trials = badtrials, eeg = NULL, markers = NULL)
    }
}

#' Compute Global Field Power
#'
#' \code{compGfp} computes Global Field Power (the standard deviation of 
#' channel values for each sampling point)
#' @param dat numeric matrix or array with named dimensions (one of which
#' must be "chan")
#' @param keep_channels logical value; if TRUE, the original channels are 
#' retained, and the GFP values are added with a channel code of "GFP"
#' (default = FALSE)
#' @export
#' @return The function returns a matrix or an array. Note that if
#' the original channels are not retained, the channel dimension is dropped.
compGfp <- function(dat, keep_channels = FALSE) {
    dat <- aperm(dat, c("chan", setdiff(names(dimnames(dat)), "chan")))
    dat.mat <- matrix(dat, nrow = nrow(dat))
    tempn <- colSums(!is.na(dat.mat))
    out <- sweep(dat.mat, 2, colMeans(dat.mat, na.rm = TRUE), "-")
    out <- sqrt(colSums(out^2, na.rm = TRUE)/tempn)
    if (keep_channels) {
        out <- rbind(dat.mat, out)
        dims <- dim(dat)
        dims[1] <- dims[1]+1
        dims.n <- dimnames(dat)
        dims.n[[1]] <- c(dims.n[[1]], "GFP")
    } else {
        dims <- dim(dat)[-1]
        dims.n <- dimnames(dat)
        dims.n[[1]] <- NULL
    }
    arrayIP(out, dims, dims.n)
    # return
    out
}

#' Scale channels
#' 
#' \code{scaleChan} normalizes the data across channels so that for
#' each sampling point, the mean of the channel amplitudes is zero and
#' the standard deviation is one.
#' @param dat numeric matrix or array with names dimensions, one of which
#' must be "chan"
#' @param keep_dimorder logical value; if TRUE (default), the order of the 
#' dimensions is kept intact, otherwise the channel dimension will be 
#' the first dimension in the resulting matrix or array
#' @export
#' @return A numeric matrix or array
scaleChan <- function(dat, keep_dimorder = TRUE) {
    dimnames.orig <- names(dimnames(dat))
    dat <- aperm(dat, c("chan", setdiff(names(dimnames(dat)), "chan")))
    dim.dat <- dim(dat)
    dimnames.dat <- dimnames(dat)
    dat <- matrix(dat, nrow = nrow(dat))
    tempn <- colSums(!is.na(dat))
    dat <- sweep(dat, 2, colMeans(dat, na.rm = TRUE), "-")
    dat <- sweep(dat, 2, sqrt(colSums(dat^2, na.rm = TRUE)/tempn), "/")
    dat <- array(dat, dim.dat, dimnames = dimnames.dat)
    if (keep_dimorder) dat <- aperm(dat, dimnames.orig)
    # return
    dat
}

#' Compute centroids
#'
#' \code{centroid} computes the centroids (separately for the negative and
#' positive values)
#' @param dat numeric vector of amplitudes at a given sampling point
#' @param ch_pos matrix or data.frame containing channel positions (in the same
#' order as dat)
#' @param proj2map logical value; if TRUE (default), the centroids are projected
#' onto a 2D plane
#' @param proj_unitsphere logical value; if TRUE, and proj2map is also TRUE, the 
#' projection assumes unit radius
#' @param ... additional parameters to \code{\link{project3dMap}}
#' @export
#' @return A data.frame containing the positions of the negative and positive 
#' centroids 
centroid <- function(dat, ch_pos, proj2map = TRUE, proj_unitsphere = FALSE, ...) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("x", "y", "z") %in% colnames(ch_pos))) {
        ch_pos <- sph2cart(ch_pos)
    }
    ind.neg <- dat < 0
    ind.pos <- dat >= 0
    centroids <- as.data.frame(
        rbind(colMeans(ch_pos[ind.neg, c("x", "y", "z")]),
              colMeans(ch_pos[ind.pos, c("x", "y", "z")]))
    )
    rownames(centroids) <- c("negative", "positive")
    if (proj2map) {
        pos <- cart2sph(centroids)
        if (proj_unitsphere) {
            pos$r <- 1
        }
        out <- project3dMap(pos, ...)
    }
    # return
    out
}

#' Average single-trials
#' 
#' \code{avgTrials} performs single-trial averaging
#' @param dat numeric array, the segmented ERPs
#' @param markers data.frame containing marker definitions
#' @param which_factors numeric or character vector indicating which columns
#' from markers should be used for averaging. If NULL (default), all columns
#' are used.
#' @export
#' @return numeric array
avgTrials <- function(dat, markers, which_factors = NULL) {
    stopifnot(is.numeric(dat))
    stopifnot(!is.null(markers))
    message("\n****\nAverage single-trials ... ", appendLF = FALSE)
    if (is.character(which_factors)) {
        if (!all(which_factors %in% colnames(markers))) {
            stop("markers and which_factors do not match")
        }
    } else if (is.numeric(which_factors)) {
        if (!all(which_factors %in% (1:ncol(markers)))) {
            stop("markers and which_factors do not match")
        }
    }
    if (!is.null(which_factors)) markers <- markers[, which_factors]
    groups <- interaction(markers, drop = TRUE, sep = "|")
    out <- fnDims(dat, "trial", 
                  function(x, g) sweep(rowsum(x, g), 1, table(g), "/"), 
                  list(g = groups), vectorized = TRUE, 
                  newdims = list(factor_level = levels(groups)))
    tempfac <- strsplit(dimnames(out)$factor_level, "\\|")
    setattr(out, "processing_steps",
            c(attr(out, "processing_steps"),
              list(list(
                  what = "averaging",
                  call = match.call(),
                  factors = as.data.frame(
                      matrix(unlist(tempfac, use.names = FALSE), 
                             nrow = length(tempfac), 
                             ncol = length(tempfac[[1]]), 
                             byrow = TRUE,
                             dimnames = list(1:length(tempfac), 
                                             colnames(markers))))
                  ))
            )
    )
    message("Done")
    # return
    out
}

# TANOVA functions =========== 

#' Topographical ANOVA (TANOVA) and related methods
#' 
#' \code{tanova} performs point-to-point topographical ANOVA on arrays. Related
#' methods are GFP-analysis (which is based on the intensity of the signal) and
#' DISS-analysis (focusing only on the global dissimilarity by comparing 
#' normalized topographies). See References for further details. 
#' @param arraydat a numeric array with named dimnames containing the EEG (or 
#' other) data. Missing values are not allowed.
#' @param factordef a named list of factor definitions, containing the following 
#' elements:
#' \itemize{
#' \item{"between"}{"character vector of between-subject factors (default: NULL)"}
#' \item{"within"}{"character vector of within-subject factors (default: NULL)"}
#' \item{"w_id"}{"name of the dimension which identifies the subjects 
#' (default: "id")}
#' }
#' @param bwdat a data.frame which contains the identification codes 
#' (factordef$w_id) and the group memberships (factordef$between) of the 
#' subjects. Missing values are not allowed.
#' @param type a character value of "tanova" (default), "dissimilarity", or 
#' "gfp"
#' @param verbose logical value indicating whether means for each factor level 
#' should be also returned
#' @param nperm integer value giving the number of permutations (default: 999L)
#' @param useparallel logical value; if TRUE (default), computations are done
#' in parallel
#' @param ncores integer value corresponding to the number of cores; 
#' if NULL (default), it is set to the maximum number of cores available
#' @param par_method parallelization method; can be "snow" (default) or "mc"
#' (multicore)
#' @param cl a cluster definition for snow-type parallelization; if NULL 
#' (default), it is set up automatically
#' @param iaterms_last logival variable whether all interaction terms should 
#' follow all main effect terms (default: TRUE)
#' @param seed an integer value which specifies a seed (default: NULL)
#' @param pcrit the significance level for duration and global count p-value 
#' correction (default: 0.05)
#' @details The function assumes that the input array contains at least two 
#' named dimensions: chan (corresponding to the channels [electrodes]) and time 
#' (corresponding to time points). All dimensions which are not listed as 
#' within-subject factors are treated in a similar way as chan and time, that is 
#' separate TANOVA-s are computed for each level of those dimensions.
#' @note The function computes type I p-values - this is correct if the design
#' is fully balanced and orthogonal (if the number of between-subject 
#' factors is one, it may have slightly unequal group sizes).
#' @export
#' @return A list object with effect statistics, uncorrected p-values, and 
#' two types of corrected p-values: length correction or global correction
#' @references Koenig (2011) RAGU
tanova <- function(arraydat, factordef, bwdat = NULL, 
                   type = c("tanova", "dissimilarity", "gfp"), 
                   verbose = TRUE, nperm = 999L, useparallel = FALSE, ncores = NULL,
                   par_method = c("snow", "mc"), cl = NULL,
                   iaterms_last = TRUE, seed = NULL, pcrit = 0.05) {
    # some checks
    stopifnot(is.array(arraydat))
    stopifnot(is.list(factordef))
    if (!is.null(factordef$between) && is.null(bwdat)) {
        stop("No between-participant data provided")
    }
    if (useparallel) {
        if (is.null(ncores)) ncores <- parallel::detectCores()
    }
    #
    out <- list(call = match.call())
    #
    effSize <- function(means) {
        datfn <- if (type == "gfp") abs else compGfp
        out <- lapply(means, function(x) avgDims(datfn(x), "modelterm"))
        out <- rearrangeList(out, "modelterm")
        out <- drop(out)
        return( out )
    }
    # 
    par_method <- match.arg(par_method)
    #
    type <- match.arg(type)
    if (type == "gfp") {
        if (verbose) {
            arraydat_orig <- arraydat
            arraydat_orig <- preAnova(arraydat, factordef, bwdat, 
                                      useparallel, ncores, par_method)$arraydat
        }
        arraydat <- compGfp(arraydat)
    }
    if (type == "dissimilarity") arraydat <- scaleChan(arraydat)
    #
    temp <- preAnova(arraydat, factordef, bwdat, useparallel, ncores, par_method)
    # assign variables to this environment, and potentially overwrite existing ones
    assignList(temp, verbose = FALSE)
    rm(temp)
    # formula of anova
    aov_formula <- as.formula(paste(
        as.character(quote(arraydat)), "~", 
        sub("^\\*", "", paste(
            paste(as.character(factordef$between), collapse = "*"), 
            paste(as.character(factordef$within), collapse = "*"),
            sep = "*"))
    )
    )
    # compute marginal means
    origdims <- vapply(origdimnames, length, 0L)
    modeldims <- c(factordef$w_id, factordef$within)
    keepdims <- setdiff(names(origdimnames), modeldims)
    term_means <- marginalMeans(aov_formula, dat, arraydat, 
                                origdimnames[keepdims], 
                                keep_term_order = !iaterms_last, 
                                residualmean = TRUE)
    es_obs <- effSize(term_means)
    out <- c(out, list(effect = es_obs, residual_means = term_means))
    names(out)[1] <- paste(names(out)[1], type, sep = "_")
    if (verbose) {
        if (type == "gfp") {
            factor_means <- marginalMeans(aov_formula, dat, arraydat_orig, 
                                          origdimnames[keepdims], 
                                          keep_term_order = !iaterms_last, 
                                          residualmean = FALSE)
            
        } else {
            factor_means <- marginalMeans(aov_formula, dat, arraydat, 
                                          origdimnames[keepdims], 
                                          keep_term_order = !iaterms_last, 
                                          residualmean = FALSE)
        }
        out <- c(out, list(factor_means = factor_means))
    }
    if (nperm > 1) {
        #
        permfn <- function(i) {
            x <- marginalMeans(aov_formula, dat, arraydat[randind[i,], ], 
                               origdimnames[keepdims], 
                               keep_term_order = !iaterms_last, 
                               residualmean = TRUE) 
            return( effSize(x) )
        }
        #
        # generate random orders (dim(randind) = nperm X nrow(dat))
        if (!is.null(seed)) set.seed(seed)
        randind <- shuffleSet(nrow(dat), nperm, 
                              how(within = Within(type = "free"), 
                                  plots = Plots(strata = dat[,factordef$w_id], 
                                              type = "free")))
        #
        if (!useparallel) {
            es_perm <- lapply(1:nperm, permfn)
        } else if (par_method == "snow") {
            if (is.null(cl)) cl <- makePSOCKcluster(par_params$ncores, outfile = "")
            clusterExport(cl, 
                          varlist = par_params$varlist2snow,
                          envir = environment())
            es_perm <- parLapply(cl, 1:nperm, permfn)
            stopCluster(cl)
            rm(cl)
        } else if (par_method == "mc") {
            es_perm <- mclapply(1:nperm, permfn, mc.cores = par_params$ncores)
        } 
        # permuted F values
        es_perm <- c(list(es_obs), es_perm)
        names(es_perm) <- seq.int(nperm + 1)
        es_perm <- rearrangeList(es_perm, "perm")
        #if (!is.null(perm_fname)) save(es_perm, file = perm_fname)
        # p-values
        pvalues_perm <- fnDims(es_perm, "perm", colRanks, 
                               arg_list = list(preserveShape = TRUE),
                               vectorized = TRUE, keep_dimorder = TRUE)
        pvalues_perm <- (nperm + 2 - pvalues_perm) / (nperm + 1)
        pvalues <- avgDims(subsetArray(pvalues_perm, 
                                       list(perm = 1), drop = FALSE),
                           "perm")
        pvalues_perm <- (pvalues_perm < pcrit)
        # consecutive sign. criterion
        pvalues_consec <- fnDims(pvalues, "time", I)
        pvalues_maxconsec <- mergeDims(pvalues_perm, 
                                       setdiff(names(dimnames(pvalues_perm)), 
                                               c("time", "perm")))
        pvalues_maxconsec <- apply(pvalues_maxconsec, 1, function(x) {
            res <- fnDims(x, "time", matrixRle, vectorized = TRUE)
            res <- quantile(res$lengths[res$values > 0], 1 - pcrit)
            return( res )
        })
        temp <- matrixRle(
            array2mat(pvalues, "time", return_attributes = FALSE,
                      keep_dimnames = FALSE) < pcrit)
        ind <- ( (temp$lengths < pvalues_maxconsec[temp$matrixcolumn]) &
                     (temp$values == 1) )
        temp$values[ind] <- 0
        temp <- inverse.matrixRle(temp)
        pvalues_consec[temp == 0] <- 1
        pvalues_consec <- aperm(pvalues_consec, names(dimnames(pvalues))) 
        # number of sign. time points criterion
        pvalues_global <- fnDims(pvalues_perm, "time", colSums, vectorized = TRUE)
        tempfn <- function(x) colSums(sweep(x, 2, x[1, ], ">="))
        pvalues_global <- fnDims(pvalues_global, "perm", tempfn, vectorized = TRUE)
        pvalues_global <- pvalues_global / (nperm + 1)
        # return
        out <- c(out, list(perm_pvalues = pvalues,
                           perm_pvalues_consec = pvalues_consec,
                           perm_pvalues_global = pvalues_global))
    } 
    # return
    out
}

# Peak Anova functions =========== 

#' ANOVA on individual peaks
#' 
#' \code{peakAnova} performs automatic peak detection and runs traditional 
#' ANOVA on them
#' @param arraydat a numeric array with named dimnames containing the EEG (or 
#' other) data. Missing values are not allowed.
#' @param factordef a named list of factor definitions, containing the following 
#' elements:
#' \itemize{
#' \item{"between"}{"character vector of between-subject factors (default: NULL)"}
#' \item{"within"}{"character vector of within-subject factors (default: NULL)"}
#' \item{"w_id"}{"name of the dimension which identifies the subjects 
#' (default: "id")}
#' }
#' @param peakdef a named list of peak definitions. A peak definition is a 
#' three-element integer vector which specifies the time window (start and 
#' end point) and the polarity (1 for positivity and -1 for negativity).
#' @param bwdat a data.frame which contains the identification codes 
#' (factordef$w_id) and the group memberships (factordef$between) of the 
#' subjects. Missing values are not allowed.
#' @param type a character value of "tanova" (default), "dissimilarity", or 
#' "gfp"
#' @param verbose logical value indicating whether means for each factor level 
#' should be also returned
#' @param avg_around_peak integer value which specifies the number of data points
#' which are considered +- in the "mean around the peak" measure
#' @param nperm integer value giving the number of permutations (default: 999L)
#' @param useparallel logical value; if TRUE (default), computations are done
#' in parallel
#' @param ncores integer value corresponding to the number of cores; 
#' if NULL (default), it is set to the maximum number of cores available
#' @param par_method parallelization method; can be "snow" (default) or "mc"
#' (multicore)
#' @param cl a cluster definition for snow-type parallelization; if NULL 
#' (default), it is set up automatically
#' @param iaterms_last logival variable whether all interaction terms should 
#' follow all main effect terms (default: TRUE)
#' @param seed an integer value which specifies a seed (default: NULL)
#' @details The function assumes that the input array contains at least two 
#' named dimensions: chan (corresponding to the channels [electrodes]) and time 
#' (corresponding to time points). All dimensions which are not listed as 
#' within-subject factors are treated in a similar way as chan and time, that is 
#' separate TANOVA-s are computed for each level of those dimensions.
#' @note The function computes type I p-values - this is correct if the design
#' is fully balanced and orthogonal (if the number of between-subject 
#' factors is one, it may have slightly unequal group sizes).
#' @export
#' @return A list object with effect statistics, uncorrected p-values, and 
#' two types of corrected p-values: length correction or global correction
#' @references Koenig (2011) RAGU
peakAnova <- function(arraydat, factordef, peakdef, bwdat = NULL,
                      verbose = TRUE, avg_around_peak = 2, 
                      nperm = 1, useparallel = FALSE, ncores = NULL,
                      par_method = c("snow", "mc"), cl = NULL,
                      iaterms_last = TRUE, seed = NULL) {
    # some checks
    stopifnot(is.array(arraydat) | missing(arraydat))
    stopifnot(is.list(factordef) | missing(factordef))
    stopifnot(length(peakdef) != 2 | missing(peakdef))
    if (!is.null(factordef$between) && is.null(bwdat)) {
        stop("No between-participant data provided")
    }
    if (useparallel) {
        if (is.null(ncores)) ncores <- parallel::detectCores()
    }
    #
    out <- list(call = match.call())
    #
    sumSq <- function(arraydat, f_dat, aov_form, labels = FALSE) {
        out <- summary(aov(as.formula(aov_form), data = f_dat))
        if (!labels) {
            out <- matrixIP(unlist(lapply(out, "[", "Mean Sq"), 
                                   use.names = FALSE),
                            nrow(out[[1]]), ncol(arraydat))
        } else {
            out <- matrixIP(unlist(lapply(out, "[", "Mean Sq"), 
                                   use.names = FALSE),
                            nrow(out[[1]]), ncol(arraydat),
                            list(term = gsub(" ","",rownames(out[[1]])), 
                                 peak = seq_along(out)))
        }         
        return(out)
    }
    # 
    par_method <- match.arg(par_method)
    # prepare data
    temp <- preAnova(arraydat, factordef, bwdat, useparallel, ncores, par_method)
    # assign variables to this environment, and potentially overwrite existing ones
    assignList(temp, verbose = FALSE)
    rm(temp)
    timedim <- which(names(dimnames(arraydat)) == "time")
    origtimedim <- which(names(origdimnames) == "time")
    origdimnames_peaks <- origdimnames
    arraydat_peaks <- subsetArray(arraydat, 
                                  list(time = dimnames(arraydat)$time[seq_along(peakdef)]))
    origdimnames_peaks$time <- 
        dimnames(arraydat_peaks)$time <- seq_along(peakdef)
    names(origdimnames_peaks)[origtimedim] <- 
        names(dimnames(arraydat_peaks))[timedim] <- "peak"
    # formula of anova
    aov_formula_char <- 
        if (is.null(factordef$within)) {
            paste(
                as.character(quote(arraydat)), " ~ ", 
                paste(as.character(factordef$between), collapse = "*"),
                sep = "")
        } else if (is.null(factordef$between)) {
            paste(
                as.character(quote(arraydat)), " ~ ", 
                paste(as.character(factordef$within), collapse = "*"),
                sep = "")
        } else {
            paste(
                as.character(quote(arraydat)), " ~ ", 
                paste(
                    paste(as.character(factordef$between), collapse = "*"), 
                    paste(as.character(factordef$within), collapse = "*"),
                    sep = "*"),
                sep = "")
        }
    aov_formula <- as.formula(aov_formula_char)
    mean_formula <- as.formula(gsub("\\*", ":", aov_formula_char))
    # compute marginal means
    origdims <- vapply(origdimnames, length, 0L)
    modeldims <- c(factordef$w_id, factordef$within)
    keepdims <- setdiff(names(origdimnames), modeldims)
    # find indices corresponding to time windows
    tpoints <- as.numeric(origdimnames$time)
    peakdef <- lapply(peakdef, function(x) 
        c(which(origdimnames$time %in% as.character(x[1:2])), x[3]))
    # find peaks
    fmeans <- marginalMeans(mean_formula, dat, arraydat, 
                            origdimnames[keepdims], 
                            keep_term_order = !iaterms_last, residualmean = FALSE)
    datrowind <- match(
        interaction(dat[, strsplit(names(fmeans), ":")[[1]] ]), 
        rownames(fmeans[[1]]) )
    peakind_facs <- strsplit(rownames(fmeans[[1]]), "\\.")
    peakind_facs <- setNames(data.frame(
        matrix(unlist(peakind_facs, use.names = FALSE),
               length(peakind_facs), length(peakind_facs[[1]]), TRUE)),
        strsplit(names(fmeans), ":")[[1]])
    peakind <- sapply(peakdef, function(x) {
        minmax <- if (x[3] < 0) which.min else which.max
        x[1] -1 + apply(fmeans[[1]][,seq(x[1],x[2])], 1, minmax)})
    #
    # peak amplitudes
    arraydat_peaks[] <- sapply(1:ncol(peakind), function(i) {
        out <- matrixIP(0, nrow(arraydat), avg_around_peak*2 + 1)
        avgind <- outer(peakind[,i], c(-avg_around_peak:avg_around_peak), "+")
        out[] <- arraydat[
            cbind(seq_along(datrowind), as.integer(avgind[datrowind,]))]
        return(rowMeans(out))
    })
    #
    # Compute ANOVA on peak amplitudes
    Fvals_obs <- arrayAnovaSub(arraydat_peaks, 
                               factordef, origdimnames_peaks, dat, verbose)
    out <- c(out, list(F_obs = Fvals_obs))
    # if verbose, save factor means
    if (verbose) {
        factor_means <- data.frame(
            group = peakind_facs[rep(1:nrow(peakind_facs), 
                                   length(peakdef)), , drop = F],
            peak = factor(rep(seq_along(peakdef), each = nrow(fmeans[[1]]))),
            ampl = c(marginalMeans(mean_formula, dat, arraydat_peaks, 
                                 origdimnames_peaks["peak"], 
                                 keep_term_order = !iaterms_last, 
                                 residualmean = FALSE)[[1]]),
            lat = tpoints[c(peakind)])
        out <- c(out, list(factor_means = factor_means))
    }
    #
    # Compute permutation statistics for peak latencies
    #
    # observed sum of squares
    lateff_obs <- sumSq(peakind, peakind_facs, aov_formula_char, labels = TRUE)
    # permutations
    if (nperm > 1) {
        #
        permfn <- function(i) {
            facmeans <- marginalMeans(mean_formula, dat, arraydat[randind[i,], ], 
                                      origdimnames[keepdims], 
                                      keep_term_order = !iaterms_last, 
                                      residualmean = FALSE)[[1]]
            peaki <- sapply(peakdef, function(x) {
                minmax <- if (x[3] < 0) which.min else which.max
                x[1] - 1 + apply(facmeans[, seq(x[1], x[2])], 1, minmax)})
            out <- sumSq(peaki, peakind_facs, aov_formula_char)
            return(out)
        }
        #
        # generate random orders (dim(randind) = nperm X nrow(dat))
        if (!is.null(seed)) set.seed(seed)
        randind <- mclapply(1:nperm, function(i) 
            shuffle(nrow(dat),  
                    how(within = Within(type = "free"), 
                        plots = Plots(strata = dat[,factordef$w_id], 
                                    type = "free"))),
            mc.cores = par_params$ncores)
        randind <- matrix(unlist(randind, use.names = FALSE), 
                          nperm, nrow(dat), TRUE)
        #
        if (!useparallel) {
            lateff_perm <- lapply(1:nperm, permfn)
        } else if (par_method == "snow") {
            if (is.null(cl)) cl <- makePSOCKcluster(par_params$ncores, outfile = "")
            clusterExport(cl, 
                          varlist = par_params$varlist2snow,
                          envir = environment())
            lateff_perm <- parLapply(cl, 1:nperm, permfn)
            #             stopCluster(cl)
            #             rm(cl)
        } else if (par_method == "mc") {
            lateff_perm <- mclapply(1:nperm, permfn, mc.cores = par_params$ncores)
        } 
        lateff_perm <- cbind(c(lateff_obs), 
                             matrixIP(unlist(lateff_perm, use.names = FALSE),
                                      ncol = nperm))
        lateff_perm <- sweep(lateff_perm[, -1, drop = F], 1, 
                             lateff_perm[, 1], "-")
        pvalues <- arrayIP(
            (rowSums(lateff_perm >= 0) + 1) / (nperm + 1),
            dim(lateff_obs), dimnames(lateff_obs))
        out <- c(out, list(lat_pvalues = pvalues))
    } 
    # return
    out
}

#
# <<< plotting functions >>> -----------
#

#' Plot ERP curves
#' 
#' \code{plotERParray} is a generalization of \code{matplot} onto array inputs.
#' @param dat a numeric array with named dimensions
#' @param xdim character value; the name of the dimension of dat which defines 
#' the x axis (default: "time")
#' @param sepdim character value; the name of the dimension of dat which 
#' separates the lines (default: "chan")
#' @param title character value; the title of the plot
#' @param subtitle.col the colour of the subtitles
#' @param gfp_plot logical value; if TRUE (default), the GFP curves are also 
#' plotted
#' @param gfp_col the colour of the GFP curves (if plotted)
#' @param gfp_lwd the thickness of the GFP curves (if plotted)
#' @param minus_up logical value; if set to TRUE, minus values are plotted 
#' upwards. If set to FALSE, minus values are plotted downwards. If NULL 
#' (default), it is set to TRUE or FALSE depending on the gfp_plot parameter. 
#' NULL
#' @param grid_labels character vector; the name of the x and y axes
#' @param grid_dim integer vector giving the number of rows and columns
#' in which the plots are arranged; if set to NULL (default), the arrangement 
#' of the plot is set up automatically
#' @param ... additional parameters to be passed to \code{matlines}
#' @export
plotERParray <- function(dat, xdim = "time", sepdim = "chan", 
                         title = "", subtitle.col = "black",
                         gfp_plot = TRUE, gfp_col = "black", gfp_lwd = 1.3, 
                         minus_up = NULL, grid_labels = c("time", "ampl"), 
                         grid_dim = NULL, ...) {
    emptyplot <- function() {
        plot(0, 0, xlim = c(-1,1), type = "n", axes = FALSE, 
             frame.plot = FALSE, xlab = "", ylab = "")
    }
    if (length(dim(dat)) > 3) {
        dat <- mergeDims(dat, setdiff(names(dimnames(dat)), c(xdim, sepdim)))
    }
    wrapdim <- setdiff(names(dimnames(dat)), c(xdim, sepdim))
    dat <- aperm(dat, c(xdim, sepdim, wrapdim))
    subtitle.col <- rep(subtitle.col, length_out = dim(dat)[3])
    if (gfp_plot) gfpdat <- compGfp(dat)
    x <- as.numeric(as.character(dimnames(dat)[[1]]))
    xrange <- if (is.null(list(...)$xlim)) range(x) else list(...)$xlim 
    if (is.null(list(...)$ylim)) {
        yr <- range(dat)
        yrange <- mean(yr) + c(-1, 1)*(yr[2]-mean(yr))*1.02
        if (is.null(minus_up)) {
            minus_up <- if (!gfp_plot) TRUE else FALSE
        }
        if (minus_up) yrange <- -yrange
    } else {
        yrange <- list(...)$ylim
    } 
    if (is.null(grid_dim)) grid_dim = rep(ceiling(sqrt(dim(dat)[3])), 2)
    layoutmat <- cbind(
        c(0, rep(2, grid_dim[1]), 0, 0),
        rbind(rep(1, grid_dim[2]),
              matrix(1:prod(grid_dim), grid_dim[1], grid_dim[2], TRUE) + 3,
              rep(0, grid_dim[2]),
              rep(3, grid_dim[2])),
        c(rep(0, grid_dim[1] + 3)))
    layout(layoutmat, widths = c(0.3, rep(1, grid_dim[2]), 0.3),
           heights = c(0.3, rep(1, grid_dim[1]), 0.3, 0.3))
    par(mar = rep(0, 4))
    emptyplot(); text(0, 0, title, cex = 1.5)
    emptyplot(); text(0, 0, grid_labels[2], cex = 1.3, srt = 90)
    emptyplot(); text(0, 0, grid_labels[1], cex = 1.3)
    par(mar = c(0, 0, 2, 0))
    for (i in 1:dim(dat)[3]) {
        matplot(x, dat[,,i], xlim = xrange, ylim = yrange, yaxs = "i",  
                axes = FALSE, frame.plot = TRUE, type = "n")
        grid()
        mtext(dimnames(dat)[[3]][i], 3, cex = 0.7, col = subtitle.col[i])
        matlines(x, dat[,,i], axes = FALSE, ...)
        if (i == max(dim(dat)[3])) {
            axis(1)
            axis(4)
        }
        if (gfp_plot) lines(x, gfpdat[, i], col = gfp_col, lwd = gfp_lwd)
    }
}


#' Plot the topography of the signals
#' 
#' \code{plot2dview} creates 2D topoplots
#' @param dat a numeric vector or matrix. If a matrix is provided, it must 
#' contain the participants in rows and the channels in columns
#' @param ch_pos a data.frame with the channel coordinates
#' @param r a numeric value; the radius of the head
#' @param timepoint an integer value
#' @param ampl_range numeric vector, the minimum and maximum value of the 
#' amplitudes (default: c(-5, 5))
#' @param resol integer value, the resolution of the projection (default: 100L)
#' @param resolcol integer value, the resolution of the colours (default: 1000L)
#' @param projection character value (default = "laea"). See 
#' \link{http://www.remotesensing.org/geotiff/proj_list/} for common projections
#' @param projref projection reference (pole [ = default] or equator)
#' @param gfp a numeric vector of the GFP values in the whole time window
#' @importFrom sgeostat in.polygon
#' @importFrom gplots colorpanel redgreen greenred bluered redblue
#' @export
plot2dview <- function(dat, ch_pos, r = 1, timepoint = NULL, ampl_range = c(-5, 5), 
                       resol = 100L, resolcol = 1000L,
                       projection = "laea", projref = c("pole", "equator"),
                       origo = c(lat = ifelse(projref == "pole", 90, 0),
                               long = ifelse(projref == "pole", 270, 0)),
                       plot_centroid = TRUE, centroid_circle = 0.5, centroid_unitsphere = FALSE,
                       plot_bar = TRUE, plot_ch = TRUE, plot_chnames = TRUE, 
                       title = NULL,
                       gfp = NULL, gfp_max = NULL, ...) {
    #
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("theta", "phi") %in% colnames(ch_pos))) {
        ch_pos <- cart2sph(ch_pos)
    }
    if (!is.null(ch_pos$r)) r <- ch_pos$r
    if (is.matrix(dat)) {
        if (all(c("chan", "id") %in% names(dimnames(dat))))
            dat <- aperm(dat, c("id", "chan"))
        if (ncol(dat) != nrow(ch_pos)) 
            stop("Wrong data format! Should be: Rows = participants, Columns = channels")
        subjdat <- dat
        dat <- colMeans(dat, na.rm = TRUE)
    } else {
        centroid_circle <- NA
    }
    #
    projref <- match.arg(projref)
    posgeo <- sph2geo(ch_pos)
    #
    boundarypos <- 
        if (projref == "pole") {
            data.frame(
                long = seq(-180, 180, length.out = 180),
                lat = rep(min(posgeo$lat), 180))
        } else {
            data.frame(
                long = c(rep(border*max(abs(posgeo$long)), 90),
                       rep(-border*max(abs(posgeo$long)), 90)),
                lat = c(seq(-90, 90, length.out = 90), 
                        seq(90, -90, length.out = 90)))
        }
    boundarypos <- unique(project3dMap(boundarypos, projection = "laea", 
                                       projref = projref, origo = origo))
    #
    xlim <- range(boundarypos$x) * 1.1
    ylim <- range(boundarypos$y) * 1.1
    ylim.cm <- ylim.cp <- 0.1
    if (plot_bar) {
        ypos.bar <- c(0, ylim[1] * (1 + ylim.cm))
        ylim.cm <- ylim.cm + 0.1
        ypos.bar[1] <- 0.75 * ypos.bar[2] + 0.25 * ylim[1] * (1 + ylim.cm)    
    }
    if (!is.null(gfp)) {
        ypos.gfp <- c(ylim[2] * (1 + ylim.cp * 1.1), 0)
        ylim.cp <- ylim.cp + 0.3
        ypos.gfp[2] <- 0.1 * ypos.gfp[1] + 0.9 * ylim[2] * (1 + ylim.cp)
        if (is.null(gfp_max)) gfp_max <- max(gfp) * 1.05
        gfp <- ypos.gfp[1] + gfp * diff(ypos.gfp) / gfp_max
    } else {
        ypos.gfp <- c(ylim[2], ylim[2])
    }
    if (!is.null(title) | !is.null(timepoint)) {
        ylim.cp <- ylim.cp + 0.2
        ypos.title <- 0.3 * ypos.gfp[2] + 0.7 * ylim[2] * (1+ylim.cp)
    }
    ylim[1] <- ylim[1] * (1 + ylim.cm)
    ylim[2] <- ylim[2] * (1 + ylim.cp)
    #
    gridx <- seq(min(boundarypos$x), max(boundarypos$x), length.out = resol)
    gridy <- seq(min(boundarypos$y), max(boundarypos$y), length.out = resol)
    gridpos <- expand.grid(x = gridx, y = gridy)
    ind <- in.polygon(gridpos$x, gridpos$y, 
                      boundarypos$x * 1.1, boundarypos$y * 1.1)
    gridgeo <- project3dMap(gridpos[ind,], projection = "laea", 
                            projref = projref, origo = origo, inverse = TRUE)
    gridcart <- geo2cart(gridgeo)
    z <- matrixIP(NA_real_, resol, resol)
    z[ind] <- chanInterp(dat, ch_pos, gridcart)
    par(mar = rep(0, 4))
    image(gridx, gridy, z, useRaster = TRUE, col = bluered(resolcol),
          xlim = xlim, ylim = ylim,
          zlim = ampl_range, xlab = "", ylab = "", axes = FALSE)
    gridx <- seq(min(boundarypos$x), max(boundarypos$x), 
                 length.out = 1000) * 1.1
    gridy <- seq(min(boundarypos$y), max(boundarypos$y), 
                 length.out = 1000) * 1.1
    gridpos <- expand.grid(x = gridx, y = gridy)
    ind <- !in.polygon(gridpos$x, gridpos$y, 
                       boundarypos$x*1, boundarypos$y*1)
    z <- matrixIP(NA_integer_, 1000, 1000)
    z[ind] <- 1L
    image(gridx, gridy, z, useRaster = TRUE, col = "white", add = TRUE)
    lines(boundarypos, col = "white", lwd = 1)
    #
    if (plot_ch) {
        ch_pos_xy <- project3dMap(ch_pos)
        points(ch_pos_xy,, pch = 20)
        if (plot_chnames) {
            indleft <- ch_pos_xy$x <= 0
            indright <- ch_pos_xy$x > 0
            text(ch_pos_xy[indleft,], , rownames(ch_pos[indleft, ]), pos = 4)
            text(ch_pos_xy[indright,], , rownames(ch_pos[indright, ]), pos = 2)
        }
    }
    #
    if (plot_bar) {
        xpos.bar <- seq(min(boundarypos$x / 3), max(boundarypos$x / 3), 
                        length.out = resolcol)
        image(xpos.bar, ypos.bar,
              matrixIP(seq(ampl_range[1], ampl_range[2], length.out = resolcol), 
                       resolcol, 2),
              zlim = ampl_range, col = bluered(resolcol), add = TRUE)
        text(min(xpos.bar), ypos.bar[1], 
             substitute(paste(k, " ", mu, "V", sep = ""), 
                        list(k = ampl_range[1])), pos = 2)
        text(max(xpos.bar), ypos.bar[1], 
             substitute(paste(k, " ", mu, "V", sep = ""), 
                        list(k = ampl_range[2])), pos = 4)
    }
    if (!is.null(title) | !is.null(timepoint)) {
        if (is.null(title)) {
            title <- paste("Time:", timepoint, "ms")
        }
        text(0, ypos.title, title, cex = 1.1)
    }
    if (!is.null(gfp)) {
        lgfp <- length(gfp)
        gfp.x <- seq(xlim[1] * 0.8, xlim[2] * 0.8, length.out = lgfp)
        rect(min(gfp.x), ypos.gfp[1], max(gfp.x), ypos.gfp[2], 
             border = NA, col = "grey60")
        axis(1, at = gfp.x[seq(1, lgfp, 25)], 
             labels = names(gfp)[seq(1, lgfp, 25)], 
             pos = ypos.gfp[1], cex.axis = 0.7, padj = -1)
        lines(gfp.x, gfp, col = "white")
        smpl <- which(as.numeric(names(gfp)) == timepoint)
        if (length(smpl) == 1) {
            points(gfp.x[smpl], gfp[smpl], pch = 16, col = "red", cex = 1.1)
        }
    }
    if (plot_centroid) {
        c1 <- centroid(dat, ch_pos, proj_unitsphere = centroid_unitsphere)
        if (!is.na(centroid_circle)) {
            if (reqFn("plotrix")) {
                cc <- apply(subjdat, 1, centroid, ch_pos, 
                            proj_unitsphere = centroid_unitsphere)
                temp <- matrixIP(unlist(lapply(cc, function(x) abs(x-c1))), 
                                 ncol = length(cc))
                rads <- apply(temp, 2, quantile, 
                              probs = centroid_circle, na.rm = TRUE)
                plotrix::draw.ellipse(x = c1$x, y = c1$y, 
                             a = rads[c(1, 3)], b = rads[c(2, 4)], 
                             col = c(rgb(0, 0, 0.5, 0.15), rgb(0.5, 0, 0, 0.15)), 
                             border = NA)
            }
        }
        points(c1["neg", "x"], c1["neg", "y"], col = "white", 
               pch = 15, cex = 2.5)
        points(c1["neg", "x"], c1["neg", "y"], col = "blue", 
               pch = 15, cex = 2)
        points(c1["pos", "x"], c1["pos", "y"], col = "white", 
               pch = 15, cex = 2.5)
        points(c1["pos", "x"], c1["pos", "y"], col = "red", 
               pch = 15, cex = 2)
    }
}


#' Plot topography in complex layout
#' 
#' \code{complexplot2dview} calls plot2dview for the given groups/conditions
#' @param dat a list of numeric vectors or matrices, which can be plotted by
#' \code{\link{plot2dview}}
#' @param ch_pos channel position matrix
#' @param datgrid arrangement of the grid
#' @seealso \code{\link{plot2dview}}
#' @export
complexplot2dview <- function(dat, ch_pos, timepoint, 
                              datgrid = NULL, layout_matrix = NULL, 
                              heights = NULL, widths = NULL, 
                              title_row = NULL, title_col = NULL, 
                              title_main = paste("Time:", timepoint, "ms"), 
                              plot_bar = TRUE, ampl_range = NULL, resolcol = 1000,
                              gfp = NULL, ...) {
    #
    emptyplot <- function() plot(0, 0, xlim = c(-1, 1), type = "n", 
                                 axes = FALSE, frame.plot = FALSE, xlab = "", ylab = "")
    if (is.null(layout_matrix)) {
        if (is.null(datgrid)) {
            datgrid <- c(floor(sqrt(length(dat))), 
                         ceiling(sqrt(length(dat))))
        }
        rownum <- datgrid[1]
        colnum <- datgrid[2]
        datnum <- prod(datgrid)
        layout_matrix <- rbind(
            matrix(c(
                0, rep(1, colnum), 
                0, 2:(colnum + 1)), 2, colnum + 1, TRUE), 
            matrix(c(
                (colnum + 2):(colnum + rownum + 1), 
                (colnum + rownum + 1) + 1:datnum), rownum, colnum + 1))
        layout_matrix <- rbind(layout_matrix, 
                               c(0, rep(max(layout_matrix) + 1, colnum)))
        if (is.null(heights)) 
            heights <- c(0.3, 0.2, rep(1, rownum), 0.2)
        if (is.null(widths)) 
            widths <- c(0.2, rep(1, colnum))
    } else {
        if (is.null(heights)) heights <- rep(1, nrow(layout_matrix))
        if (is.null(widths)) widths <- rep(1, ncol(layout_matrix))
    }
    layout(layout_matrix, widths = widths, heights = heights)
    par(mar = rep(0, 4))
    emptyplot(); text(0, 0, title_main, cex = 1.5)
    for (i in 1:length(unique(layout_matrix[2, -1]))) {
        emptyplot()
        text(0, 0, title_col[i], cex = 1.3)
    }
    for (i in 1:length(unique(layout_matrix[-c(1, 2, nrow(layout_matrix)),1]))) {
        emptyplot()
        text(-1, 0, title_row[i], cex = 1.3, pos = 4)
    }    
    if (is.null(ampl_range)) {
        ampl_range <- c(min(unlist(dat, use.names = FALSE)), 
                        max(unlist(dat, use.names = FALSE)))
        ampl_range <- round(c((1 - 0.2 * sign(ampl_range)[1]) * ampl_range[1],
                              (1 + 0.2 * sign(ampl_range)[2]) * ampl_range[2]), 
                            1)
    }
    if (length(unlist(gfp, use.names = FALSE)) == 1 && gfp == TRUE) {
        gfp <- lapply(dat, compGfp)    
    }
    gfp_max <- if (is.null(gfp)) NULL else max(unlist(gfp, use.names = FALSE))
    for (i in 1:length(dat)) {
        plot2dview(dat[[i]], ch_pos = ch_pos, timepoint = timepoint, 
                   gfp = gfp[[i]], gfp_max = gfp_max,
                   plot_bar = FALSE, ampl_range = ampl_range, title = "", ...)
    }
    if (plot_bar) {
        xpos.bar <- seq(0.4, 0.6, length.out = resolcol)
        ypos.bar <- c(0.45, 0.55)
        image(xpos.bar, ypos.bar,
              matrixIP(seq(ampl_range[1], ampl_range[2], length.out = resolcol), 
                       resolcol, 2),
              xlim = c(0, 1), ylim = c(0, 1), zlim = ampl_range, 
              useRaster = TRUE, col = bluered(resolcol), 
              xlab = "", ylab = "", axes = FALSE 
        )
        text(min(xpos.bar), ypos.bar[1], 
             substitute(paste(k, " ", mu, "V", sep = ""), 
                        list(k = ampl_range[1])), pos = 2)
        text(max(xpos.bar), ypos.bar[1], 
             substitute(paste(k, " ", mu, "V", sep = ""), 
                        list(k = ampl_range[2])), pos = 4)
    }
}

# plot p-values with across-country similarity matrices 
reprPlot <- function(p, e, title = "Reproducibility", plim = c(0.05, 0.01), 
                     plim_shading = TRUE, plim_outline = 0.05,
                     corr_range = c(-0.99, 0.99), 
                     corr_method = c("pearson", "spearman")) {
    #
    emptyplot <- function() {
        plot(0, 0, xlim = c(-1, 1), type = "n", 
             axes = FALSE, frame.plot = FALSE, xlab = "", ylab = "")
    }
    pPolygon <- function(x, y, direction = "h") {
        plim <- sort(unique(c(0, plim, 1)))
        y_cat <- findInterval(y, plim)
        y_phase <- c(1, cumsum(abs(diff(y_cat))) + 1)
        yy <- 0.04 + -log(y) / (-log(p_min)) * 0.15
        cols <- grey(c(seq(0.4, 0.6, 
                           length.out = (length(plim)-2)), 0.95))
        if (direction == "h") {
            for (i in unique(y_phase)) {
                ind <- y_phase == i
                polygon(c(x[ind], rev(x[ind])), 
                        c(-yy[ind], rep(-0.04, sum(ind))),
                        col = cols[y_cat[ind][1]], border = "grey20")
                if (plim_shading && y_cat[ind]<(length(plim)-1)) {
                    acol <- col2rgb(cols[y_cat[ind][1]], TRUE) / 255
                    acol[4] <- 0.2
                    rect(x[min(which(ind))], 0,
                         x[max(which(ind))], 1,
                         col = rgb(acol[1], acol[2], acol[3], acol[4]), 
                         border = rgb(0, 1, 0, acol[4]), lty = "1F")
                }
            }
            if (!is.na(plim_outline)) {
                sig <- rep(NA_integer_, length(x))
                sig[y <= plim_outline] <- 1L
                lines(x, sig, col = "green", lwd = 3)
            }
        } else {
            for (i in unique(y_phase)) {
                ind <- y_phase == i
                polygon(c(-yy[ind], rep(-0.04, sum(ind))), 
                        c(x[ind], rev(x[ind])),
                        col = cols[y_cat[ind][1]], border = "grey20")
                if (plim_shading && y_cat[ind]<(length(plim)-1)) {
                    acol <- col2rgb(cols[y_cat[ind][1]], TRUE) / 255
                    acol[4] <- 0.2
                    rect(0, x[min(which(ind))],
                         1, x[max(which(ind))],
                         col = rgb(acol[1], acol[2], acol[3], acol[4]), 
                         border = rgb(0, 1, 0, acol[4]), lty = "1F")
                }
            }
            if (!is.na(plim_outline)) {
                sig <- rep(NA_integer_, length(x))
                sig[y <= plim_outline] <- 1L
                lines(sig, x, col = "green", lwd = 3)
            }
        }
    }                     
    implot <- function(corrmat, pp) {
        zlim <- atanh(corr_range)
        mat <- atanh(corrmat)
        mat[mat < zlim[1]] <- zlim[1]
        mat[mat > zlim[2]] <- zlim[2]
        image(mat, zlim = zlim, xlim = c(-0.2, 1), ylim = c(-0.2, 1),
              axes = FALSE, xlab = "", ylab = "", col = bluered(1000))
        lines(c(0, 1), c(0, 1), col = "grey70")
        tpoints <- as.numeric(rownames(corrmat))
        seqtpoints <- tpoints[seq(1, length(tpoints), 25)]
        at <- (seqtpoints - min(tpoints)) / (max(tpoints) - min(tpoints))
        x <- (tpoints - min(tpoints)) / (max(tpoints) - min(tpoints))
        for (ii in 1:2) {
            axis(ii + 2, at = at, labels = tpoints[tpoints %in% seqtpoints],
                 cex.axis = 0.8)
            y <- pp[, ii]
            pPolygon(x, y, c("h", "v")[ii])
        }
    }
    plim <- sort(unique(c(0, plim, 1)))
    corr_method <- match.arg(corr_method)
    e <- aperm(e, c("chan", "time", "nation"))
    p <- aperm(p, c("time", "nation"))
    p_min <- min(p)
    natnum <- ncol(p)
    nations <- colnames(p)
    imlayout <- matrixIP((2 * natnum + 1), natnum, natnum)
    diagpos <- diag(matrixIP(1:natnum^2, natnum, natnum))
    for (i in 1:length(imlayout)) {
        imlayout[i] <- 
            if (i %in% diagpos) 0 
        else max(imlayout)+1
    }
    layout_matrix <- 
        rbind(c(0, rep(1, natnum)),
              c(0, rep(2:(natnum + 1))),
              cbind((natnum + 2):(2 * natnum + 1),
                    imlayout))
    layout(layout_matrix, 
           widths = c(0.4, rep(1, natnum)),
           heights = c(0.2, 0.1, rep(1, natnum)))
    par(mar = rep(0, 4))
    emptyplot(); text(0, 0, title, cex = 1.5)
    for (i in 1:(2*natnum)) {
        emptyplot()
        if (i < 5) { 
            text(0, 0, nations[i], cex = 1.2)
        } else {
            text(-1, 0, nations[i-4], cex = 1.2, pos = 4)
        }
    }
    par(mar = c(1, 1, 3, 3))
    for (nat in 1:natnum) {
        for (i in 1:natnum) {
            if (nat != i) {
                corrs <- cor(e[,,nat], e[,,i], method = corr_method)
                implot(corrs, p[, c(nat, i)])
            }
        }
    }
}

#' Image plot of p-values
#' 
#' \code{imagePvalues} creates an image plot from a matrix or array of p-values
#' @param pvalues numeric matrix or array of p-values with named dimensions (at
#' least chan and time)
#' @param pcrit numeric vector of significancy limits 
#' (default: 0.001, 0.01, 0.05)
#' @param grid character vector or formula defining the layout of panels
#' @param wrap character vector or formula defining the dimension which 
#' separates panels (only considered if grid is NULL)
#' @export
#' @return A ggplot object
imagePvalues <- function(pvalues, pcrit = c(0.001, 0.01, 0.05),
                         grid = NULL, wrap = NULL) {
    require(ggplot2)
    chans <- dimnames(pvalues)$chan
    timebreaks <- as.numeric( as.character( dimnames(pvalues)$time ) )
    timebreaks <- range( timebreaks %/% 100)
    timebreaks <- seq(timebreaks[1], timebreaks[2]) * 100
    pvalues_df <- array2df(pvalues, response_name = "p", 
                           dim_types = list(time = "numeric"))
    pvalues_df$pcrit <- findInterval(pvalues_df$p, pcrit)
    for (i in colnames(pvalues_df)) {
        if (is.character(pvalues_df[[i]])) {
            pvalues_df[[i]] <- factor(pvalues_df[[i]], 
                                      levels = dimnames(pvalues)[[i]])
        }
    }
    pp <- ggplot(pvalues_df, aes(x = time, y = chan)) + 
        geom_tile(aes(fill = pcrit), size = 0) + 
        scale_fill_gradient(guide = "legend", 
                            high = "white", 
                            low = "#a50f15",
                            name = "p-value",
                            breaks = 0:2,
                            labels = paste("< ", pcrit, sep = ""),
                            limits = c(0, length(pcrit))) + 
        scale_y_discrete(limits = rev(chans), expand = c(0.05, 0),
                         name = "channels") + 
        scale_x_continuous(breaks = timebreaks, expand = c(0.05, 0.1))
    if (!is.null(grid)) {
        pp <- pp + facet_grid( as.formula(grid) )
    } else if (!is.null(wrap)) {
        pp <- pp + facet_wrap( as.formula(wrap) )
    }
    # return
    pp
}


#' Image plot of channel values
#' 
#' \code{imageValues} creates an image plot from a matrix or array of channel x 
#' time values
#' @param dat numeric matrix or array of values with named dimensions (at least
#' chan and time)
#' @param grid character vector or formula defining the layout of panels
#' @param wrap character vector or formula defining the dimension which 
#' separates panels (only considered if grid is NULL)
#' @export
#' @return A ggplot object
imageValues <- function(dat, grid = NULL, wrap = NULL) {
    require(ggplot2)
    chans <- dimnames(dat)$chan
    timebreaks <- as.numeric( as.character( dimnames(dat)$time ) )
    timebreaks <- range( timebreaks %/% 100)
    timebreaks <- seq(timebreaks[1], timebreaks[2]) * 100
    dat_df <- array2df(dat, response_name = "effect", 
                       dim_types = list(time = "numeric"))
    for (i in colnames(dat_df)) {
        if (is.character(dat_df[[i]])) {
            dat_df[[i]] <- factor(dat_df[[i]], 
                                  levels = dimnames(dat)[[i]])
        }
    }
    pp <- ggplot(dat_df, aes(x = time, y = chan)) + 
        geom_tile(aes(fill = effect), size = 0) + 
        scale_y_discrete(limits = rev(chans), expand = c(0.05, 0),
                         name = "channels") + 
        scale_x_continuous(breaks = timebreaks, expand = c(0.05, 0.1))
    if (!is.null(grid)) {
        pp <- pp + facet_grid( as.formula(grid) )
    } else if (!is.null(wrap)) {
        pp <- pp + facet_wrap( as.formula(wrap) )
    }
    # return
    pp
}


# simple function to plot TFCE effects
tfce.plot <- function(arraydat, breaks = c(0, 0.001, 0.01, 0.05), 
                      colors = rev(brewer.pal(length(breaks), "Reds")[-1]), title = "",
                      gridlines_step = 50) {
    extradim <- setdiff(names(dimnames(arraydat)), 
                        c("chan", "time", "modelterm", "nation"))
    if (!is.null(extradim)) arraydat <- avgDims(arraydat, extradim)
    arraydat <- aperm(arraydat, c("chan", "time", "modelterm", "nation"))
    dimnms <- dimnames(arraydat)
    tpoints <- as.numeric(dimnms$time)
    gridlines <- seq(
        min(tpoints) %/% gridlines_step * gridlines_step,
        max(tpoints) %/% gridlines_step * gridlines_step,
        gridlines_step)
    dims <- vapply(dimnms, length, 0L)
    emptyplot <- function() {
        plot(0, 0, xlim = c(-1, 1), 
             type = "n", axes = FALSE, frame.plot = FALSE, xlab = "", ylab = "")
    }
    layoutmat <- rbind(
        c(0, rep(1, dims[3])),
        c(0, (1:dims[3]) + dims[4] + 1),
        cbind(1:dims[4] + 1,
              matrix((1:prod(dims[3:4])) + sum(dims[3:4]) + 1, 
                     dims[4], dims[3], TRUE)))
    layout(layoutmat, 
           widths = c(0.3, rep(1, dims[3])),
           heights = c(0.5, 0.3, rep(1, dims[4])))
    par(mar = c(0, 0, 0, 0))
    # row and column labels
    emptyplot()
    text(0, 0, title, cex = 1.3)
    for (i in 1:dims[4]) {
        emptyplot()
        text(0, 0, dimnms[[4]][i], srt = 90)
    }
    for (i in 1:dims[3]) {
        emptyplot()
        text(0, 0, dimnms[[3]][i])
    }
    par(mar = c(2, 2, 0, 0.5))
    # plot images
    for (i in 1:dims[4]) {
        for (ii in 1:dims[3]) {
            temp <- aperm(arraydat[,,ii,i], c("time", "chan"))
            plot(0, 0,
                 xlim = range(tpoints), ylim = range(1:dims["chan"]),  
                 type = "n", 
                 xlab = "", ylab = "", yaxt = "n")
            abline(v = gridlines, lty = 1, col = "grey95")
            image(tpoints, 1:dims["chan"], temp, 
                  breaks = breaks, col = colors, add = TRUE,
                  xlab = "", ylab = "", yaxt = "n")
            axis(2, at = 1:dims["chan"], labels = FALSE, tick = FALSE)
            text(par("usr")[1] - 21, 1:dims["chan"], cex = 0.6, pos = 2,
                 labels = dimnms$chan, xpd = TRUE)
            abline(v = 0, col = "grey40")
        }
    }
}

###

plot.tanova <- function(results, plot_title = "", only_p = FALSE) {
    reshapefn <- function(slot, headername) {
        x <- results[[slot]]
        x <- array2df(x, response_name = headername, 
                      dim_types = list(time = "numeric"))
        return(x)
    }
    #
    pcrit <- as.list(results$call)$pcrit
    if (is.null(pcrit)) pcrit <- formals(tanova)$pcrit
    #
    dat <- array2df(results$effect, "effect", 
                    dim_types = list(time = "numeric"))
    dat$pvalue <- -log(c(results$perm_pvalues))
    dat$pvalue_consec <- -log(c(results$perm_pvalues_consec))
    dat$pcrit <- factor(dat$pvalue_consec > -log(pcrit))
    legendtitle <- paste("pvalue <", substitute(pcrit), sep = " ")
    #
    if (only_p) {
        qp <- ggplot(dat[order(dat$time),], 
                     aes(x = time, y = pvalue, col = pcrit, group = NA)) + 
            geom_hline(yintercept = -log(pcrit), lty = 3) + 
            ylab("-log(P-value)")
    } else {
        qp <- ggplot(dat[order(dat$time),], 
                     aes(x = time, y = effect, col = pcrit, group = NA)) + 
            ylab("Effect")
    }
    qp <- qp + geom_line() + facet_wrap(~modelterm) +
        ggtitle(plot_title) + 
        scale_colour_manual(name = legendtitle,
                            values = c("FALSE"="grey60","TRUE"="red")) + 
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(color = "grey95"),
              panel.grid.minor = element_line(color = "grey97"),
              panel.border = element_rect(color = "grey70", fill = NA))
    print(qp)
}


# plot neurodys tanova results
fastplot_tanova <- function(results, plot_title = "", pcrit = 0.05, 
                            only_p = FALSE) {
    reshapefn <- function(slot, headername) {
        x <- lapply(results, "[[", slot)
        #x <- rearrangeList(x, "nation")
        if (length(x) > 1) {
            temp <- matrix(unlist(strsplit(names(x), "-")), 
                           nrow = length(x), byrow = TRUE)
            names(x) <- temp[, 2]
            x <- rearrangeList(x, temp[1, 1])
        } else {
            x <- x[[1]]
        }
        x <- as.data.frame.table(x, responseName = headername)
        x$time <- as.numeric(as.character(x$time))
        return(x)
    }
    dimn <- names(results[[1]])
    dat <- reshapefn(dimn[grep("effect_", dimn)], "effect_size")
    dat$pvalue <- -log(reshapefn("perm_pvalues", "pvalue")$pvalue)
    dat$pvalue_consec <- reshapefn("perm_pvalues_consec", "pvalue")$pvalue
    dat$pcrit <- factor(dat$pvalue_consec<pcrit)
    legendtitle <- paste("pvalue <", substitute(pcrit), sep = " ")
    #
    if (only_p) {
        qp <- ggplot(dat[order(dat$time),], 
                     aes(x = time, y = pvalue, col = pcrit, group = NA)) + 
            geom_hline(yintercept = -log(pcrit), lty = 3) + 
            ylab("-log(P-value)")
    } else {
        qp <- ggplot(dat[order(dat$time),], 
                     aes(x = time, y = effect_size, col = pcrit, group = NA)) + 
            ylab("Effect size")
    }
    qp <- qp + geom_line() + facet_grid(nation~modelterm) +
        ggtitle(plot_title) + 
        scale_colour_manual(name = legendtitle,
                            values = c("FALSE"="grey70","TRUE"="red"))
    print(qp)
}


# plot peak results
fastplot_peak <- function(results, plot_title = "", 
                          pcrit = 0.05, only_p = FALSE) {
    reshapefn <- function(slot, headername, attr_slot = NULL) {
        if (is.null(attr_slot)) {
            x <- lapply(results, "[[", slot)
        } else {
            x <- lapply(results, 
                        function(x) attr(x[[slot]], attr_slot))
        }
        x <- rearrangeList(x, "nation")
        x <- as.data.frame.table(x, responseName = headername)
        x$peak <- factor(x$peak, levels = rev(levels(x$peak)))
        return(x)
    }
    dat_ampl <- reshapefn("F_obs", "pvalue", attr_slot = "pvalues")
    dat_ampl <- dat_ampl[with(dat_ampl, order(nation, peak, modelterm)),
                         c("nation", "peak", "modelterm", "pvalue")]
    dat_ampl$measure <- "amplitude"
    dat_lat <- reshapefn("lat_pvalues", "pvalue")
    dat_lat <- dat_lat[with(dat_lat, order(nation, peak, term)),
                       c("nation", "peak", "term", "pvalue")]
    colnames(dat_lat)[3] <- "modelterm"
    dat_lat$measure <- "latency"
    dat <- rbind(dat_ampl, dat_lat)
    dat$pcrit <- factor(dat$pvalue < pcrit)
    legendtitle <- paste("pvalue <", substitute(pcrit), sep = " ")
    qp <- ggplot(dat, aes(x = nation, y = peak, fill = pcrit)) + 
        geom_tile(col = "white") + facet_grid(modelterm ~ measure) + 
        ggtitle(plot_title) + 
        scale_fill_manual(name = legendtitle,
                          values = c("FALSE"="grey60","TRUE"="red3"))
    print(qp)
}

#
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow = 2, byrow = TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
    require(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots == 1) {
        print(plots[[1]])   
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))        
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

#
# <<< complex functions >>> -----------
#

#' Data preparation mainly aimed at facilitating plotting in lattice or ggplot2
#' 
#' \code{prepare2plot} provides several options to transform the eeg data array 
#' to a data file which enables direct plotting in lattice or ggplot2 afterwards.
#' It can also be used for analyses purposes where compact code is desirable.
#' @param dat an array of ERP data. Must have named dimnames, one of which must
#' be id (corresponding to participants' identification codes)
#' @param datid data.frame containing information on participants. Must have an
#' id column (participants' identification codes)
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
#
prepare2plot <- function(dat, datid, 
                         bwFac = NULL, wiFac = NULL, 
                         collFac = NULL, diffFac = NULL, 
                         compGFP = TRUE, keep_channels = !compGFP,
                         sc = NULL, 
                         datfr = TRUE, iaFac = NULL, ...) {
    bool_collFac <- !is.null(collFac) && collFac%in%names(dimnames(dat))
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


#
# <<< Rcpp functions >>> -----------
#

#' Low-level TFCE function which calls C++ code
#' 
#' \code{tfceFn} performs TFCE correction. This function is not 
#' intended for direct use.
#' @param x numeric matrix or array
#' @param chn channel neighbourhood matrix
#' @param eh numeric vector of E and H parameters
#' @param nr_steps number of threshold steps (default: 50L)
#' @export
#' @keywords internal
#' @return numeric matrix or array of the same dimensions as x
tfceFn <- function(x, chn, eh, nr_steps = 50L, channel_dim=1L) {
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
                    ChN = ChN, EH = eh, numSteps = nr_steps)
    }
    if (any(ind.neg)) {
        sdat <- -x
        sdat[ind.pos] <- 0
        out <- out - tfce(inData = sdat, chan_dim = chan_dim, 
                          ChN = ChN, EH = eh, numSteps = nr_steps)
    }
    # return
    out
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
    if (length(dim(x)) == 1) {
        col <- T
        x <- matrix(x, ncol = 1)
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
