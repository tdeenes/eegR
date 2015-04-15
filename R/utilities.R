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

#' Replace elements of a vector
#' 
#' \code{Replace} replaces the values in \code{x} which are equal to 
#' \code{from} to the values given in \code{to}. The values in \code{to} are
#' recycled to match the length of \code{from}.
#' @param x a raw, integer, character or numeric (double) vector
#' @param from vector of elements to replace
#' @param to the elements
#' @param digits integer indicating the number of decimal places to be used in
#' the comparison of double values (ignored if \code{is.double(x)} is FALSE)
#' @return a vector with replaced values
#' @export
#' @examples
#' x <- c(NA, c("a", "a", "z", "e", "q"))
#' 
#' # note that 'from' might contain values which are not present in 'x'
#' xr <- Replace(x, from = c("w", "a", "e"), to = c("ww", "aa", "ee"))
#' 
#' # only 'a' and 'e' were replaced
#' xr
#' stopifnot(identical(length(x), length(xr)))
#' stopifnot(identical(xr[c(2, 3, 5)], c("aa", "aa", "ee")))
#' 
#' # the missing value is not affected
#' stopifnot(is.na(xr[1]))
Replace <- function(x, from, to, digits = 10L) {
    if (!identical(typeof(x), typeof(from)) || 
        !identical(typeof(x), typeof(to)))
        stop("The type of 'x', 'from', and 'to' must be identical")
    checkVector(x, strict = TRUE)
    checkVector(from, strict = TRUE, any.missing = FALSE)
    checkVector(to, strict = TRUE, any.missing = FALSE)
    if (is.double(x)) {
        back <- TRUE
        x <- format(x, digits = digits)
        from <- format(from, digits = digits)
        to <- format(to, digits = digits)
    } else {
        back <- FALSE
    }
    un <- if (uniqueN(x) == length(x)) TRUE else FALSE
    from <- unique(from)
    to <- rep_len(to, length(from))
    if (un) {
        ind <- match(from, x)
        discard <- is.na(ind)
        x[ind[!discard]] <- to[!discard]
        if (back) x <- as.double(x)
        # return
        x
    } else {
        dt <- data.table(ind = seq_along(x), x = x)
        dt2 <- data.table(x = from, to = to)
        dt <- merge(dt, dt2, by = "x", all.x = TRUE)
        setkey(dt, to)
        dt[is.na(to), to := x]
        setkey(dt, ind)
        out <- dt$to
        if (back) out <- as.double(out)
        setattr(out, "names", names(x))
        # return
        out
    }
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

#' Create list with substituted names
#' 
#' \code{listS} creates a named list where names are substituted with the 
#' content of the referenced variable.
#' @param ... objects; if not named, listS is equilent to \code{\link{list}}. 
#' Names which should be substituted should start with a dot (.) or INDICES has 
#' to be provided. 
#' @param indices_ character or numeric vector indicating the position of those 
#' list elements whose name should be substituted. If provided, dotted names
#' are treated as original names and not substituted.
#' @export
#' @return A list with substituted names.
listS <- function(..., indices_ = NULL) {
    call_env <- parent.frame()
    subst <- function(x) {
        vapply(x, 
               function(xx) as.character(eval(parse(text = xx), 
                                              call_env)), 
               character(1))
    }
    #
    list_def <- list(...)
    if (is.null(onames <- names(list_def))) {
        return( list_def )
    }
    if (is.null(indices_)) {
        ind <- grep("^[.]", onames)
        onames[ind] <- subst(sub("^[.]", "", onames[ind]))
    } else {
        onames[indices_] <- subst(onames[indices_])
    }
    names(list_def) <- onames
    # return
    list_def
}