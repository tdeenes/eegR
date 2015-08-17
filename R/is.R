#
# <<< convenience functions to test logical statements >>> --------
#


# < helper functions > ------

#' Coerce to numeric
#' 
#' \code{asNumeric} is similar to \code{\link{as.numeric}} except that 1) 
#' it results in an error for character strings for which as.numeric would 
#' send a warning, 2) keeps integer, Date and POSIXt objects as they are (i.e. 
#' does not convert them to double), 3) keeps 'dim' and 'dimnames' 
#' attributes if requested, and 4) always copies the original object, even if no
#' actual coercion occurs. \code{asNumeric_} returns the original object if 
#' possible and does not alter the shape of multidimensional inputs.
#' @param x the object to coerce
#' @param keep_dim logical value whether \code{dim} and \code{dimnames}
#' attributes should be preserved
#' @return \code{asNumeric} returns the original object (at the same memory 
#' address) if it is numeric (that is, is.numeric(x) returns TRUE) or inherits 
#' "Date" or "POSIXt" classes. Otherwise, it returns a double vector if 
#' \code{as.numeric} is successful and stops with an error if 
#' \code{as.numeric} results in a warning or error.
#' @export
#' @examples
#' # if 'x' is double, asNumeric() simply returns it, just like as.numeric()
#' x <- rnorm(10)
#' stopifnot(identical(as.numeric(x), asNumeric(x)))
#' 
#' # if 'x' is integer, asNumeric() simply returns it, 
#' # but as.numeric() converts to double
#' x <- 1:10
#' stopifnot(identical(x, asNumeric(x)))
#' stopifnot(all.equal(as.numeric(x), asNumeric(x))) # not fully identical!
#' 
#' # if 'x' is Date or POSIXt, asNumeric() simply returns it
#' x1 <- as.Date("2013-03-14")
#' x2 <- as.POSIXlt(x1)
#' stopifnot(identical(x1, asNumeric(x1)))
#' stopifnot(identical(x2, asNumeric(x2)))
#' 
#' # if 'x' contains only numeric-like characters (or NAs),
#' # asNumeric() returns the same as as.numeric() 
#' x <- c(NA, as.character(1:10))
#' stopifnot(identical(as.numeric(x), asNumeric(x)))
#' 
#' # if 'x' contains any non-numeric characters (e.g. letters),
#' # asNumeric() returns an error
#' x <- c("a", as.character(1:10))
#' x_num <- try(asNumeric(x), silent = TRUE)
#' stopifnot(inherits(x_num, "try-error"))
#' 
#' # if 'x' is a matrix or array, asNumeric() converts it to a vector by 
#' # default, just like as.numeric()
#' ( y <- x <- matrix(rnorm(10), 5, 2) )
#' ( x_num <- asNumeric(x) )
#' stopifnot(is.null(attr(x_num, "dim")))
#' 
#' # however, you can ask to keep the original shape
#' ( x_num2 <- asNumeric(x, TRUE) )
#' 
#' # if you use asNumeric_() instead, it keeps the shape
#' ( x_num3 <- asNumeric_(x) )
#' stopifnot(identical(x_num3, x))
#' 
#' # note that since 'x' is numeric, no coercion is needed, so 
#' # so asNumeric_() does not create a copy, but asNumeric() does
#' address(x)
#' address(x_num2)
#' address(x_num3)
#' stopifnot(identical(address(x), address(x_num3)))
#' 
asNumeric <- function(x, keep_dim = FALSE) {
    dims <- dim(x)
    dimn <- dimnames(x)
    x <- asNumeric_(copy(x))
    if (!keep_dim) {
        setattr(x, "dim", NULL)
        setattr(x, "dimnames", NULL)
    } else if (is.null(dim(x)) && !is.null(dims)) {
        setattr(x, "dim", dims)
        setattr(x, "dimnames", dimn)
    }
    x
}

#' @describeIn asNumeric
#' @export
asNumeric_ <- function(x) {
    if (is.numeric(x) || inherits(x, c("POSIXt", "Date"))) {
        x
    } else {
        if (!is.numeric(x)) {
            ._ok <- TRUE
            x <- tryCatch(as.numeric(x), 
                          warning = function(w) ._ok <<- FALSE,
                          error = function(e) ._ok <<- FALSE)
            if (!._ok) {
                stop(paste0("if 'x' is not numeric or date, it must be ",
                            "numeric-like character, e.g. x = c('1', '2') ", 
                            "or any other type for which as.numeric() returns ",
                            "numeric values"))
            }
        }
        x
    }
}

#' Match the length of an object to an other object with warning
#' 
#' \code{repLen} is a wrapper around \code{rep_len} to recycle or crop a
#' vector so that its length matches a pre-defined length. It throws a warning
#' in such cases.
#' @param x the object
#' @param len the desired length
#' @keywords internal
repLen <- function(x, len, x_name = NULL) {
    len_x <- length(x)
    if (is.null(x_name)) x_name <- deparse(substitute(x))
    if (len_x > 1L) {
        if (len_x < len) {
            warning(sprintf("'%s' recycled to be of length %i", x_name, len),
                    call. = FALSE)
            rep_len(x, len)
        } else if (len_x > len) {
            warning(sprintf("'%s' cropped to be of length %i", x_name, len),
                    call. = FALSE)
            rep_len(x, len)
        } else {
            x
        }
    } else {
        x
    }
}


# < main functions > ------

#' Test logical statements on data points
#' 
#' \code{is*} are functions which produce functions to test logical statements 
#' on each data point of an atomic object (vector, matrix, or array). See 
#' Details for the general idea, Functions for the short function-specific
#' descriptions and the section Use cases for short examples.
#' @name is
#' @param strict. logical value whether the tested condition is obligatory
#' (TRUE, the default) or optional (FALSE). Only used for combining multiple
#' tests (not available at the moment).
#' @param ... further arguments passed to the wrapped function (see Functions
#' for details)
#' @details All \code{is*} functions return a function which has only one
#' argument, \code{x} (the data object to test). A second major rule is that 
#' all functions returned by an \code{is*} function return a logical object
#' of the same shape as the input object. The third rule is that the function 
#' returned by any \code{is*} function has an attribute called \code{.MUST}.
#' This affects the way how the results are combined for joint logical tests
#' (more on this later). 
#' \cr
#' @section Use cases:
#' \subsection{\code{is*} as simple pre-defined test}{
#' A very basic use case is to pre-define a logical rule and apply it to several
#' objects. For example one can define\cr
#' \code{rangeTest <- isBetween(200, 300)}\cr
#' and use this rule to check for whatever numeric or numeric-like character 
#' vector, matrix or array whether the values are between 200 and 300. If 
#' \code{x} denotes the object, \code{rangeTest(x)} returns a logical object 
#' of the same length and shape as \code{x} with TRUEs for values between 200 
#' and 300 and FALSE otherwise.
#' }
#' \subsection{\code{is*} as function argument}{
#' \code{is*} functions come in especially handy if they are used as function 
#' arguments.\cr 
#' For example the code \code{subsetArray(erps, time = isBetween(200, 300))} is
#' a very compact and readable way of subsetting the \code{erps} data array
#' on its \code{time} dimension, selecting only that part of the array where 
#' \code{time >= 200 & time <= 300}. Note that this allows the definition
#' of a subsetting rule without knowing in advance how the object which is to
#' be subsetted looks like. For example you can define 
#' \code{sub_def <- isBetween(200, 300)} beforehand and use \code{sub_def} as
#' an argument in \code{subsetArray} for all data arrays if all of them have a
#' \code{time} dimension measured in the same unit (note that the resolution and 
#' the time range may be different).
#' }
#' 
NULL

#' @describeIn is produces a function to test whether the values in an atomic 
#' object match a pre-specified character pattern. See \code{\link{grepl}} for 
#' further details and additional arguments it accepts.
#' @param pattern. a character regexp pattern, see the 'pattern' argument in 
#' \code{\link{grepl}}
#' @export
#' @examples
#' # example for isPattern
#' # note how we pass the 'ignore.case' argument to the internally used
#' # grepl() function (see ?grepl)
#' check_start_a <- isPattern("^a", ignore.case = TRUE)
#' stopifnot(identical(
#'     check_start_a(c("a", "A", "ba")), 
#'     structure(c(TRUE, TRUE, FALSE), .MUST = TRUE)
#' ))
#' 
isPattern <- function(pattern., strict. = TRUE, ...) {
    assertString(pattern., na.ok = FALSE, .var.name = "pattern")
    options <- list(...)
    assertList(options, names = "unique", .var.name = "arguments in '...'")
    function(x) {
        assertCharacter(x, .var.name = "x")
        out <- do("grepl", pattern., as.vector(x), arg_list = options)
        setattr(out, "dim", dim(x))
        setattr(out, "dimnames", dimnames(x))
        setattr(out, ".MUST", strict.)
        out
    }
}

#' @describeIn is produces a function to test whether the values in an atomic 
#' object are equal to a (vector of) reference value(s) with a given tolerance.
#' The returned function accepts objects of any type for which \code{asNumeric} 
#' does not fail (see the documentation of \code{\link{asNumeric}}).
#' @param ref. the reference value or a vector (or matrix/array) of reference 
#' values. 
#' @param tol. tolerance
#' @export
#' @examples
#' #
#' # example for isEqual
#' check_zero <- isEqual(0, tol. = 1e-4)
#' stopifnot(identical(
#'     check_zero(c("0", "0.000001", "1")), 
#'     structure(c(TRUE, TRUE, FALSE), .MUST = TRUE)
#' ))
#' 
#' # note that just like for other is* functions, the reference and the 
#' # tolerance can be vectors (or even matrices or arrays), not only scalars,
#' # and its values are recycled or cropped if necessarily to match the length
#' # of the 'x' object
#' check_values <- isEqual(c(0, 1, 2), strict. = FALSE)
#' stopifnot(identical(
#'     check_values(c(0, 1.1, 2, 1e-9)), 
#'     structure(c(TRUE, FALSE, TRUE, TRUE), .MUST = FALSE)
#' ))
#' 
isEqual <- function(ref., tol. = .Machine$double.eps^0.5, 
                     strict. = TRUE, ...) {
    ref. <- asNumeric_(ref.)
    tol. <- asNumeric_(tol.)
    function(x) {
        len_x <- length(x)
        ref. <- repLen(ref., len_x, "ref.")
        tol. <- repLen(tol., len_x, "tol.")
        x <- asNumeric_(x)
        out <- abs(x - ref.) <= tol.
        setattr(out, ".MUST", strict.)
        out
    }
}

#' @describeIn is produces a function to test whether
#' the values in an atomic object are between pre-specified limits.
#' The returned function accepts objects of any type for which \code{asNumeric} 
#' does not fail (see the documentation of \code{\link{asNumeric}}).
#' @param lwr.,upr. numeric values referring to the lower and upper limit
#' (the defaults are -Inf and Inf, respectively). Both lwr. and upr. can be
#' scalars and vectors as well.
#' @param closed. the limit which is closed; can be "both" (the default), 
#' "none", lwr" or "upr"
#' @export
#' @examples
#' #
#' # example for isBetween();
#' # note that lwr_ and upr_ might be given in reversed order
#' check_0_100 <- isBetween(100, 0)
#' stopifnot(identical(
#'     check_0_100(c(-1, 1, 101)), 
#'     structure(c(FALSE, TRUE, FALSE), .MUST = TRUE)
#' ))
#' 
isBetween <- function(lwr. = -Inf, upr. = Inf, 
                       closed. = c("both", "none", "lwr", "upr"),
                       strict. = TRUE, ...) {
    lwr. <- asNumeric_(lwr.)
    upr. <- asNumeric_(upr.)
    closed. <- match.arg(closed.)
    assertLogical(strict., any.missing = FALSE, len = 1L, .var.name = "strict")
    #
    lwr <- lwr.; upr <- upr.
    lwr. <- pmin(lwr, upr)
    upr. <- pmax(lwr, upr)
    lwr <- NULL; upr <- NULL
    #
    function(x) {
        x <- asNumeric_(x)
        len_x <- length(x)
        out <- rep_len(TRUE, len_x)
        if (!identical(lwr., -Inf)) {
            lwr. <- repLen(lwr., len_x, "lwr.")
            if (closed. %in% c("both", "lwr")) `>` <- `>=`
            out <- out & x > lwr.
        }
        if (!identical(upr., -Inf)) {
            upr. <- repLen(upr., len_x, "upr.")
            if (closed. %in% c("both", "upr")) `<` <- `<=`
            out <- out & x < upr.
        }
        setattr(out, "dim", dim(x))
        setattr(out, "dimnames", dimnames(x))
        setattr(out, ".MUST", strict.)
        out
    }
}

#' @describeIn is produces a function to test whether the values in an atomic 
#' object are negative. The returned function accepts 
#' objects of any type for which \code{asNumeric} 
#' does not fail (see the documentation of \code{\link{asNumeric}}).
#' @export
isNegative <- function(strict. = TRUE, ...) {
    isBetween(-Inf, -.Machine$double.eps^0.5, strict. = strict., ...)
}

#' @describeIn is produces a function to test whether the values in an atomic 
#' object are positive. The returned function accepts 
#' objects of any type for which \code{asNumeric} 
#' does not fail (see the documentation of \code{\link{asNumeric}}).
#' @export
isPositive <- function(strict. = TRUE, ...) {
    isBetween(.Machine$double.eps^0.5, Inf, strict. = strict., ...)
}

#' @describeIn is produces a function to test whether the the values in an 
#' atomic object are local or global extrema. See \code{\link{findExtrema}}
#' for additional details and further arguments. The returned function accepts 
#' objects of any type for which \code{asNumeric} 
#' does not fail (see the documentation of \code{\link{asNumeric}}).
#' @param what. specifies whether minima ("min"), maxima ("max") or
#' extrema ("both", the default) should be tested
#' @export
#' @examples
#' #
#' # example for isExtremum
#' # note how we can pass the 'global' and 'tail' arguments to the internal 
#' # findExtrema() function (see ?findExtrema)
#' check_global_extr <- isExtremum(global = TRUE, tail = "do_not_care")
#' stopifnot(identical(
#'     check_global_extr(c(-1, 1, 0, 100, 50)), 
#'     structure(c(TRUE, FALSE, FALSE, TRUE, FALSE), .MUST = TRUE)
#' ))
#' 
isExtremum <- function(what. = c("both", "min", "max"), strict. = TRUE, ...) {
    what. <- match.arg(what.)
    assertLogical(strict., any.missing = FALSE, len = 1L, .var.name = "strict")
    options <- list(...)
    assertList(options, names = "unique", .var.name = "arguments in '...'")
    # return
    function(x) {
        out <- do("findExtrema", x, arg_list = options)
        out <- switch(what., 
                      max = out == 1L, 
                      min = out == -1L, 
                      both = (out == 1L) | (out == -1L))
        setattr(out, ".MUST", strict.)
    }
}

#' @describeIn is a shorthand for \code{isExtremum("max", strict., ...)}. See 
#' \code{\link{findExtrema}} for additional details and further arguments. 
#' @export
isLocalMaximum <- function(strict. = TRUE, ...) {
    isExtremum("max", strict., ...)
}

#' @describeIn is a shorthand for \code{isExtremum("min", strict., ...)}. See 
#' \code{\link{findExtrema}} for additional details and further arguments. 
#' @export
isLocalMinimum <- function(strict. = TRUE, ...) {
    isExtremum("min", strict., ...)
}

#' @describeIn is a shorthand for 
#' \code{isExtremum("max", strict., global = TRUE, ...)}. See 
#' \code{\link{findExtrema}} for additional details and further arguments. 
#' @export
isMaximum <- function(strict. = TRUE, ...) {
    isExtremum("max", strict., global = TRUE, ...)
}

#' @describeIn is a shorthand for 
#' \code{isExtremum("min", strict., global = TRUE, ...)}. See 
#' \code{\link{findExtrema}} for additional details and further arguments. 
#' @export
isMinimum <- function(strict. = TRUE, ...) {
    isExtremum("min", strict., global = TRUE, ...)
}
