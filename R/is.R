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

#' Negation for the is* family of functions
#' 
#' \code{`!.IsFunction`} implements logical negation for \code{is*} function. 
#' The same effect can be achieved be setting the 'negate.' argument of the 
#' corresponding function.
#' @name negate
#' @aliases !
#' @param fn a function of class IsFunction
#' @export
#' @keywords internal
"!.IsFunction" <- function(fn) {
    function(...) !fn(...)
}

#' Internal function which combines is* functions
#' 
#' \code{isCombineOperator} is the workhorse function behind logical operators
#' for objects of class IsFunction.
#' @param lhs the function (of class "IsFunction") on the left-hand side
#' @param rhs the function (of class "IsFunction") on the right-hand side
#' @param op the operator (\code{`&`}, \code{`|`} or \code{xor})
#' @keywords internal
isCombineOperator <- function(lhs, rhs, op) {
    strict_l <- isStrict(lhs)
    strict_r <- isStrict(rhs)
    structure(
        function(x) {
            out_l <- lhs(x)
            out_r <- rhs(x)
            out <- op(out_l, out_r)
            if ((strict_l & strict_r) || (any(out))) {
                out
            } else if (strict_l) {
                out_l
            } else if (strict_r) {
                out_r
            } else {
                out[] <- TRUE
                out
            }
        },
        .MUST = if (strict_l | strict_r) TRUE else FALSE,
        class = "IsFunction"
    )
}


#' Combination of is* conditions
#' 
#' Logical tests as defined by \code{is*} functions can be combined by the 
#' standard logical operators \code{\link{&}}, \code{\link{|}}, and
#' \code{\link{xor}}.
#' @name combine
#' @param x,y two logical statements (\code{\link{is}} functions) which should
#' be combined
#' @aliases & &.IsFunction | |.IsFunction xor xor.IsFunction
#' @details If both tests are strict tests (that is, they were called by 
#' setting 'strict.' to TRUE [the default], which can be checked by 
#' \code{\link{isStrict}}), the tests are combined in the standard fashion.\cr
#' However, if at least one of the tests is optional (strict. = FALSE), the 
#' combined result is going to be returned by the new test only if it contains 
#' at least one TRUE value. Otherwise, it will return the result of the strict 
#' test alone, or if both tests are optional, only FALSE values will be 
#' returned.\cr
#' Note that \code{xor} has been made generic and a new \code{xor.default}
#' method was added besides \code{xor.IsFunction}. The new default method
#' is much faster than \code{base::xor} and is authored by Jens Oehlschlagel
#' (see \code{?xor} in package \bold{bit}).
#' @examples
#' # define two strict constraints and combine them
#' check_above_10 <- isBetween(lwr. = 10, open. = "lwr")
#' check_local_max <- isLocalMaximum()
#' check_local_max_above_10 <- check_local_max & check_above_10 
#' 
#' # call the combined test on the following vector
#' vec <- c(0, 1, 0, 11, 0, 21, 0)
#' ( res <- check_local_max_above_10(vec) )
#' stopifnot(identical(
#'     c(res),
#'     c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE)
#' ))
#'
#' # make the above_10 constraint optional
#' check_maybe_above_10 <- isBetween(lwr. = 10, open. = "lwr", strict. = FALSE)
#' check_local_max_maybe_above_10 <- check_local_max & check_maybe_above_10
#' 
#' # call the combined test on the vector
#' ( res2 <- check_local_max_maybe_above_10(vec) )
#' 
#' # the results are identical, because the combined test is TRUE for at least
#' # one data point
#' stopifnot(identical(
#'     c(res), c(res2)
#' ))
#' 
#' # now change the test vector with no local maximum above 10, 
#' # and run the test
#' vec3 <- c(0, 1, 0)
#' ( res3 <- check_local_max_maybe_above_10(vec3) )
#' 
#' # here, the optional test is ignored
#' stopifnot(identical(
#'     res3, check_local_max(vec3)
#' ))
#' 
#' 
NULL

#' @rdname combine
#' @export
`&.IsFunction` <- function(x, y) {
    isCombineOperator(x, y, `&`)
}

#' @rdname combine
#' @export
`|.IsFunction` <- function(x, y) {
    isCombineOperator(x, y, `|`)
}

#' @rdname combine
#' @export
xor.IsFunction <- function(x, y) {
    isCombineOperator(x, y, xor)
}

#' @rdname combine
#' @export
xor <- function(x, y) UseMethod(x, y)

#' @rdname combine
#' @export
xor.default <- function(x, y) as.logical(x) != as.logical(y)



#' Check if an is* function is strict or optional
#' 
#' \code{isStrict} is a helper function to check whether a function of class
#' \code{IsFunction} (see \code{\link[eegR]{is}}) is strict or optional. This
#' attribute is important if the tests are combined, see \code{\link{combine}}.
#' @param ... the functions to be tested. Can be a single list of functions as
#' well.
#' @seealso
#' \code{\link[eegR]{is}}, \code{\link{combine}}
#' @export
isStrict <- function(...) {
    fn_list <- list(...)
    if (is.list(fn_list[[1L]])) {
        if (length(fn_list) > 1L) {
            stop(paste0("if the functions are provided as a list of ",
                        "functions, only one such list is allowed"))
        }
        fn_list <- fn_list[[1L]]
    }
    # return
    vapply(fn_list,
           function(x) {
               if (!inherits(x, "IsFunction"))
                   stop("all input objects must be of class IsFunction")
               attr(x, ".MUST")
           },
           logical(1L))
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
            if (!identical(len_x %% len, 0L)) {
                warning(paste0("and the new length is not a multiple of the ",
                               "original length"),
                        call. = FALSE)
            }
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


#' Expansion in is* functions, with a fallback to repLen
#' 
#' \code{isExpandInto} is a wrapper around \code{expandInto} which is used
#' in the is* functions to transform the arguments like ref., tol., etc. 
#' to match the input array.
#' @keywords internal
isExpandInto <- function(from, to, from_name) {
    out <- try(expandInto(from, to, fill = FALSE), silent = TRUE)
    if (inherits(out, "try-error")) {
        warning(sprintf(
            paste0("The expansion of the argument '%s' ",
                   "as an array was not successful. Probably ",
                   "the dimension names or identifiers do not match ",
                   "between '%s' and 'x'. Instead of treating '%s' ",
                   "as an array, it is regarded as a vector in the ",
                   "followings. Note that this can lead to painful ",
                   "bugs later, so it is better to check the arguments. "),
            from_name, from_name, from_name))
        out <- repLen(from, length(to), from_name)
    }
    # return
    out
}


# < main functions > ------

#' Test logical statements on data points
#' 
#' \code{is*} are functions which produce functions of class \code{IsFunction} 
#' to test logical statements on each data point of an atomic object (vector, 
#' matrix, or array). See Details for the general idea, Functions for the short 
#' function-specific descriptions and the section Use cases for short examples.
#' @name is
#' @param strict. logical value whether the tested condition is obligatory
#' (TRUE, the default) or optional (FALSE). Only used for combining multiple
#' tests (not available at the moment).
#' @param negate. logical value whether the return value of the function should
#' be negated (default: FALSE). You can also use the standard \code{\link{!}} 
#' operator for negation which has a specific method for the \code{IsFunction} 
#' class (see Examples).
#' @param ... further arguments passed to the wrapped function (see Functions
#' for details)
#' @details All \code{is*} functions return a function of the class 
#' \code{IsFunction} which has only one argument, '\code{x}' (the data object 
#' to test). '\code{x}' can be an atomic vector, matrix, or array; all other 
#' object types result in error.\cr
#' A second major rule is that all functions returned by an \code{is*}
#' function return a logical object of the same shape as the input object. 
#' The third rule is that the function returned by any \code{is*} function is
#' a hard (strict) or soft (optional) constraint. This affects the way how the 
#' results are combined for joint logical tests (see \code{\link{combine}}). 
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
#' @seealso
#' \code{\link{combine}} how to combine is* tests, and \code{\link{isStrict}}
#' how to check afterwards whether a test is strict or optional
NULL

#' @describeIn is produces a function to test whether the values in an atomic 
#' character object match a pre-specified character pattern. See 
#' \code{\link{grepl}} for further details and additional arguments.
#' @param pattern. a character regexp pattern, see the 'pattern' argument in 
#' \code{\link{grepl}}
#' @export
#' @examples
#' # example for isPattern
#' # note how we pass the 'ignore.case' argument to the internally used
#' # grepl() function (see ?grepl)
#' check_start_a <- isPattern("^a", ignore.case = TRUE)
#' 
#' # check the results -> note that the class and the attribute has to be 
#' # removed by calling 'c' (or 'as.vector') to make the results identical
#' stopifnot(identical(
#'     c(check_start_a(c("a", "A", "ba"))), 
#'     c(TRUE, TRUE, FALSE)
#' ))
#' 
isPattern <- function(pattern., strict. = TRUE, negate. = FALSE, 
                      subset. = list(), ...) {
    assertString(pattern., na.ok = FALSE, .var.name = "pattern.")
    options <- list(...)
    assertList(options, names = "unique", .var.name = "arguments in '...'")
    structure(
        function(x) {
            tempfn <- function(dat, pattern., options) {
                out <- do("grepl", pattern., as.vector(dat), 
                          arg_list = options)
                setattr(out, "dim", dim(dat))
                setattr(out, "dimnames", dimnames(dat))
            }
            #
            assertCharacter(x, .var.name = "x")
            out <- 
                if (length(dim(x)) > 0L && length(subset.) > 0L) {
                    expandInto(
                        tempfn(subsetArray(x, subset.), pattern., options),
                        x,
                        fill = FALSE
                    )
                } else {
                    tempfn(x, pattern., options)    
                }
            if (!negate.) out else !out
        }, .MUST = strict., class = "IsFunction")
}

#' @describeIn is produces a function to test whether the values in an atomic 
#' object are equal to a (vector of) reference value(s) with a given tolerance.
#' The returned function accepts objects of any type for which \code{asNumeric} 
#' does not fail (see the documentation of \code{\link{asNumeric}}).
#' @param ref. the reference value or a vector (or matrix/array) of reference 
#' values. 
#' @param tol. tolerance value, the maximum difference which is still tolerated.
#' Can be a vector (or matrix/array) as well.
#' @export
#' @examples
#' #
#' # example for isEqual
#' check_zero <- isEqual(0, tol. = 1e-4)
#' stopifnot(identical(
#'     c(check_zero(c("0", "0.000001", "1"))), 
#'     c(TRUE, TRUE, FALSE)
#' ))
#' 
#' # note that just like for other is* functions, the reference and the 
#' # tolerance can be vectors (or even matrices or arrays), not only scalars,
#' # and its values are recycled or cropped if necessarily to match the length
#' # of the 'x' object
#' check_values <- isEqual(c(0, 1, 2))
#' stopifnot(identical(
#'     c(check_values(c(0, 1.1, 2, 1e-9))), 
#'     c(TRUE, FALSE, TRUE, TRUE)
#' ))
#' 
#' # it is also possible to base the check only on a subset of x and expand 
#' # the results back to the original array
#' check_width <- isEqual(1, tol. = 0.1, subset. = list(measure = "width"))
#' ( x <- matrix(c(1.001, 1.6, 1.2, 1.003), 2, 2, 
#'               dimnames = list(observation = c("a", "b"), 
#'                               measure = c("width", "height"))) )
#' ( res <- check_width(x) )
#' stopifnot(identical(dim(res), dim(x)))
#' stopifnot(identical(
#'     c(res),
#'     c(TRUE, FALSE, TRUE, FALSE)
#' ))
#' 
#' # isEqual, just like all other is* functions, can be negated:
#' check_not_zero <- isEqual(0, negate. = TRUE)
#' stopifnot(check_not_zero(1))
#' 
#' # the same can be achieved by using ! for negation
#' check_not_zero2 <- !isEqual(0)
#' stopifnot(identical(check_not_zero(1), check_not_zero2(1)))
#' 
isEqual <- function(ref., tol. = .Machine$double.eps^0.5, 
                    strict. = TRUE, negate. = FALSE, subset. = list(), ...) {
    ref. <- asNumeric_(ref.)
    tol. <- asNumeric_(tol.)
    assertLogical(strict., len = 1L, any.missing = FALSE, 
                  .var.name = "strict.")
    assertLogical(negate., len = 1L, any.missing = FALSE, 
                  .var.name = "negate.")
    structure(
        function(x) {
            tempfn <- function(dat, ref., tol.) {
                len_dat <- length(dat)
                ref. <-
                    if (is.array(ref.)) {
                        isExpandInto(ref., dat, "ref.")
                    } else {
                        repLen(ref., len_dat, "ref.")
                    }
                tol. <-
                    if (is.array(tol.)) {
                        isExpandInto(tol., dat, "tol.")
                    } else {
                        repLen(tol., len_dat, "tol.")
                    }    
                abs(asNumeric_(dat) - ref.) <= tol.
            }
            #
            out <- 
                if (length(dim(x)) > 0L && length(subset.) > 0L) {
                    expandInto(
                        tempfn(subsetArray(x, subset.), ref., tol.),
                        x,
                        fill = FALSE
                    )
                } else {
                    tempfn(x, ref., tol.)
                }
            if (!negate.) out else !out
        }, .MUST = strict., class = "IsFunction")
}


#' @describeIn is produces a function to test whether the values in an atomic 
#' object are the same as a (vector of) reference value(s). It is a wrapper
#' around \code{\link{==}}. If you want to test numeric values, \code{isEqual}
#' is a much better alternative.
#' @export
#' @examples
#' #
#' # example for isSame;
#' check_not_A <- isSame("A", negate. = TRUE)
#' stopifnot(identical(
#'     c(check_not_A(c("a", "A"))), 
#'     c(TRUE, FALSE)
#' ))
#' 
isSame <- function(ref., strict. = TRUE, negate. = FALSE, 
                   subset. = list(), ...) {
    structure(
        function(x) {
            tempfn <- function(dat, ref., negate.) {
                ref. <-
                    if (is.array(ref.)) {
                        isExpandInto(ref., dat, "ref.")
                    } else {
                        repLen(ref., length(dat), "ref.")
                    }
                if (!negate.) dat == ref. else dat != ref.
            }
            #
            if (length(dim(x)) > 0L && length(subset.) > 0L) {
                expandInto(
                    tempfn(subsetArray(x, subset.), ref., negate.),
                    x,
                    fill = FALSE
                )
            } else {
                tempfn(x, ref., negate.)
            }
        }, .MUST = strict., class = "IsFunction")
}


#' @describeIn is produces a function to test whether
#' the values in an atomic object are between pre-specified limits.
#' The returned function accepts objects of any type for which \code{asNumeric} 
#' does not fail (see the documentation of \code{\link{asNumeric}}).
#' @param lwr.,upr. numeric values referring to the lower and upper limit
#' (the defaults are -Inf and Inf, respectively). Both lwr. and upr. can be
#' scalars and vectors as well.
#' @param open. the side of the interval which is open (that is, it does not
#' contain the given endpoint); can be "none" (the default), "both", lwr" or 
#' "upr"
#' @export
#' @examples
#' #
#' # example for isBetween();
#' # note that lwr. and upr. might be given in reversed order
#' check_0_100 <- isBetween(100, 0)
#' stopifnot(identical(
#'     c(check_0_100(c(-1, 1, 101))), 
#'     c(FALSE, TRUE, FALSE)
#' ))
#' 
isBetween <- function(lwr. = -Inf, upr. = Inf, 
                      open. = c("none", "both", "lwr", "upr"),
                      strict. = TRUE, negate. = FALSE, subset. = list(), ...) {
    lwr. <- asNumeric_(lwr.)
    upr. <- asNumeric_(upr.)
    open. <- match.arg(open.)
    assertLogical(strict., any.missing = FALSE, len = 1L, .var.name = "strict")
    #
    lwr <- lwr.; upr <- upr.
    lwr. <- pmin(lwr, upr)
    upr. <- pmax(lwr, upr)
    lwr <- NULL; upr <- NULL
    #
    structure(
        function(x) {
            tempfn <- function(dat, lwr., upr., open.) {
                len_dat <- length(dat)
                out <- rep_len(TRUE, len_dat)
                if (!identical(lwr., -Inf)) {
                    lwr. <-
                        if (is.array(lwr.)) {
                            isExpandInto(lwr., dat, "lwr.")
                        } else {
                            repLen(lwr., len_dat, "lwr.")
                        }
                    if (open. %in% c("none", "upr")) `>` <- `>=`
                    out <- out & dat > lwr.
                }
                if (!identical(upr., -Inf)) {
                    upr. <-
                        if (is.array(upr.)) {
                            isExpandInto(upr., dat, "upr.")
                        } else {
                            repLen(upr., len_dat, "lwr.")
                        }
                    if (open. %in% c("none", "lwr")) `<` <- `<=`
                    out <- out & dat < upr.
                }
                setattr(out, "dim", dim(dat))
                setattr(out, "dimnames", dimnames(dat))
            }
            #
            x <- asNumeric_(x)
            out <- 
                if (length(dim(x)) > 0L && length(subset.) > 0L) {
                    expandInto(
                        tempfn(subsetArray(x, subset.), lwr., upr., open.),
                        x,
                        fill = FALSE
                    )
                } else {
                    tempfn(x, lwr., upr., open.)
                }
            if (!negate.) out else !out
        }, .MUST = strict., class = "IsFunction")
}

#' @describeIn is produces a function to test whether the values in an atomic 
#' object are negative. The returned function accepts 
#' objects of any type for which \code{asNumeric} 
#' does not fail (see the documentation of \code{\link{asNumeric}}).
#' @export
isNegative <- function(strict. = TRUE, negate. = FALSE, ...) {
    isBetween(-Inf, 0, open. = "upr", strict. = strict., 
              negate. = FALSE, ...)
}

#' @describeIn is produces a function to test whether the values in an atomic 
#' object are positive. The returned function accepts 
#' objects of any type for which \code{asNumeric} 
#' does not fail (see the documentation of \code{\link{asNumeric}}).
#' @export
isPositive <- function(strict. = TRUE, ...) {
    isBetween(0, Inf, open. = "lwr", strict. = strict.,
              negate. = FALSE, ...)
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
#'     c(check_global_extr(c(-1, 1, 0, 100, 50))), 
#'     c(TRUE, FALSE, FALSE, TRUE, FALSE)
#' ))
#' 
isExtremum <- function(what. = c("both", "min", "max"), strict. = TRUE, 
                       negate. = FALSE, ...) {
    what. <- match.arg(what.)
    assertLogical(strict., any.missing = FALSE, len = 1L, .var.name = "strict")
    options <- list(...)
    assertList(options, names = "unique", .var.name = "arguments in '...'")
    # return
    structure(
        function(x) {
            out <- do("findExtrema", x, arg_list = options)
            out <- switch(what., 
                          max = out == 1L, 
                          min = out == -1L, 
                          both = (out == 1L) | (out == -1L))
            if (!negate.) out else !out
        }, .MUST = strict., class = "IsFunction")
}

#' @describeIn is a shorthand for \code{isExtremum("max", strict., ...)}. See 
#' \code{\link{findExtrema}} for additional details and further arguments. 
#' @export
isLocalMaximum <- function(strict. = TRUE, negate. = FALSE, ...) {
    isExtremum("max", strict., negate., ...)
}

#' @describeIn is a shorthand for \code{isExtremum("min", strict., ...)}. See 
#' \code{\link{findExtrema}} for additional details and further arguments. 
#' @export
isLocalMinimum <- function(strict. = TRUE, negate. = FALSE, ...) {
    isExtremum("min", strict., negate., ...)
}

#' @describeIn is a shorthand for 
#' \code{isExtremum("max", strict., global = TRUE, ...)}. See 
#' \code{\link{findExtrema}} for additional details and further arguments. 
#' @export
isMaximum <- function(strict. = TRUE, negate. = FALSE, ...) {
    isExtremum("max", strict., negate., global = TRUE, ...)
}

#' @describeIn is a shorthand for 
#' \code{isExtremum("min", strict., global = TRUE, ...)}. See 
#' \code{\link{findExtrema}} for additional details and further arguments. 
#' @export
isMinimum <- function(strict. = TRUE, negate. = FALSE, ...) {
    isExtremum("min", strict., negate., global = TRUE, ...)
}




