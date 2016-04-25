#
# <<< convenience functions to test logical statements >>> --------
#


# < helper functions > ------

#' Coerce to numeric
#'
#' \code{as.Numeric} is similar to \code{\link{as.numeric}} except that 1)
#' it results in an error for character strings for which as.numeric would
#' send a warning, 2) keeps integer, Date and POSIXt objects as they are (i.e.
#' does not convert them to double), 3) coerces logicals to integer, not double,
#' 4) keeps 'dim' and 'dimnames' attributes if requested, and 5) always copies
#' the original object, even if no actual coercion occurs. \code{as.Numeric_}
#' returns the original object if possible and does not alter the shape of
#' multidimensional inputs.
#' @param x the object to coerce
#' @param keep_dim logical value whether \code{dim} and \code{dimnames}
#' attributes should be preserved. The default is FALSE. (But note that
#' \code{as.Numeric_} does not have this argument, because it always preserves
#' both attributes.)
#' @return \code{as.Numeric} returns the original object (at the same memory
#' address) if it is numeric (that is, is.numeric(x) returns TRUE) or inherits
#' "Date" or "POSIXt" classes. Otherwise, it returns a double vector if
#' \code{as.numeric} is successful and stops with an error if
#' \code{as.numeric} results in a warning or error.
#' @export
#' @examples
#' # if 'x' is double, as.Numeric() simply returns it, just like as.numeric()
#' x <- rnorm(10)
#' stopifnot(identical(as.numeric(x), as.Numeric(x)))
#'
#' # if 'x' is integer, as.Numeric() simply returns it,
#' # but as.numeric() converts to double
#' x <- 1:10
#' stopifnot(identical(x, as.Numeric(x)))
#' stopifnot(all.equal(as.numeric(x), as.Numeric(x))) # not fully identical!
#'
#' # if 'x' is Date or POSIXt, as.Numeric() simply returns it
#' x1 <- as.Date("2013-03-14")
#' x2 <- as.POSIXlt(x1)
#' stopifnot(identical(x1, as.Numeric(x1)))
#' stopifnot(identical(x2, as.Numeric(x2)))
#'
#' # if 'x' contains only numeric-like characters (or NAs),
#' # as.Numeric() returns the same as as.numeric()
#' x <- c(NA, as.character(1:10))
#' stopifnot(identical(as.numeric(x), as.Numeric(x)))
#'
#' # if 'x' contains any non-numeric characters (e.g. letters),
#' # as.Numeric() returns an error
#' x <- c("a", as.character(1:10))
#' x_num <- try(as.Numeric(x), silent = TRUE)
#' stopifnot(inherits(x_num, "try-error"))
#'
#' # if 'x' is a matrix or array, as.Numeric() converts it to a vector by
#' # default, just like as.numeric()
#' ( y <- x <- matrix(rnorm(10), 5, 2) )
#' ( x_num <- as.Numeric(x) )
#' stopifnot(is.null(attr(x_num, "dim")))
#'
#' # however, you can ask to keep the original shape
#' ( x_num2 <- as.Numeric(x, TRUE) )
#'
#' # if you use as.Numeric_() instead, it keeps the shape
#' ( x_num3 <- as.Numeric_(x) )
#' stopifnot(identical(x_num3, x))
#'
#' # note that since 'x' is numeric, no coercion is needed, so
#' # so as.Numeric_() does not create a copy, but as.Numeric() does
#' address(x)
#' address(x_num2)
#' address(x_num3)
#' stopifnot(identical(address(x), address(x_num3)))
#'
as.Numeric <- function(x, keep_dim = FALSE) {
    dims <- dim(x)
    dimn <- dimnames(x)
    x <- as.Numeric_(copy(x))
    if (!keep_dim) {
        setattr(x, "dim", NULL)
        setattr(x, "dimnames", NULL)
    } else if (is.null(dim(x)) && !is.null(dims)) {
        setattr(x, "dim", dims)
        setattr(x, "dimnames", dimn)
    }
    x
}

#' @describeIn as.Numeric Modify by reference
#' @export
as.Numeric_ <- function(x) {
    if (is.numeric(x) || inherits(x, c("POSIXt", "Date"))) {
        x
    } else if (is.logical(x)) {
        storage.mode(x) <- "integer"
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
        rep_len(x, len)
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

#' Get 'dim_index' of the subset. and expand. arguments in the is* functions
#'
#' \code{checkIsSubsetExpand} is a tiny helper function to check the
#' 'dim_index' attributes of the 'subset.' and 'expand.' arguments in
#' \code{is*} functions. It also checks for potential dimension overlaps.
#' @keywords internal
checkSubsetExpand <- function(subset., expand.) {
    # helper function
    checkFn <- function(x, var_name) {
        if (length(x) > 0) {
            if (!is.list(x)) {
                stop(paste0(
                    "'", var_name, "' must be a list (or NULL)"
                    ), call. = FALSE)
            }
            if (is.null(names(x))) {
                out <- attr(x, "dim_index")
                if (is.null(out)) {
                    stop(paste0(
                        "if '", var_name, "' is not named, it must ",
                        "have a non-NULL 'dim_index' attribute"
                        ), call. = FALSE)
                }
                out
            } else {
                names(x)
            }
        } else {
            NULL
        }
    }
    # main
    args <- list(subset. = subset., expand. = expand.)
    dim_index <- setNames(
        lapply(names(args), function(x) checkFn(args[[x]], x)),
        names(args))
    if (any(vapply(dim_index, length, 1L) == 0L)) {
        dim_index
    } else if (!do.call("identical", unname(lapply(dim_index, typeof)))) {
        stop(paste0(
            "both 'subset.' and 'expand.' must be either named lists ",
            "or unnamed lists with 'dim_index' attribute"
            ), call. = FALSE)
    } else if (anyDuplicated(unlist(dim_index, use.names = FALSE))) {
        stop("dimensions in 'subset.' and 'expand.' may not overlap",
             call. = FALSE)
    } else {
        lapply(dim_index, function(x) if (is.character(x)) NULL else x)
    }
}

#' Call the internal function in is* functions on subsets and expand
#'
#' \code{isMainCompFn} calls the internal function in the given \code{is*}
#' function on the input vector or on the appropriate slices of the input
#' array.
#' @param FUN. the internal function whose first two arguments are the data (x)
#' and the options (options.). Additional arguments must be provided in dots
#' (...).
#' @param subset.,expand.,which.dims.,options.,negate. obligatory arguments, see
#' \code{is}
#' @param ... arguments passed to FUN.
#' @keywords internal
isMainCompFn <- function(FUN., x, subset., expand.,
                         which_dims., options., negate., ...) {
    do_sub <- length(subset.) > 0L
    do_exp <- length(expand.) > 0L
    out <-
        if (length(dim(x)) > 0L && (do_sub | do_exp)) {
            if (do_sub) {
                out <- array_(FALSE, dim(x), dimnames(x))
                x_sub <- subsetArray(x, subset., which_dims.$subset.,
                                     drop. = FALSE)
                subsetArray(out, subset., which_dims.$subset., drop. = FALSE) <-
                    if (length(expand.) > 0L) {
                        expandInto(
                            FUN.(subsetArray(x_sub,
                                             expand.,
                                             which_dims.$expand.,
                                             drop. = FALSE),
                                 options., negate., ...),
                            x_sub,
                            fill = FALSE
                        )
                    } else {
                        FUN.(x_sub, options., negate., ...)
                    }
                out
            } else if (do_exp) {
                expandInto(
                    FUN.(subsetArray(x,
                                     expand.,
                                     which_dims.$expand.,
                                     drop. = FALSE),
                         options., negate., ...),
                    x,
                    fill = FALSE
                )
            }
        } else {
            FUN.(x, options., negate., ...)
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
#' @param subset.,expand. a named list of character, numeric, or logical
#' vectors or a subsetting function (of class IsFunction) indicating which
#' levels of which dimensions to subset or expand on (see Details). If an
#' unnamed list, it must have an attribute 'dim_index', referring to the numeric
#' indices of the dimensions to subset or expand on. Use it only if you are
#' \emph{absolutely sure} about the datasets you want to test later on!
#' @param options. a list of further arguments passed to the wrapped function 
#' (see Functions for details)
#' @details All \code{is*} functions return a function of the class
#' \code{IsFunction} which has only one argument, '\code{x}' (the data object
#' to test). '\code{x}' can be an atomic vector, matrix, or array; all other
#' object types result in error.\cr
#' A second major rule is that all functions returned by an \code{is*}
#' function return a logical object of the same shape as the input object.
#' The third rule is that the function returned by any \code{is*} function is
#' a hard (strict) or soft (optional) constraint, controled by the 'strict.'
#' argument. This affects the way how the results are combined for joint
#' logical tests (see \code{\link{combine}}).
#' \cr
#' Additionally, all \code{is*} functions has 'negate.', 'subset.', and
#' 'expand.' arguments. By setting 'negate.' to TRUE, the function can be
#' conceived of as a \code{not*} function. The same can be achieved by using
#' the standard negation operator (see Examples).\cr
#' The 'subset.' and 'expand.' arguments are only considered if the returned
#' function will be called on matrices or arrays. If 'subset.' is provided, the
#' logical condition is tested only on the subsetted slice of the array, and
#' all other data points in the array become FALSE (or TRUE, if negate. is
#' TRUE). The 'expand.' argument does the same except that it expands the
#' subsetted results to the original array instead of fixing the outer-subset
#' data points to FALSE (or TRUE). See Examples.
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
#' #
#' # example for isPattern #
#' #
#' # note how we pass the 'ignore.case' argument to the internally used
#' # grepl() function (see ?grepl)
#' check_start_a <- isPattern("^a", options. = list(ignore.case = TRUE))
#' check_start_a(c("a", "A", "ba"))
#' \dontshow{
#' # check the results -> note that the class and the attribute has to be
#' # removed by calling 'c' (or 'as.vector') to make the results identical
#' stopifnot(identical(
#'     c(check_start_a(c("a", "A", "ba"))),
#'     c(TRUE, TRUE, FALSE)
#' ))
#' }
isPattern <- function(pattern., strict. = TRUE, negate. = FALSE,
                      subset. = list(), expand. = list(), options. = list()) {
    assertString(pattern., na.ok = FALSE, .var.name = "pattern.")
    assertFlag(strict., .var.name = "strict.")
    assertFlag(negate., .var.name = "negate.")
    which_dims. <- checkSubsetExpand(subset., expand.)
    assertList(options., names = "unique", .var.name = "options.")
    structure(
        function(x) {
            tempfn <- function(dat, options., negate., pattern.) {
                out <- do("grepl", pattern., as.vector(dat),
                          arg_list = options.)
                setattr(out, "dim", dim(dat))
                setattr(out, "dimnames", dimnames(dat))
                if (!negate.) out else !out
            }
            #
            assertCharacter(x, .var.name = "x")
            out <- isMainCompFn(tempfn, x, subset., expand., which_dims.,
                                options., negate.,
                                pattern. = pattern.)
        }, .MUST = strict., class = "IsFunction")
}

#' @describeIn is produces a function to test whether the values in an atomic
#' object are equal to a (vector of) reference value(s) with a given tolerance.
#' The returned function accepts objects of any type for which \code{as.Numeric}
#' does not fail (see the documentation of \code{\link{as.Numeric}}).
#' @param ref. the reference value or a vector (or matrix/array) of reference
#' values.
#' @param tol. tolerance value, the maximum difference which is still tolerated.
#' Can be a vector (or matrix/array) as well.
#' @export
#' @examples
#' #
#' # example for isEqual #
#' #
#' check_zero <- isEqual(0, tol. = 1e-4)
#' check_zero(c("0", "0.000001", "1"))
#' \dontshow{
#' stopifnot(identical(
#'     c(check_zero(c("0", "0.000001", "1"))),
#'     c(TRUE, TRUE, FALSE)
#' ))
#' }
#'
#' # note that just like for other is* functions, the reference and the
#' # tolerance can be vectors (or even matrices or arrays), not only scalars,
#' # and its values are recycled or cropped if necessarily to match the length
#' # of the 'x' object
#' check_values <- isEqual(c(0, 1, 2))
#' check_values(c(0, 1.1, 2, 1e-9))
#'
#' \dontshow{
#' stopifnot(identical(
#'     c(check_values(c(0, 1.1, 2, 1e-9))),
#'     c(TRUE, FALSE, TRUE, TRUE)
#' ))
#' }
#'
#' # it is also possible to base the check only on a subset of x and also
#' # on a reference slice of x which is then expanded back to the original
#' # array;
#' # here we want to test that in the "pre-test" phase of an experiment which
#' # observations had a width of 1 with 0.1 tolerance (expanded also to height)
#' check_width <- isEqual(1, tol. = 0.1,
#'                        subset. = list(time = "pre"),
#'                        expand. = list(measure = "width")
#'                        )
#' ( x <- array(c(1.001, 1.6, 1.2, 1.003, 1.01, 1, 4, 3.5),
#'              dim = c(2, 2, 2),
#'              dimnames = list(observation = c("a", "b"),
#'                              measure = c("width", "height"),
#'                              time = c("pre", "post"))) )
#' ( res <- check_width(x) )
#'
#' # note that all "time = post" data points are FALSE, and also that
#' # all height measures are TRUE where width is TRUE
#' array(paste(x, res, sep = ": "), dim(x), dimnames(x))
#' \dontshow{
#' stopifnot(identical(dim(res), dim(x)))
#' stopifnot(identical(
#'     c(res),
#'     c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
#' ))
#' #
#' subs <- list("pre")
#' attr(subs, "dim_index") <- 3L
#' exp <- list("width")
#' attr(exp, "dim_index") <- 2L
#' check_width <- isEqual(1, tol. = 0.1, subset. = subs, expand. = exp)
#' res2 <- check_width(x)
#' stopifnot(identical(res, res2))
#' }
#'
#' # isEqual, just like all other is* functions, can be negated:
#' check_not_zero <- isEqual(0, negate. = TRUE)
#' check_not_zero(1)
#' \dontshow{
#' stopifnot(check_not_zero(1))
#' }
#'
#' # the same can be achieved by using ! for negation
#' check_not_zero2 <- !isEqual(0)
#' check_not_zero2(1)
#' \dontshow{
#' stopifnot(identical(check_not_zero(1), check_not_zero2(1)))
#' }
isEqual <- function(ref., tol. = .Machine$double.eps^0.5,
                    strict. = TRUE, negate. = FALSE,
                    subset. = list(), expand. = list(), 
                    options. = list()) {
    ref. <- as.Numeric_(ref.)
    tol. <- as.Numeric_(tol.)
    assertFlag(strict., .var.name = "strict.")
    assertFlag(negate., .var.name = "negate.")
    which_dims. <- checkSubsetExpand(subset., expand.)
    assertList(options., names = "unique", .var.name = "options.")
    structure(
        function(x) {
            tempfn <- function(dat, options., negate., ref., tol.) {
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
                # return
                if (!negate.) {
                    abs(dat - ref.) <= tol.
                } else {
                    abs(dat - ref.) > tol.
                }
            }
            #
            x <- as.Numeric_(x)
            isMainCompFn(tempfn, x, subset., expand., which_dims.,
                         options., negate.,
                         ref. = ref., tol. = tol.)
        }, .MUST = strict., class = "IsFunction")
}


#' @describeIn is produces a function to test whether the values in an atomic
#' object are the same as a (vector of) reference value(s). It is a wrapper
#' around \code{\link{==}}. If you want to test numeric values, \code{isEqual}
#' is a much better alternative.
#' @export
#' @examples
#' #
#' # example for isSame #
#' #
#' check_not_A <- isSame("A", negate. = TRUE)
#' check_not_A(c("a", "A"))
#' \dontshow{
#' stopifnot(identical(
#'     c(check_not_A(c("a", "A"))),
#'     c(TRUE, FALSE)
#' ))
#' }
isSame <- function(ref., strict. = TRUE, negate. = FALSE,
                   subset. = list(), expand. = list(), options. = list()) {
    assertFlag(strict., .var.name = "strict.")
    assertFlag(negate., .var.name = "negate.")
    which_dims. <- checkSubsetExpand(subset., expand.)
    assertList(options., names = "unique", .var.name = "options.")
    structure(
        function(x) {
            tempfn <- function(dat, options., negate., ref.) {
                ref. <-
                    if (is.array(ref.)) {
                        isExpandInto(ref., dat, "ref.")
                    } else {
                        repLen(ref., length(dat), "ref.")
                    }
                if (!negate.) dat == ref. else dat != ref.
            }
            # return
            isMainCompFn(tempfn, x, subset., expand., which_dims.,
                         options., negate.,
                         ref. = ref.)
        }, .MUST = strict., class = "IsFunction")
}


#' @describeIn is produces a function to test whether
#' the values in an atomic object are between pre-specified limits.
#' The returned function accepts objects of any type for which \code{as.Numeric}
#' does not fail (see the documentation of \code{\link{as.Numeric}}).
#' @param lwr.,upr. numeric values referring to the lower and upper limit
#' (the defaults are -Inf and Inf, respectively). Both lwr. and upr. can be
#' scalars and vectors as well.
#' @param open. the side of the interval which is open (that is, it does not
#' contain the given endpoint); can be "none" (the default), "both", lwr" or
#' "upr"
#' @export
#' @examples
#' #
#' # example for isBetween() #
#' #
#' # note that lwr. and upr. might be given in reversed order
#' check_0_100 <- isBetween(100, 0)
#' check_0_100(c(-1, 1, 101))
#' \dontshow{
#' stopifnot(identical(
#'     c(check_0_100(c(-1, 1, 101))),
#'     c(FALSE, TRUE, FALSE)
#' ))
#' }
isBetween <- function(lwr. = -Inf, upr. = Inf,
                      open. = c("none", "both", "lwr", "upr"),
                      strict. = TRUE, negate. = FALSE,
                      subset. = list(), expand. = list(), 
                      options. = list()) {
    lwr. <- as.Numeric_(lwr.)
    upr. <- as.Numeric_(upr.)
    open. <- match.arg(open.)
    assertFlag(strict., .var.name = "strict.")
    assertFlag(negate., .var.name = "negate.")
    which_dims. <- checkSubsetExpand(subset., expand.)
    assertList(options., names = "unique", .var.name = "options.")
    #
    lwr <- lwr.; upr <- upr.
    lwr. <- pmin(lwr, upr)
    upr. <- pmax(lwr, upr)
    lwr <- NULL; upr <- NULL
    #
    structure(
        function(x) {
            tempfn <- function(dat, options., negate., lwr., upr., open.) {
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
                if (!negate.) out else !out
            }
            #
            x <- as.Numeric_(x)
            # return
            isMainCompFn(tempfn, x, subset., expand., which_dims.,
                         options., negate.,
                         lwr. = lwr., upr. = upr., open. = open.)
        }, .MUST = strict., class = "IsFunction")
}

#' @describeIn is produces a function to test whether the values in an atomic
#' object are negative. The returned function accepts
#' objects of any type for which \code{as.Numeric}
#' does not fail (see the documentation of \code{\link{as.Numeric}}).
#' @export
isNegative <- function(strict. = TRUE, negate. = FALSE,
                       subset. = list(), expand. = list(), options. = list()) {
    isBetween(-Inf, 0, "upr", strict., negate., subset., expand., options.)
}

#' @describeIn is produces a function to test whether the values in an atomic
#' object are positive. The returned function accepts
#' objects of any type for which \code{as.Numeric}
#' does not fail (see the documentation of \code{\link{as.Numeric}}).
#' @export
isPositive <- function(strict. = TRUE, negate. = FALSE,
                       subset. = list(), expand. = list(), options. = list()) {
    isBetween(0, Inf, "lwr", strict., negate., subset., expand., options.)
}

#' @describeIn is produces a function to test whether the the values in an
#' atomic object are local or global extrema. See \code{\link{findExtrema}}
#' for additional details and further arguments. The returned function accepts
#' objects of any type for which \code{as.Numeric}
#' does not fail (see the documentation of \code{\link{as.Numeric}}).
#' @param what. specifies whether minima ("min"), maxima ("max") or
#' extrema ("both", the default) should be tested
#' @export
#' @examples
#' #
#' # example for isExtremum #
#' #
#' # note how we can pass the 'global' and 'tail' arguments to the internal
#' # findExtrema() function (see ?findExtrema)
#' check_global_extr <- isExtremum(options. = list(global = TRUE, 
#'                                                 tail = "do_not_care"))
#' check_global_extr(c(-1, 1, 0, 100, 50))
#' \dontshow{
#' stopifnot(identical(
#'     c(check_global_extr(c(-1, 1, 0, 100, 50))),
#'     c(TRUE, FALSE, FALSE, TRUE, FALSE)
#' ))
#' }
isExtremum <- function(what. = c("both", "min", "max"), strict. = TRUE,
                       negate. = FALSE,
                       subset. = list(), expand. = list(), options. = list()) {
    what. <- match.arg(what.)
    assertFlag(strict., .var.name = "strict.")
    assertFlag(negate., .var.name = "negate.")
    which_dims. <- checkSubsetExpand(subset., expand.)
    assertList(options., names = "unique", .var.name = "options.")
    if (is.null(options.$constant)) {
        options.$constant <- switch(what.,
                                    both = 3L,
                                    min = 1L,
                                    max = 2L)
    }
    # return
    structure(
        function(x) {
            tempfn <- function(x, options., negate., what.) {
                out <- do("findExtrema", x, arg_list = options.)
                out <- switch(what.,
                              both = out > 0L,
                              min = out == 1L,
                              max = out == 2L)
                if (!negate.) out else !out
            }
            #
            x <- as.Numeric_(x)
            # return
            isMainCompFn(tempfn, x, subset., expand., which_dims.,
                         options., negate.,
                         what. = what.)
        }, .MUST = strict., class = "IsFunction")
}

#' @describeIn is a shorthand for \code{isExtremum("max", ...)}. See
#' \code{\link{findExtrema}} for additional details and further arguments.
#' @export
isLocalMaximum <- function(strict. = TRUE, negate. = FALSE,
                           subset. = list(), expand. = list(), 
                           options. = list()) {
    isExtremum("max", strict., negate., subset., expand., options.)
}

#' @describeIn is a shorthand for \code{isExtremum("min", ...)}. See
#' \code{\link{findExtrema}} for additional details and further arguments.
#' @export
isLocalMinimum <- function(strict. = TRUE, negate. = FALSE,
                           subset. = list(), expand. = list(), 
                           options. = list()) {
    isExtremum("min", strict., negate., subset., expand., options.)
}

#' @describeIn is a shorthand for
#' \code{isExtremum("max", global = TRUE, ...)}. See
#' \code{\link{findExtrema}} for additional details and further arguments.
#' @export
isMaximum <- function(strict. = TRUE, negate. = FALSE,
                      subset. = list(), expand. = list(), options. = list()) {
    assertList(options., names = "unique", .var.name = "options.")
    options.$global <- TRUE
    isExtremum("max", strict., negate., subset., expand., options.)
}

#' @describeIn is a shorthand for
#' \code{isExtremum("min", global = TRUE, ...)}. See
#' \code{\link{findExtrema}} for additional details and further arguments.
#' @export
isMinimum <- function(strict. = TRUE, negate. = FALSE,
                      subset. = list(), expand. = list(), options. = list()) {
    assertList(options., names = "unique", .var.name = "options.")
    options.$global <- TRUE
    isExtremum("min", strict., negate., subset., expand., options.)
}




