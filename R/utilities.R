#
# <<< simple utility functions >>> --------
#

#' Address in RAM of a variable
#' 
#' \code{address} returns the pointer address of its argument.
#' @param x anything
#' @details Sometimes useful in determining whether a value has been copied or 
#' not, programatically.
#' @return A character vector length 1.
#' @note This is just an imported and re-exported function (without any 
#' modification) from \bold{data.table}.
#' @export
address <- data.table::address

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
    assertVector(x, strict = TRUE, .var.name = "x")
    assertVector(from, strict = TRUE, any.missing = FALSE, .var.name = "x")
    assertVector(to, strict = TRUE, any.missing = FALSE, .var.name = "x")
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

#' Execute a function call not unlike \code{do.call}.
#'
#' This function serves as an efficient replacement for 
#' \code{\link[base]{do.call}}; arguments can be passed via \code{...} to avoid 
#' any copying of potentially large objects.
#' @param what either a function or a non-empty character string naming the 
#' function to be called
#' @param ... arguments to \code{what}, usually specified as \code{key = value}
#' pairs
#' @param arg_list a list of arguments to the function call. The names attribute
#' of arg_list gives the argument names.
#' @return The result of the (evaluated) function call.
#' @note This function was inspired by \code{do.call2} in package 
#' \pkg{BBmisc}
#' @export
#' @examples
#' # create a largish data.frame
#' x <- data.frame(a = seq_len(1e7), b = seq_len(1e7)/10)
#' 
#' # check that do() and do.call() returns the same;
#' # suppose we want to call head() to display the first 10 rows
#' stopifnot(identical(head(x, n = 10L), 
#'                     do("head", arg_list = list(x, n = 10L))))
#' stopifnot(identical(do.call("head", list(x, n = 10L)), 
#'                     do("head", x, n = 10L)))
#' 
#' #
#' # speed comparisons
#' # 
#' 
#' # a little helper function (do not use for serious measurements)
#' test <- function(expr) {
#'     gc(reset = TRUE)
#'     cat("CPU time:\n")
#'     print(system.time(expr, gcFirst = FALSE))
#'     cat("\nRAM usage:\n")
#'     print(gc())
#' }
#' 
#' # a direct call for comparison
#' test(head(x, n = 10L))
#' 
#' # do.call() can be substantially slower because it might make a copy
#' test(do.call("head", list(x, n = 10L)))
#' 
#' # do() is almost as fast as a direct call in this case
#' test(do("head", x, n = 10L))
#' 
#' # try to avoid using the 'arg_list' argument for passing large objects
#' test(do("head", n = 10L, arg_list = list(x = x)))
#' 
do <- function(what, ..., arg_list = list()) {
    mc <- match.call(expand.dots = FALSE)[["..."]]
    to_call <- 
        if (is.function(what)) {
            c(list(what), mc, arg_list)
        } else if (is.character(what)) {
            c(list(as.name(what[[1L]])), mc, arg_list)
        } else {
            stop("'what' must be either a function or a non-empty character string naming the function to be called")
        }
    expr <- as.call(to_call)
    eval(expr, parent.frame())
}

#' Pass \code{arg = .(key1 = value1, key2 = value2)} function arguments
#' 
#' \code{argumentDeparser} passes \code{arg = .(key1 = value1, key2 = value2)} 
#' function arguments to the appropriate function
#' @param arg an unevaluated call, see Details
#' @param replace_dot a character string which defines the function to which
#' \code{arg} should be passed to. If not provided, the name of \code{arg} is
#' guessed using \code{deparse} and \code{substitute}, and "Params" is attached
#' to call the corresponding parameter setter function.
#' @param transform_logical logical; if TRUE (default), a single logical 
#' \code{arg} argument is treated in a special way (see Details)
#' @param null_params a list of parameters which is passed to the parameter
#' setter function (see 'replace_dot') if \code{argumentDeparser} would return
#' NULL
#' @details This function is not intended for direct use. It allows  
#' \code{method = .(key = value)} argument definition in high-level functions 
#' by substituting \code{.} to the appriopriate \code{methodParams}
#' function. If \code{transform_logical} is TRUE, the call 
#' \code{method = TRUE} is transformed to \code{method = methodParams()}, 
#' thereby it returns the default parameter setting. Similarily, 
#' \code{method = FALSE} is treated as \code{method = NULL}.
#' @export
#' @keywords internal
#' @examples
#' mymethodParams <- function(x = 3, y = 4) {
#'     list(x = x, y = y)
#' }
#' tempfn <- function(mymethod = NULL, ...) {
#'     argumentDeparser(substitute(mymethod), "mymethodParams", ...)
#' }
#' stopifnot(is.null(tempfn()))
#' stopifnot(identical(tempfn(null_params = list(x = 1)),
#'                     list(x = 1, y = 4)))
#' stopifnot(identical(FALSE, 
#'                     tempfn(mymethod = FALSE, transform_logical = FALSE)))
#' stopifnot(identical(tempfn(mymethod = TRUE),
#'                     tempfn(mymethod = mymethodParams())))
#' new_y = 1:5
#' stopifnot(identical(tempfn(mymethod = .(y = new_y)),
#'                     tempfn(mymethod = mymethodParams(y = new_y))))
argumentDeparser <- function(arg, replace_dot, 
                             transform_logical = TRUE,
                             null_params = NULL) {
    if (missing(replace_dot)) {
        argname <- deparse(substitute(replace_dot))
        if (is.null(argname)) 
            stop("Provide the 'replace_dot' argument, its name could not be figured out automagically")
        replace_dot <- paste0(argname, "Params")
    }
    if (is.symbol(arg)) arg <- eval(arg, parent.frame())
    out <- 
        if (transform_logical && identical(arg, TRUE)) {
            do.call(replace_dot, list())
        } else if (transform_logical && identical(arg, FALSE)) {
            NULL
        } else {
            if (identical(arg[[1L]], as.name("."))) {
                arg[[1L]] <- as.name(replace_dot)
            }
            eval(arg, parent.frame())
        }
    if (is.null(out) && !is.null(null_params)) {
        assertList(null_params, .var.name = "null_params")
        out <- do.call(replace_dot, null_params)
    }
    # return
    out
}


#' Set seed
#' 
#' \code{setSeed} is called for its side effect, namely it specifies a random 
#' seed by calling \code{\link{set.seed}}.
#' @param seed either NULL, in which case \code{setSeed} does nothing, or an 
#' integer value, in which case it is interpreted as the 'seed' argument in
#' \code{\link{set.seed}}, or a list of arguments passed to 
#' \code{\link{set.seed}}
setSeed <- function(seed) {
    if (!is.null(seed)) {
        if (is.list(seed)) {
            do.call("set.seed", seed)
        } else {
            set.seed(seed = seed)
        }
    }
}

#' Argument verification if there might be an "all" option
#' 
#' \code{matchArg} is a wrapper around \code{\link[base]{match.arg}} to check
#' arguments which have an "all" option in their choices. Note that the defaults
#' are not the same as in \code{\link[base]{match.arg}}.
#' @param arg a character vector (of length one unless several_ok is TRUE) or 
#' NULL
#' @param choices a character vector of candidate values (if NULL [default], it 
#' is deduced from the formals of the calling function)
#' @param several_ok logical specifying whether arg should be allowed to have 
#' more than one element (default: TRUE)
#' @return If 'arg' has an element 'all' (or 'ALL', 'All', etc.), all other 
#' choices are returned. Otherwise see \code{\link[base]{match.arg}}.
#' @keywords internal 
matchArg <- function(arg, choices = NULL, several_ok = TRUE) {
    arg <- tolower(arg)
    if (is.null(choices)) {
        formal_args <- formals(sys.function(sys.parent()))
        choices <- eval(formal_args[[deparse(match.call()$arg)]])
    }
    choices <- setdiff(tolower(choices), "all")
    # return choices if "all" and matching choice(s) otherwise
    if (length(arg) > 0L && any(arg == "all")) {
        if (!several_ok) stop("'arg' may not be 'all' if 'several_ok' is FALSE")
        choices
    } else {
        match.arg(arg, choices, several.ok = several_ok)
    }
}

#' Attach dimension and type attributes to test statistics
#' 
#' Used internally in the workhorse functions for arrayAnova, arrayTtest, etc.
#' to set the dimension, dimension names and verbose type of the given test
#' statistic
#' @param x the vector or matrix of the test statistic
#' @param obj if not NULL (default), a list object returned by 
#' \code{\link{preTtest}} or \code{\link{preAnova}}
#' @param type if not NULL (default), a character string with a verbose 
#' description of the test statistic (e.g. "Traditional F-value", "Welch t", 
#' "Generalized Eta Squared", etc.)
#' @param dimorder if not NULL (default), a numeric or character vector 
#' indicating the order (or subset) of the dimensions
#' @keywords internal
setattributes <- function(x, obj = NULL, type = NULL, dimorder = NULL) {
    if (!is.null(obj)) {
        if (is.null(dimorder)) dimorder <- seq_along(obj$teststat_dimid)
        array_(x, obj$full_dims[obj$teststat_dimid][dimorder], 
               arg_check = FALSE)    
    }
    if (!is.null(type))
        setattr(x, "type", type)
    # return
    invisible(x)
}
