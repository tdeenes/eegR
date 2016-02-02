#
# <<< print and summary >>>----
#

#' Extract slots from arrayAnova, arrayTtest, or tanova object
#'
#' \code{extract} is a generic function for extracting slots from
#' \code{\link{arrayTtest}}, \code{\link{arrayAnova}}, or \code{\link{tanova}}
#' objects, each produced by the corresponding function.
#' @param object an object of class \code{arrayTtest}, \code{arrayAnova}, or
#' \code{tanova}
#' @param what a character vector indicating which slot(s) to extract. If
#' 'what' is 'all' (default), all relevant slots is extracted (see the Usage
#' section for the relevant slot names in each particular method). Otherwise,
#' the slots to extract are found by partial matching to the following choices:
#' \itemize{
#' \item{stat: }{Test statistic (e.g. F values in \code{arrayAnova})}
#' \item{stat_corr: }{Corrected test statistic (e.g. TFCE values in
#' \code{arrayAnova})}
#' \item{p: }{Uncorrected P-values}
#' \item{p_corr: }{Corrected P-values (e.g. TFCE-corrected P-values in
#' \code{arrayAnova}, or minimum duration correction in \code{tanova})}
#' \item{es: }{Effect size statistic (e.g. Generalized Eta Squared in
#' \code{arrayAnova})}
#' \item{means: }{Marginal means for each model term} (TODO! Temporarily removed...)
#' }
#' @param time_window numeric vector of length two indicating the time window
#' which should be extracted. If NULL (default), all time points are returned.
#' @param term character vector indicating the name of the model terms which
#' should be extracted. If NULL (default), all model terms are returned.
#' @param chan character or integer indices of channels to extract. If NULL
#' (default), all channels are returned.
#' @param drop logical value; if TRUE, singleton dimensions are dropped from
#' the arrays after subsetting on time, term, and/or chan (default: FALSE)
#' @param ... not used yet
#' @export
#' @return The function returns an array of the requested component or a named
#' list of arrays if multiple components were requested.
#' @seealso See the examples of \code{\link{arrayAnova}}
extract <- function(...) UseMethod("extract")

#' Default extractor
#'
#' Default extractor which is behind \code{\link{extract.arrayTtest}},
#' \code{\link{arrayAnova}}, and \code{\link{tanova}}, and does nothing
#' (at the moment) for objects of other classes.
#' @export
#' @keywords internal
extract.default <- function(object, what,
                            time_window = NULL, term = NULL, chan = NULL,
                            drop = FALSE, ...) {
    # helper function
    extractor <- function(x, params) {
        out <- copy(x)
        if (is.null(out)) return(out)
        if (identical(unname(dim(out)), unname(params$dims)) )
            setattr(out, "dimnames", params$dimn)
        subs <- params$subs[intersect(names(params$subs),
                                      names(dimnames(out)))]
        if (length(subs) > 0L) {
            out <- subsetArray(out, subset. = subs, drop. = params$drop,
                               keep_attributes. = FALSE)
            setattr(out, "label", attr(x, "label"))
        }
        # return
        out
    }
    #
    if (!inherits(object, c("arrayTtest", "arrayAnova", "tanova")))
        stop("No extractor for this class")
    exparams <- list(
        dims = object$dim, dimn = object$dimnames,
        drop = drop
    )
    if (is.null(term)) term <- exparams$dimn$modelterm
    if (is.null(chan)) chan <- exparams$dimn$chan
    time <-
        if (is.null(time_window)) {
            dimnames(object$stat)$time
        } else if (length(time_window) != 2L) {
            stop("time_window must contain two values")
        } else if (anyNA(time_window)) {
            stop("time_window has missing value(s)")
        } else {
            time_window <- sort(time_window)
            timep <- as.numeric(exparams$dimn$time)
            time <- as.character(timep[timep >= time_window[1] &
                                           timep <= time_window[2]])
        }
    subs <- list(chan = chan, time = time, modelterm = term)
    exparams$subs <- subs[!vapply(subs, is.null, FALSE)]
    out <- setNames(vector("list", length(what)), what)
    if ("stat" %in% what) {
        out[["stat"]] <- extractor(object[["stat"]], exparams)
    }
    if ("p" %in% what) {
        out[["p"]] <- extractor(attr(object[["stat"]], "p_value"), exparams)
    }
    if ("es" %in% what) {
        out[["es"]] <- extractor(attr(object[["stat"]], "effect_size"), exparams)
    }
    if ("means" %in% what) {
        m <- attr(object$stat, "marginal_means")
        if (!is.null(m)) {
            m <- lapply(m, subsetArray, subset. = sub[c("chan", "time")],
                        drop. = drop)
            out[["means"]] <-
                if (drop && length(term) == 1L) m[[term]] else m[term]
        } else {
            out[["means"]] <- NULL
        }
    }
    if ("stat_corr" %in% what) {
        out[["stat_corr"]] <- extractor(object[["stat_corr"]], exparams)
    }
    if ("p_corr" %in% what) {
        out[["p_corr"]] <- extractor(object[["p_corr"]], exparams)
    }
    if ("df" %in% what) {
        out[["df"]] <- extractor(attr(object[["stat"]], "Df"), exparams)
    }
    what_null <- vapply(out[what], is.null, FALSE)
    if (any(what_null)) {
        warning(paste0(paste(what[what_null], collapse = ", "),
                       " values are not available"))
    }
    if (length(what) == 1L) out <- out[[1L]]
    # return
    out
}

#' @describeIn extract
extract.arrayTtest <- function(object,
                               what = c("all", "stat", "p", "es",
                                        "stat_corr", "p_corr",
                                        "means", "df"),
                               time_window = NULL, term = NULL, chan = NULL,
                               drop = FALSE,
                               ...) {
    what <- matchArg(what)
    NextMethod(what = what, time_window = time_window, term = term,
               chan = chan, drop = drop, ...)
}

#' @describeIn extract
extract.arrayAnova <- function(object,
                               what = c("all", "stat", "p", "es",
                                        "stat_corr", "p_corr",
                                        "means", "df"),
                               time_window = NULL, term = NULL, chan = NULL,
                               drop = FALSE,
                               ...) {
    what <- matchArg(what)
    NextMethod(what = what, time_window = time_window, term = term,
               chan = chan, drop = drop, ...)
}

#' @describeIn extract
extract.tanova <- function(object,
                           what = c("all", "stat", "p", 
                                    "stat_corr", "p_corr",
                                    "means"),
                           time_window = NULL, term = NULL,
                           drop = FALSE,
                           ...) {
    what <- matchArg(what)
    NextMethod(what = what, time_window = time_window, term = term,
               chan = NULL, drop = drop, ...)
}


#' Compute summary statistics for \code{arrayTtest}, \code{arrayAnova}, or
#' \code{tanova} objects
#'
#' \code{summary.arrayTtest}, \code{summary.arrayAnova}, and
#' \code{summary.tanova} compute the minimum, maximum and median of
#' criterion-referenced test statistics, effect sizes and P-values for
#' particular subviews (e.g. channel X time subviews).
#' @name summary
#' @param object an object of class \code{\link{arrayTtest}},
#' \code{\link{arrayAnova}}, or \code{\link{tanova}}
#' @param basis character value indicating the criterion node (default:
#' "p_corr"), see Details
#' @param basis_dim integer or character vector indicating the dimensions
#' which form relevant subparts for which summary statistics should be
#' computed (default: c("chan", "time"))
#' @param crit numeric scalar, the level of criterion (see Details)
#' @param ... arguments passed to \code{\link{extract}}, e.g.
#' \code{time_window = c(100, 200)}
#' @details If 'basis' is "p_corr" (the default) or "p", denoting corrected and
#' uncorrected P-values, respectively, argument 'crit' refers to the desired
#' alpha-level (set to 0.05 by default). In this case the basis values should be
#' \emph{below} or equal to the criterion. If basis is "stat" (test statistic),
#' "stat_corr" (corrected test statistic) or "es" (effect size), absolute basis
#' values \emph{above} or equal to the 'crit' criterion are selected.\cr
#' The summary statistics are computed for those data points which correspond
#' to the selected positions in the 'basis_dim' subviews. For example,
#' supposing the object is the result of an uncorrected ANOVA calculation (see
#' \code{\link{arrayAnova}}), and separate analyses were conducted for each
#' channel X time data points but no other dimension, \code{summary} returns
#' the list of \code{min}, \code{max}, and \code{median} of F-values and effect
#' sizes (generalized eta squares) taking into consideration only those
#' channel X time values for each model term, which were significant at the
#' level of 0.05.\cr
#' The returned list also contains the descriptive statistics of the selected
#' P-values together with the number and ratio of significant values. Note that
#' if no data points meet the given criterion in a particular subview, no
#' summaries can be computed, resulting in \code{NA} values in the
#' corresponding cells.
#' @export
#' @return The methods return a list of summaries which can be printed with
#' the corresponding \code{print} method.
summary.arrayTtest <- function(object,
                               basis = c("p_corr", "p",
                                         "es", "stat", "stat_corr"),
                               basis_dim = c("chan", "time"),
                               crit = NULL,
                               ...) {
    basis <- tolower(basis)
    basis <- match.arg(basis)
    if (is.null(crit) && grepl("p", basis)) crit <- 0.05
    out <- summary.eegr(object, basis, basis_dim, crit, ...)
    setattr(out, "class", c("summary.arrayTtest", "summary.eegr"))
    # return
    out
}

#' @rdname summary
#' @export
summary.arrayAnova <- function(object,
                               basis = c("p_corr", "p",
                                         "es", "stat", "stat_corr"),
                               basis_dim = c("chan", "time"),
                               crit = NULL,
                               ...) {
    basis <- tolower(basis)
    basis <- match.arg(basis)
    if (is.null(crit) && grepl("p", basis)) crit <- 0.05
    out <- summary.eegr(object, basis, basis_dim, crit, ...)
    setattr(out, "class", c("summary.arrayAnova", "summary.eegr"))
    # return
    out
}

#' @rdname summary
#' @export
summary.tanova <- function(object,
                           basis = c("p_corr", "p", "stat", "stat_corr"),
                           basis_dim = "time",
                           crit = NULL,
                           ...) {
    basis <- tolower(basis)
    basis <- match.arg(basis)
    if (is.null(crit) && grepl("p", basis)) crit <- 0.05
    out <- summary.eegr(object, basis, basis_dim, crit, ...)
    setattr(out, "class", c("summary.tanova", "summary.eegr"))
    # return
    out
}

#' Workhorse function for eegR-summary methods
#'
#' \code{summary.eegr} is the main function behind eegR-summary methods
#' (\code{\link{summary.arrayTtest}}, \code{\link{summary.arrayAnova}},
#' \code{\link{summary.arrayAnova}}).
#' @param object an object of class \code{\link{arrayTtest}},
#' \code{\link{arrayAnova}}, or \code{\link{tanova}}
#' @param basis character value indicating the criterion node
#' @param basis_dim integer or character vector indicating the dimensions
#' which form relevant subparts for which summary statistics should be
#' computed
#' @param crit numeric scalar, the level of criterion
#' @param ... arguments passed to \code{\link{extract}}, e.g.
#' \code{time_window = c(100, 200)}
#' @keywords internal
summary.eegr <- function(object, basis, basis_dim, crit, ...) {
    # helper functions
    applyFn <- function(x, size = FALSE) {
        out <- t(colQuantiles(x, probs = c(0, 1, 0.5),
                              na.rm = TRUE, drop = FALSE))
        if (size) {
            size <- colSums(!is.na(x))
            ratio <- size/nrow(x)
            out <- rbind(out, size, ratio)
        }
        out
    }
    Aperm <- function(x) {
        ndimn <- names(dimnames(x))
        dims <-
            if ("modelterm" %in% ndimn) {
                c("modelterm", "descriptives")
            } else {
                "descriptives"
            }
        # return
        if (length(dim(x)) > 1L) {
            apermArray(x, first = dims, keep_attributes. = TRUE)
        } else {
            x
        }
    }
    #
    # basis must be extracted
    assertString(basis, .var.name = "basis")
    args <- list(...)
    args$what <-
        if (is.null(args$what)) {
            "all"
        } else {
            unique(c(tolower(args$what), basis))
        }
    #
    # extract nodes
    suppressWarnings(out <- do("extract", object, arg_list = args))
    if (!is.list(out)) out <- setNames(list(out), args$what)
    #
    # basis must be present
    if (!basis %in% names(out))
        stop(sprintf("The node corresponding to basis ('%s') is not available",
                     basis))
    #
    # use absolute values of the test statistics
    abs_ind <- grep("stat", names(out))
    if (length(abs_ind) > 0L)
        out[abs_ind] <- lapply(out[abs_ind], abs)
    #
    # find values above (or below) the criteria
    assertNumber(crit, .var.name = "crit")
    ind <-
        if (basis %in% c("p_corr", "p")) {
            out[[basis]] <= crit
        } else {
            out[[basis]] >= crit
        }
    #
    # compute summary statistics
    excl <- "means"
    if (!is.null(dimnames(attr(object$stat, "Df"))))
        excl  <- c(excl, "df")
    for (i in setdiff(names(out), excl)) {
        if (i != basis) {
            arg_size <- FALSE
            descr <- c("min", "max", "median")
        } else {
            arg_size <- TRUE
            descr <- c("min", "max", "median", "size", "ratio")
        }
        temp <- out[[i]]
        temp[!ind] <- NA
        temp <- fnDims(
            temp, target_dim = basis_dim, applyFn,
            arg_list = list(size = arg_size),
            newdims = list(descriptives = descr),
            vectorized = TRUE)
        setattr(temp, "label", attr(out[[i]], "label"))
        out[[i]] <- Aperm(temp)
    }
    #
    # some decoration
    out[["call"]] <- object[["call"]]
    out[["basis"]] <- list(basis = basis, crit = crit, dim = basis_dim)
    #
    # return
    out
}

#' Print summaries of \code{\link{arrayTtest}}, \code{\link{arrayAnova}}, or
#' \code{\link{tanova}} objects
#'
#' \code{print.summary.eegr} prints summaries of \code{\link{arrayTtest}},
#' \code{\link{arrayAnova}}, or \code{\link{tanova}} objects
#' @name print
#' @param x an object of class \code{summary.arrayTtest},
#' \code{summary.arrayAnova}, or \code{summary.tanova}
#' @param ... further arguments passed to other methods (not used yet)
#' @param digits the minimum number of digits to be printed (default: 3)
#' @param quote logical value indicating whether or not strings should be
#' printed with surrounding quotes (default: FALSE [no quotes printed])
#' @param fill_na a character string which is used to indicate NA values in
#' printed output, that is descriptive statistics which could not be computed
#' because no data points met the given criterion (default: "-")
#' @export
#' @return The function invisibly returns the input object.
print.summary.eegr <- function(x, digits = 3L, quote = FALSE, fill_na = "-",
                               ...) {
    # function to remove "descriptives" header
    rD <- function(x) {
        ind <- match("descriptives", names(dimnames(x)))
        names(dimnames(x))[ind] <- ""
        attr(x, "label") <- NULL
        x
    }
    # check arguments
    assertCount(digits, .var.name = "digits")
    assertFlag(quote, .var.name = "quote")
    # main
    cat("\n----\nCall\n----\n",
        paste(deparse(x$call), sep = "\n", collapse = "\n"))
    cat("\n\n----\nDesriptive statistics\n----\n")
    rel <- if (grepl("p", x$basis$basis)) "<" else ">"
    basis <- attr(x[[x$basis$basis]], "label")
    cat("\n( Basis:", paste(basis, rel, x$basis$crit, ")\n", sep = " "))
    for (i in seq_along(x[setdiff(names(x), c("call", "basis"))])) {
        cat("\n<<<", attr(x[[i]], "label"),  ">>>\n")
        print(rD(x[[i]]), digits = digits, quote = quote, na.print = fill_na,
              ...)
    }
    # return
    invisible(x)
}
