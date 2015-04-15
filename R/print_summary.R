#
# <<< print and summary >>>----
#

#' Extract slots from arrayAnova, arrayTtest, or tanova object
#' 
#' \code{extract} extracts slots from arrayAnova, arrayTtest, or tanova objects,
#' each produced by the corresponding function.
#' @param obj an object as returned by \code{\link{arrayTtest}},
#' \code{\link{arrayAnova}}, or \code{\link{tanova}}
#' @param ... parameters to the correspondig extractor method
#' (e.g. \code{\link{extract.arrayAnova}})
#' @note Only \code{arrayAnova} is implemented at the moment.
#' @export
#' @return an array of the requested component or a named list of arrays if
#' multiple components were requested
#' @seealso See the examples of \code{\link{arrayAnova}}
extract <- function(obj, ...) UseMethod("extract")

#' Default extractor
extract.default <- function(obj, ...) {
    stop("No extractor for this class")
}

#' Extract slots from an \code{arrayAnova} object
#' 
#' Extractor for \code{\link{arrayAnova}} objects
#' @param obj an object of class \code{arrayAnova}
#' @param what a character vector indicating the slots to extract
#' @param time_window numeric vector of length two indicating the time window
#' which should be extracted. If NULL (default), all time points are returned.
#' @param term character vector indicating the name of the model terms which
#' should be extracted. If NULL (default), all model terms are returned.
#' @param chan character or integer indices of channels to extract. If NULL 
#' (default), all channels are returned.
#' @param drop logical value; if TRUE, singleton dimensions are dropped from
#' the arrays after subsetting on time, term, and/or chan (default: FALSE)
#' @export
extract.arrayAnova <- function(obj, 
                               what = c("stat", "p", "es", 
                                        "stat_corr", "p_corr",
                                        "means", "df"),
                               time_window = NULL, term = NULL, chan = NULL, 
                               drop = FALSE) {
    what <- tolower(what)
    what <- match.arg(what, several.ok = TRUE)
    dimn <- dimnames(obj$effect_F_obs)
    if (is.null(term)) term <- dimn$modelterm
    if (is.null(chan)) chan <- dimn$chan
    time <- 
        if (is.null(time_window)) {
            dimnames(obj$effect_F_obs)$time
        } else if (length(time_window) != 2L) {
            stop("time_window must contain two values")
        } else if (anyNA(time_window)) {
            stop("time_window has missing value(s)")
        } else {
            time_window <- sort(time_window)
            timep <- as.numeric(dimn$time)
            time <- as.character(timep[timep > time_window[1] & 
                                           timep < time_window[2]])
        }
    sub <- list(chan = chan, time = time, modelterm = term)
    out <- setNames(vector("list", length(what)), what)
    if ("stat" %in% what) {
        out[["stat"]] <- subsetArray(obj$effect_F_obs, sub, drop = drop)
        setattr(out[["stat"]], "Df.term", NULL)
        setattr(out[["stat"]], "Df.resid", NULL)
        setattr(out[["stat"]], "pvalues", NULL)
        setattr(out[["stat"]], "ges", NULL)
        setattr(out[["stat"]], "factor_means", NULL)
    }
    if ("p" %in% what) {
        out[["p"]] <- subsetArray(attr(obj$effect_F_obs, "pvalues"), sub, 
                                  drop = drop)
    }
    if ("es" %in% what) {
        out[["es"]] <- subsetArray(attr(obj$effect_F_obs, "ges"), sub, 
                                   drop = drop)
    }
    if ("means" %in% what) {
        m <- attr(obj$effect_F_obs, "factor_means")
        m <- lapply(m, subsetArray, subsets = sub[c("chan", "time")], 
                    drop = drop)
        out[["means"]] <- 
            if (drop && length(term) == 1L) m[[term]] else m[term]
    }
    if ("stat_corr" %in% what) {
        out[["stat_corr"]] <- subsetArray(obj$effect_tfce_obs, sub, drop = drop)
    }
    if ("p_corr" %in% what) {
        out[["p_corr"]] <- subsetArray(obj$perm_pvalues, sub, drop = drop)
        setattr(out[["p_corr"]], "type", attr(obj$perm_pvalues, "type"))
    }
    if ("df" %in% what) {
        out[["df"]] <- cbind(attr(obj$effect_F_obs, "Df.term"),
                             attr(obj$effect_F_obs, "Df.resid"))
        setattr(out[["df"]], "dimnames", list(modelterm = dimn$modelterm, 
                                              c("Df.term", "Df.resid")))
        out[["df"]] <- out[["df"]][term, , drop = drop]
    }
    what_null <- vapply(out[what], is.null, FALSE)
    if (any(what_null)) {
        warning(paste0(paste(what[what_null], collapse = ", "), 
                       " values are not available"))
    }
    if (length(out) == 1) out <- out[[1]]
    # return
    out
}

#' Compute summary statistics for \code{arrayAnova} object
#' 
#' \code{summary.arrayAnova} computes the minimum, maximum and median of
#' F-values, generalized eta squared effect sizes and P-values for all 
#' channel X time points which meet a certain criterion.
#' @param obj an object of class \code{\link{arrayAnova}}
#' @param basis character value indicating the criterion; if \code{p_corr} 
#' (the default), channel X time points below a specified corrected alpha-level 
#' are selected (choose \code{p} for uncorrected p-values); if \code{es} or 
#' \code{stat} (or \code{stat_corr}), the effect size or (corrected) test 
#' statistic should exceed a given limit.
#' @param crit the level of criterion
#' @param ... arguments passed to \code{\link{extract}}, e.g. 
#' \code{time_window = c(100, 200)}
#' @export
#' @return a list of summaries
summary.arrayAnova <- function(obj, 
                               basis = c("p_corr", "p", "es", "stat", "stat_corr"),
                               crit = 0.05,
                               ...) {
    #
    applyFn <- function(x, size = FALSE) {
        x[!ind] <- NA
        out <- t(colQuantiles(x, probs = c(0, 1, 0.5), na.rm = TRUE))
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
        aperm(x, c(dims, setdiff(ndimn, dims)))
    }
    #
    basis <- tolower(basis)
    basis <- match.arg(basis)
    what <- 
        if (any(grepl("corr", basis))) {
            c("stat", "stat_corr", "p_corr", "es", "means", "df")
        } else {
            c("stat", "p", "es", "means", "df")
        }
    out <- extract.arrayAnova(obj, what = what, ...)
    names(out)[names(out) == "p_corr"] <- "p"
    if ("p_corr" %in% what) {
        if (is.null(out[["p"]]) && basis == "p_corr") 
            stop("TFCE corrected or permuted p-values are not available. Check 'basis'")
    } else {
        if (is.null(out[["p"]])) 
            stop("Traditional P-values are not available")
    }
    ind <- 
        if (basis %in% c("p_corr", "p")) {
            out[["p"]] < crit
        } else {
            out[[basis]] > crit
        }
    aslots <- setdiff(names(out), c("p", "means", "df"))
    out[aslots] <- lapply(out[aslots], function(slot) 
        Aperm(fnDims(
            slot, target_dim = c("chan", "time"), applyFn,
            newdims = list(descriptives = c("min", "max", "median")),
            vectorized = TRUE)
        ))
    p_type <- 
        if ("p_corr" %in% what) {
            attr(out[["p"]], "type")
        } else {
            "trad"
        }
    out[["p"]] <- Aperm(
        fnDims(out[["p"]], target_dim = c("chan", "time"),
               applyFn, arg_list = list(size = TRUE),
               newdims = list(descriptives = c("min", "max", "median",
                                               "size", "ratio")),
               vectorized = TRUE))
    setattr(out[["p"]], "type", p_type)
    out[["call"]] <- obj[["call"]]
    out[["basis"]] <- data.frame(basis = basis, crit = crit, 
                                 stringsAsFactors = FALSE)
    setattr(out, "class", "summary.arrayAnova")
    # return
    out
}

#' Print summaries of \code{\link{arrayTtest}}, \code{\link{arrayAnova}}, or
#' \code{\link{tanova}} objects
#' 
#' \code{print.summary} prints summaries of \code{\link{arrayTtest}}, 
#' \code{\link{arrayAnova}}, or \code{\link{tanova}} objects
#' @param obj an object of class \code{summary.arrayTtest}, 
#' \code{summary.arrayAnova}, or \code{summary.tanova}
#' @param ... arguments passed to the specific methods
#' @export
#' @return The function invisibly returns the input object.
print.summary <- function(obj, ...) UseMethod("print.summary")

print.summary.default <- function(obj, ...) {
    print.default(obj, ...)
}

#' Print summaries of \code{\link{arrayAnova}} objects
#' 
#' \code{print.summary.arrayAnova} prints summaries of  
#' \code{\link{summary.arrayAnova}} objects
#' @param obj an object of class \code{summary.arrayAnova}
#' @param ... not used yet
#' @export
#' @return The function invisibly returns the input object.
print.summary.arrayAnova <- function(obj, ...) {
    # function to remove "descriptives" header
    rD <- function(x) {
        names(dimnames(x))[2] <- ""
        x
    }
    # main
    cat("\n----\nCall\n----\n", 
        paste(deparse(obj$call), sep = "\n", collapse = "\n"))
    cat("\n\n----\nDesriptive statistics\n----\n")
    rel <- if (grepl("p", obj$basis$basis)) "<" else ">"
    basis <- Replace(obj$basis$basis, 
                     c("stat", "stat_corr", "es", "p", "p_corr"),
                     c("F-statistic", 
                       "TFCE-corrected F-statistic",
                       "Effect size",
                       "Uncorrected P-values",
                       "Corrected P-values"))
    cat("\n( Basis:", paste(basis, rel, obj$basis$crit, ")\n", sep = " "))
    cat("\n<<< F-values >>>\n")
    print(formatC(
        rD(obj$stat), digits = 2, format = "f"), quote = FALSE)
    cat("\nDegrees of freedom:\n")
    print(obj$df)
    p_type <- Replace(attr(obj$p, "type"), 
                      c("trad", "perm", "tfce"),
                      c("\n<<< (Traditional)",
                        "\n<<< (Permuted)",
                        "\n<<< (TFCE-corrected)"))
    cat(p_type, "P-values >>>\n")
    print(
        formatC(rD(subsetArray(obj$p, 
                               list(descriptives = c("min", "max", "median")))), 
                digits = 3, format = "f"),
        quote = FALSE
    )
    cat("\nNumber and ratio of significant time points:\n")
    p <- subsetArray(obj$p, list(descriptives = c("size", "ratio")))
    pchar <- arrayIP("", dim(p), dimnames(p))
    subsetArray(pchar, list(descriptives = "size")) <- 
        formatC(subsetArray(obj$p, list(descriptives = "size")), 
                digits = 0, format = "f")
    subsetArray(pchar, list(descriptives = "ratio")) <- 
        formatC(subsetArray(obj$p, list(descriptives = "ratio")), 
                digits = 3, format = "f")
    print(rD(pchar), quote = FALSE)
    # return
    invisible(obj)
}