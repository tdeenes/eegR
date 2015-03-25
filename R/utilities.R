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