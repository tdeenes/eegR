
###############################################################################
# Functions borrowed from the NMF package to get and set doParallel backends.
# This was needed because importing from the NMF package has the unfortunate 
# consequence of polluting the search path by a whole bunch of other packages.
#
# Original Author: Renaud Gaujoux 
# Minor modifications: Denes Toth
###############################################################################

#' \code{getDoBackend} returns the internal data of the currently registered 
#' foreach \%dopar\% backend.
#' @author Renaud Gaujoux, Denes Toth
#' @export
#' @keywords internal
getDoBackend <- function() {
    fe_ns <- asNamespace('foreach')
    fe <- get('.foreachGlobals', fe_ns)
    if (!exists("fun", where = fe, inherits = FALSE)) return(NULL)
    #
    getDoPar <- get('getDoPar', fe_ns)
    # this returns the registered %dopar% function + associated data
    # -> add info function from foreach internal environment
    info <- 
        if (exists("info", where = fe, inherits = FALSE)) {
            get('info', fe, inherits=FALSE) 
        } else {
            function(data, item) NULL
        }
    cleanup <-
        if (exists("cleanup", where = fe, inherits = FALSE)) {
            get('cleanup', fe, inherits=FALSE)
        }
    # return
    c(getDoPar(), info, cleanup)
}

#' \code{setDoBackend} is identical to \code{\link[foreach]{setDoPar}}, but 
#' returns the internal of the previously registered backend.
#' 
#' @param data internal data of a foreach \%dopar\% backend.
#' @author Renaud Gaujoux, Denes Toth
#' @export
#' @keywords internal
setDoBackend <- function(data) {
    # get old backend data
    ob <- getDoBackend()
    #
    if (!is.null(data)) {
        bdata <- data
        if (is(data, "foreach_backend"))
            data <- data[!names(data) %in% c("name", "cleanup")]
        do.call("setDoPar", data)
        setBackendCleanup(bdata)
    } else {
        do.call('setDoPar', list(NULL))
        fe <- get(".foreachGlobals", asNamespace("foreach"))
        if (exists("fun", envir = fe, inherits = FALSE))
            remove("fun", envir = fe)
        setBackendCleanup(NULL)
    }
    # return old backend
    invisible(ob)
}

# setup cleanup procedure for the current backend
setBackendCleanup <- function(object, verbose = FALSE){
    fe <- get(".foreachGlobals", asNamespace("foreach"))
    name <- getDoParName()
    if (!is.null(fun <- object$cleanup)) {
        if (verbose) message("# Registering cleaning up function for '", 
                             name, "'... ", appendLF = FALSE)
        assign("cleanup", fun, fe)
        if (verbose) message("OK")
    } else if (exists("cleanup", envir = fe, inherits = FALSE)) {
        if (verbose) message("# Removing cleaning up function for '", 
                             name, "'... ", appendLF = FALSE)
        remove("cleanup", envir = fe)
        if (verbose) message("OK")
    }
    invisible(object)
}

# run cleanup procedure for a given backend object
doBackendCleanup <- function(object, ..., run = TRUE, verbose = FALSE){
    name <- object$name
    if (!is.null(fun <- object$cleanup)) {
        if (verbose) message("# Cleaning up '", name, "'... ", appendLF = FALSE)
        res <- try(fun(), silent=TRUE) 
        if (verbose) 
            message( if (is(res, "try-error")) "ERROR" else "OK")
        if (identical(TRUE, res)) object$cleanup <- NULL
        if (verbose) 
            message("OK", if (!is.null(res)) paste0(" [", res, "]"))
    }
    invisible(object)
}



