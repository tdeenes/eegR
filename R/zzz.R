#
# <<< initialize >>> --------
#

#' A package to analyze EEG signals
#' 
#' The package 'eegR' has been developed to process electroencephalography (EEG)
#' signals. Beyond common EEG signal processing functionalities, eegR provides 
#' advanced tools to analyze single-trial and averaged event related potentials 
#' (ERPs).
#' IMPORTANT NOTE: This is only a pre-alpha version which will be restructured
#' substantially in short term. Do not rely on any of the functions.
#' @name eegR-package
#' @aliases eegR-package eegR
#' @author Denes Toth <\email{toth.denes@@ttk.mta.hu}>
#' @docType package
#' @useDynLib eegR 
#' @import abind data.table matrixStats ggplot2 permute parallel doParallel foreach iterators
#' @importFrom checkmate checkVector
#' @importFrom Rcpp evalCpp
#' @importFrom Kmisc factor_ counts
#' @importFrom orthopolynom polynomial.values legendre.polynomials
#' @importFrom sgeostat in.polygon
#' @importFrom gplots colorpanel redgreen greenred bluered redblue
#' @importFrom RColorBrewer brewer.pal
NULL

.onAttach <- function(lib, pkg) {
    packageStartupMessage(paste0("*** eegR ",
                                 utils::packageVersion("eegR"),
                                 " loaded ***"), 
                          appendLF = TRUE)
}

.onUnload <- function (libpath) {
    library.dynam.unload("eegR", libpath)
}
