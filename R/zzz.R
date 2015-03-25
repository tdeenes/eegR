#
# <<< initialize >>> --------
#

#' @useDynLib eegR 
#' @import abind matrixStats permute parallel doParallel foreach iterators
#' @importFrom data.table setattr
#' @importFrom Rcpp evalCpp
#' @importFrom Kmisc str_collapse
#' @importFrom orthopolynom polynomial.values legendre.polynomials
#' @importFrom sgeostat in.polygon
#' @importFrom gplots colorpanel redgreen greenred bluered redblue
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
