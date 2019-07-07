#
# <<< initialize >>> --------
#

#' A package to analyze EEG signals
#' 
#' The eegR package has been developed to process
#' electroencephalography (EEG) signals. This is only a pre-alpha version
#' which has basic support for pre-processing EEG signals and provides
#' some more advanced tools to plot and analyze averaged event related
#' potentials (ERPs).
#' @name eegR-package
#' @aliases eegR-package eegR
#' @author Denes Toth <\email{toth.denes@@ttk.mta.hu}>
#' @docType package
#' @useDynLib eegR
#' @import abind data.table matrixStats ggplot2 permute proj4
#' @import parallel doParallel foreach iterators
#' @importFrom checkmate assertAtomic assertAtomicVector assertVector assertMatrix assertArray
#' @importFrom checkmate assertDataFrame assertList assertClass assertChoice
#' @importFrom checkmate assertLogical assertCharacter assertNumeric assertIntegerish 
#' @importFrom checkmate assertString assertNumber assertScalar assertCount assertFlag assertInt
#' @importFrom checkmate assertFunction
#' @importFrom checkmate checkAtomic checkAtomicVector checkVector checkMatrix checkArray
#' @importFrom checkmate checkDataFrame checkList checkClass checkChoice
#' @importFrom checkmate checkLogical checkCharacter checkNumeric checkIntegerish 
#' @importFrom checkmate checkString checkNumber checkScalar checkCount checkFlag checkInt
#' @importFrom checkmate checkNull
#' @importFrom checkmate testVector testIntegerish testNamed
#' @importFrom checkmate assert
#' @importFrom Rcpp evalCpp
#' @importFrom orthopolynom polynomial.values legendre.polynomials
#' @importFrom sgeostat in.polygon
#' @importFrom gplots colorpanel redgreen greenred bluered redblue
#' @importFrom RColorBrewer brewer.pal
#' @importFrom fdasrvf time_warping
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
