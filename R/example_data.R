#
# <<< document example data >>> --------
#

#' Averaged ERPs in a visual word recognition experiment
#' 
#' A dataset containing averaged event-related potentials (ERPs) from 20 
#' participants. The data were collected in a same-different matching task, and
#' simplified for present purposes. Only brain responses to the target stimuli 
#' from 10 dyslexic and 10 control participants are included. There were 3 types
#' of stimulus classes (A, B, and C) and 3 types of stimulus pairs (identical,
#' substituted, and transposed). The data were downsampled to 500 Hz, and cover 
#' the following time window: -50 ms to 500 ms. 
#' @docType data
#' @keywords datasets
#' @name erps
#' @format An array with 5 dimensions: stimulus class (stimclass, 3 level) x 
#' pair type (pairtype, 3 level) x channels (chan, 33 levels) x time (time, 276 
#' levels) x participant (id, 20 levels). Additionally, the array has an 
#' attribute "id" which is a data.frame of the group memberships (dyslexic or 
#' control), and an attribute "chan" which is a data.frame of the electrode 
#' positions in spherical coordinates.
NULL