% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/erp_preprocess.R
\name{artifactRejection}
\alias{artifactRejection}
\title{Artifact rejection}
\usage{
artifactRejection(
  dat,
  markers = NULL,
  artrej_options = artrejOptions(),
  return_data = TRUE,
  return_details = TRUE,
  print_result = TRUE
)
}
\arguments{
\item{dat}{numeric array (EEG-data) with the following named dimensions 
(dimension order does not matter): chan, time, trial}

\item{markers}{if not NULL (default), a matrix or data.frame containing the 
characteristics of the trials (markers)}

\item{artrej_options}{a named list containing the parameters for the 
artifact rejection criteria. See \code{\link{artrejOptions}} for details.}

\item{return_data}{logical value, if TRUE (default), dat and markers without
rejected trials are returned}

\item{return_details}{logical value, if TRUE (default), the full array of 
results (e.g., bad trials for each channel and for each criterion) is 
returned as an attribute of bad_trials (see Values section)}

\item{print_result}{logical value, if TRUE (default), a summary of the 
results is printed to the console}
}
\value{
A named list containing bad_trials (trials identified with artifacts)
and the modified input data (dat and markers without contaminated trials)
}
\description{
\code{artifactRejection} performs artifact rejection on segmented data.
}
