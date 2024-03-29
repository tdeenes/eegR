% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova.R
\name{arrayAnova}
\alias{arrayAnova}
\title{Perform ANOVA (potentially with TFCE correction) on arrays}
\usage{
arrayAnova(
  .arraydat,
  factordef,
  bwdat = NULL,
  verbose = TRUE,
  perm = NULL,
  tfce = NULL,
  parallel = NULL,
  seed = NULL
)
}
\arguments{
\item{.arraydat}{a numeric array with named dimnames containing the EEG (or
other) data. Missing values are not allowed.}

\item{factordef}{a named list of factor definitions, containing the following
elements:
\itemize{
\item{between: }{character vector of between-subject factors (default: NULL)}
\item{within: }{character vector of within-subject factors (default: NULL)}
\item{w_id: }{name of the dimension which identifies the subjects
(default: "id")}
\item{observed: }{character vector of observed (i.e., not manipulated)
variables (default: NULL). Only used for effect-size calculations.}
}}

\item{bwdat}{a data.frame which contains the identification codes
(factordef$w_id) and all subject-level variables (usually factors) listed in
'factordef$between'. Missing values are not allowed.}

\item{verbose}{logical value indicating if p-values and effect sizes should
be computed for the traditional ANOVA results}

\item{perm}{either 1) NULL (the default) or FALSE (the same as NULL),
both of which mean no permutation, or 2) TRUE, which means
permutation with default parameters, or 3) an object as returned by
\code{\link{permParams}} with custom parameters (see Examples and also
\code{\link{permParams}}).\cr
Custom parameters can be also provided by \code{perm = .(key = value)} to
save typing (this works by calling \code{\link{permParams}} with the
given parameters).}

\item{tfce}{either 1) NULL (the default) or FALSE (the same as NULL), both of
which mean no TFCE correction, or 2) TRUE, which means TFCE correction with
default parameters, or 3) an object as returned by \code{\link{tfceParams}}
with custom TFCE parameters (see Examples and also \code{\link{tfceParams}}).\cr
Custom parameters can be also provided by \code{tfce = .(key = value)} to
save typing (this works by calling \code{\link{tfceParams}} with the given
parameters).}

\item{parallel}{either 1) NULL (the default) or FALSE (the same as NULL),
both of which mean single-core computation, or 2) TRUE, which means
parallelization with default parameters, or 3) an object as returned by
\code{\link{parallelParams}} with custom parameters (see Examples and also
\code{\link{parallelParams}}).\cr
Custom parameters can be also provided by \code{parallel = .(key = value)} to
save typing (this works by calling \code{\link{parallelParams}} with the
given parameters).}

\item{seed}{an integer value which specifies a seed (default: NULL), or a
list of arguments passed to \code{\link{set.seed}}}
}
\value{
A list object with the following numeric arrays:
\itemize{
\item{stat: }{a numeric array of F-values, with attribute 'Df' (an integer
matrix) containing the term- and residual degrees of freedom for each model
term. Optionally (if 'verbose' is TRUE), the F-value array has two extra
attributes: 'p_value' and 'effect_size', consisiting of the traditional
p-values and the generalized eta squared effect size measures
(see References), respectively.}
\item{stat_corr: }{a numeric array of TFCE-corrected F-values (if requested)}
\item{p_corr: }{a numeric array of permutation-based p-values (if requested)}
}
}
\description{
\code{arrayAnova} performs point-to-point ANOVAs on arrays. Permutation-based
p-values and Threshold-free Cluster Enhancement (TFCE) correction can be
requested.
}
\details{
The function assumes that the input array contains at least three
named dimensions: "chan" (corresponding to the channels [electrodes]), "time"
(corresponding to time points), and a subject identifier as given by
\code{factordef$w_id}. All dimensions which are not listed as
within-subject factors are treated in a similar way as chan and time, that is
separate ANOVA-s are computed for each level of those dimensions.
}
\note{
The function computes type I p-values - this is correct if the design
is fully balanced and orthogonal (if the number of between-subject
factors is one, it may have slightly unequal group sizes).
}
\examples{
# example dataset
data(erps)
dat_id <- attr(erps, "id") # to get group memberships
chan_pos <- attr(erps, "chan") # needed for TFCE correction

# make the dataset unbalanced to illustrate that if there is only one
# between-subject factor, Type 1 and Type 2 results are identical
erps <- subsetArray(erps, list(id = 2:20))
dat_id <- dat_id[2:20, ]

# average the data in each 12 ms time-bin to decrease the computational
# burden (not needed in serious analyses)
tempdat <- avgBin(erps, "time", 6)

# analyze the effect of the reading group (between-subject factor) and the
# two experimental conditions (stimclass, pairtye; within-subject factors)
# for each channel and time sample (without requesting TFCE correction); this
# means we run 1518 repeated measures ANOVAs
system.time(
    result_eegr <- arrayAnova(tempdat,
                              list(between = "group",
                                   within = c("stimclass", "pairtype"),
                                   w_id = "id",
                                   observed = "group"),
                              bwdat = dat_id)
)

# if package 'ez' is installed, you can compare the results; we take only a
# subset of the data (choosing "02" channel at time point "207") because
# ezANOVA is much slower
if (requireNamespace("ez")) {
    sub <- list(chan = "O2", time = "207")
    tempdat_ez <- transformArray(y ~ ., tempdat, subset = sub)
    tempdat_ez$group <- factor(dat_id$group[match(tempdat_ez$id, dat_id$id)])
    result_ez <- ez::ezANOVA(tempdat_ez, y, id, .(stimclass, pairtype),
                             between = group, observed = group, type = 2)
    # compare results
    ez_F <- result_ez$ANOVA$F   # F-values
    ez_p <- result_ez$ANOVA$p   # p-values
    ez_es <- result_ez$ANOVA$ges   # effect sizes
    eegr_F <- as.vector(
        subsetArray(extract(result_eegr, "stat"), sub))
    eegr_p <- as.vector(
        subsetArray(extract(result_eegr, "p"), sub))
    eegr_es <- as.vector(
        subsetArray(extract(result_eegr, "es"), sub))
    stopifnot(
        all.equal(ez_F, eegr_F),
        all.equal(ez_p, eegr_p),
        all.equal(ez_es, eegr_es)
    )
}

# the between-subject variable could be numeric, too. Let's create an 'age'
# variable and analyze the effects.
dat_id$age <- rnorm(nrow(dat_id), 10, 1)
result_eegr_c <- arrayAnova(tempdat,
                            list(between = "age",
                                 within = c("stimclass", "pairtype"),
                                 w_id = "id",
                                 observed = "age"),
                            bwdat = dat_id)

# if package 'car' is installed, you can compare the results; we take only a
# subset of the data (choosing "01" channel at time point "195")
if (requireNamespace("car")) {
    #
    # subsetting indices
    subs <- list(chan = "O1", time = "195")
    #
    # extract F values and p values from the eegR results
    eegr_F <- subsetArray(extract(result_eegr_c, "stat"), subs)
    eegr_p <- subsetArray(extract(result_eegr_c, "p"), subs)
    #
    # run the same analysis with the 'car' package (first we have to reshape
    # the data)
    tempdat_car <- subsetArray(tempdat, subs)
    tempdat_car <- mergeDims(tempdat_car,
                             list("id", c("stimclass", "pairtype")))
    tempdat_car <- data.frame(tempdat_car)
    tempdat_car$id <- dat_id$id
    tempdat_car$age <- dat_id$age
    idata <- expand.grid(dimnames(tempdat)[c("stimclass", "pairtype")])
    maov <- lm(cbind(A.ident, B.ident, C.ident,
                     A.subst, B.subst, C.subst,
                     A.transp, B.transp, C.transp) ~ age,
               data = tempdat_car)
    maov <- car::Anova(maov, idata = idata, idesign= ~stimclass*pairtype,
                       type = 2)
    maov <- summary(maov, multivariate = FALSE)$univ
    car_F <- maov[names(eegr_F), "F"]
    car_p <- maov[names(eegr_p), "Pr(>F)"]
    #
    # compare results
    stopifnot(
        all.equal(as.vector(car_F), as.vector(eegr_F)),
        all.equal(as.vector(car_p), as.vector(eegr_p))
    )
}

# in order to use TFCE correction, the channel neigbourhood matrix is needed
# (see ?tfceParams and ?chanNb)
ChN <- chanNb(chan_pos, alpha = 0.7)

# now analyze the data by collapsing the pairtypes, and apply TFCE correction
# (note: this will take a couple of seconds); use more randomization
# runs (n should be several thousand instead of 499L) in serious analyses
tempdat <- avgDims(tempdat, "pairtype")
result_tfce <- arrayAnova(tempdat,
                          list(between = "group",
                               within = "stimclass",
                               w_id = "id",
                               observed = "group"),
                          bwdat = dat_id,
                          perm = .(n = 499L),
                          tfce = .(ChN = ChN),
                          parallel = .(ncores = 2))

# plot the corrected and uncorrected results
modelplot(result_tfce)
modelplot(result_tfce, type = "unc")

# compare traditional and TFCE p-values
p_all <- extract(result_tfce, c("p", "p_corr"))
p_all <- bindArrays(trad = p_all$p, tfce = p_all$p_corr,
                    along_name = "method")

# plot p-values after -log transformation to increase discriminability;
# note how the sporadic effects disappear
p_plot <- imageValues(-log(p_all)) # returns a ggplot object
p_plot

}
\references{
The TFCE correction follows:\cr
Mensen, A. and Khatami, R. (2013)
Advanced EEG analysis using threshold-free cluster-enhancement and
non-parametric statistics. Neuroimage, 67, 111-118.
doi:10.1016/j.neuroimage.2012.10.027  \cr

The Generalized Eta Squared effect size statistic is described in:\cr
Olejnik, S., Algina, J. (2003) Generalized eta and omega squared statistics:
Measures of effect size for some common research designs. Psychological
Methods 8: pp. 434-447. doi:10.1037/1082-989X.8.4.434 \cr
Bakeman, R. (2005) Recommended effect size statistics for repeated measures
designs. Behavior Research Methods, 37 (3), 379-384.
}
\seealso{
See also the related methods to explore the results, e.g.
\code{\link{extract.arrayAnova}}, \code{\link{summary.arrayAnova}}, and the
plotting functions \code{\link{modelplot}}, or the lower-level
\code{\link{imageValues}} and \code{\link{imagePvalues}}.
}
