#
# <<< ANOVA on arrays >>> ----
#

#' Compute marginal means in an ANOVA design
#' 
#' \code{marginalMeans} computes marginal means in ANOVA designs
#' @param form formula of the model
#' @param f_dat data.frame of factors
#' @param a_dat matrix which contains the dependent variables. It must have as
#' many rows as f_dat.  Usually a_dat is the result of calling array2mat on the
#' orginal data array.
#' @param dimn the original dimension names for the array corresponding to the
#' column names of a_dat
#' @param keep_term_order logival variable whether to keep the order of the 
#' model terms as provided in form (TRUE) or all interaction terms should follow
#' all main effect terms (FALSE, default)
#' @param residualmean logical variable; if TRUE, each factor term is extracted
#' from the data after computing the given marginal mean (default: FALSE)
#' @param whichterm character, numeric or logical vector; indices of model terms
#' which should be analyzed (default: NULL, meaning all terms are included)
#' @param no_icpt logical value; if TRUE, global mean is subtracted before 
#' computing the marginal means (default: FALSE)
#' @export
#' @return A named list containing the marginal means for each model terms
# TODO: examples
marginalMeans <- function(form, f_dat, a_dat, dimn, keep_term_order = FALSE, 
                          residualmean = FALSE, whichterm = NULL, 
                          no_icpt = FALSE) {
    labels <- attr(terms(as.formula(form), keep.order = keep_term_order), 
                   "term.labels")
    termL <- if (is.null(whichterm) || is.na(whichterm)) {
        strsplit(labels, ":")
    } else if (is.character(whichterm)) {    
        strsplit(whichterm[whichterm %in% labels], ":")
    } else {
        strsplit(labels[whichterm], ":")
    }
    marg_means <- vector("list", length(termL))
    if (residualmean && !no_icpt) a_dat <- sweep(a_dat, 2, colMeans(a_dat), "-")
    for (i in 1:length(termL)) {
        groups <- factor_(interaction(f_dat[,termL[[i]]], drop = TRUE))
        groupfreq <- tabulate(groups)
        names(groupfreq) <- levels(groups)
        marg_means[[i]] <- arrayIP(rowsum(a_dat, groups)/groupfreq, 
                                   c(length(groupfreq), vapply(dimn, length, 0L)),
                                   c(list(modelterm = names(groupfreq)), dimn))
        setattr(marg_means[[i]], "freq", groupfreq)
        if (residualmean) {
            a_dat <- a_dat - 
                as.vector(subsetArray(marg_means[[i]], 
                                      list(modelterm = as.numeric(groups))))
        }
    }
    names(marg_means) <- labels
    # return
    marg_means
}



#' Compute Sum of Squares
#' @keywords internal
sumsq <- function(form, f_dat, a_dat, dimn, keep_term_order = FALSE, 
                  whichterm = NULL, no_icpt = FALSE, return_means = FALSE,
                  return_array = TRUE) {
    labels <- attr(terms(as.formula(form), keep.order = keep_term_order), 
                   "term.labels")
    termL <- if (is.null(whichterm) || is.na(whichterm)) {
        strsplit(labels, ":")
    } else if (is.character(whichterm)) {    
        strsplit(whichterm[whichterm %in% labels], ":")
    } else {
        strsplit(labels[whichterm], ":")
    }
    ssq <- matrixIP(0, length(termL), ncol(a_dat))
    df <- rep.int(0L, length(termL))
    names(df) <- labels
    if (return_means) {
        marg_means <- vector("list", length(termL))
        names(marg_means) <- labels
    }
    #
    if (!no_icpt) a_dat <- sweep(a_dat, 2, colMeans(a_dat), "-")
    #
    for (i in 1:length(termL)) {
        groups <- factor_(interaction(f_dat[,termL[[i]]], drop = TRUE))
        groupfreq <- tabulate(groups)
        names(groupfreq) <- levels(groups)
        tmeans <- rowsum(a_dat, groups)/groupfreq
        a_dat <- a_dat - tmeans[as.numeric(groups),]
        ssq[i, ] <- colSums(groupfreq*tmeans^2)
        if (length(termL[[i]]) == 1) {
            df[i] <- length(groupfreq)-1
        } else {
            df[i] <- prod(df[termL[[i]]])
        }
        if (return_means) {
            marg_means[[i]] <- array(tmeans, 
                                     c(length(groupfreq), vapply(dimn, length, 0L)),
                                     c(list(modelterm = names(groupfreq)), dimn))
        }
    }
    #
    if (return_array) {
        arrayIP(ssq, c(nrow(ssq), vapply(dimn, length, 0L)),
                c(list(modelterm = labels), dimn))
    } else {
        rownames(ssq) <- labels
        names(dimnames(ssq))[1] <- "modelterm"
    }
    setattr(ssq, "Df", df)
    if (return_means) setattr(ssq, "term_means", marg_means)
    # return
    ssq
}

#' Parameter checks and data preparation in \code{arrayAnova}
#' 
#' \code{preAnova} performs parameter checking, creates a data.frame for 
#' the ANOVA design and reshapes arraydat accordingly, and prepares parallel
#' computations.
#' @param arraydat numeric array
#' @param factordef a list of between-subject factors (\code{between}), 
#' within-subject factors (\code{within}) and subject identifier (\code{w_id})
#' @param useparallel logical; do computations in parallel or not
#' @param ncores integer; number of cores for parallel computations
#' @param par_method character; parallelization method
#' @param usetfce logical; use TFCE correction or not (FALSE, default) 
#' @keywords internal
preAnova <- function(arraydat, factordef, bwdat, 
                     useparallel, ncores, par_method, usetfce = FALSE) {
    params <- list()
    # fast check of arguments
    if (is.null(factordef$w_id)) factordef$w_id <- "id"
    if (is.data.frame(arraydat)) arraydat <- as.matrix(arraydat)
    if (!is.array(arraydat) & !is.matrix(arraydat)) {
        stop("Improper input data class (should be data.frame, matrix or array))")
    }
    if (!is.null(factordef$between)) {
        if (!all(factordef$between %in% colnames(bwdat))) {
            stop("Between-subject factors and 'bwdat' dataset do not match.")
        }
        bwdat[,factordef$between] <- lapply(
            data.frame(bwdat[,factordef$between]), factor)
    }
    if (useparallel) {
        params$ncores <- ncores
        if (par_method == "snow") {
            fns <- ls(envir = parent.env(environment()))
            fns <- fns[sapply(fns, function(x) is.function(get(x)))]
            params$varlist2snow <-  
                if (usetfce) {
                    c(fns, "randind", "ChN","EH") 
                } else {
                    c(fns, "randind")
                }
        } 
    }
    #
    op <- options(contrasts = c("contr.sum", "contr.poly"))
    origdimnames <- dimnames(arraydat)
    modeldims <- c(factordef$w_id, factordef$within)
    keepdims <- setdiff(names(origdimnames), modeldims)
    arraydat <- mergeDims(arraydat, list(modeldims, keepdims))
    dat <- data.frame(matrix(unlist(strsplit(rownames(arraydat), "[.]"), 
                                    use.names = FALSE), 
                             nrow = nrow(arraydat), byrow = TRUE))
    colnames(dat) <- modeldims
    dat[,factordef$between] <- bwdat[
        match(dat[,factordef$w_id],bwdat[,factordef$w_id]),factordef$between]
    # return
    list(dat = dat, arraydat = arraydat, factordef = factordef, 
         origdimnames = origdimnames, par_params = params)
}

#' Workhorse function of arrayAnova
#' 
#' \code{arrayAnovaSub} computes F-values in ANOVA designs, and provides 
#' additional p-values, degrees of freedom, generalized effect sizes if 
#' requested.
#' @param a_dat a numeric matrix of responses
#' @param f_def a named list of factor definitions, containing the following 
#' elements:
#' \itemize{
#' \item{between:}{character vector of between-subject factors (default: NULL)}
#' \item{within:}{character vector of within-subject factors (default: NULL)}
#' \item{w_id:}{name of the dimension which identifies the subjects 
#' (default: "id")}
#' }
#' @param d_names (list of) dimension names of the original array in 
#' \code{arrayAnova}
#' @param f_dat a data.frame which contains the identification codes 
#' (f_def$w_id) and the between- and/or within-subject factors
#' @param verbose logical value indicating if p-values and effect sizes should 
#' be computed for the traditional ANOVA results
#' @keywords internal
arrayAnovaSub <- function(a_dat, f_def, d_names, f_dat, verbose = TRUE) {
    origdims <- vapply(d_names, length, 0L)
    modeldims <- c(f_def$w_id, f_def$within)
    keepdims <- setdiff(names(d_names), modeldims)
    #
    # distinct algorithms for between-subject, within-subject, and mixed models
    if (is.null(f_def$within)) {
        # SHOULD BE IMPROVED - use it with one between-subject factor and a 
        # balanced design
        aov_formula <- as.formula(paste(
            as.character(quote(a_dat)), "~", 
            paste(as.character(f_def$between), collapse = "*")))
        results <- summary(aov(aov_formula, data = f_dat))
        results <- arrayIP(
            unlist(results, use.names = FALSE),
            c(dim(results[[1]]), length(results)),
            list(modelterm = gsub(" ", "", rownames(results[[1]])),
                 measures = colnames(results[[1]]),
                 columns = gsub("[*Response ]", "", names(results))))
        termrows <- tolower(rownames(results)) != "residuals"
        Fvals <- arrayIP(results[termrows,"F value",], 
                         c(sum(termrows), origdims[keepdims]),
                         c(list(modelterm = rownames(results)[termrows]),
                           d_names[keepdims]))
        if (verbose) {
            Df <- as.numeric(results[termrows, "Df", 1])
            rDf <- results[!termrows, "Df", 1]
            pvalues <- results[termrows,"Pr(>F)",]
            # effect size (Olejnik and Algina, 2003; Bakeman, 2005)
            results <- subsetArray(results, list(measures = "Sum Sq"), drop = FALSE)
            matrixIP(results, nrow(results), dim(results)[3], 
                     dimnames = dimnames(results)[-2])
            SSr <- results["Residuals", ]
            if (!is.null(f_def$observed)) {
                ind <- apply(sapply(f_def$observed, 
                                    function(x) 
                                        grepl(x, rownames(results))),
                             1, any)
                SS1 <- colSums(results[ind, , drop = F])
                SS2 <- sweep(results, 1, ind, "*")
                ges <- results / sweep(sweep(results - SS2, 2, SSr, "+"), 2,
                                       SS1, "+")
            } else {
                ges <- results / sweep(results, 2, SSr, "+")
            }
            ges <- ges[-nrow(ges), ]
            factor_means <- marginalMeans(aov_formula, f_dat, a_dat, 
                                          d_names[keepdims], residualmean = FALSE)
        }
    } else if (is.null(f_def$between)) {
        # compute marginal means
        if (verbose) {
            aov_formula <- as.formula(paste(
                as.character(quote(a_dat)), "~", 
                paste(as.character(f_def$within), collapse = "*")
            )
            )
            factor_means <- marginalMeans(aov_formula, f_dat, a_dat, 
                                          d_names[keepdims], residualmean = FALSE)
        }
        # model formula
        aov_formula <- as.formula(paste(
            as.character(quote(a_dat)), "~", 
            paste(
                f_def$w_id, 
                paste(as.character(f_def$within), collapse = "*"),
                sep = "*")
        )
        )
        #
        # run model
        #
        a_datsc <- sweep(a_dat, 2, colMeans(a_dat), "-")
        model <- sumsq(aov_formula, f_dat, a_datsc, d_names[keepdims], 
                       no_icpt = TRUE, return_array = FALSE)
        model_df <- attr(model, "Df")
        # compute indices to handle error terms
        modnames <- dimnames(model)$modelterm
        m_ind <- !grepl(f_def$w_id, modnames)
        err_ind <- match(
            paste(f_def$w_id, ":", modnames[m_ind], sep = ""),
            modnames)
        Df <- model_df[m_ind]
        rDf <- model_df[err_ind]
        Fvals <- (model[m_ind, , drop = F] / Df) / 
            (model[err_ind, , drop = F] / rDf)
        if (verbose) {
            pvalues <- mapply(pf, Fvals, Df, rDf, MoreArgs = list(lower.tail = FALSE))
            dim(pvalues) <- dim(Fvals)
            # effect size (Olejnik and Algina, 2003; Bakeman, 2005)
            SSr <- colSums(model[grep(f_def$w_id, 
                                      rownames(model)), ])
            ges <- model[m_ind, , drop = F] / 
                sweep(model[m_ind, , drop = F], 2, SSr, "+")
        }
        arrayIP(Fvals, c(sum(m_ind), origdims[keepdims]),
                c(list(modelterm = modnames[m_ind]), d_names[keepdims]))
    } else {   
        aov_formula1 <- as.formula(paste(
            as.character(quote(a_dat)), "~", 
            paste(
                paste(as.character(f_def$between), collapse = "*"), 
                paste(as.character(f_def$within), collapse = "*"),
                sep = "*")
        )
        )
        aov_formula2 <- as.formula(paste(
            as.character(quote(a_dat)), "~", 
            paste(
                f_def$w_id, 
                paste(as.character(f_def$within), collapse = "*"),
                sep = "*")
        )
        )
        # compute marginal means
        if (verbose) 
            factor_means <- marginalMeans(aov_formula1, f_dat, a_dat, 
                                          d_names[keepdims], residualmean = FALSE)
        #
        # run model
        #
        a_datsc <- sweep(a_dat, 2, colMeans(a_dat), "-")
        model1 <- sumsq(aov_formula1, f_dat, a_datsc, d_names[keepdims], 
                        no_icpt = TRUE, return_array = FALSE)
        model1_df <- attr(model1, "Df")
        model2 <- sumsq(aov_formula2, f_dat, a_datsc, d_names[keepdims], 
                        no_icpt = TRUE, return_array = FALSE)
        model2_df <- attr(model2, "Df")
        # compute indices to handle error terms
        modnames <- dimnames(model2)$modelterm
        m2indices <- grep(f_def$w_id, modnames)
        modnames <- modnames[m2indices]
        bwnames <- unlist(lapply(seq_along(f_def$between), 
                                 function(i) apply(combn(f_def$between, i), 2, 
                                                   paste, collapse = ":")))
        sumindices <- sapply(bwnames, function(i) 
            match(sub(f_def$w_id, i, modnames), dimnames(model1)$modelterm))
        modnames <- sub(paste(f_def$w_id,"\\:?", sep = ""), "", modnames)
        fillindices <- lapply(1:nrow(sumindices), function(i)
            c(sumindices[i,], 
              na.omit(match(modnames[i], 
                            dimnames(model1)$modelterm))))
        # compute F values
        model1_err <- arrayIP(0, dim(model1))
        model1_err_df <- rep.int(0L, length(model1_df))
        for (i in 1:nrow(sumindices)) {
            errs <- model2[m2indices[i],] - colSums(model1[sumindices[i,], , drop = F])
            model1_err[fillindices[[i]], ] <- 
                rep(c(errs), each = length(fillindices[[i]])) 
            model1_err_df[fillindices[[i]]] <- model2_df[m2indices[i]] - 
                sum(model1_df[sumindices[i, ]])
        }
        Fvals <- model1/model1_df*model1_err_df/model1_err
        if (verbose) {
            Df <- as.numeric(model1_df)
            rDf <- model1_err_df
            pvalues <- mapply(pf, Fvals, Df, rDf, MoreArgs = list(lower.tail = FALSE))
            # effect size (Olejnik and Algina, 2003; Bakeman, 2005)
            SSr <- colSums(model1_err[grep(f_def$w_id, 
                                           rownames(model2)), ])
            if (!is.null(f_def$observed)) {
                ind <- apply(sapply(f_def$observed, 
                                    function(x) grepl(x, rownames(model1))),
                             1, any)
                SS1 <- colSums(model1[ind, ])
                SS2 <- sweep(model1, 1, ind, "*")
                ges <- model1 / sweep(sweep(model1 - SS2, 2, SSr, "+"), 2,
                                      SS1, "+")
            } else {
                ges <- model1 / sweep(model1, 2, SSr, "+")
            }
        }
        dimnms <- c(list(modelterm = rownames(Fvals)), d_names[keepdims])
        arrayIP(Fvals, c(nrow(Fvals), origdims[keepdims]), dimnms)
    }
    # rearrange dimensions to be compatible with arrayTtest output
    if (verbose) {
        tempfn <- function(x) aperm(array(x, dim(Fvals), dimnames(Fvals)), 
                                    c(keepdims, "modelterm"))
        pvalues <- tempfn(pvalues)
        ges <- tempfn(ges)
    }
    Fvals <- aperm(Fvals, c(keepdims, "modelterm"))
    if (verbose) {
        setattr(Fvals, "Df.term", Df)
        setattr(Fvals, "Df.resid", rDf)
        setattr(Fvals, "pvalues", pvalues)
        setattr(Fvals, "ges", ges)
        setattr(Fvals, "factor_means", factor_means)
    }
    # return
    Fvals
}

#' Perform ANOVA (potentially with TFCE correction) on arrays
#' 
#' \code{arrayAnova} performs point-to-point ANOVAs on arrays. Permutation-based
#' p-values and Threshold-free Cluster Enhancement (TFCE) correction can be requested.
#' @param arraydat a numeric array with named dimnames containing the EEG (or 
#' other) data. Missing values are not allowed.
#' @param factordef a named list of factor definitions, containing the following 
#' elements:
#' \itemize{
#' \item{between:}{character vector of between-subject factors (default: NULL)}
#' \item{within:}{character vector of within-subject factors (default: NULL)}
#' \item{w_id:}{name of the dimension which identifies the subjects 
#' (default: "id")}
#' }
#' @param bwdat a data.frame which contains the identification codes 
#' (factordef$w_id) and the group memberships (factordef$between) of the 
#' subjects. Missing values are not allowed.
#' @param verbose logical value indicating if p-values and effect sizes should 
#' be computed for the traditional ANOVA results
#' @param nperm integer value giving the number of permutations (default: 999L)
#' @param useparallel logical value; if TRUE (default), computations are done
#' in parallel
#' @param ncores integer value corresponding to the number of cores; 
#' if NULL (default), it is set to the maximum number of cores available
#' @param par_method parallelization method; can be "snow" (default) or "mc"
#' (multicore)
#' @param cl a cluster definition for snow-type parallelization; if NULL 
#' (default), it is set up automatically
#' @param usetfce logical value whether TFCE (threshold-free cluster enhancement)
#' correction should also be computed (default: TRUE)
#' @param tfce_options a named list containing the channel neighbourhood matrix
#' (named ChN) and the vector of the E/H parameters (named EH)
#' @param seed an integer value which specifies a seed (default: NULL)
#' @details The function assumes that the input array contains at least two 
#' named dimensions: chan (corresponding to the channels [electrodes]) and time 
#' (corresponding to time points). All dimensions which are not listed as 
#' within-subject factors are treated in a similar way as chan and time, that is 
#' separate ANOVA-s are computed for each level of those dimensions.
#' @note The function computes type I p-values - this is correct if the design
#' is fully balanced and orthogonal (if the number of between-subject 
#' factors is one, it may have slightly unequal group sizes).
#' @export
#' @return A list object with F values, TFCE-corrected F-values and
#' permutation-based p-values (if requested)
#' @references The TFCE correction follows Mensen, A. and Khatami, R. (2013):
#' Advanced EEG analysis using threshold-free cluster-enhancement and 
#' non-parametric statistics. Neuroimage, 67, 111-118. 
#' doi:10.1016/j.neuroimage.2012.10.027 
#' @examples
#' # example dataset
#' data(erps)
#' dat_id <- attr(erps, "id") # to get group memberships
#' chan_pos <- attr(erps, "chan") # needed for TFCE correction
#' 
#' # average the data in each 12 ms time-bin to decrease the computational 
#' # burden (not needed in serious analyses)
#' tempdat <- avgBin(erps, "time", 6)
#' 
#' # analyze the effect of the reading group (between-subject factor) and the
#' # two experimental conditions (stimclass, pairtye; within-subject factors) 
#' # for each channel and time sample (without requesting TFCE correction); this
#' # means we run 1518 repeated measures ANOVAs
#' system.time(
#'     result_eegr <- arrayAnova(tempdat, 
#'                               list(between = "group",
#'                                    within = c("stimclass", "pairtype"),
#'                                    w_id = "id"),
#'                               bwdat = dat_id,
#'                               nperm = 0L, usetfce = FALSE)
#' )
#' 
#' # if package 'ez' is installed, you can compare the results; we take only a 
#' # subset of the data because ezANOVA is much slower
#' if (require(ez)) {
#'     sub <- list(chan = "O2", time = "207")
#'     tempdat_ez <- transformArray(y ~ ., tempdat, subset = sub)
#'     tempdat_ez$group <- factor(dat_id$group[match(tempdat_ez$id, dat_id$id)])
#'     result_ez <- ezANOVA(tempdat_ez, y, id, .(stimclass, pairtype),
#'                          between = group)
#'     # compare results
#'     ez_F <- result_ez$ANOVA$F   # F values
#'     ez_p <- result_ez$ANOVA$p   # p-values
#'     ez_ges <- result_ez$ANOVA$ges   # generalized effect sizes
#'     eegr_F <- as.vector(
#'         subsetArray(result_eegr$effect_F_obs, sub))
#'     eegr_p <- as.vector(
#'         subsetArray(attr(result_eegr$effect_F_obs, "pvalues"), sub))
#'     eegr_ges <- as.vector(
#'         subsetArray(attr(result_eegr$effect_F_obs, "ges"), sub))
#'     stopifnot(
#'         all.equal(ez_F, eegr_F),
#'         all.equal(ez_p, eegr_p),
#'         all.equal(ez_ges, eegr_ges)
#'     )
#' }
#' 
#' # in order to use TFCE correction, the channel neigbourhood matrix is needed
#' # (see ?chanNb)
#' ChN <- chanNb(chan_pos, alpha = 0.7)
#' 
#' # now analyze the data by collapsing the pairtypes, and apply TFCE correction 
#' # (note: this will take a couple of seconds); use more randomization 
#' # runs (nperm should be several thousand) in serious analyses
#' tempdat <- avgDims(tempdat, "pairtype")
#' result_tfce <- arrayAnova(tempdat, 
#'                           list(between = "group",
#'                                within = "stimclass",
#'                                w_id = "id"),
#'                           bwdat = dat_id,
#'                           nperm = 499L, useparallel = TRUE, ncores = 2,
#'                           tfce_options = list(ChN = ChN))
#' 
#' # compare traditional and TFCE p-values
#' p_all <- extract(result_tfce, c("p", "p_corr"))
#' p_all <- bindArrays(trad = p_all$p, tfce = p_all$p_corr, 
#'                     along_name = "method")
#' 
#' # plot p-values after -log transformation to increase discriminability;
#' # note how the sporadic effects disappear
#' p_plot <- imageValues(-log(p_all)) # returns a ggplot object
#' p_plot
#' 
#' # you can also use the imagePvalues function to plot discretized p-values
#' imagePvalues(p_all, pcrit = c(0.01, 0.05, 0.1))
#'
arrayAnova <- function(arraydat, factordef, bwdat = NULL, verbose = TRUE, 
                       nperm = 999L, useparallel = FALSE, ncores = NULL, 
                       par_method = "snow", cl = NULL,
                       usetfce = TRUE, tfce_options = NULL, seed = NULL) {
    # some checks
    stopifnot(is.array(arraydat))
    stopifnot(is.list(factordef))
    nperm <- as.integer(nperm)
    if (!is.null(factordef$between) && is.null(bwdat)) {
        stop("No between-participant data provided")
    } 
    if (usetfce) {
        if (is.null(tfce_options)) stop("Provide at least tfce_options$ChN")
        ChN <- tfce_options$ChN
        EH <- if (is.null(tfce_options$EH)) c(0.66, 2) else tfce_options$EH
    }
    if (useparallel) {
        if (is.null(ncores)) ncores <- detectCores()
    }
    #
    out <- list(call = match.call())
    #
    #
    temp <- preAnova(arraydat, factordef, bwdat, 
                     useparallel, ncores, par_method, usetfce)
    dat <- temp$dat
    arraydat <- temp$arraydat
    factordef <- temp$factordef
    origdimnames <- temp$origdimnames
    par_params <- temp$par_params
    rm(temp)
    #
    Fvals_obs <- arrayAnovaSub(arraydat, factordef, origdimnames, 
                               dat, verbose)
    mergedimnames <- setdiff(names(dimnames(Fvals_obs)), c("chan", "time"))
    if (usetfce) {
        nrchan <- length(origdimnames$chan)
        nrtime <- length(origdimnames$time)
        tfce_obs <- fnDims(Fvals_obs, c("chan", "time"), 
               function(x) tfceFn(matrixIP(x, nrchan, nrtime), 
                                  chn = ChN, eh = EH), 
               keep_dimorder = TRUE)
    }
    if (nperm > 1L) {
        #
        permfn <- function(i) {
            x <- arrayAnovaSub(arraydat, 
                               factordef, origdimnames, 
                               dat[randind[i,], ], FALSE)
            if (usetfce) {
                fnDims(x, c("chan", "time"), 
                       function(xx) max(tfceFn(matrixIP(xx, nrchan, nrtime), 
                                              chn = ChN, eh = EH)), 
                       keep_dimorder = TRUE)
            } else {
                fnDims(x, c("chan", "time"), max)
            }
        }
        #
        # generate random orders (dim(randind) = nperm X nrow(dat))
        if (!is.null(seed)) set.seed(seed)
        randind <- shuffleSet(nrow(dat), nperm, 
                              how(within = Within(type = "free"), 
                                  plots = Plots(strata = dat[,factordef$w_id], 
                                                type = "free")))
        #
        if (!useparallel) {
            maxperm <- lapply(1:nperm, permfn)
        } else if (par_method == "snow") {
            if (is.null(cl)) cl <- makePSOCKcluster(par_params$ncores, 
                                                    outfile = "")
            clusterExport(cl, 
                          varlist = par_params$varlist2snow,
                          envir = environment())
            maxperm <- parLapply(cl, 1:nperm, permfn)
            stopCluster(cl)
            rm(cl)
        } else if (par_method == "mc") {
            maxperm <- mclapply(1:nperm, permfn, mc.cores = par_params$ncores)
        } 
        maxperm <- matrixIP(unlist(maxperm, use.names = FALSE), ncol = nperm)
        sig <- if (usetfce) tfce_obs else Fvals_obs
        sig <- mergeDims(sig, mergedimnames)
        for (i in 1:nrow(sig)) {
            sig[i,,] <- (rowSums(outer(c(sig[i,,]), maxperm[i,], "<="))+1) / 
                (nperm + 1)
        }
        sig <- revMergeDims(sig)
        type <- if (usetfce) "tfce" else "perm"
        setattr(sig, "type", type)
    } 
    out$effect_F_obs <- Fvals_obs
    if (usetfce) out$effect_tfce_obs <- tfce_obs
    if (nperm > 1L) out$perm_pvalues <- sig
    setattr(out, "class", "arrayAnova")
    # return
    out
}

#' Extract interaction
#' 
#' \code{extractInteraction} computes the highest-order (interaction) effect 
#' (i.e. the difference between differences) from means of each cell of the 
#' design matrix
#' @param dat a named list with array elemens, see details below
#' @param sep the character which separates the name of factor levels in the 
#' model terms (default: ".")
#' @param sep_fixed logical indicating as sep should be treated "as is" while 
#' splitting the model term names
#' @details This function is not intended for direct use, and it only works for 
#' specifically formatted data. The input data is supposed to be the return 
#' value of \code{\link{marginalMeans}}; a named list which contains the 
#' marginal means for each effect of an ANOVA. The names are the names of the 
#' factor(s) - for interaction effects the factor names are separated by ":". 
#' The list elements are three-dimensional arrays with the following dimension 
#' names: modelterm, chan and time. For pure main effects, this function 
#' computes all pairwise differences between the levels of the factor. For 
#' interaction effects, the function computes all pairwise differences at the 
#' level of the highest order interaction. That means for the A:B interaction 
#' where factor A has levels Aa & Ab and factor B has levels Ba, Bb and Bc, 
#' the function computes the following differences: (Aa-Ab) - (Ba-Bb), 
#' (Aa-Ab) - (Ba-Bc), (Aa-Ab) - (Bb-Bc).
#' @export
#' @return A named list is returned with the same names as the input, but the 
#' list elements are ia_level x chan x time arrays.  
extractInteraction <- function(dat, sep = ".", sep_fixed = TRUE) {
    # function to compute the weights (+1 or -1) of the factor levels (f)
    iasign <- function(f) {
        combin <- lapply(f, function(x) combn(levels(x), 2))
        combin_ind <- expand.grid(lapply(combin, 
                                         function(x) seq.int(ncol(x))))
        combin_names <- expand.grid(lapply(f, function(x) 
            combn(levels(x), 2, FUN = paste, collapse = "-")))
        combin_names <- apply(combin_names, 1, paste, collapse = ".")
        out <- sapply(1:nrow(combin_ind), function(ii) {
            out <- rep(0, nrow(f))
            ind <- sapply(names(combin), function(x) 
                f[, x] %in% combin[[x]][, combin_ind[ii, x]])
            ind <- apply(ind, 1, all)
            if (sum(ind) > 2) {
                tempf <- sapply(droplevels(f[ind, , drop = F]), as.numeric)
                tempf <- sign(1.5 - tempf)
                out[which(ind)] <- apply(tempf, 1, prod)
            } else {
                out[which(ind)] <- c(1, -1)
            }
            return(out)
        })
        colnames(out) <- combin_names
        return(out)
    }
    # function to compute the interaction effects
    compia <- function(x, modterm_name) {
        modterm <- dimnames(x)$modelterm
        facs <- strsplit(modterm, sep, fixed = sep_fixed)
        facs <- data.frame(matrix(unlist(facs), 
                                  ncol = length(facs[[1]]), 
                                  byrow = TRUE))
        colnames(facs) <- unlist( strsplit(modterm_name, ":") )
        signs <- iasign(facs)
        rownames(signs) <- modterm
        avg_ia <- fnDims(x, "modelterm", 
                         function(xx, y) t( crossprod(xx, y) ), 
                         list(y = signs),
                         newdims = list(ia_level = colnames(signs)),
                         vectorized = TRUE)
        return( avg_ia )
    }
    # ------
    # run computations
    out <- mapply(compia, dat, names(dat), SIMPLIFY = FALSE)
    # return
    out
}

# TANOVA functions =========== 

#' Topographical ANOVA (TANOVA) and related methods
#' 
#' \code{tanova} performs point-to-point topographical ANOVA on arrays. Related
#' methods are GFP-analysis (which is based on the intensity of the signal) and
#' DISS-analysis (focusing only on the global dissimilarity by comparing 
#' normalized topographies). See References for further details. 
#' @param arraydat a numeric array with named dimnames containing the EEG (or 
#' other) data. Missing values are not allowed.
#' @param factordef a named list of factor definitions, containing the following 
#' elements:
#' \itemize{
#' \item{"between"}{"character vector of between-subject factors (default: NULL)"}
#' \item{"within"}{"character vector of within-subject factors (default: NULL)"}
#' \item{"w_id"}{"name of the dimension which identifies the subjects 
#' (default: "id")}
#' }
#' @param bwdat a data.frame which contains the identification codes 
#' (factordef$w_id) and the group memberships (factordef$between) of the 
#' subjects. Missing values are not allowed.
#' @param type a character value of "tanova" (default), "dissimilarity", or 
#' "gfp"
#' @param verbose logical value indicating whether means for each factor level 
#' should be also returned
#' @param nperm integer value giving the number of permutations (default: 999L)
#' @param useparallel logical value; if TRUE (default), computations are done
#' in parallel
#' @param ncores integer value corresponding to the number of cores; 
#' if NULL (default), it is set to the maximum number of cores available
#' @param par_method parallelization method; can be "snow" (default) or "mc"
#' (multicore)
#' @param cl a cluster definition for snow-type parallelization; if NULL 
#' (default), it is set up automatically
#' @param iaterms_last logival variable whether all interaction terms should 
#' follow all main effect terms (default: TRUE)
#' @param seed an integer value which specifies a seed (default: NULL)
#' @param pcrit the significance level (or a vector of multiple levels) for the 
#' consecutive length correction (default: 0.05)
#' @details The function assumes that the input array contains at least two 
#' named dimensions: chan (corresponding to the channels [electrodes]) and time 
#' (corresponding to time points). All dimensions which are not listed as 
#' within-subject factors are treated in a similar way as chan and time, that is 
#' separate TANOVA-s are computed for each level of those dimensions.
#' @note The function computes type I p-values - this is correct if the design
#' is fully balanced and orthogonal (if the number of between-subject 
#' factors is one, it may have slightly unequal group sizes).
#' @export
#' @return A list object with effect statistics, uncorrected p-values, and 
#' two types of corrected p-values: length correction or global correction
#' @seealso \code{\link{plotTanova}} plots the result of \code{tanova}.
#' @references Ported from the MATLAB toolbox by Koenig T, Kottlow M, Stein M, 
#' Melie García L (2011) Ragu: A free tool for the analysis of EEG and MEG 
#' event-related scalp field data using global randomization statistics. 
#' Computational Intelligence and Neuroscience, 2011:938925
#' @examples
#' # example dataset
#' data(erps)
#' dat_id <- attr(erps, "id") # to get group memberships
#' 
#' # average the data in each 12 ms time-bin to decrease the computational 
#' # burden (not needed in serious analyses)
#' tempdat <- avgBin(erps, "time", 6)
#' 
#' # Analyze the effect of the reading group (between-subject factor) and the
#' # two experimental conditions (stimclass, pairtye; within-subject factors) 
#' # for each channel and time sample, taking into account both intensity and 
#' # distributional differences (the default: type = "tanova")
#' # Note that the number of permutations should be increased for serious 
#' # purposes.
#' result_tanova <- tanova(tempdat, 
#'                         list(between = "group",
#'                              within = c("stimclass", "pairtype"),
#'                              w_id = "id"),
#'                         bwdat = dat_id,
#'                         nperm = 499L,
#'                         useparallel = TRUE, ncores = 2)
#' 
#' # plot results (for now, only p-values)
#' plotTanova(result_tanova, only_p = TRUE)
tanova <- function(arraydat, factordef, bwdat = NULL, 
                   type = c("tanova", "dissimilarity", "gfp"), 
                   verbose = TRUE, nperm = 999L, useparallel = FALSE, ncores = NULL,
                   par_method = c("snow", "mc"), cl = NULL,
                   iaterms_last = TRUE, seed = NULL, pcrit = 0.05) {
    # some checks
    stopifnot(is.array(arraydat))
    stopifnot(is.list(factordef))
    nperm <- as.integer(nperm)
    if (!is.null(factordef$between) && is.null(bwdat)) {
        stop("No between-participant data provided")
    }
    if (useparallel) {
        if (is.null(ncores)) ncores <- detectCores()
    }
    pcrit <- pcrit[pcrit > (1 / (nperm + 1))]
    if (length(pcrit) < 1) stop("The number of permutations is too low for this pcrit")
    #
    out <- list(call = match.call())
    #
    effSize <- function(means) {
        datfn <- if (type == "gfp") abs else compGfp
        out <- lapply(means, function(x) avgDims(datfn(x), "modelterm"))
        out <- rearrangeList(out, "modelterm")
        out <- drop(out)
        return( out )
    }
    # 
    consecLimit <- function(sigmat, alpha) {
        out <- matrixRle(sigmat)
        quantile(out$lengths[out$values > 0], 1 - alpha)
    }
    #
    par_method <- match.arg(par_method)
    #
    type <- match.arg(type)
    if (type == "gfp") {
        if (verbose) {
            arraydat_orig <- arraydat
            arraydat_orig <- preAnova(arraydat, factordef, bwdat, 
                                      useparallel, ncores, par_method)$arraydat
        }
        arraydat <- compGfp(arraydat)
    }
    if (type == "dissimilarity") arraydat <- scaleChan(arraydat)
    #
    temp <- preAnova(arraydat, factordef, bwdat, useparallel, ncores, par_method)
    # assign variables to this environment, and potentially overwrite existing ones
    #assignList(temp, verbose = FALSE)
    dat <- temp$dat
    arraydat <- temp$arraydat
    factordef <- temp$factordef
    origdimnames <- temp$origdimnames
    par_params <- temp$par_params
    rm(temp)
    # formula of anova
    aov_formula <- as.formula(paste(
        as.character(quote(arraydat)), "~", 
        sub("^\\*|\\*$", "", paste(
            paste(as.character(factordef$between), collapse = "*"), 
            paste(as.character(factordef$within), collapse = "*"),
            sep = "*"))
    )
    )
    # compute marginal means
    origdims <- vapply(origdimnames, length, 0L)
    modeldims <- c(factordef$w_id, factordef$within)
    keepdims <- setdiff(names(origdimnames), modeldims)
    term_means <- marginalMeans(aov_formula, dat, arraydat, 
                                origdimnames[keepdims], 
                                keep_term_order = !iaterms_last, 
                                residualmean = TRUE)
    es_obs <- effSize(term_means)
    out <- c(out, list(effect = es_obs, residual_means = term_means))
    names(out)[1] <- paste(names(out)[1], type, sep = "_")
    if (verbose) {
        if (type == "gfp") {
            factor_means <- marginalMeans(aov_formula, dat, arraydat_orig, 
                                          origdimnames[keepdims], 
                                          keep_term_order = !iaterms_last, 
                                          residualmean = FALSE)
            
        } else {
            factor_means <- marginalMeans(aov_formula, dat, arraydat, 
                                          origdimnames[keepdims], 
                                          keep_term_order = !iaterms_last, 
                                          residualmean = FALSE)
        }
        out <- c(out, list(factor_means = factor_means))
    }
    if (nperm > 1L) {
        #
        permfn <- function(i) {
            x <- marginalMeans(aov_formula, 
                               dat[randind[i, ], ], 
                               arraydat,
                               origdimnames[keepdims], 
                               keep_term_order = !iaterms_last, 
                               residualmean = TRUE)
            # return
            effSize(x)
        }
        #
        # generate random orders (dim(randind) = nperm X nrow(dat))
        if (!is.null(seed)) set.seed(seed)
        randind <- shuffleSet(nrow(dat), nperm, 
                              how(within = Within(type = "free"), 
                                  plots = Plots(strata = dat[,factordef$w_id], 
                                                type = "free")))
        #
        if (!useparallel) {
            es_perm <- lapply(1:nperm, permfn)
        } else if (par_method == "snow") {
            if (is.null(cl)) cl <- makePSOCKcluster(par_params$ncores, outfile = "")
            clusterExport(cl, 
                          varlist = par_params$varlist2snow,
                          envir = environment())
            es_perm <- parLapply(cl, 1:nperm, permfn)
            stopCluster(cl)
            rm(cl)
        } else if (par_method == "mc") {
            es_perm <- mclapply(1:nperm, permfn, mc.cores = par_params$ncores)
        } 
        # permuted F values
        es_perm <- c(list(es_obs), es_perm)
        names(es_perm) <- seq.int(nperm + 1)
        es_perm <- rearrangeList(es_perm, "perm")
        #if (!is.null(perm_fname)) save(es_perm, file = perm_fname)
        # p-values
        pvalues_perm <- fnDims(es_perm, "perm", colRanks, 
                               arg_list = list(preserveShape = TRUE),
                               vectorized = TRUE, keep_dimorder = TRUE)
        pvalues_perm <- (nperm + 2 - pvalues_perm) / (nperm + 1)
        pvalues <- avgDims(subsetArray(pvalues_perm, 
                                       list(perm = 1), drop = FALSE),
                           "perm")
        # consecutive sign. criterion
        pvalues_perm <- aperm(pvalues_perm, 
                              c("time", "perm", 
                                setdiff(names(dimnames(pvalues_perm)),
                                        c("time", "perm"))))
        temp <- dimnames(pvalues_perm)[-2]
        pvalues_consec <- arrayIP(1, vapply(temp, length, 0L), temp)
        apply_dim <- setdiff(names(dimnames(pvalues_perm)), c("time", "perm"))
        for (alpha in sort(pcrit, TRUE)) {
            consec_limit <- apply(pvalues_perm < alpha, 
                                       apply_dim, 
                                       consecLimit, alpha = alpha)
            temp <- matrixRle(
                array2mat(pvalues < alpha, "time", 
                          return_attributes = FALSE,
                          keep_dimnames = FALSE))
            ind <- temp$lengths < consec_limit[temp$matrixcolumn]
            temp$values[ind] <- 0
            temp <- inverse.matrixRle(temp)
            pvalues_consec[temp > 0] <- alpha
        }
        pvalues_consec <- aperm(pvalues_consec, names(dimnames(pvalues))) 
        ## number of sign. time points criterion (deprecated)
        # pvalues_global <- fnDims(pvalues_perm < 0.05, "time", colSums, 
        #                          vectorized = TRUE)
        # tempfn <- function(x) colSums(sweep(x, 2, x[1, ], ">="))
        # pvalues_global <- fnDims(pvalues_global, "perm", tempfn, 
        #                          vectorized = TRUE)
        # pvalues_global <- pvalues_global / (nperm + 1)
        ## return
        out <- c(out, list(perm_pvalues = pvalues,
                           perm_pvalues_consec = pvalues_consec))
    } 
    # return
    out
}


# Peak Anova functions =========== 

#' ANOVA on individual peaks
#' 
#' \code{peakAnova} performs automatic peak detection and runs traditional 
#' ANOVA on them
#' @param arraydat a numeric array with named dimnames containing the EEG (or 
#' other) data. Missing values are not allowed.
#' @param factordef a named list of factor definitions, containing the following 
#' elements:
#' \itemize{
#' \item{"between"}{"character vector of between-subject factors (default: NULL)"}
#' \item{"within"}{"character vector of within-subject factors (default: NULL)"}
#' \item{"w_id"}{"name of the dimension which identifies the subjects 
#' (default: "id")}
#' }
#' @param peakdef a named list of peak definitions. A peak definition is a 
#' three-element integer vector which specifies the time window (start and 
#' end point) and the polarity (1 for positivity and -1 for negativity).
#' @param bwdat a data.frame which contains the identification codes 
#' (factordef$w_id) and the group memberships (factordef$between) of the 
#' subjects. Missing values are not allowed.
#' @param verbose logical value indicating whether means for each factor level 
#' should be also returned
#' @param avg_around_peak integer value which specifies the number of data points
#' which are considered +- in the "mean around the peak" measure
#' @param nperm integer value giving the number of permutations (default: 999L)
#' @param useparallel logical value; if TRUE (default), computations are done
#' in parallel
#' @param ncores integer value corresponding to the number of cores; 
#' if NULL (default), it is set to the maximum number of cores available
#' @param par_method parallelization method; can be "snow" (default) or "mc"
#' (multicore)
#' @param cl a cluster definition for snow-type parallelization; if NULL 
#' (default), it is set up automatically
#' @param iaterms_last logival variable whether all interaction terms should 
#' follow all main effect terms (default: TRUE)
#' @param seed an integer value which specifies a seed (default: NULL)
#' @details The function assumes that the input array contains at least two 
#' named dimensions: chan (corresponding to the channels [electrodes]) and time 
#' (corresponding to time points). All dimensions which are not listed as 
#' within-subject factors are treated in a similar way as chan and time, that is 
#' separate TANOVA-s are computed for each level of those dimensions.
#' @note The function computes type I p-values - this is correct if the design
#' is fully balanced and orthogonal (if the number of between-subject 
#' factors is one, it may have slightly unequal group sizes).
#' @export
#' @return A list object with effect statistics, uncorrected p-values, and 
#' two types of corrected p-values: length correction or global correction
#' @references Koenig (2011) RAGU
peakAnova <- function(arraydat, factordef, peakdef, bwdat = NULL,
                      verbose = TRUE, avg_around_peak = 2, 
                      nperm = 1, useparallel = FALSE, ncores = NULL,
                      par_method = c("snow", "mc"), cl = NULL,
                      iaterms_last = TRUE, seed = NULL) {
    # some checks
    stopifnot(is.array(arraydat) | missing(arraydat))
    stopifnot(is.list(factordef) | missing(factordef))
    stopifnot(length(peakdef) != 2 | missing(peakdef))
    nperm <- as.integer(nperm)
    if (!is.null(factordef$between) && is.null(bwdat)) {
        stop("No between-participant data provided")
    }
    if (useparallel) {
        if (is.null(ncores)) ncores <- parallel::detectCores()
    }
    #
    out <- list(call = match.call())
    #
    sumSq <- function(arraydat, f_dat, aov_form, labels = FALSE) {
        out <- summary(aov(as.formula(aov_form), data = f_dat))
        if (!labels) {
            out <- matrixIP(unlist(lapply(out, "[", "Mean Sq"), 
                                   use.names = FALSE),
                            nrow(out[[1]]), ncol(arraydat))
        } else {
            out <- matrixIP(unlist(lapply(out, "[", "Mean Sq"), 
                                   use.names = FALSE),
                            nrow(out[[1]]), ncol(arraydat),
                            list(term = gsub(" ","",rownames(out[[1]])), 
                                 peak = seq_along(out)))
        }         
        return(out)
    }
    # 
    par_method <- match.arg(par_method)
    # prepare data
    temp <- preAnova(arraydat, factordef, bwdat, useparallel, ncores, par_method)
    # assign variables to this environment, and potentially overwrite existing ones
    #assignList(temp, verbose = FALSE)
    dat <- temp$dat
    arraydat <- temp$arraydat
    factordef <- temp$factordef
    origdimnames <- temp$origdimnames
    par_params <- temp$par_params
    rm(temp)
    #
    timedim <- which(names(dimnames(arraydat)) == "time")
    origtimedim <- which(names(origdimnames) == "time")
    origdimnames_peaks <- origdimnames
    arraydat_peaks <- subsetArray(arraydat, 
                                  list(time = dimnames(arraydat)$time[seq_along(peakdef)]))
    origdimnames_peaks$time <- 
        dimnames(arraydat_peaks)$time <- seq_along(peakdef)
    names(origdimnames_peaks)[origtimedim] <- 
        names(dimnames(arraydat_peaks))[timedim] <- "peak"
    # formula of anova
    aov_formula_char <- 
        if (is.null(factordef$within)) {
            paste(
                as.character(quote(arraydat)), " ~ ", 
                paste(as.character(factordef$between), collapse = "*"),
                sep = "")
        } else if (is.null(factordef$between)) {
            paste(
                as.character(quote(arraydat)), " ~ ", 
                paste(as.character(factordef$within), collapse = "*"),
                sep = "")
        } else {
            paste(
                as.character(quote(arraydat)), " ~ ", 
                paste(
                    paste(as.character(factordef$between), collapse = "*"), 
                    paste(as.character(factordef$within), collapse = "*"),
                    sep = "*"),
                sep = "")
        }
    aov_formula <- as.formula(aov_formula_char)
    mean_formula <- as.formula(gsub("\\*", ":", aov_formula_char))
    # compute marginal means
    origdims <- vapply(origdimnames, length, 0L)
    modeldims <- c(factordef$w_id, factordef$within)
    keepdims <- setdiff(names(origdimnames), modeldims)
    # find indices corresponding to time windows
    tpoints <- as.numeric(origdimnames$time)
    peakdef <- lapply(peakdef, function(x) 
        c(which(origdimnames$time %in% as.character(x[1:2])), x[3]))
    # find peaks
    fmeans <- marginalMeans(mean_formula, dat, arraydat, 
                            origdimnames[keepdims], 
                            keep_term_order = !iaterms_last, residualmean = FALSE)
    datrowind <- match(
        interaction(dat[, strsplit(names(fmeans), ":")[[1]] ]), 
        rownames(fmeans[[1]]) )
    peakind_facs <- strsplit(rownames(fmeans[[1]]), "\\.")
    peakind_facs <- setNames(data.frame(
        matrix(unlist(peakind_facs, use.names = FALSE),
               length(peakind_facs), length(peakind_facs[[1]]), TRUE)),
        strsplit(names(fmeans), ":")[[1]])
    peakind <- sapply(peakdef, function(x) {
        minmax <- if (x[3] < 0) which.min else which.max
        x[1] -1 + apply(fmeans[[1]][,seq(x[1],x[2])], 1, minmax)})
    #
    # peak amplitudes
    arraydat_peaks[] <- sapply(1:ncol(peakind), function(i) {
        out <- matrixIP(0, nrow(arraydat), avg_around_peak*2 + 1)
        avgind <- outer(peakind[,i], c(-avg_around_peak:avg_around_peak), "+")
        out[] <- arraydat[
            cbind(seq_along(datrowind), as.integer(avgind[datrowind,]))]
        return(rowMeans(out))
    })
    #
    # Compute ANOVA on peak amplitudes
    Fvals_obs <- arrayAnovaSub(arraydat_peaks, 
                               factordef, origdimnames_peaks, dat, verbose)
    out <- c(out, list(F_obs = Fvals_obs))
    # if verbose, save factor means
    if (verbose) {
        factor_means <- data.frame(
            group = peakind_facs[rep(1:nrow(peakind_facs), 
                                     length(peakdef)), , drop = F],
            peak = factor_(rep(seq_along(peakdef), each = nrow(fmeans[[1]]))),
            ampl = c(marginalMeans(mean_formula, dat, arraydat_peaks, 
                                   origdimnames_peaks["peak"], 
                                   keep_term_order = !iaterms_last, 
                                   residualmean = FALSE)[[1]]),
            lat = tpoints[c(peakind)])
        out <- c(out, list(factor_means = factor_means))
    }
    #
    # Compute permutation statistics for peak latencies
    #
    # observed sum of squares
    lateff_obs <- sumSq(peakind, peakind_facs, aov_formula_char, labels = TRUE)
    # permutations
    if (nperm > 1L) {
        #
        permfn <- function(i) {
            facmeans <- marginalMeans(mean_formula, dat[randind[i,], ], 
                                      arraydat, 
                                      origdimnames[keepdims], 
                                      keep_term_order = !iaterms_last, 
                                      residualmean = FALSE)[[1]]
            peaki <- sapply(peakdef, function(x) {
                minmax <- if (x[3] < 0) which.min else which.max
                x[1] - 1 + apply(facmeans[, seq(x[1], x[2])], 1, minmax)})
            out <- sumSq(peaki, peakind_facs, aov_formula_char)
            return(out)
        }
        #
        # generate random orders (dim(randind) = nperm X nrow(dat))
        if (!is.null(seed)) set.seed(seed)
        randind <- mclapply(1:nperm, function(i) 
            shuffle(nrow(dat),  
                    how(within = Within(type = "free"), 
                        plots = Plots(strata = dat[,factordef$w_id], 
                                      type = "free"))),
            mc.cores = par_params$ncores)
        randind <- matrix(unlist(randind, use.names = FALSE), 
                          nperm, nrow(dat), TRUE)
        #
        if (!useparallel) {
            lateff_perm <- lapply(1:nperm, permfn)
        } else if (par_method == "snow") {
            if (is.null(cl)) cl <- makePSOCKcluster(par_params$ncores, outfile = "")
            clusterExport(cl, 
                          varlist = par_params$varlist2snow,
                          envir = environment())
            lateff_perm <- parLapply(cl, 1:nperm, permfn)
            #             stopCluster(cl)
            #             rm(cl)
        } else if (par_method == "mc") {
            lateff_perm <- mclapply(1:nperm, permfn, mc.cores = par_params$ncores)
        } 
        lateff_perm <- cbind(c(lateff_obs), 
                             matrixIP(unlist(lateff_perm, use.names = FALSE),
                                      ncol = nperm))
        lateff_perm <- sweep(lateff_perm[, -1, drop = F], 1, 
                             lateff_perm[, 1], "-")
        pvalues <- arrayIP(
            (rowSums(lateff_perm >= 0) + 1) / (nperm + 1),
            dim(lateff_obs), dimnames(lateff_obs))
        out <- c(out, list(lat_pvalues = pvalues))
    } 
    # return
    out
}

#' Correct p-values in time-series
#' 
#' \code{pvalueConsec} sets a minimal consecutive length criterion on the 
#' p-values in a time series
#' @param dat a numeric matrix or array containing p-values, or a list of such 
#' matrices/arrays. dat must have a named time dimension.
#' @param sig_level numeric value (default: 0.05), the level of significance
#' @param min_length numeric value (default: 10), the minimum number of 
#' consecutive significant time points
#' @export
#' @return An object of the same dimension as the input data
pvalueConsec <- function(dat, sig_level = 0.05, min_length = 10) {
    pCorr <- function(x) {
        temp <- rle(x)
        ind <- (temp$values == 1) & (temp$lengths < min_length)
        temp$values[ind] <- 0
        return( inverse.rle(temp) )
    }
    pCorrMat <- function(datarr) {
        datarr <- array2mat(datarr, "time")
        datarr[apply(datarr <= sig_level, 2, pCorr) == 0] <- 1
        datarr <- mat2array(datarr)
        return(datarr)
    }
    if (is.list(dat)) {
        out <- lapply(dat, pCorrMat)
    } else {
        out <- pCorrMat(dat)
    }
    # return
    out
}