#
# <<< ANOVA on arrays >>> ----
#

#' TFCE-correction of the test statistic
#' 
#' \code{anovaTfce} is called internally by \code{arrayAnova} and 
#' \code{arrayTtest} and performs TFCE-correction of the given test statistic
#' @param x numeric matrix or array of test statistic 
#' (chan_time X other_dimensions or chan X time X other_dimenions)
#' @param target_dim the dimension(s) which hold chan and time values (usually 
#' 1L if x is a matrix and 1:2L if x is an array)
#' @param has_neg logical scalar whether x might have negative values (TRUE for
#' a t-statistic and FALSE for an F-statistic)
#' @param nr_chan number of channels
#' @param nr_time number of time points
#' @param tfce a list of TFCE parameters
#' @keywords internal
anovaTfce <- function(x, target_dim, has_neg, nr_chan, nr_time, tfce) {
    fnDims(x, target_dim, 
           function(y) tfceFn(
               matrix_(y, nr_chan, nr_time, arg_check = FALSE), 
               chn = tfce$ChN, eh = tfce$EH,
               nr_steps = tfce$steps, has_neg = has_neg),
           arg_check = FALSE, parallel = FALSE)
}

#' Permutation of the test statistic
#' 
#' \code{permPvalues} is called internally by \code{arrayAnova} and 
#' \code{arrayTtest} and performs permutation analysis.
#' @param obj a list returned by \code{preAnova}
#' @param nperm an integer value indicating the number of permutations
#' @param seed an integer value which specifies a seed (default: NULL), or a 
#' list of arguments passed to \code{\link{set.seed}}
#' @param stat_obs a numeric array of the observed test statistic
#' @param tfce NULL (default) or a list of the TFCE parameters 
#' @return \code{\link{permPvalues}} returns a logical array whose elements 
#' indicate whether the absolute observed test statistic in the given chanXtime 
#' cell is smaller than the absolute randomized test statistic
#' @keywords internal
permPvalues <- function(obj, nperm, seed, stat_obs, tfce = NULL) {
    # helper function
    permfn <- function(obj, perm_vec, tfce, 
                       stat_fun, abs_stat_obs, has_neg) {
        stat_perm <- stat_fun(obj, new = as.vector(perm_vec))
        setattr(stat_perm, "dim", 
                c(obj$nr_chanXtime, obj$otherdims$size))
        if (!is.null(tfce)) 
            stat_perm <- anovaTfce(stat_perm, 1L, has_neg, 
                                   obj$nr_chan, obj$nr_time, tfce)
        if (has_neg) stat_perm <- abs(stat_perm)
        stat_perm <- colMaxs(stat_perm)
        # return
        if (length(stat_perm) > 1L) {
            sweep(abs_stat_obs, obj$otherdims$index, stat_perm, "<")
        } else {
            abs_stat_obs < stat_perm
        }
    }
    # TFCE yes / no
    use_tfce <- if (is.null(tfce)) FALSE else TRUE
    # type: ANOVA / t-test
    if (grepl("between|within|mixed", obj$type)) {
        rand_sample <- anovaRandomIndices(obj, nperm, seed)
        stat_fun <- compF
        has_neg <- FALSE
    } else if (grepl("independent", obj$type)) {
        rand_sample <- isTtestRandomGroups(obj, nperm, seed)
        stat_fun <- isTtest
        has_neg <- TRUE
    } else if (grepl("one|paired", obj$type)) {
        rand_sample <- osTtestRandomSigns(obj, nperm, seed)
        stat_fun <- osTtest
        has_neg <- TRUE
    } else {
        stop("Wrong 'type' node in 'obj'")
    }
    # take the absolute values of the observed statistic
    if (has_neg) stat_obs <- abs(stat_obs)
    # create larger chunks
    chunks <- min(10L, ceiling(nperm/max(1L, getDoParWorkers())))
    # avoid code-style warning
    xi <- NULL; rm(xi)
    # compute sig
    sig <- 
        foreach(x = iter(rand_sample, by = "col", chunksize = chunks), 
                .combine = "+", .inorder = FALSE) %dopar% 
            {
                foreach(xi = iter(x, by = "col", chunksize = 1L), 
                        .combine = "+", .inorder = FALSE) %do% 
                    permfn(obj, xi, tfce, stat_fun, stat_obs, has_neg)
            }
    sig <- (sig + 1) / (nperm + 1)
    # set type
    attr_type <- 
        if (use_tfce) {
            "Permuted P-value (with TFCE correction)"
        } else {
            "Permuted P-value (with max(Channel X Time) correction)"
        }
    setattr(sig, "type", attr_type)
}

#' Compute marginal means in an ANOVA design
#' 
#' \code{marginalMeans} computes marginal means in ANOVA designs
#' @param form formula of the model
#' @param f_dat data.frame of factors
#' @param a_dat matrix which contains the dependent variables. It must have as
#' many rows as f_dat. Usually a_dat is the result of calling 
#' \code{\link{array2mat}} on the orginal data array.
#' @param dimn the original dimension names for the array corresponding to the
#' column names of \code{a_dat}
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
#' @keywords internal
#' @return A named list containing the marginal means for each model term
# TODO: examples, argument checking, handle numeric variable
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
        marg_means[[i]] <- array_(rowsum(a_dat, groups)/groupfreq, 
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

#' Compute adjusted or unadjusted cell means in ANOVA designs
#' 
#' \code{modelMeans} is a generic function for computing adjusted or unadjusted 
#' cell means in ANOVA designs.
#' @param model an object returned by \code{\link{arrayAnova}} or 
#' \code{\link{tanova}}
#' @param .arraydat a numeric array with named dimnames containing the EEG (or 
#' other) data. Missing values are not allowed.
#' @param factordef a named list of factor definitions, containing the following 
#' elements:
#' \itemize{
#' \item{between: }{character vector of between-subject factors (default: NULL)}
#' \item{within: }{character vector of within-subject factors (default: NULL)}
#' \item{w_id: }{name of the dimension which identifies the subjects 
#' (default: "id")}
#' }
#' @param bwdat a data.frame which contains the identification codes 
#' (factordef$w_id) and all subject-level variables (usually factors) listed in
#' 'factordef$between'. Missing values are not allowed.
#' @param term a character vector of model terms; a separate array of cell means
#' is returned for each term. The default is NULL, in which case all terms
#' (that is, the means corresponding to all main effects and interactions) are
#' returned. For custom terms, note that interaction terms must be given as 
#' "factorA*factorB" and not as "factorA:factorB".
#' @param adjusted a logical value if the cell means should be adjusted for 
#' unbalanced data (default: FALSE). See Details.
#' @param ... not used yet
#' @details If 'adjusted' is set to FALSE (the default), \code{modelMeans} fits 
#' separate models for each model term and computes the model-based predictions 
#' for a reference grid which contains the mean and the 1SD value for continuous
#' covariates, and all levels of the factor variables in the design. If 
#' 'adjusted' is set to TRUE, the means are simply the marginal means of the
#' predicted means of the highest-order interaction term. This is equivalent to 
#' the so-called LS means (least-squares means) in the SAS terminology. 
#' This distinction is only relevant if the design is unbalanced, that is, the 
#' sample sizes in a between-subject design are unequal (note that the input 
#' array may not contain missing values, therefore the within subject part of 
#' the design is always balanced). For such datasets the lower order 
#' interaction or main effect means can substantially differ from the averages 
#' of the corresponding higher-order means if 'adjusted' is FALSE. 
#' @export
#' @return The function returns a named list of arrays; the names are the model 
#' terms.
#' @examples 
#' # example data
#' data(erps)
#' dat_id <- attr(erps, "id") # to get group memberships
#' 
#' # make the dataset unbalanced to illustrate the difference between adjusted
#' # and unadjusted means
#' erps <- subsetArray(erps, list(id = 2:20))
#' dat_id <- dat_id[2:20, ]
#' 
#' # compute the means for the group*stimclass*pairtype interaction and all
#' # other lower-order terms
#' fdef <- list(between = "group", 
#'              within = c("stimclass", "pairtype"))
#' means <- modelMeans(erps, fdef, bwdat = dat_id)
#' 
#' # the same for the least-squares means
#' means_adj <- modelMeans(erps, fdef, bwdat = dat_id, adjusted = TRUE)
#' 
#' # this is an unbalanced dataset, so the means for the highest order
#' # interaction are identical, but the marginal means are not
#' stopifnot(identical(
#'     TRUE,
#'     all.equal(means$`group*stimclass*pairtype`,
#'               means_adj$`group*stimclass*pairtype`))
#' )
#' stopifnot(!identical(
#'     TRUE,
#'     all.equal(means$`stimclass*pairtype`,
#'               means_adj$`stimclass*pairtype`))
#' )
#' 
#' # simple means for the Fz channel at time 200, for identical pairs in the
#' # "A" stimulus class, separately for the two reading groups
#' simple_means <- sapply(c("control", "dl"), function(g)
#'     mean(subsetArray(erps, 
#'                      list(chan = "Fz", time = "200", 
#'                           stimclass = "A", pairtype = "ident",
#'                           id = dat_id$group == g)))
#' )
#' stopifnot(identical(
#'     TRUE,
#'     all.equal(means$`group*stimclass*pairtype`["Fz", "200", , "A", "ident"], 
#'               simple_means))
#' )
#' 
#' # the same ignoring the groups
#' simple_mean_nogroup <- mean(subsetArray(erps, 
#'                                         list(chan = "Fz", time = "200", 
#'                                         stimclass = "A", pairtype = "ident")))
#' stopifnot(identical(
#'     TRUE,
#'     all.equal(means$`stimclass*pairtype`["Fz", "200", "A", "ident"], 
#'               simple_mean_nogroup))
#' )
#' stopifnot(identical(
#'     TRUE,
#'     all.equal(means_adj$`stimclass*pairtype`["Fz", "200", "A", "ident"], 
#'               mean(simple_means)))
#' )
#' 
modelMeans <- function(...) UseMethod("modelMeans")

#' @export
#' @describeIn modelMeans
modelMeans.default <- function(.arraydat, factordef, 
                               bwdat = NULL, term = NULL, 
                               adjusted = FALSE, ...) {
    # workhorse function
    compute <- function(tm, dat, y, factordef, ydims, 
                        tm0 = NULL, model0 = NULL) {
        vars <- strsplit(tm, "\\*")[[1L]]
        if (is.null(tm0)) {
            form <- as.formula(paste("~", tm, sep = ""))
            mod <- preSumSq(form, dat, y = y, between = character(), 
                            adjusted = FALSE)
            coeff <- mod$coeff
            pred_levels <- lapply(dat[vars], function(x) {
                if (is.numeric(x)) c(0, 1) else levels(x)
            })
            pred_grid <- expand.grid(pred_levels, KEEP.OUT.ATTRS = FALSE)
            pred_mm <- model.matrix(form, pred_grid) 
            pred <- pred_mm %*% coeff[colnames(pred_mm), , drop = FALSE]
            setattr(pred, "dim", c(vapply(pred_levels, length, 0L), ydims$dim))
            setattr(pred, "dimnames", c(pred_levels, ydims$dimn))
            pred <- aperm(pred, c(ydims$orig_dimid, vars))
        } else {
            vars0 <- strsplit(tm0, "\\*")[[1L]]
            pred <- 
                if (identical(vars, vars0)) {
                    model0
                } else {
                    avgDims(model0, setdiff(vars0, vars))
                }
        }
        # return
        pred
    }
    # set contrasts
    opcons <- options("contrasts")
    options(contrasts = c("contr.helmert", "contr.poly"))
    on.exit(options(opcons))
    #
    input <- preAnova(
        .arraydat = .arraydat, factordef = factordef, bwdat = bwdat, 
        verbose = FALSE, tfce = FALSE, perm = permParams(n = 0L))
    dat <- input$dat
    y <- input$.arraydat
    if (is.null(term)) {
        term <- gsub(":", "*", input$full_dimnames$modelterm)
    } 
    #
    origpos <- attr(input$teststat_dimid, "origpos")
    names(origpos) <- input$teststat_dimid
    ydimid <- setdiff(input$teststat_dimid, "modelterm")
    orig_ydimid <- ydimid[order(origpos[ydimid])]
    ydimn <- input$full_dimnames[ydimid]
    ydim <- vapply(ydimn, length, 0L)
    y_dims <- list(dim = ydim, dimn = ydimn, orig_dimid = orig_ydimid)
    #
    out <- 
        if (!adjusted) {
            lapply(term, compute, dat = dat, y = y, factordef = factordef,
                   ydims = y_dims)
        } else {
            term0 <- input$full_dimnames$modelterm
            term0 <- term0[length(term0)]
            term0 <- gsub(":", "*", term0)
            model0 <- compute(term0, dat = dat, y = y, factordef = factordef,
                              ydims = y_dims)
            lapply(term, compute, dat = dat, y = y, factordef = factordef,
                   ydims = y_dims, tm0 = term0, model0 = model0)
        }
    # return
    setattr(out, "names", term)
    out
}

#' @export
#' @describeIn modelMeans
modelMeans.arrayAnova <- function(model = NULL, term = NULL, adjusted = FALSE,
                                  ...) {
    mcall <- model$call
    mcall[[1]] <- as.name("modelMeans.default")
    mcall$term <- term
    mcall$adjusted <- adjusted
    eval(mcall)
}

#' @export
#' @describeIn modelMeans
modelMeans.tanova <- modelMeans.arrayAnova


#' Force orthogonality in a model matrix
#' 
#' \code{forceOrthogonal} checks if a model matrix is orthogonal and makes
#' adjustments if it is not.
#' @param model_matrix the model matrix
#' @param between the names of between-subject factors
#' @keywords internal
forceOrthogonal <- function(model_matrix, between = character()) {
    if (ncol(model_matrix) > 2L) {
        orth <- cor(model_matrix[, -1L])
        orth <- orth[lower.tri(orth)]
        if (max(abs(orth)) > .Machine$double.eps) {
            if (length(between) > 1L)
                warning("The design is unbalanced, and this function calculates 
sequential (type-I) tests. Consider re-running the test
with reordered between-subject factors (see ?arrayAnova).")
            model_matrix[, 2L] <- model_matrix[, 2L] - mean(model_matrix[, 2L])
            for (i in 3:ncol(model_matrix)) {
                mmi <- model_matrix[,1:(i-1L)]
                coeff <- solve(crossprod(mmi), 
                               crossprod(mmi, model_matrix[, i, drop=FALSE]))
                model_matrix[,i] <- model_matrix[,i]- mmi %*% coeff
            }    
        }
    }
    # return
    model_matrix
}

#' Prepare sum-of-squares computations
#' 
#' \code{preSumSq} prepares the objects which are needed by 
#' \code{\link{compSumSq}}. Thus, it returns the model.matrix (X), the X'X 
#' matrix, and all unique rows of Xs, and also the residuals if residuals are
#' permuted instead of raw observations.
#' @param form model formula (or character string to \code{as.formula})
#' @param form_data a data.frame containing all variables of the formula
#' @param keep_term_order a logical value indicating whether the terms should 
#' keep their positions. If FALSE the terms are reordered so that main effects 
#' come first, followed by the interactions, all second-order, all third-order 
#' and so on. Effects of a given order are kept in the order specified.
#' @param scaled a logical value indicating if numeric variables in
#' 'form_data' are already scaled (TRUE) or not (FALSE, the default).
#' @param y a numeric matrix of response data if residuals should be computed
#' (usually coming from \code{preAnova})
#' @param between the names of between-subject factors. Only used for warning
#' message in case of unbalanced data.
#' @param adjusted adjustment for unbalanced data (default: TRUE)
#' @return A list with the following components: model_matrix, 
#' model_matrix_labels, xpx, unique_contrasts, residuals (optional).
#' @keywords internal
preSumSq <- function(form, form_data, 
                     keep_term_order = FALSE,
                     scaled = FALSE, y = NULL, between = NULL,
                     adjusted = TRUE) {
    term_object <- terms(as.formula(form), keep.order = keep_term_order)
    term_labels <- attr(term_object, "term.labels")
    term_factors <- attr(term_object, "factors")
    vars <- rownames(term_factors)
    if (!scaled) {
        scfn <- function(x) {
            if (is.numeric(x) && is.null(attr(x, "scaled:scale"))) {
                scale(x) 
            } else {
                x
            }
        }
        form_data[vars] <- lapply(form_data[vars], scfn)
    }
    mm <- model.matrix(term_object, data = form_data)
    # check orthogonality
    if (adjusted) mm <- forceOrthogonal(mm, between)
    #
    xpx <- crossprod(mm)
    mm_labels <- c("(Intercept)", term_labels[attr(mm, "assign")])
    unique_contrasts <- 
        lapply(term_labels, function(term) {
            ind <- which(mm_labels == term)
            fastUnique(mm[, ind, drop = FALSE],
                       freq = TRUE)
        })
    setattr(unique_contrasts, "names", term_labels)
    mm_labels <- colnames(mm)
    setattr(mm, "dimnames", NULL)
    setattr(mm, "assign", NULL)
    setattr(mm, "contrasts", NULL)
    #
    out <- list(
        model_matrix = mm, model_matrix_labels = mm_labels, xpx = xpx, 
        unique_contrasts = unique_contrasts)
    #
    if (!is.null(y)) {
        out$coeff <- solve(xpx, crossprod(mm, y))
        pred <- mm %*% out$coeff
        out$residuals <- y - pred
    }
    # return
    out
}

#' Compute sum of squares
#' 
#' \code{compSumSq} computes sum of squares which are needed by 
#' \code{\link{compF}} and \code{\link{permPvalues}}. 
#' @param obj a list object as returned by \code{\link{preAnova}}
#' @param new_indices if not NULL (default), it must be a vector of numeric
#' indices to be passed to the model matrix. This parameter is used for 
#' permutation-based analyses.
#' @return a numeric matrix of the sum-of-squares with an extra attribute for 
#' the degrees of freedom ('Df'), or a list of such matrices
#' @keywords internal
compSumSq <- function(obj, new_indices = NULL) {
    ssqFn <- function(x, y, r, new_indices) {
        mm <- x$model_matrix
        if (!is.null(new_indices)) {
            if (!is.null(r)) {
                y <- r[new_indices, ]
            } else {
                mm <- mm[new_indices, ]
            }
        } 
        mm_labels <- x$model_matrix_labels
        xpx <- x$xpx
        unc <- x$unique_contrasts
        ssq <- matrix_(0, ncol(y), length(unc),
                       dimnames = list(NULL, modelterm = names(unc)))
        df <- rep.int(0L, length(unc))
        setattr(df, "names", names(unc))
        coeff <- solve(xpx, crossprod(mm, y))
        setattr(coeff, "dimnames", list(mm_labels, NULL)) 
        for (i in seq_along(unc)) {
            freq <- attr(unc[[i]], "freq")
            fmeans <- unc[[i]] %*% coeff[colnames(unc[[i]]), , drop = FALSE]   
            ssq[, i] <- colSums(fmeans^2 * freq)
            df[i] <- ncol(unc[[i]])
        }
        setattr(ssq, "Df", df)
        #
        ssq
    }
    #
    input <- obj$pre_sumsq
    y <- obj$.arraydat
    residuals <- input[[1]]$residuals
    # return
    lapply(input, ssqFn, y = y, r = residuals, new_indices = new_indices)
}

#' Random permutations for \code{\link{arrayAnova}} and \code{\link{tanova}}
#' 
#' \code{anovaRandomIndices} creates a matrix of \code{nperm} random
#' permutations 
#' @param input a list returned by \code{\link{preAnova}}
#' @param nperm integer value; the number of permutations
#' @param seed an integer value which specifies a seed (default: NULL), or a 
#' list of arguments passed to \code{\link{set.seed}}
#' @return A matrix with as many rows as the design matrix and 'nperm' columns
#' @seealso \code{\link{arrayAnova}} and \code{\link{tanova}}
#' @keywords internal
anovaRandomIndices <- function(input, nperm, seed) {
    setSeed(seed)
    out <- shuffleSet(
        nrow(input$dat), nperm, 
        how(within = Within(type = "free"), 
            plots = Plots(strata = input$dat[[input$factordef$w_id]], 
                          type = "free")))
    # return
    t(out)
}

#' Compute F-values
#' 
#' \code{compF} computes F-values. This is an internal function called by
#' the high-level function \code{\link{arrayAnova}}.
#' @param obj a list object as returned by \code{\link{preAnova}}
#' @param new_indices if not NULL (default), it must be a vector of numeric
#' indices to be passed to the model matrix. This parameter is used for 
#' permutation-based analyses.
#' @param attribs a logical value if dimension and verbose attributes should
#' be attached to the resulting array (default: FALSE)
#' @param verbose a logical value if p-values, effect sizes, and degrees of
#' freedom should be added to the returned array (as 'p_value', 'effect_size',
#' and 'Df' attributes)
#' @return This function returns a numeric matrix of F-values if 'attribs' is 
#' FALSE, and an array of F-values if 'attribs' is TRUE. If the verbose variant 
#' was requested, the array has the following extra attributes: 'p_value' and 
#' 'effect_size', having the same dimension as the array of F-values, and
#' an integer matrix named 'Df' with term- and residual degrees of freedom for
#' each model term.
#' @keywords internal
compF <- function(obj, new_indices = NULL, attribs = FALSE, verbose = FALSE) {
    if (!attribs) verbose <- FALSE
    type <- obj$type
    f_def <- obj$factordef
    ssq <- compSumSq(obj, new_indices = new_indices)
    Df <- lapply(ssq, attr, "Df")
    if (type == "between") {
        model_ssq <- ssq[[1]]
        model_df <- Df[[1]]
        resid_ssq <- obj$pre_sumsq[[1]]$ssq_total - rowSums(model_ssq)
        resid_df <- attr(resid_ssq, "Df") - sum(model_df)
        Fvals <- sweepMatrix(model_ssq, 2, model_df/resid_df, "/", 
                             has_NA = FALSE)
        Fvals <- Fvals/resid_ssq
        if (verbose) {
            pvalues <- t(pf(t(Fvals), model_df, resid_df, lower.tail = FALSE))
            if (!is.null(f_def$observed)) {
                ind <- grepl(paste(f_def$observed, collapse = "|"), 
                             colnames(model_ssq))
                SS1 <- rowSums(model_ssq[, ind, drop = FALSE])
                SS2 <- model_ssq; SS2[, !ind] <- 0
                ges <- model_ssq / (model_ssq - SS2 + resid_ssq + SS1)
            } else {
                ges <- model_ssq / (model_ssq + resid_ssq)
            }
        }
    } else if (type == "within") {
        model_ssq <- ssq[[2]]
        model_df <- Df[[2]]
        modnames <- dimnames(model_ssq)$modelterm
        m_ind <- !grepl(f_def$w_id, modnames)
        err_ind <- match(
            paste(f_def$w_id, ":", modnames[m_ind], sep = ""),
            modnames)
        resid_ssq <- model_ssq[, err_ind, drop = FALSE]
        resid_df <- model_df[err_ind]
        model_ssq <- model_ssq[, m_ind, drop = FALSE]
        model_df <- model_df[m_ind]
        Fvals <- (model_ssq/resid_ssq) * (resid_df/model_df)
        if (verbose) {
            pvalues <- t(pf(t(Fvals), model_df, resid_df, lower.tail = FALSE))
            SSr <- rowSums(resid_ssq)
            ges <- model_ssq / (model_ssq + SSr)
        }
    } else if (type == "mixed") {
        model_ssq <- ssq[[1]]
        model_df <- Df[[1]]
        bwind <- grepl(paste(f_def$between, collapse = "|"),
                       colnames(model_ssq))
        residnames <- paste(f_def$w_id, 
                            gsub("^:|:$", "",
                                 gsub(paste(f_def$between, collapse = "|"), 
                                      "", colnames(model_ssq))), 
                            sep = ":")
        residnames <- gsub(":+", ":", gsub("^:|:$", "", residnames))
        resid_ssq <- ssq[[2]][, unique(residnames), drop = FALSE]
        resid_df <- Df[[2]][colnames(resid_ssq)]
        for (i in colnames(resid_ssq)) {
            ind <- which(residnames == i & bwind)
            resid_ssq[, i] <- resid_ssq[, i] - 
                rowSums(model_ssq[, ind, drop = FALSE])
            resid_df[i] <- resid_df[i] - sum(model_df[ind])
        }
        resid_ssq <- resid_ssq[, residnames, drop = FALSE]
        resid_df <- resid_df[residnames]
        Fvals <- model_ssq / resid_ssq
        Fvals <- sweepMatrix(Fvals, 2, resid_df/model_df, "*", has_NA = FALSE)
        if (verbose) {
            pvalues <- t(pf(t(Fvals), model_df, resid_df, lower.tail = FALSE))
            SSr <- rowSums(resid_ssq[, !duplicated(residnames)])
            if (!is.null(f_def$observed)) {
                ind <- grepl(paste(f_def$observed, collapse = "|"), 
                             colnames(model_ssq))
                SS1 <- rowSums(model_ssq[, ind, drop = FALSE])
                SS2 <- model_ssq; SS2[, !ind] <- 0
                ges <- model_ssq / (model_ssq - SS2 + SSr + SS1)
            } else {
                ges <- model_ssq / (model_ssq + SSr)
            }
        }
    }
    #
    if (attribs) {
        setattributes(Fvals, obj, "Traditional F statistic")
        if (verbose) {
            setattributes(pvalues, obj, "Traditional P-value")
            setattributes(ges, obj, "Generalized Eta Squared")
            setattr(ges, "Df", NULL)
            Df <- cbind(Df_model = model_df, 
                        Df_residual = rep_len(resid_df, length(model_df)))
            names(dimnames(Df)) <- c("modelterm", "")
            setattr(Df, "type", "Degrees of freedom")
            setattr(Fvals, "p_value", pvalues)
            setattr(Fvals, "effect_size", ges)
            setattr(Fvals, "Df", Df)
        }
    }
    # return
    Fvals
}


#' Parameter checks and data preparation in \code{arrayAnova}
#' 
#' \code{preAnova} performs parameter checking, creates a data.frame for 
#' the ANOVA design and reshapes .arraydat accordingly, prepares the computation
#' of sum-of-squares. It is called by \code{\link{arrayAnova}} and 
#' \code{\link{tanova}}.
#' @param .arraydat numeric array
#' @param factordef a list of between-subject factors (named 'between'), 
#' within-subject factors ('within'), subject identifier ('w_id'), and a 
#' character vector of observed (i.e., not manipulated) variables ('observed').
#' The last element ('observed') is only used in effect size computations.
#' @param bwdat a data.frame with the subject-level variables (usually factors)
#' @param tfce logical value whether TFCE is requested
#' @param perm logical value whether permutation is requested
#' @param tanova_type type of tanova ("tanova", "dissimilarity", "gfp") if this 
#' function is called from \code{tanova}
#' @keywords internal
#' @seealso \code{\link{arrayAnova}} and \code{\link{tanova}} are the 
#' higher-level function exported to the user, \code{\link{preSumSq}}, 
#' \code{\link{compF}} and \code{\link{compTanovaEffect}} are internal 
#' functions relying on the returned object of \code{preAnova}.
preAnova <- function(.arraydat, factordef, bwdat, verbose, tfce, perm, 
                     tanova_type = NULL) {
    # helper function
    checkObligDimnames <- function(method, dimn, full_dimid) {
        print_dimn <- sub("$", "'", sub("^", "'", dimn))
        print_dim <- if (length(dimn) > 1L) "dimensions" else "dimension"
        print_fac <- if (length(dimn) > 1L) "factors" else "factor"
        print_allow <- if (length(dimn) > 1L) "are" else "is"
        if (!all(dimn %in% full_dimid)) {
            stop(sprintf("For %s, the input array must have %s %s.", 
                         method, print_dimn, print_dim))
        }
        if (any(dimn %in% factordef$within)) {
            stop(sprintf("For %s, %s %s not allowed as within-subject %s.", 
                         method, print_dimn, print_allow, print_fac))
        }
    }
    #
    # fast check of arguments
    #
    # check if .arraydat is an array
    if (is.data.frame(.arraydat)) .arraydat <- as.matrix(.arraydat)
    assertArray(.arraydat, mode = "numeric", min.d = 2L,
                .var.name = ".arraydat")
    has_NA <- anyNA(.arraydat)
    if (has_NA) {
        stop("Assertion on '.arraydat' failed: Contains missing values")
    }
    # check if factordef is a list
    assertList(factordef, types = "character", names = "unique",
               .var.name = "factordef")
    # check if w_id is provided
    if (is.null(factordef$w_id)) factordef$w_id <- "id"
    # save dimension names
    pre_dimnames <- dimnames(.arraydat)
    full_dimnames <- fillMissingDimnames(pre_dimnames,
                                         dim(.arraydat))
    full_dimid_origord <- names(full_dimnames)
    if (!identical(pre_dimnames, full_dimnames)) {
        dimnames(.arraydat) <- full_dimnames
    }
    full_dimid <- names(full_dimnames)
    # check dimension names
    if (!factordef$w_id %in% full_dimid) {
        stop(sprintf(".arraydat has no '%s' dimension identifier", 
                     factordef$w_id))
    }
    if (anyNA(full_dimnames) || anyDuplicated(full_dimnames)) {
        stop("Provide array with unique named dimensions")
    }
    # if TFCE is requested, "chan" and "time" dimensions are obligatory
    nr_chanXtime <- nr_chan <- nr_time <- NULL
    if (tfce) {
        checkObligDimnames("TFCE", c("chan", "time"), full_dimid)
        nr_chan <- length(full_dimnames$chan)
        nr_time <- length(full_dimnames$time)
        nr_chanXtime <- nr_chan * nr_time
    } else if (perm$n > 1L) {
        nr_chanXtime <- c(length(full_dimnames$chan),
                          length(full_dimnames$time))
        nr_chanXtime[nr_chanXtime == 0L] <- 1L
        nr_chanXtime <- prod(nr_chanXtime)
    } 
    # if tanova is requested, "chan" dimension is obligatory
    # create also a logical variable for tanova
    if (!is.null(tanova_type)) {
        checkObligDimnames(toupper(tanova_type), "chan", full_dimid)
        nr_chan <- length(full_dimnames$chan)
        tanova <- TRUE
    } else {
        tanova <- FALSE
    }
    # check type of anova
    type <- 
        if (!is.null(factordef$between)) {
            if (is.null(factordef$within)) {
                "between"
            } else {
                "mixed"
            }
        } else if (!is.null(factordef$within)) {
            "within"
        } else {
            stop("Provide at least one between- and/or within-subject factor")
        }
    # check between subject factors
    if (!is.null(factordef$between)) {
        #  check bwdat
        assertDataFrame(bwdat, .var.name = "bwdat")
        if (!all(c(factordef$between, factordef$w_id) %in% colnames(bwdat))) {
            stop("Between-subject factors and 'bwdat' dataset do not match.")
        }
        # force to data.frame, keep only relevant variables
        bwdat <- droplevels(
            as.data.frame(bwdat[, c(factordef$w_id, factordef$between)]))
        # check for missing values
        if (anyNA(bwdat)) stop("Missing values in bwdat are not allowed")
        # align .arraydat and bwdat on w_id
        id <- intersect(full_dimnames[[factordef$w_id]],
                        bwdat[, factordef$w_id])
        if (length(id) < length(full_dimnames[[factordef$w_id]])) {
            .arraydat <- subsetArray(.arraydat, 
                                     setNames(list(id), factordef$w_id))
            full_dimnames <- dimnames(.arraydat)
            bwdat <- bwdat[match(id, bwdat[[factordef$w_id]]), ]
        }
    }
    # check if .arraydat has proper dimension identifiers
    if (!is.null(factordef$within)) {
        if (!all(factordef$within %in% names(full_dimnames)))
            stop("Missing within-subject dimension identifiers in .arraydat")
    }
    #
    # reshape data
    #
    # if gfp or dissimilarity, modify the input array
    if (tanova) {
        if (tanova_type == "gfp") {
            .arraydat <- compGfp(.arraydat)
            full_dimnames <- dimnames(.arraydat)
        }
        if (tanova_type == "dissimilarity") 
            .arraydat <- scaleChan(.arraydat)
    }
    # set contrasts
    opcons <- options("contrasts")
    options(contrasts = c("contr.helmert", "contr.poly"))
    on.exit(options(opcons))
    # model dimensions should be in rows
    modeldims <- c(factordef$w_id, factordef$within)
    # all other dimensions become columns; "chan" and "time" come first
    keepdims <- setdiff(names(full_dimnames), modeldims)
    if ("time" %in% keepdims) keepdims <- unique(c("time", keepdims))
    if ("chan" %in% keepdims) keepdims <- unique(c("chan", keepdims))
    .arraydat <- mergeDims(.arraydat, list(modeldims, keepdims),
                           return_attributes = FALSE,
                           keep_dimnames = FALSE)
    # model data
    dat <- expand.grid(full_dimnames[modeldims], KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)
    colnames(dat) <- modeldims
    dat[,factordef$between] <- bwdat[,factordef$between]
    # transform characters to factors, and scale numeric variables
    dat[] <- lapply(
        names(dat), function(n) {
            x <- dat[[n]]
            if (storage.mode(x) == "character") {
                level <- 
                    if (n %in% factordef$between) {
                        levels(bwdat[[n]])
                    } else {
                        full_dimnames[[n]]
                    }
                factor_(x, levels = level)
            } else if (is.numeric(x)) {
                scale(x)
            } else {
                x
            }
        })
    #
    # formulas
    #
    mean_formula <- as.formula(paste(
        "~",
        gsub("^\\*|\\*$", "", 
             paste(
                 paste(factordef$between, collapse = "*"), 
                 paste(factordef$within, collapse = "*"),
                 sep = "*")
        )
    ))
    within_formula <- as.formula(paste(
        "~",
        gsub("^\\*|\\*$", "", 
             paste(
                 paste(factordef$w_id, collapse = "*"), 
                 paste(factordef$within, collapse = "*"),
                 sep = "*")
        )
    ))
    model_formula <- 
        if (type == "between" | tanova) {
            list(mean_formula)
        } else {
            list(mean_formula, within_formula)
        }
    #
    # model matrix and others
    #
    pre_sumsq <- vector("list", length(model_formula))
    for (i in seq_along(model_formula)) {
        tempy <- if (i == 1L && perm$type == "residuals") .arraydat else NULL
        pre_sumsq[[i]] <- preSumSq(model_formula[[i]], form_data = dat, 
                                   scaled = TRUE, y = tempy, 
                                   between = character())
    }
    if (type == "between" && !tanova) {
        pre_sumsq[[1]]$ssq_total <- colSums(
            sweepMatrix(.arraydat, 2, colMeans(.arraydat), "-", has_NA = FALSE)^2)
        setattr(pre_sumsq[[1]]$ssq_total, "Df", 
                length(levels(dat[[factordef$w_id]])) - 1L)
    }
    #
    # dimensions for later use
    #
    # test statistic
    pre_dimnames$modelterm <- full_dimnames$modelterm <- 
        attr(terms(mean_formula), "term.labels")
    full_dimid <- names(full_dimnames)
    full_dims <- vapply(full_dimnames, length, integer(1L))
    setattr(full_dims, "names", full_dimid)
    teststat_dimid <- c(keepdims, "modelterm")
    setattr(teststat_dimid, "origpos", 
            match(teststat_dimid, c(full_dimid_origord, "modelterm")))
    # not chan or time
    otherdims <- 
        if (!is.null(nr_chanXtime)) {
            other_dimid <- setdiff(teststat_dimid, c("chan", "time"))
            list(dimid = other_dimid,
                 index = match(other_dimid, teststat_dimid),
                 size = prod(full_dims[other_dimid]))
        } else {
            list(NULL)
        }
    # 
    #
    # return
    #
    list(.arraydat = .arraydat, dat = dat, factordef = factordef, 
         pre_sumsq = pre_sumsq, type = type, has_NA = has_NA,
         mean_formula = mean_formula,
         full_dimnames = full_dimnames, full_dims = full_dims,
         pre_dimnames = pre_dimnames,
         teststat_dimid = teststat_dimid,
         nr_chan = nr_chan, nr_time = nr_time, nr_chanXtime = nr_chanXtime,
         otherdims = otherdims)
}



#' Perform ANOVA (potentially with TFCE correction) on arrays
#' 
#' \code{arrayAnova} performs point-to-point ANOVAs on arrays. Permutation-based
#' p-values and Threshold-free Cluster Enhancement (TFCE) correction can be 
#' requested.
#' @param .arraydat a numeric array with named dimnames containing the EEG (or 
#' other) data. Missing values are not allowed.
#' @param factordef a named list of factor definitions, containing the following 
#' elements:
#' \itemize{
#' \item{between: }{character vector of between-subject factors (default: NULL)}
#' \item{within: }{character vector of within-subject factors (default: NULL)}
#' \item{w_id: }{name of the dimension which identifies the subjects 
#' (default: "id")}
#' \item{observed: }{character vector of observed (i.e., not manipulated) 
#' variables (default: NULL). Only used for effect-size calculations.}
#' }
#' @param bwdat a data.frame which contains the identification codes 
#' (factordef$w_id) and all subject-level variables (usually factors) listed in
#' 'factordef$between'. Missing values are not allowed.
#' @param verbose logical value indicating if p-values and effect sizes should 
#' be computed for the traditional ANOVA results
#' @param perm either 1) NULL (the default) or FALSE (the same as NULL), 
#' both of which mean no permutation, or 2) TRUE, which means 
#' permutation with default parameters, or 3) an object as returned by 
#' \code{\link{permParams}} with custom parameters (see Examples and also 
#' \code{\link{permParams}}).\cr 
#' Custom parameters can be also provided by \code{perm = .(key = value)} to 
#' save typing (this works by calling \code{\link{permParams}} with the 
#' given parameters).
#' @param tfce either 1) NULL (the default) or FALSE (the same as NULL), both of 
#' which mean no TFCE correction, or 2) TRUE, which means TFCE correction with 
#' default parameters, or 3) an object as returned by \code{\link{tfceParams}} 
#' with custom TFCE parameters (see Examples and also \code{\link{tfceParams}}).\cr 
#' Custom parameters can be also provided by \code{tfce = .(key = value)} to 
#' save typing (this works by calling \code{\link{tfceParams}} with the given 
#' parameters).
#' @param parallel either 1) NULL (the default) or FALSE (the same as NULL), 
#' both of which mean single-core computation, or 2) TRUE, which means 
#' parallelization with default parameters, or 3) an object as returned by 
#' \code{\link{parallelParams}} with custom parameters (see Examples and also 
#' \code{\link{parallelParams}}).\cr 
#' Custom parameters can be also provided by \code{parallel = .(key = value)} to 
#' save typing (this works by calling \code{\link{parallelParams}} with the 
#' given parameters).
#' @param seed an integer value which specifies a seed (default: NULL), or a 
#' list of arguments passed to \code{\link{set.seed}}
#' @details The function assumes that the input array contains at least three 
#' named dimensions: "chan" (corresponding to the channels [electrodes]), "time" 
#' (corresponding to time points), and a subject identifier as given by 
#' \code{factordef$w_id}. All dimensions which are not listed as 
#' within-subject factors are treated in a similar way as chan and time, that is 
#' separate ANOVA-s are computed for each level of those dimensions. 
#' @note The function computes type I p-values - this is correct if the design
#' is fully balanced and orthogonal (if the number of between-subject 
#' factors is one, it may have slightly unequal group sizes).
#' @export
#' @return A list object with the following numeric arrays:
#' \itemize{
#' \item{stat: }{a numeric array of F-values, with attribute 'Df' (an integer 
#' matrix) containing the term- and residual degrees of freedom for each model 
#' term. Optionally (if 'verbose' is TRUE), the F-value array has two extra 
#' attributes: 'p_value' and 'effect_size', consisiting of the traditional
#' p-values and the generalized eta squared effect size measures 
#' (see References), respectively.}
#' \item{stat_corr: }{a numeric array of TFCE-corrected F-values (if requested)}
#' \item{p_corr: }{a numeric array of permutation-based p-values (if requested)}
#' }
#' @seealso See also the related methods to explore the results, e.g. 
#' \code{\link{extract.arrayAnova}}, \code{\link{summary.arrayAnova}}, and the
#' plotting functions \code{\link{modelplot}}, or the lower-level 
#' \code{\link{imageValues}} and \code{\link{imagePvalues}}.
#' @references The TFCE correction follows:\cr
#' Mensen, A. and Khatami, R. (2013)
#' Advanced EEG analysis using threshold-free cluster-enhancement and 
#' non-parametric statistics. Neuroimage, 67, 111-118. 
#' doi:10.1016/j.neuroimage.2012.10.027  \cr
#' 
#' The Generalized Eta Squared effect size statistic is described in:\cr
#' Olejnik, S., Algina, J. (2003) Generalized eta and omega squared statistics: 
#' Measures of effect size for some common research designs. Psychological 
#' Methods 8: pp. 434-447. doi:10.1037/1082-989X.8.4.434 \cr
#' Bakeman, R. (2005) Recommended effect size statistics for repeated measures 
#' designs. Behavior Research Methods, 37 (3), 379-384.
#' @examples
#' # example dataset
#' data(erps)
#' dat_id <- attr(erps, "id") # to get group memberships
#' chan_pos <- attr(erps, "chan") # needed for TFCE correction
#' 
#' # make the dataset unbalanced to illustrate that if there is only one 
#' # between-subject factor, Type 1 and Type 2 results are identical
#' erps <- subsetArray(erps, list(id = 2:20))
#' dat_id <- dat_id[2:20, ]
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
#'                                    w_id = "id",
#'                                    observed = "group"),
#'                               bwdat = dat_id)
#' )
#' 
#' # if package 'ez' is installed, you can compare the results; we take only a 
#' # subset of the data (choosing "02" channel at time point "207") because 
#' # ezANOVA is much slower
#' if (requireNamespace("ez")) {
#'     sub <- list(chan = "O2", time = "207")
#'     tempdat_ez <- transformArray(y ~ ., tempdat, subset = sub)
#'     tempdat_ez$group <- factor(dat_id$group[match(tempdat_ez$id, dat_id$id)])
#'     result_ez <- ez::ezANOVA(tempdat_ez, y, id, .(stimclass, pairtype),
#'                              between = group, observed = group, type = 2)
#'     # compare results
#'     ez_F <- result_ez$ANOVA$F   # F-values
#'     ez_p <- result_ez$ANOVA$p   # p-values
#'     ez_es <- result_ez$ANOVA$ges   # effect sizes
#'     eegr_F <- as.vector(
#'         subsetArray(extract(result_eegr, "stat"), sub))
#'     eegr_p <- as.vector(
#'         subsetArray(extract(result_eegr, "p"), sub))
#'     eegr_es <- as.vector(
#'         subsetArray(extract(result_eegr, "es"), sub))
#'     stopifnot(
#'         all.equal(ez_F, eegr_F),
#'         all.equal(ez_p, eegr_p),
#'         all.equal(ez_es, eegr_es)
#'     )
#' }
#' 
#' # the between-subject variable could be numeric, too. Let's create an 'age'
#' # variable and analyze the effects.
#' dat_id$age <- rnorm(nrow(dat_id), 10, 1)
#' result_eegr_c <- arrayAnova(tempdat, 
#'                             list(between = "age",
#'                                  within = c("stimclass", "pairtype"),
#'                                  w_id = "id",
#'                                  observed = "age"),
#'                             bwdat = dat_id)
#'                            
#' # if package 'car' is installed, you can compare the results; we take only a 
#' # subset of the data (choosing "01" channel at time point "195")
#' if (requireNamespace("car")) {
#'     #
#'     # subsetting indices
#'     subs <- list(chan = "O1", time = "195") 
#'     #
#'     # extract F values and p values from the eegR results
#'     eegr_F <- subsetArray(extract(result_eegr_c, "stat"), subs)
#'     eegr_p <- subsetArray(extract(result_eegr_c, "p"), subs)
#'     #
#'     # run the same analysis with the 'car' package (first we have to reshape
#'     # the data)
#'     tempdat_car <- subsetArray(tempdat, subs)
#'     tempdat_car <- mergeDims(tempdat_car, 
#'                              list("id", c("stimclass", "pairtype")))
#'     tempdat_car <- data.frame(tempdat_car)
#'     tempdat_car$id <- dat_id$id
#'     tempdat_car$age <- dat_id$age
#'     idata <- expand.grid(dimnames(tempdat)[c("stimclass", "pairtype")])
#'     maov <- lm(cbind(A.ident, B.ident, C.ident, 
#'                      A.subst, B.subst, C.subst,
#'                      A.transp, B.transp, C.transp) ~ age,
#'                data = tempdat_car)
#'     maov <- car::Anova(maov, idata = idata, idesign= ~stimclass*pairtype,
#'                        type = 2)
#'     maov <- summary(maov, multivariate = FALSE)$univ
#'     car_F <- maov[names(eegr_F), "F"]
#'     car_p <- maov[names(eegr_p), "Pr(>F)"]
#'     #
#'     # compare results
#'     stopifnot(
#'         all.equal(as.vector(car_F), as.vector(eegr_F)),
#'         all.equal(as.vector(car_p), as.vector(eegr_p))
#'     )
#' }
#' 
#' # in order to use TFCE correction, the channel neigbourhood matrix is needed
#' # (see ?tfceParams and ?chanNb)
#' ChN <- chanNb(chan_pos, alpha = 0.7)
#' 
#' # now analyze the data by collapsing the pairtypes, and apply TFCE correction 
#' # (note: this will take a couple of seconds); use more randomization 
#' # runs (n should be several thousand instead of 499L) in serious analyses
#' tempdat <- avgDims(tempdat, "pairtype")
#' result_tfce <- arrayAnova(tempdat, 
#'                           list(between = "group",
#'                                within = "stimclass",
#'                                w_id = "id",
#'                                observed = "group"),
#'                           bwdat = dat_id,
#'                           perm = .(n = 499L), 
#'                           tfce = .(ChN = ChN),
#'                           parallel = .(ncores = 2))
#' 
#' # plot the corrected and uncorrected results
#' modelplot(result_tfce)
#' modelplot(result_tfce, type = "unc")
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
arrayAnova <- function(.arraydat, factordef, bwdat = NULL, verbose = TRUE, 
                       perm = NULL, tfce = NULL, parallel = NULL, 
                       seed = NULL) {
    # deparse tfce and parallel
    mcall <- match.call() 
    perm <- argumentDeparser(substitute(perm), "permParams",
                             null_params = list(n = 0L))
    tfce <- argumentDeparser(substitute(tfce), "tfceParams")
    ob <- getDoBackend()
    parallel <- argumentDeparser(substitute(parallel), "parallelParams",
                                 null_params = list(ncores = 0L))
    if (parallel$cl_new) {
        on.exit(stopCluster(parallel$cl))
    }
    on.exit(setDoBackend(ob), add = TRUE)
    #
    # compute statistics
    out <- arrayTtestAnova(
        "ANOVA", .arraydat, 
        factordef = factordef, bwdat = bwdat, verbose = verbose,
        perm = perm, tfce = tfce, parallel = parallel, seed = seed)
    #
    # replace call to the original
    out$call <- mcall
    # set class
    setattr(out, "class", "arrayAnova")
    # return
    out
}

#' United function for \code{\link{arrayTtest}} and \code{\link{arrayAnova}}
#' @keywords internal
arrayTtestAnova <- function(test, 
                            .arraydat, .arraydat2 = NULL, 
                            factordef = NULL, bwdat = NULL, 
                            paired = FALSE, groups = NULL, 
                            mu = 0, var_equal = FALSE, id_dim = "id", 
                            verbose = TRUE, perm = permParams(n = 0L), 
                            tfce = NULL, parallel = NULL, seed = NULL) {
    # helper function for back-transform to original
    backFn <- function(x, obj) {
        if (is.null(dim(x)) || !is.null(dimnames(x))) return(x)
        origpos <- order(attr(obj$teststat_dimid, "origpos"))
        apermArray(x, origpos, keep_attributes = TRUE)
    }
    # workaround to avoid CRAN warnings
    x <- NULL; i <- NULL
    rm(x, i)
    # prepare arguments
    nperm <- as.integer(perm$n)
    use_tfce <- !is.null(tfce)
    #
    if (test == "ANOVA") {
        #
        # prepare data
        input <- preAnova(.arraydat, factordef, bwdat, verbose = verbose,
                          tfce = use_tfce, perm = perm)
        rm(.arraydat)
        #
        # observed test statistic
        stat_obs <- compF(input, attribs = TRUE, verbose = verbose)
    } else if (test == "T_TEST") {
        #
        # prepare data
        input <- preTtest(.arraydat, .arraydat2, paired, groups, mu, var_equal, 
                          id_dim, verbose, tfce = use_tfce, perm = nperm > 1L)
        rm(.arraydat, .arraydat2)
        #
        # observed test statistic
        stat_obs <- 
            if (input$type == "independent_samples") {
                isTtest(input, attribs = TRUE)
            } else {
                osTtest(input, attribs = TRUE)
            }
    }
    #
    # TFCE correction
    if (use_tfce) {
        has_neg <- if (test == "ANOVA") FALSE else TRUE
        verbose_type <- sprintf(
            "TFCE-corrected %s statistic",
            if (test == "ANOVA") "F" else "t")
        tfce_obs <- anovaTfce(stat_obs, 1:2, 
                              has_neg = has_neg,
                              nr_chan = input$nr_chan, 
                              nr_time = input$nr_time, tfce = tfce)
        setattr(tfce_obs, "type", verbose_type)
    }
    #
    # permutations
    if (nperm > 1L) {
        obs <- if (use_tfce) tfce_obs else stat_obs
        perm_pvalues <- permPvalues(input, nperm, seed, obs, tfce)
    }
    #
    # prepare output
    out <- list(call = match.call())
    #
    # back-transform to original shape
    if (verbose) {
        setattr(stat_obs, "Df", backFn(attr(stat_obs, "Df"), input))
        setattr(stat_obs, "p_value", backFn(attr(stat_obs, "p_value"), input))
        setattr(stat_obs, "effect_size", 
                backFn(attr(stat_obs, "effect_size"), input))
    }
    stat_obs <- backFn(stat_obs, input)
    out$stat <- stat_obs
    if (use_tfce) {
        tfce_obs <- backFn(tfce_obs, input)
        out$stat_corr <- tfce_obs
    }
    if (nperm > 1L) {
        perm_pvalues <- backFn(perm_pvalues, input)
        out$p_corr <- perm_pvalues
    }
    #
    # attach dimnames
    dimind <- sort(attr(input$teststat_dimid, "origpos"))
    out$dimnames <- input$pre_dimnames[dimind]
    out$dim <- dim(out$stat)
    setattr(out$dim, "names", names(out$dimnames))
    #
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


#' Compute generalized GFP-effects
#' 
#' \code{compTanovaEffect} computes the generalized effect size measure as
#' given in Koenig & Melie-Garcia, 2009, p177. This is an internal function 
#' called by the high-level function \code{\link{tanova}}.
#' @param obj a list object as returned by \code{\link{preAnova}}
#' @param new_indices if not NULL (default), it must be a vector of numeric
#' indices to be passed to the model matrix. This parameter is used for 
#' permutation-based analyses.
#' @param attribs a logical value if dimension attributes should be attached to 
#' the resulting array (default: FALSE)
#' @return This function returns a numeric matrix of effects if 'attribs' is 
#' FALSE, and an array of effects if 'attribs' is TRUE.
#' @keywords internal
compTanovaEffect <- function(obj, new_indices = NULL, attribs = FALSE) {
    model <- obj$pre_sumsq[[1]]
    mm <- model$model_matrix
    if (!is.null(new_indices)) {
        if (!is.null(model$residuals)) {
            y <- model$residuals[new_indices,]
        } else {
            mm <- mm[new_indices, ]
            y <- obj$.arraydat
        }
    } else {
        y <- obj$.arraydat
    }
    mm_labels <- model$model_matrix_labels
    xpx <- model$xpx
    unc <- model$unique_contrasts
    #
    gfp <- !"chan" %in% obj$teststat_dimid
    #
    outdimn <- setdiff(obj$teststat_dimid, "chan")
    out <- matrix_(0, prod(obj$full_dims[setdiff(outdimn, "modelterm")]), 
                   obj$full_dims["modelterm"])
    #
    coeff <- solve(xpx, crossprod(mm, y))
    setattr(coeff, "dimnames", list(mm_labels, NULL))
    if (gfp) {
        for (i in seq_along(unc)) {
            fmeans <- unc[[i]] %*% coeff[colnames(unc[[i]]), , drop = FALSE]
            out[, i] <- colMeans(abs(fmeans))
        }
    } else {
        for (i in seq_along(unc)) {
            fmeans <- unc[[i]] %*% coeff[colnames(unc[[i]]), , drop = FALSE]
            setattr(fmeans, "dim", c(nrow(fmeans), obj$nr_chan, nrow(out)))
            out[, i] <- colMeans(compGfp(fmeans, channel_dim = 2L))
        }
    }
    if (attribs) {
        setattributes(out, obj, "Effect (generalized GFP)", outdimn)
        dimn <-  setNames(vector("list", length(outdimn)), outdimn)
        setattr(out, "dimnames", dimn)
    }
    #
    out
}

#' Compute p-values and length-corrected p-values in TANOVA
#' 
#' \code{compPvalueTanova} is a helper function to compute p-values and 
#' length-corrected p-values in the high-level \code{\link{tanova}} function.
#' @param effect_perm an array of permuted effects (including the observed 
#' effect)
#' @param pcrit a vector of p-value limits for consecutive length correction
#' @return a list of two arrays: the p-values and the corrected p-values
#' @keywords internal
compPvalueTanova <- function(effect_perm, pcrit, obj) {
    # helper function
    consecLimit <- function(sigvec, .dim, alpha) {
        out <- matrixRle(matrix_(sigvec, .dim))
        quantile(out$lengths[out$values > 0], 1 - alpha)
    }
    # p-values
    dims <- dim(effect_perm)
    names(dims) <- names(dimnames(effect_perm))
    nperm <- dims["perm"]
    pvalues_perm <- fnDims(effect_perm, "perm", colRanks, 
                           arg_list = list(preserveShape = TRUE),
                           vectorized = TRUE, keep_dimorder = FALSE)
    pvalues_perm <- (nperm + 1 - pvalues_perm) / nperm
    pvalues <- subsetArray(pvalues_perm, list(perm = 1L))
    setattr(pvalues, "dimnames", dimnames(pvalues_perm)[-1L])
    setattr(pvalues, "type", "Permuted P-value")
    out <- list(p = pvalues)
    # consecutive sign. criterion
    if ("time" %in% names(dims)) {
        pvalues_perm <- apermArray(pvalues_perm, first = c("time", "perm"))
        pvalues_consec <- array_(1, dim(pvalues_perm)[-2L], 
                                 dimnames(pvalues_perm)[-2L])
        for (alpha in sort(pcrit, TRUE)) {
            consec_limit <- fnDims(pvalues_perm < alpha, 
                                   1:2, consecLimit, 
                                   arg_list = list(alpha = alpha))
            temp <- matrixRle(
                array2mat(pvalues < alpha, "time", 
                          return_attributes = FALSE,
                          keep_dimnames = FALSE))
            ind <- temp$lengths < consec_limit[temp$matrixcolumn]
            temp$values[ind] <- 0
            temp <- inverse.matrixRle(temp)
            pvalues_consec[temp > 0] <- alpha
        }
        pvalues_consec <- apermArray(pvalues_consec, names(dimnames(pvalues)))
        setattr(pvalues_consec, "type", 
                "Permuted P-value (with min. length correction)")
        out[["p_corr"]] <- pvalues_consec
    }
    # return
    out
}


#' Topographical ANOVA (TANOVA) and related methods
#' 
#' \code{tanova} performs point-to-point topographical ANOVA on arrays. Related
#' methods are GFP-analysis (which is based on the intensity of the signal) and
#' DISS-analysis (focusing only on the global dissimilarity by comparing 
#' normalized topographies). See References for further details. 
#' @param .arraydat a numeric array with named dimnames containing the EEG (or 
#' other) data. Missing values are not allowed.
#' @param factordef a named list of factor definitions, containing the following 
#' elements:
#' \itemize{
#' \item{between: }{character vector of between-subject factors (default: NULL)}
#' \item{within: }{character vector of within-subject factors (default: NULL)}
#' \item{w_id: }{name of the dimension which identifies the subjects 
#' (default: "id")}
#' }
#' @param bwdat a data.frame which contains the identification codes 
#' (factordef$w_id) and all subject-level variables (usually factors) listed in
#' 'factordef$between'. Missing values are not allowed.
#' @param type a character value of "tanova" (default), "dissimilarity", or 
#' "gfp"
#' @param verbose logical value indicating if p-values should be computed
#' @param perm either 1) NULL or FALSE (the same as NULL), 
#' both of which mean no permutation, or 2) TRUE, which means 
#' permutation with default parameters (the default), or 3) an object as 
#' returned by \code{\link{permParams}} with custom parameters (see Examples 
#' and also \code{\link{permParams}}).\cr 
#' Custom parameters can be also provided by \code{perm = .(key = value)} to 
#' save typing (this works by calling \code{\link{permParams}} with the 
#' given parameters).
#' @param parallel either 1) NULL (the default) or FALSE (the same as NULL), 
#' both of which mean single-core computation, or 2) TRUE, which means 
#' parallelization with default parameters, or 3) an object as returned by 
#' \code{\link{parallelParams}} with custom parameters (see Examples and also 
#' \code{\link{parallelParams}}).\cr 
#' Custom parameters can be also provided by \code{parallel = .(key = value)} to 
#' save typing (this works by calling \code{\link{parallelParams}} with the 
#' given parameters).
#' @param seed an integer value which specifies a seed (default: NULL), or a 
#' list of arguments passed to \code{\link{set.seed}}
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
#' @seealso See also the related methods to explore the results, e.g. 
#' \code{\link{extract.tanova}}, \code{\link{summary.tanova}}, and the plotting
#' function \code{\link{modelplot.tanova}}.
#' @references Koenig, T., Melie-Garcia, L. (2009) Statistical analysis of 
#' multichannel scalp field data. In: Michel, C.M., Koenig, T., Brandeis, D.,
#' Gianotti, L.R.R., Wackermann, J. (eds) Electrical neuroimaging. Cambridge
#' University Press\cr
#' See also the MATLAB toolbox by Koenig T, Kottlow M, Stein M, 
#' Melie-Garcia L (2011) Ragu: A free tool for the analysis of EEG and MEG 
#' event-related scalp field data using global randomization statistics. 
#' Computational Intelligence and Neuroscience, 2011:938925
#' See also 
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
#' # purposes (here we use the default: 999L).
#' result_tanova <- tanova(tempdat, 
#'                         list(between = "group",
#'                              within = c("stimclass", "pairtype"),
#'                              w_id = "id"),
#'                         bwdat = dat_id,
#'                         parallel = .(ncores = 2))
#' 
#' # plot results (for now, only p-values and only for time > 0)
#' modelplot(result_tanova, what = "p", time_window = c(0, Inf))
tanova <- function(.arraydat, factordef, bwdat = NULL, 
                   type = c("tanova", "dissimilarity", "gfp"), 
                   verbose = TRUE, perm = TRUE, parallel = NULL,
                   seed = NULL, pcrit = 0.05) {
    # helper function for back-transform to original
    backFn <- function(x, obj) {
        origpos <- attr(obj$teststat_dimid, "origpos")
        origpos <- origpos[obj$teststat_dimid != "chan"]
        x <- apermArray(x, order(origpos), keep_attributes = TRUE)
        setattr(x, "dimnames", NULL)
    }
    # CRAN check
    x <- i <- NULL
    rm(x, i)
    # some arguments must be pre-checked (e.g. deparse parallel)
    type <- match.arg(type)
    mcall <- match.call() 
    ob <- getDoBackend()
    perm <- argumentDeparser(substitute(perm), "permParams",
                             null_params = list(n = 0L))
    parallel <- argumentDeparser(substitute(parallel), "parallelParams",
                                 null_params = list(ncores = 0L))
    if (parallel$cl_new) {
        on.exit(stopCluster(parallel$cl))
    }
    on.exit(setDoBackend(ob), add = TRUE)
    #
    # prepare data
    input <- preAnova(.arraydat, factordef, bwdat, verbose,
                      tfce = FALSE, perm = perm,
                      tanova_type = type)
    rm(.arraydat)
    #
    # prepare output
    out <- list(call = mcall)
    #
    # calculate effect
    es_obs <- compTanovaEffect(input, attribs = TRUE)
#     if (verbose) {
#         factor_means <- marginalMeans(input$mean_formula, 
#                                       input$dat, input$.arraydat, 
#                                       dimnames(es_obs)[-length(dim(es_obs))], 
#                                       keep_term_order = !iaterms_last, 
#                                       residualmean = FALSE)
#         out$factor_means <- factor_means
#     }
    if (perm$n > 1L) {
        # generate random orders (dim(randind) = nrow(dat) x nperm)
        randind <- anovaRandomIndices(input, perm$n, seed)
        # run calculations
        chunks <- min(10L, ceiling(perm$n/max(1L, length(parallel$cl))))
        es_perm <- 
            foreach(x = iter(randind, by = "col", chunksize = chunks), 
                    .combine = c, .inorder = FALSE,
                    .options.snow = parallel$snow_options,
                    .options.multicore = parallel$mc_options) %dopar% 
                    {
                        foreach(i = 1:ncol(x)) %do% 
                            compTanovaEffect(input, x[,i])
                    }
        # reshape permuted tanova values
        es_perm <- c(es_obs, unlist(es_perm, use.names = FALSE))
        setattr(es_perm, "dim", c(dim(es_obs), perm$n + 1L))
        dimn <-  c(dimnames(es_obs), list(perm = NULL))
        setattr(es_perm, "dimnames", dimn)
        # p-values
        pvals <- compPvalueTanova(es_perm, pcrit, input)
    }
    # back-transform to original dimorder
    out$stat <- backFn(es_obs, input)
    rm(es_obs)
    if (perm$n > 1L) {
        pvals <- lapply(pvals, backFn, obj = input)
        setattr(out$stat, "p_value", pvals[["p"]])
        if (!is.null(pvals[["p_corr"]])) 
            out$p_corr <- pvals[["p_corr"]]
    }
    #
    # attach dimnames
    notchan <- input$teststat_dimid != "chan"
    dimind <- sort(attr(input$teststat_dimid, "origpos")[notchan])
    out$dimnames <- input$pre_dimnames[dimind]
    out$dim <- dim(out$stat)
    setattr(out$dim, "names", names(out$dimnames))
    # return
    setattr(out, "class", "tanova")
    out
}


# Peak Anova functions =========== 

#' ANOVA on individual peaks
#' 
#' \code{peakAnova} performs automatic peak detection and runs traditional 
#' ANOVA on them
#' @param .arraydat a numeric array with named dimnames containing the EEG (or 
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
#' @details OLD VERSION, DO NOT USE IT! updates will follow...
#' @note The function computes type I p-values - this is correct if the design
#' is fully balanced and orthogonal (if the number of between-subject 
#' factors is one, it may have slightly unequal group sizes).
#' @export
#' @return A list object.
peakAnova <- function(.arraydat, factordef, peakdef, bwdat = NULL,
                      verbose = TRUE, avg_around_peak = 2, 
                      nperm = 1, useparallel = FALSE, ncores = NULL,
                      par_method = c("snow", "mc"), cl = NULL,
                      iaterms_last = TRUE, seed = NULL) {
    # some checks
    stopifnot(is.array(.arraydat) | missing(.arraydat))
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
    sumSq <- function(.arraydat, f_dat, aov_form, labels = FALSE) {
        out <- summary(aov(as.formula(aov_form), data = f_dat))
        if (!labels) {
            out <- matrix_(unlist(lapply(out, "[", "Mean Sq"), 
                                  use.names = FALSE),
                           nrow(out[[1]]), ncol(.arraydat))
        } else {
            out <- matrix_(unlist(lapply(out, "[", "Mean Sq"), 
                                  use.names = FALSE),
                           nrow(out[[1]]), ncol(.arraydat),
                           list(term = gsub(" ","",rownames(out[[1]])), 
                                peak = seq_along(out)))
        }         
        return(out)
    }
    # 
    par_method <- match.arg(par_method)
    # prepare data
    temp <- preAnova(.arraydat, factordef, bwdat, useparallel, ncores, par_method)
    # assign variables to this environment, and potentially overwrite existing ones
    #assignList(temp, verbose = FALSE)
    dat <- temp$dat
    .arraydat <- temp$.arraydat
    factordef <- temp$factordef
    origdimnames <- temp$origdimnames
    par_params <- temp$par_params
    rm(temp)
    #
    timedim <- which(names(dimnames(.arraydat)) == "time")
    origtimedim <- which(names(origdimnames) == "time")
    origdimnames_peaks <- origdimnames
    .arraydat_peaks <- subsetArray(.arraydat, 
                                  list(time = dimnames(.arraydat)$time[seq_along(peakdef)]))
    origdimnames_peaks$time <- 
        dimnames(.arraydat_peaks)$time <- seq_along(peakdef)
    names(origdimnames_peaks)[origtimedim] <- 
        names(dimnames(.arraydat_peaks))[timedim] <- "peak"
    # formula of anova
    aov_formula_char <- 
        if (is.null(factordef$within)) {
            paste(
                as.character(quote(.arraydat)), " ~ ", 
                paste(as.character(factordef$between), collapse = "*"),
                sep = "")
        } else if (is.null(factordef$between)) {
            paste(
                as.character(quote(.arraydat)), " ~ ", 
                paste(as.character(factordef$within), collapse = "*"),
                sep = "")
        } else {
            paste(
                as.character(quote(.arraydat)), " ~ ", 
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
    fmeans <- marginalMeans(mean_formula, dat, .arraydat, 
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
    .arraydat_peaks[] <- sapply(1:ncol(peakind), function(i) {
        out <- matrix_(0, nrow(.arraydat), avg_around_peak*2 + 1)
        avgind <- outer(peakind[,i], c(-avg_around_peak:avg_around_peak), "+")
        out[] <- .arraydat[
            cbind(seq_along(datrowind), as.integer(avgind[datrowind,]))]
        return(rowMeans(out))
    })
    #
    # Compute ANOVA on peak amplitudes
    Fvals_obs <- arrayAnovaSub(.arraydat_peaks, 
                               factordef, origdimnames_peaks, dat, verbose)
    out <- c(out, list(F_obs = Fvals_obs))
    # if verbose, save factor means
    if (verbose) {
        factor_means <- data.frame(
            group = peakind_facs[rep(1:nrow(peakind_facs), 
                                     length(peakdef)), , drop = F],
            peak = factor_(rep(seq_along(peakdef), each = nrow(fmeans[[1]]))),
            ampl = c(marginalMeans(mean_formula, dat, .arraydat_peaks, 
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
                                      .arraydat, 
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
                             matrix_(unlist(lateff_perm, use.names = FALSE),
                                     ncol = nperm))
        lateff_perm <- sweep(lateff_perm[, -1, drop = F], 1, 
                             lateff_perm[, 1], "-")
        pvalues <- array_(
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