#
# <<< t-tests on arrays >>>----
#

#' Workhorse functions of arrayTest
#' 
#' \code{isTtest} and \code{osTtest} compute t-values and additional p-values & 
#' degrees of freedom if requested for independent samples (\code{isTtest}) and
#' one-sample or paired-samples (\code{osTest}).
#' @name Ttest
#' @param dat numeric matrix or array with named dimension names
#' @param groups integer vector with values 1 and 2, corresponding to the group
#' membership of the subjects (see \code{id_dim})
#' @param id_dim character; the dimension of \code{dat} which refers to the 
#' subjects (default: "id")
#' @param mu a number indicating the true value of the mean (or difference in 
#' means if one is performing a two sample test)
#' @param var_equal a logical variable indicating whether to treat the two 
#' variances as being equal. If TRUE then the pooled variance is used to 
#' estimate the variance otherwise the Welch (or Satterthwaite) approximation 
#' to the degrees of freedom is used.
#' @param verbose logical value; if TRUE, degrees of freedom and p-values are
#' also returned (default: FALSE)
#' @param has_NA logical indicating whether there is any missing value in 
#' \code{dat}. If NULL (default), it is explicitly checked by calling 
#' \code{\link{anyNA}} on \code{dat}.
#' @keywords internal
#' @return An array of the same dimensions as \code{dat} except for the
#' \code{id_dim} dimension. Auxiliary results (degrees of freedom and p-values)
#' are added as attributes (if requested).
NULL

#' @rdname Ttest
# independent samples t-test
isTtest <- function(dat, groups, id_dim = "id", mu = 0, var_equal = TRUE, 
                    verbose = FALSE, has_NA = NULL) {
    if (is.null(has_NA)) has_NA <- anyNA(dat)
    res <- lapply(1:2, function(g) {
        datx <- subsetArray(dat, listS(.id_dim = groups==g), drop = FALSE)
        mm <- avgDims(datx, id_dim, na_rm = has_NA)
        vv <- fnDims(datx, id_dim, colVars, arg_list = list(na.rm = has_NA), 
                     vectorized = TRUE, columnwise = TRUE)
        if (has_NA) {
            nn <- fnDims(!is.na(datx), id_dim, colSums, 
                         arg_list = list(na.rm = FALSE), 
                         vectorized = TRUE, columnwise = TRUE)
            nn[nn < 1] <- NA
        } else {
            nn <- sum(groups == g)
        }
        list(mm, vv, nn)
    })
    m1 <- res[[1]][[1]]; v1 <- res[[1]][[2]]; n1 <- res[[1]][[3]]
    m2 <- res[[2]][[1]]; v2 <- res[[2]][[2]]; n2 <- res[[2]][[3]]
    out <- 
        if (var_equal) {
            sd12 <- sqrt(
                ((n1-1)*v1 + (n2-1)*v2)/(n1 + n2 - 2)
            )
            (m1 - m2)/(sd12*sqrt(1/n1 + 1/n2))
        } else {
            sd12 <- sqrt(v1/n1 + v2/n2)
            (m1 - m2)/sd12
        }
    if (verbose) {
        df <- 
            if (var_equal) {
                n1 + n2 - 2 
            } else {
                sd12^4/( (v1/n1)^2/(n1-1) + (v2/n2)^2/(n2-1) )
            }
        p <- 2 * pt(-abs(out), df)
        setattr(out, "Df", df)
        setattr(out, "pvalues", p)
    }
    out
}

#' @rdname Ttest
# one-sample or paired samples t-test
osTtest <- function(dat, id_dim = "id", mu = 0, verbose = FALSE, has_NA = NULL) {
    if (is.null(has_NA)) has_NA <- anyNA(dat)
    m1 <- avgDims(dat, id_dim, na_rm = has_NA)
    sd1 <- fnDims(dat, id_dim, colSds, 
                  arg_list = list(na.rm = has_NA), vectorized = TRUE, 
                  columnwise = TRUE)
    if (has_NA) {
        n1 <- fnDims(!is.na(dat), id_dim, colSums, 
                     arg_list = list(na.rm = FALSE), vectorized = TRUE, 
                     columnwise = TRUE)
        n1[n1 < 1] <- NA
    } else {
        if (is.character(id_dim)) id_dim <- match(id_dim, names(dimnames(dat)))
        n1 <- dim(dat)[id_dim]
    }
    out <- (m1 - mu)/(sd1/sqrt(n1))
    if (verbose) {
        pvalues <- 2 * pt(-abs(out), n1 - 1)
        setattr(out, "Df", n1 - 1)
        setattr(out, "pvalues", pvalues)
    }
    out
}


#' Point-to-point t-tests (potentially with TFCE correction) on arrays
#' 
#' \code{arrayTtest} performs point-to-point t-tests on arrays. 
#' Permutation-based p-values and Threshold-free Cluster Enhancement (TFCE) 
#' correction can be requested.
#' @param arraydat a numeric array with named dimnames containing EEG (or 
#' other) data. Missing values are not allowed. Must have at least three
#' dimensions with names "chan", "time", and "id" (see also id_dim)
#' @param arraydat2 a numeric array with named dimnames containing EEG (or 
#' other) data. If provided, see the parameter \code{paired} for running 
#' dependent or independent-samples t-tests
#' @param paired logical scalar, only used if arraydat2 is provided. If paired 
#' is FALSE (default), the function computes independent samples t-tests, 
#' otherwise paired samples t-tests are performed
#' @param groups provides an alternative (and more efficient) way to perform 
#' independent samples t-tests; a character, factor, or integer vector which 
#' defines group membership. Groups is ignored if arraydat2 is not missing. 
#' NA values code subjects to drop.
#' @param mu a numeric scalar indicating the true value of the mean (or 
#' difference between means, if two-sample tests are performed)
#' @param var_equal a logical scalar whether the variances are equal (only 
#' relevant for independent-samples t-tests). If TRUE (default), the pooled 
#' variance is used to estimate the variance, otherwise the Welch (or 
#' Satterthwaite) approximation to the degrees of freedom is used.
#' @param id_dim name of the dimension which identifies the subjects 
#' (default: "id")
#' @param verbose logical value indicating if p-values should be computed for 
#' the traditional t-test results
#' @param nperm integer value giving the number of permutations (default: 999L)
#' @param useparallel logical value; if TRUE (default), computations are done
#' in parallel (only if nperm > 1)
#' @param par_method parallelization method; can be set explicitly to "snow" or
#' "multicore" (ignored on Windows OS), or chosen automatically ("default"). 
#' Ignored if cl is provided.
#' @param ncores integer value corresponding to the number of cores; 
#' if NULL (default), it is set to the maximum number of cores available. 
#' Ignored if cl is provided.
#' @param cl a cluster definition; if NULL (default), it is set up automatically
#' @param usetfce logical value whether TFCE (threshold-free cluster enhancement)
#' correction should also be computed (default: TRUE)
#' @param tfce_options a named list containing the channel neighbourhood matrix
#' (named ChN, no default) and the vector of the E/H parameters (named EH, 
#' defaults to c(0.66, 2))
#' @param seed an integer value which specifies a seed (default: NULL)
#' @details The function assumes that the input array contains at least three 
#' named dimensions: chan (corresponding to the channels [electrodes]) and time 
#' (corresponding to time points), and \code{id_dim} (corresponding to subjects). 
#' All other dimensions are treated in a similar way as chan and time, that is 
#' separate t-tests are computed for each level of those dimensions.
#' @export
#' @return A list object with t-values, TFCE-corrected t-values and
#' permutation-based p-values (if requested)
#' @examples
#' # example dataset
#' data(erps)
#' dat_id <- attr(erps, "id") # to get group memberships
#' chan_pos <- attr(erps, "chan") # needed for TFCE correction
#' 
#' # compare controls and dyslexics with traditional indep. samples t-test;
#' # comparison is performed for each time sample, channel and experimental 
#' # condition separately (altogether 81972 tests - note the speed)
#' system.time(
#'     result_eegr <- arrayTtest(erps, groups = dat_id$group, nperm = 0L,
#'                               var_equal = FALSE, usetfce = FALSE)
#' )
#' 
#' # the built-in R function (t.test) provides the same result, but running it on
#' # the full dataset would take much more time; 
#' # here we take a subsample of the data
#' sub <- list(chan = "F4", time = "200", stimclass = "B", pairtype = "ident")
#' result_ttest <- t.test(subsetArray(erps, sub) ~ dat_id$group)
#' 
#' # to check that they are equivalent, we have to remove the attributes
#' eegr_t <- as.vector(subsetArray(result_eegr$effect_t_obs, sub))
#' eegr_p <- as.vector(subsetArray(attr(result_eegr$effect_t_obs, "pvalues"), sub))
#' stopifnot(
#'     all.equal(as.vector(result_ttest$statistic), eegr_t),
#'     all.equal(as.vector(result_ttest$p.value), eegr_p)
#' )
#' 
#' # Now let's use TFCE correction; to do that, one needs a channel neigbourhood
#' # matrix and has to use randomization. We will simplify a bit the data to 
#' # decrease the computational burden.
#' # 1) get channel neighbourhoods (type ?chanNb)
#' ChN <- chanNb(chan_pos, alpha = 0.7)
#' 
#' # 2) analyze only stimclass "B" and pairtype "ident"
#' tempdat <- subsetArray(erps, list(stimclass = "B", pairtype = "ident"))
#' 
#' # 3) run computations (now with only 499 permutations, and using 2 CPU-cores)
#' result_tfce <- arrayTtest(tempdat, groups = dat_id$group, var_equal = FALSE, 
#'                           nperm = 499L, ncores = 2L, usetfce = TRUE,
#'                           tfce_options = list(ChN = ChN))
#' 
#' # 4) compare the traditional and TFCE-corrected p-values
#' p_trad <- attr(result_tfce$effect_t_obs, "pvalues")
#' p_tfce <- result_tfce$perm_pvalues
#' p_all <- bindArrays(trad = p_trad, tfce = p_tfce, along_name = "method")
#' 
#' # 5) plot p-values after -log transform for better discriminability
#' # note how the sporadic effects disappear after TFCE correction
#' p_plot <- imageValues(-log(p_all))  # returns a ggplot object
#' p_plot
#
# TODO: compute also effect sizes + make it general for any arrays without chan
# and time dimensions + clear redundancies in the code for the two types of 
# t-tests
arrayTtest <- function(arraydat, arraydat2, paired = FALSE, groups = NULL,
                       mu = 0, var_equal = TRUE, id_dim = "id", verbose = TRUE, 
                       nperm = 999L, useparallel = TRUE, 
                       par_method = c("default", "snow", "multicore"), 
                       ncores = NULL, cl = NULL, usetfce = TRUE, 
                       tfce_options = NULL, seed = NULL) {
    # pass rcmd check
    x <- NULL; rm(x)
    # helper function for back-transform to original
    backFn <- function(x, origdimnames) {
        if (is.null(dim(x))) return(x)
        if ("TEMPX" %in% (dimn <- names(dimnames(x))))
            x <- subsetArray(x, list(TEMPX=1))
        out <- aperm(x, origdimnames[origdimnames %in% dimn])
        attribs <- attributes(x)
        attribs <- attribs[-match(c("dim", "dimnames"), names(attribs))]
        for (i in names(attribs)) setattr(out, i, attribs[[i]])
        # return
        out
    }
    # 
    if (nperm < 2) useparallel <- FALSE
    #
    if (useparallel) {
        par_method <- match.arg(par_method)
        if (newcl <- is.null(cl)) {
            if (is.null(ncores)) ncores <- detectCores()
            cl <- 
                if (.Platform$OS.type=="windows" || par_method == "snow") {
                    makePSOCKcluster(ncores)
                } else {
                    makeForkCluster(ncores)
                }
        }
        registerDoParallel(cl)
        on.exit(if (newcl) stopCluster(cl))
    } 
    else {
        registerDoSEQ()
    }
    if (usetfce) {
        if (is.null(tfce_options)) stop("Provide at least tfce_options$ChN")
        ChN <- tfce_options$ChN
        EH <- if (is.null(tfce_options$EH)) c(0.66, 2) else tfce_options$EH
    }
    if (missing(arraydat2)) {
        if (!is.null(groups)) {
            if (anyNA(groups)) 
                arraydat <- subsetArray(arraydat, listS(.id_dim=!is.na(groups)))
            groups <- as.integer(as.factor(groups))
            if (max(groups) != 2L) 
                stop("groups must contain two distinct values (e.g. 1 or 2, 'a' or 'b')")
        }
    } 
    else {
        if (paired) {
            if (!identical(dimnames(arraydat), dimnames(arraydat2)) ||
                    !identical(dim(arraydat), dim(arraydat2))) {
                stop("Input array dimensions must be identical")
            }
            arraydat <- arraydat2 - arraydat
            groups <- NULL
        } else {
            groups <- rep.int(1:2, 
                              c(length(dimnames(arraydat)[[id_dim]]),
                                length(dimnames(arraydat2)[[id_dim]])))
            dimnms <- names(dimnames(arraydat))
            arraydat <- abind(arraydat, arraydat2, 
                              along=which(names(dimnames(arraydat))==id_dim))
            setattr(dimnames(arraydat), "names", dimnms)
        }
    }
    has_NA <- anyNA(arraydat)
    origdimnames <- names(dimnames(arraydat))
    if (anyNA(origdimnames) || 
            !all(c(id_dim, "chan", "time") %in% origdimnames) ||
            any(duplicated(origdimnames))) {
        stop("Provide array(s) with named dimensions (at least id_dim, 'chan', 'time')")
    }
    arraydat <- aperm(arraydat, c(id_dim, "chan", "time",
                                  setdiff(names(dimnames(arraydat)),
                                          c(id_dim, "chan", "time"))))
    if (length(dim(arraydat)) == 3L) {
        temp <- dimnames(arraydat)
        arrayIP(arraydat, c(dim(arraydat), 1),
                c(temp, list(TEMPX="1")))
        rm(temp)
    }
    # two independent samples
    if (!is.null(groups)) {
        t_obs <- isTtest(arraydat, groups, id_dim, mu, 
                         var_equal, verbose, has_NA)
        if (nperm > 1) {
            grouplen <- length(groups)
            group2len <- sum(groups == 2L)
            maxperm <- choose(grouplen, group2len)
            if (nperm > maxperm) 
                stop("nperm is too high for this sample size (not enough unique combinations exist)")
            if (!is.null(seed)) set.seed(seed)
            if (nperm > maxperm/3) {
                randind <- combn(grouplen, group2len)
                randgroups <- matrix(1L, grouplen, ncol(randind))
                for (i in 1:ncol(randind)) randgroups[randind[,i], i] <- 2L
                randgroups <- randgroups[, sample.int(nrow(randgroups), nperm)]
            } 
            else {
                randgroups <- matrix(groups, grouplen, 1)
                repeat{
                    nc <- ncol(randgroups) - 1
                    if (nc >= nperm) {
                        randgroups <- randgroups[, 2:(nperm + 1)]
                        break
                    } 
                    else {
                        rg <- matrix(1L, grouplen, (nperm - nc) * 2)
                        for (i in 1:ncol(rg)) {
                            rg[sample.int(grouplen, group2len), i] <- 2L
                        }
                        randgroups <- cbind(randgroups, rg)
                        randgroups <- unique(randgroups, MARGIN = 2L)
                    }
                }
            }
            applydims <- seq(3, length(dim(t_obs)))
            if (usetfce) {
                tfce_obs <- arrayIP(0, dim(t_obs), dimnames(t_obs))
                tfce_obs[] <- apply(t_obs, applydims, tfceFn, 
                                    chn = ChN, eh = EH)
                stat_obs <- abs(tfce_obs)
            } 
            else {
                stat_obs <- abs(t_obs)
            }
            chunks <- min(10L, ceiling(nperm/max(1L, length(cl))))
            permres <- 
                foreach(x = iter(randgroups, by = "column", 
                                 chunksize = chunks), 
                        .combine = "+", .inorder=FALSE) %dopar% {
                            
                            out <- 
                                foreach(i = 1:ncol(x), .combine = "+") %do% {
                                    t_perm <- isTtest(arraydat, x[,i], 
                                                      id_dim, mu, var_equal, 
                                                      verbose = FALSE, has_NA)
                                    if (usetfce) {
                                        t_perm[] <- 
                                            apply(t_perm, applydims, 
                                                  tfceFn, chn = ChN, eh = EH)
                                    }
                                    t_perm <- colMaxs(
                                        mergeDims(abs(t_perm), 
                                                  list(1:2, seq_along(dim(t_perm))[-(1:2)]), 
                                                  keep_dimnames = FALSE, 
                                                  return_attributes = FALSE)
                                    )
                                    sweep(stat_obs, applydims, t_perm, "<")
                                }
                        }
            permres <- (permres + 1) / (nperm + 1)
        }
    } 
    # one-sample (can be difference of two paired samples as well)
    else {
        t_obs <- osTtest(arraydat, id_dim, mu, verbose, has_NA)
        if (nperm > 1) {
            idlen <- length(dimnames(arraydat)[[id_dim]])
            maxperm <- 2^idlen
            if (nperm > maxperm) 
                stop("nperm is too high for this sample size (not enough unique combinations exist)")
            if (!is.null(seed)) set.seed(seed)
            if (nperm > maxperm/3) {
                randi <- matrix(unlist(
                    expand.grid(rep(list(c(-1L, 1L)), idlen), 
                                KEEP.OUT.ATTRS = FALSE),
                    use.names = FALSE), idlen, maxperm, TRUE)
                randi <- randi[, sample(ncol(randi), nperm)]
            } 
            else {
                randi <- matrix(rep.int(1L, idlen))
                repeat{
                    nc <- ncol(randi) - 1
                    if (nc >= nperm) {
                        randi <- randi[, 2:(nperm + 1)]
                        break
                    } 
                    else {
                        ri <- matrixIP(sample(c(-1L, 1L), 
                                              temp <- idlen * (nperm - nc) * 2,
                                              TRUE),
                                       idlen, temp)
                        randi <- unique(cbind(randi, ri), MARGIN = 2L)
                    }
                }
            }
            applydims <- seq(3, length(dim(t_obs)))
            if (usetfce) {
                tfce_obs <- arrayIP(0, dim(t_obs), dimnames(t_obs))
                tfce_obs[] <- apply(t_obs, applydims, tfceFn, 
                                    chn = ChN, eh = EH)
                stat_obs <- abs(tfce_obs)
            } 
            else {
                stat_obs <- abs(t_obs)
            }
            chunks <- min(10L, ceiling(nperm/max(1L, length(cl))))
            permres <- 
                foreach(x = iter(randi, by = "column", 
                                 chunksize = chunks), 
                        .combine = "+", .inorder=FALSE) %dopar% {
                            out <- 
                                foreach(i = 1:ncol(x), .combine = "+") %do% {
                                    t_perm <- osTtest(arraydat * x[,i], 
                                                      id_dim, mu, 
                                                      verbose = FALSE, has_NA)
                                    if (usetfce) {
                                        t_perm[] <- 
                                            apply(t_perm, applydims, 
                                                  tfceFn, chn = ChN, eh = EH)
                                    }
                                    t_perm <- colMaxs(
                                        mergeDims(abs(t_perm), 
                                                  list(1:2, seq_along(dim(t_perm))[-(1:2)]), 
                                                  keep_dimnames=FALSE, 
                                                  return_attributes=FALSE)
                                    )
                                    sweep(stat_obs, applydims, t_perm, "<")
                                }
                        }
            permres <- (permres + 1) / (nperm + 1)
        }
    }
    # back-transform to original shape
    if (verbose) {
        setattr(t_obs, "Df", backFn(attr(t_obs, "Df"), origdimnames))
        setattr(t_obs, "pvalues", backFn(attr(t_obs, "pvalues"), origdimnames))
    }
    t_obs <- backFn(t_obs, origdimnames)
    out <- list(call = match.call(), effect_t_obs = t_obs)
    if (usetfce) {
        tfce_obs <- backFn(tfce_obs, origdimnames)
        out <- c(out, list(effect_tfce_obs = tfce_obs))
    }
    if (nperm > 0) {
        permres <- backFn(permres, origdimnames)
        out <- c(out, list(perm_pvalues = permres))
    }
    # return
    out
}