#
# <<< t-tests on arrays >>>----
#

#' Workhorse functions of arrayTest
#' 
#' \code{isTtest} and \code{osTtest} compute t-values and additional p-values & 
#' degrees of freedom if requested for independent samples (\code{isTtest}) and
#' one-sample or paired-samples (\code{osTest}).
#' @name Ttest
#' @param obj a list object returned by \code{\link{preTtest}}
#' @param new_group integer vector with values 1 and 2, corresponding to the 
#' group membership of the subjects. See \code{\link{isTtestRandomGroups}}.
#' @param new_sign integer vector with values -1L and 1L, corresponding to the 
#' sign of the value (one-sample t test) or difference of values (paired
#' samples t test). See \code{\link{osTtestRandomSigns}}.
#' @param attribs logical value whether dimension and type attributes shall be
#' attached to the returned statistic (default: FALSE)
#' @keywords internal
#' @return The functions return a vector (or array, if \code{attribs = TRUE}) 
#' of the same length as the number of columns in the input matrix. Auxiliary 
#' results ('Df', degrees of freedom; 'p_value', p-values) are added as 
#' attributes (if \code{obj$verbose == TRUE}).
NULL

#' @rdname Ttest
# independent samples t-test
isTtest <- function(obj, new_group = NULL, attribs = FALSE) {
    verbose <- if (!attribs) FALSE else obj$verbose 
    groups <- if (is.null(new_group)) obj$groups else new_group
    res <- lapply(1:2, function(g) {
        datx <- obj$.arraydat[groups == g, , drop = FALSE]
        mm <- colMeans(datx, na.rm = TRUE)
        vv <- .Call("rowVars", datx, as.integer(dim(datx)), 
                    TRUE, obj$has_NA, FALSE, 
                    PACKAGE = "matrixStats")
        if (obj$has_NA) {
            nn <- nrow(datx) - colCounts(datx, value = NA_real_)
            nn[nn < 1L] <- NA
        } else {
            nn <- nrow(datx)
        }
        list(mm, vv, nn)
    })
    m1 <- res[[1]][[1]]; v1 <- res[[1]][[2]]; n1 <- res[[1]][[3]]
    m2 <- res[[2]][[1]]; v2 <- res[[2]][[2]]; n2 <- res[[2]][[3]]
    mdiff <- m1 - m2
    out <- 
        if (!obj$var_equal) {
            sd12 <- sqrt(v1/n1 + v2/n2)
            mdiff/sd12
        } else {
            sd12 <- sqrt(
                ((n1-1)*v1 + (n2-1)*v2)/(n1 + n2 - 2)
            )
            es <- mdiff/sd12
            es/sqrt(1/n1 + 1/n2)
        }
    if (attribs) {
        setattributes(out, obj, 
                      if (obj$var_equal) "Student's t statistic" 
                      else "Welch's t statistic")
        if (verbose) {
            # degress of freedom
            Df <- 
                if (obj$var_equal) {
                    n1 + n2 - 2 
                } else {
                    sd12^4/( (v1/n1)^2/(n1-1) + (v2/n2)^2/(n2-1) )
                }
            setattributes(Df, 
                          if (length(Df) > 1L) obj else NULL, 
                          "Degrees of freedom")
            # p-value
            p <- 2 * pt(-abs(out), Df)
            setattributes(p, obj, "Traditional P-value")
            # Cohen's D
            es <- if (!obj$var_equal) abs(mdiff)/sqrt((v1 + v2)/2) else abs(es)
            setattributes(es, obj, "Cohen's D")
            # attach to out
            setattr(out, "p_value", p)
            setattr(out, "effect_size", es)
            setattr(out, "Df", Df)
        }
    }
    out
}

#' @rdname Ttest
# one-sample or paired samples t-test
osTtest <- function(obj, new_sign = NULL, attribs = FALSE) {
    verbose <- if (!attribs) FALSE else obj$verbose
    dat <- if (is.null(new_sign)) obj$.arraydat else obj$.arraydat * new_sign
    m1 <- colMeans(dat, na.rm = TRUE)
    sd1 <- .Call("rowVars", dat, as.integer(dim(dat)), 
                 TRUE, obj$has_NA, FALSE, 
                 PACKAGE = "matrixStats")
    sd1 <- sqrt(sd1)
    if (obj$has_NA) {
        n1 <- nrow(dat) - colCounts(dat, value = NA_real_)
        n1[n1 < 1L] <- NA
    } else {
        n1 <- nrow(dat)
    }
    es <- (m1 - obj$mu)/sd1
    out <- es/sqrt(n1)
    if (attribs) {
        setattributes(
            out, obj$teststat_dims, 
            if (obj$type == "paired_samples") "Paired samples t-statistic"
            else "One sample t-statistic")
        if (verbose) {
            Df <- n1 - 1
            setattributes(Df,
                          if (length(Df) > 1L) obj else NULL,
                          "Degrees of freedom")
            pvalues <- 2 * pt(-abs(out), Df)
            setattributes(pvalues, obj$teststat_dims, "Traditional P-value")
            es <- abs(es)
            setattributes(es, obj$teststat_dims, "Cohen's D")
            setattr(out, "p_value", pvalues)
            setattr(out, "effect_size", es)
            setattr(out, "Df", Df)
        }
    }
    out
}

#' Create random group memberships for independent samples t-tests
#' 
#' This function creates a random permutation of the original group
#' memberships for independent samples t-tests. Care is taken so that all 
#' permutations are unique.
#' @param obj a list object returned by \code{\link{preTtest}}; 
#' the \code{obj$groups} element is relevant indicating the original group
#' memberships
#' @param nperm integer indicating the number of permutations
#' @param seed an integer value which specifies a seed (default: NULL), or a 
#' list of arguments passed to \code{\link{set.seed}}
#' @return \code{isTtestRandomGroups} returns a matrix with as many rows as
#' the length of 'group' and as many columns as 'nperm'.
#' @keywords internal
#' @examples
#' \dontrun{
#' # this should fail: too many permutations on a small sample
#' temp <- try(isTtestRandomGroups(list(groups = rep(1:2, each = 4L), 1e3)), 
#'             silent = TRUE)
#' stopifnot(inherits(temp, "try-error"))
#' 
#' # all permutations should be unique
#' temp <- isTtestRandomGroups(list(groups = rep(1:2, c(6, 4))), 500)
#' stopifnot(!anyDuplicated(temp, MARGIN = 2L))
#' temp <- isTtestRandomGroups(list(groups = rep(1:2, each = 20L)), 1e4)
#' stopifnot(!anyDuplicated(temp, MARGIN = 2L))
#' }
isTtestRandomGroups <- function(obj, nperm, seed = NULL) {
    setSeed(seed)
    groups <- obj$groups
    grouplen <- length(groups)
    group2len <- sum(groups == 2L)
    maxperm <- choose(grouplen, group2len)
    if (nperm > maxperm) 
        stop("nperm is too high for this sample size (not enough unique combinations exist)")
    if (nperm > maxperm/3L) {
        randind <- combn(grouplen, group2len)
        randgroups <- matrix(1L, grouplen, ncol(randind))
        for (i in 1:ncol(randind)) randgroups[randind[,i], i] <- 2L
        randgroups <- randgroups[, sample.int(nrow(randgroups), nperm)]
    } else {
        randgroups <- matrix(groups, grouplen, 1)
        repeat{
            nc <- ncol(randgroups) - 1
            if (nc >= nperm) {
                randgroups <- randgroups[, 2:(nperm + 1)]
                break
            } else {
                rg <- matrix(1L, grouplen, (nperm - nc) * 2)
                for (i in 1:ncol(rg)) {
                    rg[sample.int(grouplen, group2len), i] <- 2L
                }
                randgroups <- cbind(randgroups, rg)
                randgroups <- fastUnique(randgroups, margin = 2L)
            }
        }
    }
    randgroups
}


#' Create random signs for one-sample and paired samples t-tests
#' 
#' This function creates a random sample of -1L and 1L values which can be
#' used in randomization runs of one-sample and paired samples t-tests. 
#' Care is taken so that all samples are unique.
#' @param obj a list object returned by \code{\link{preTtest}}; the number of
#' rows (let's denote here 'id_length') of \code{obj$.arraydat} is relevant
#' @param nperm integer indicating the number of permutations
#' @param seed an integer value which specifies a seed (default: NULL), or a 
#' list of arguments passed to \code{\link{set.seed}}
#' @return \code{osTtestRandomSigns} returns an integer matrix (containing 
#' only -1L and 1L values) with 'id_length' rows and 'nperm' columns.
#' @examples
#' \dontrun{
#' # this should fail: too many permutations for this sample size
#' obj <- list(.arraydat = matrix(8))
#' temp <- try(osTtestRandomSigns(obj, 1e4), silent = TRUE)
#' stopifnot(inherits(temp, "try-error"))
#' 
#' # expect an integer matrix containing unique series of -1L and 1L
#' obj <- list(.arraydat = matrix(10))
#' temp <- osTtestRandomSigns(obj, 500)
#' stopifnot(range(temp), c(-1L, 1L))
#' stopifnot(length(unique(as.vector(temp))) == 2L)
#' stopifnot(!anyDuplicated(temp, MARGIN = 2L))
#' 
#' obj <- list(.arraydat = matrix(20))
#' temp <- osTtestRandomSigns(obj, 1e4)
#' stopifnot(!anyDuplicated(temp, MARGIN = 2L))
#' }
osTtestRandomSigns <- function(obj, nperm, seed) {
    setSeed(seed)
    idlen <- nrow(obj$.arraydat)
    maxperm <- 2^idlen
    if (nperm > maxperm) 
        stop("nperm is too high for this sample size (not enough unique combinations exist)")
    if (nperm > maxperm/3) {
        randi <- matrix(unlist(
            expand.grid(rep(list(c(-1L, 1L)), idlen), 
                        KEEP.OUT.ATTRS = FALSE),
            use.names = FALSE), idlen, maxperm, TRUE)
        randi <- randi[, sample(ncol(randi), nperm)]
    } else {
        randi <- matrix(rep.int(1L, idlen))
        repeat{
            nc <- ncol(randi) - 1
            if (nc >= nperm) {
                randi <- randi[, 2:(nperm + 1)]
                break
            } else {
                ri <- matrix_(sample(c(-1L, 1L), 
                                     temp <- idlen * (nperm - nc) * 2,
                                     TRUE),
                              idlen, temp)
                randi <- fastUnique(cbind(randi, ri), margin = 2L)
            }
        }
    }
    randi
}


#' Parameter checks and data preparation in \code{arrayTtest}
#' 
#' \code{preTtest} performs parameter checking, combines .arraydat and 
#' .arraydat2 (if not NULL), reshapes .arraydat to be a matrix with subjects in 
#' rows, and returns a list with the data and all parameters. It is called by 
#' \code{\link{arrayTtest}}.
#' @inheritParams arrayTtest
#' @param tfce logical value if TFCE is requested in \code{arrayTtest} 
#' @param perm logical value if permutation is requested in \code{arrayTtest}
#' @keywords internal
#' @seealso \code{\link{arrayTtest}} is the higher-level function exported to 
#' the user. \code{\link{isTtest}} and \code{\link{osTtest}} are internal 
#' functions relying on the returned object of \code{preTtest}.
preTtest <- function(.arraydat, .arraydat2, paired, groups, 
                     mu, var_equal, id_dim, verbose, tfce, perm) {
    # check if .arraydat is an array
    if (is.data.frame(.arraydat)) .arraydat <- as.matrix(.arraydat)
    # TODO: should be possible to allow missing values, but what about
    # the permutation?
    assertArray(.arraydat, mode = "numeric", min.d = 2L, any.missing = FALSE,
                .var.name = ".arraydat")
    # dimnames
    pre_dimnames <- dimnames(.arraydat)
    full_dimnames <- fillMissingDimnames(dimnames(.arraydat), dim(.arraydat))
    full_dimid_origord <- names(full_dimnames)
    if (!identical(pre_dimnames, full_dimnames))
        dimnames(.arraydat) <- full_dimnames
    # check id_dim
    if (!id_dim %in% names(full_dimnames)) {
        stop(sprintf(".arraydat has no '%s' dimension identifier", id_dim))
    }
    # if TFCE is requested, "chan" and "time" dimensions are obligatory
    if (tfce && !all(c("chan", "time") %in% names(full_dimnames)))
        stop("If TFCE is requested, '.arraydat' must have 'chan' and 'time' dimensions")
    # check for .arraydat2
    if (is.null(.arraydat2)) {
        if (!is.null(groups)) {
            # set type of t-test
            type <- "independent_samples"
            # check 'groups'
            assertAtomic(groups, .var.name = "groups")
            # subset arraydat for is.na(groups)
            groups_keep <- !is.na(groups)
            if (any(!groups_keep)) {
                .arraydat <- subsetArray(.arraydat, 
                                         listS(.id_dim = groups_keep),
                                         drop = FALSE)
                groups <- groups[groups_keep]
            }
            # transform to 1s and 2s
            groups <- as.integer(factor_(groups))
            if (max(groups) != 2L) 
                stop("groups must contain two distinct values (e.g. 1 or 2, 'a' or 'b')")
        } else {
            # set type of t-test
            type <- "one_sample"
        }
    } else {
        # check format
        if (is.data.frame(.arraydat2)) .arraydat2 <- as.matrix(.arraydat2)
        assertArray(.arraydat2, mode = "numeric", any.missing = FALSE, 
                    min.d = 2L, .var.name = ".arraydat2")
        # check dimension names
        full_dimnames2 <- fillMissingDimnames(dimnames(.arraydat2),
                                              dim(.arraydat2))
        # check 'paired'
        assertLogical(paired, any.missing = FALSE, len = 1L, 
                      .var.name = "paired")
        # if paired == TRUE, compute the difference; 
        # if paired == FALSE, bind the two arrays and store the group membership
        # in 'groups'
        if (paired) {
            if (!identical(full_dimnames, full_dimnames2) ||
                !identical(dim(.arraydat), dim(.arraydat2))) {
                stop("Input array dimensions must be identical")
            }
            .arraydat <- .arraydat2 - .arraydat
            groups <- NULL
            type <- "paired_samples"
        } else {
            if (!identical(sort(names(full_dimnames)),
                           sort(names(full_dimnames2)))) {
                stop(".arraydat and .arraydat2 must have the same dimension identifiers")
            }
            groups <- rep.int(1:2, 
                              c(length(full_dimnames[[id_dim]]),
                                length(full_dimnames2[[id_dim]])))
            .arraydat <- bindArrays(.arraydat, .arraydat2, along_name = id_dim)
        }
    }
    # check 'mu'
    if (length(mu == 1L)) {
        assertNumber(mu, finite = TRUE, .var.name = "mu")
    } else {
        assertArray(mu, mode = storage.mode(.arraydat), any.missing = FALSE,
                    .var.name = "mu")
        expect_mu_dimnames <- full_dimnames[-match(id_dim, names(full_dimnames))]
        expect_mu_dims <- vapply(expect_mu_dimnames, length, integer(1L))
        if (!identical(dim(mu), expect_mu_dims)) {
            stop(sprintf("If 'mu' is not a scalar, it must be an array with the same dimensions (%s) as '.arraydat' without the 'id_dim' dimension",
                 paste(expect_mu_dims, collapse = "x")))
        }
        if (!is.null(dimnames(mu)) && 
            !identical(fillMissingDimnames(dimnames(mu), dim(mu)), 
                       expect_mu_dimnames)) {
            stop("If 'mu' has dimension names, they must be the same as the dimension names of '.arraydat' without the 'id_dim' dimension")
        }
    }
    #
    # reshape .arraydat (first permute dimensions if needed, finally to matrix)
    #
    # id_dim, chan, time dimension order
    nr_chanXtime <- nr_chan <- nr_time <- NULL
    if (tfce || perm) {
        firstdims <- union(c(id_dim, "chan", "time"), names(full_dimnames))
        .arraydat <- apermArray(.arraydat, first = firstdims)
        if (length(mu) > 1L) {
            mu <- apermArray(.arraydat, first = setdiff(firstdims, id_dim))    
        }
        if (tfce) {
            nr_chan <- length(full_dimnames$chan)
            nr_time <- length(full_dimnames$time)
            nr_chanXtime <- nr_chan * nr_time
        } else {
            nr_chanXtime <- c(length(full_dimnames$chan),
                              length(full_dimnames$time))
            nr_chanXtime[nr_chanXtime == 0L] <- 1L
            nr_chanXtime <- prod(nr_chanXtime)
        }
    } else {
        .arraydat <- apermArray(.arraydat, first = id_dim)
    }
    full_dimnames <- dimnames(.arraydat)
    full_dimid <- names(full_dimnames)
    full_dims <- dim(.arraydat)
    setattr(full_dims, "names", full_dimid)
    teststat_dimid <- full_dimid[full_dimid != id_dim]
    setattr(teststat_dimid, "origpos", 
            match(teststat_dimid, full_dimid_origord))
    has_NA <- anyNA(.arraydat)
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
    matrix_(.arraydat, nrow = nrow(.arraydat))
    #
    # return
    #
    list(.arraydat = .arraydat, groups = groups,  
         has_NA = has_NA,
         full_dimnames = full_dimnames, full_dims = full_dims,
         pre_dimnames = pre_dimnames,
         teststat_dimid = teststat_dimid,
         type = type, mu = mu, var_equal = var_equal, verbose = verbose,
         nr_chan = nr_chan, nr_time = nr_time, nr_chanXtime = nr_chanXtime,
         otherdims = otherdims)
}


#' Point-to-point t-tests (potentially with TFCE correction) on arrays
#' 
#' \code{arrayTtest} performs point-to-point t-tests on arrays. 
#' Permutation-based p-values and Threshold-free Cluster Enhancement (TFCE) 
#' correction can be requested.
#' @param .arraydat a numeric array with named dimnames containing EEG (or 
#' other) data. Missing values are not allowed. Must have at least three
#' dimensions with names "chan", "time", and "id" (see also id_dim)
#' @param .arraydat2 a numeric array with named dimnames containing EEG (or 
#' other) data (default: NULL). If provided, see the parameter \code{paired} 
#' for running paired or independent-samples t-tests
#' @param paired logical scalar, only used if .arraydat2 is provided. If paired 
#' is FALSE (default), the function computes independent samples t-tests, 
#' otherwise paired samples t-tests are performed
#' @param groups provides an alternative (and more efficient) way to perform 
#' independent samples t-tests; a character, factor, or integer vector which 
#' defines group membership. Groups is ignored if .arraydat2 is not missing. 
#' NA values code subjects to drop.
#' @param mu a numeric scalar indicating the true value of the mean (or 
#' difference between means, if two-sample tests are performed)
#' @param var_equal a logical scalar whether the variances are equal (only 
#' relevant for independent-samples t-tests). If TRUE, the pooled 
#' variance is used to estimate the variance. If FALSE (default), the Welch (or 
#' Satterthwaite) approximation to the degrees of freedom is used.
#' @param id_dim name of the dimension which identifies the subjects 
#' (default: "id")
#' @param verbose logical value indicating if p-values should be computed for 
#' the traditional t-test results (default: TRUE)
#' @param nperm integer value giving the number of permutations (default: 0L).
#' If \code{nperm < 2L}, no permutation is performed.
#' @param tfce either 1) NULL (the default) or FALSE (the same as NULL), both of 
#' which mean no TFCE correction, or 2) TRUE, which means TFCE correction with 
#' default parameters, or 3) an object as returned by \code{\link{tfceParams}} 
#' with custom TFCE parameters (see Examples and also \code{\link{tfceParams}}).\cr 
#' Custom parameters can be also provided by \code{tfce = .(key = value)} to 
#' save typing (this works by calling \code{\link{tfceParams}} with the given 
#' parameters).
#' @param parallel either 1) FALSE (the default), which results in single-core 
#' computation, or 2) TRUE, which means parallelization with default parameters, 
#' or 3) an object as returned by \code{\link{parallelParams}} with custom 
#' parameters (see Examples and also \code{\link{parallelParams}}), or 4) NULL,
#' which means that the registered backend (if there is any) shall be used.\cr 
#' Custom parameters can be also provided by \code{parallel = .(key = value)} to 
#' save typing (this works by calling \code{\link{parallelParams}} with the 
#' given parameters).
#' @param seed an integer value which specifies a seed (default: NULL), or a 
#' list of arguments passed to \code{\link{set.seed}}
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
#'     result_eegr <- arrayTtest(erps, groups = dat_id$group)
#' )
#' 
#' # the built-in R function (t.test) provides the same result, but running it on
#' # the full dataset would take much more time; 
#' # here we take a subsample of the data
#' sub <- list(chan = "F4", time = "200", stimclass = "B", pairtype = "ident")
#' result_ttest <- t.test(subsetArray(erps, sub) ~ dat_id$group)
#' 
#' # to check that they are equivalent, we have to remove the attributes
#' eegr_t <- as.vector(subsetArray(extract(result_eegr, "stat"), sub))
#' eegr_p <- as.vector(subsetArray(extract(result_eegr, "p"), sub))
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
#' result_tfce <- arrayTtest(tempdat, groups = dat_id$group, 
#'                           nperm = 499L, 
#'                           parallel = .(ncores = 2L),
#'                           tfce = .(ChN = ChN))
#' 
#' # 4) compare the traditional and TFCE-corrected p-values
#' p_trad <- extract(result_tfce, "p")
#' p_tfce <- extract(result_tfce, "p_corr")
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
arrayTtest <- function(.arraydat, .arraydat2 = NULL, paired = FALSE, 
                       groups = NULL, mu = 0, var_equal = FALSE, id_dim = "id", 
                       verbose = TRUE, nperm = 0L, tfce = NULL, parallel = NULL, 
                       seed = NULL) {
    # deparse tfce and parallel
    mcall <- match.call() 
    tfce <- argumentDeparser(substitute(tfce), "tfceParams")
    ob <- getDoBackend()
    parallel <- argumentDeparser(substitute(parallel), "parallelParams")
    if (is.logical(parallel) && !parallel) {
        registerDoSEQ()
    } else if (inherits(parallel, "parallelParams") && parallel$cl_new) {
        on.exit(stopCluster(parallel$cl))
    }
    on.exit(setDoBackend(ob), add = TRUE)
    #
    # compute statistics
    out <- arrayTtestAnova(
        "T_TEST", .arraydat, .arraydat2,
        paired = paired, groups = groups, 
        mu = mu, var_equal = var_equal, id_dim = id_dim, 
        verbose = verbose, nperm = nperm, 
        tfce = tfce, parallel = parallel, seed = seed)
    #
    # replace call to the original
    out$call <- mcall
    # set class
    setattr(out, "class", "arrayTtest")
    # return
    out
}

