#
# <<< ERP preprocessing functions >>> --------
#

# import functions =========== 

#' Import options for importBVdat
#' 
#' \code{importOptions} allows to set the raw data importation parameters
#' @param eeg_ext character value, the extension of the segmented EEG dataset 
#' file (default: "dat")
#' @param marker_ext character value, the extension of the raw marker file
#' (default: "vmrk")
#' @param info_ext character value, the extension of the information file 
#' (default: "vhdr")
#' @param marker_skip numeric value; number of rows to be skipped while 
#' importing the marker file (default: 14)
#' @param marker_segment character vector; name(s) of stimuli defining the
#' the segments (default: "Target")
#' @param marker_badcat character string(s) identifying bad segments in the 1st 
#' column of the marker file (default: "Bad Interval")
#' @param marker_badstim character string(s) identifying bad stimuli or 
#' responses in the 2nd column of the marker file (default: "resp_false")
#' @param marker_keepstim character vector identifying which stimuli or responses 
#' should be kept by checking the 2nd column of the marker file (default: "*")
#' @param marker_badchan character vector of length 1 or 2; markers identifying 
#' bad intervals of channels. If marker_badchan has only one element, "_start" 
#' and "_stop" are automatically appended to the end of the string.
#' @param marker_time0 character string identifying the marker name at time 0 
#' in the imported EEG dataset
#' @param segment_dpoints numeric vector of data point indices defining a 
#' segment (default: -100:1023)
#' @param marker_regexp logical value; should marker_segment/badcat/.../time0 
#' strings be handled as regular expressions (default) or treated as they are
#' @param marker_ignorecase logical value; should the case of marker_segment/
#' .../time0 definitions be ignored (default)
#' @param marker_header,marker_fill,marker_asis,marker_sep parameters 
#' to be passed to \code{\link{read.table}} while importing the marker file
#' @export
#' @return A list with named parameters
importOptions <- function(eeg_ext = "dat", marker_ext = "vmrk", info_ext = "vhdr",
                          marker_skip = 14, marker_segment = "Target", 
                          marker_badcat = "Bad Interval", 
                          marker_badstim = "resp_false",
                          marker_keepstim = "*", 
                          marker_badchan = c("badchan_start",
                                             "badchan_stop"),
                          marker_time0 = "Time 0",
                          segment_dpoints = (-100):1023,
                          marker_regexp = TRUE,
                          marker_ignorecase = TRUE,  
                          marker_header = FALSE, marker_fill = TRUE, 
                          marker_asis = TRUE, marker_sep = ",") {
    #
    if (length(marker_badchan) == 1) {
        marker_badchan <- paste(marker_badchan, c("start", "stop"), sep = "_")
    }
    # return
    mget(ls())
}

#' Import binary file exported from BrainVision
#' 
#' \code{importBVdat} imports a binary file exported from the BrainVision 
#' software
#' @param file_name character string; the name of the input files without
#' extensions
#' @param file_path character string; the path to the files if they are not in 
#' the working directory (default)
#' @param id character string denoting the identification code of the 
#' participant
#' @param import_options a list, which should be given by calling 
#' \code{\link{importOptions}}
#' @note This is a custom function tailored for the special datasets collected 
#' in our lab. Use it with extra care for general purposes!
#' @export
#' @return A list with three named elements: eeg (array), markers (data.frame), 
#' channels (data.frame)
importBVdat <- function(file_name, file_path = getwd(), id = "",
                        import_options = importOptions()) {
    message(paste("\nImport ", file_name, "...", sep = ""))
    mygrepl <- function(patterns, ...) {
        if (is.null(patterns)) {
            rep(FALSE, length(list(...)$x))
        } else {
            rowSums(sapply(patterns, 
                           grepl, 
                           ignore.case = import_options$marker_ignorecase, 
                           fixed = !import_options$marker_regexp,
                           ...)) > 0
        }
    }
    extractInfo <- function(x, type = "char") {
        out <- tolower(strsplit(info[grep(x, info)], "=")[[1]][2])
        if (type == "num") out <- as.numeric(out)
        return(out)
    }
    chanInfo <- function() {
        chan <- strsplit(info[grep("Ch.{1,3}=", info)], "=")
        chan_names <- sapply(chan[1:(length(chan)/2)],
                             function(x) strsplit(x[[2]], ",")[[1]][1])
        chan_pos <- t(sapply(chan[-(1:(length(chan)/2))],
                             function(x) as.numeric(strsplit(x[[2]], ",")[[1]])))
        colnames(chan_pos) <- c("r", "theta", "phi")
        return( data.frame(chan_pos, row.names = chan_names) )
    }
    #
    #assignList(import_options, verbose = FALSE)
    #
    # info from vhdr file
    con1 <- file(file.path(file_path, paste(file_name, 
                                            import_options$info_ext, 
                                            sep = ".")))
    info <- readLines(con1)
    close(con1)
    eeg_orientation <- extractInfo("DataOrientation=")
    nr_chan <- extractInfo("NumberOfChannels=", "num")
    nr_tpoints <- extractInfo("SegmentDataPoints=", "num")
    eeg_length <- extractInfo("DataPoints=", "num") * nr_chan
    eeg_Hz <- 1e6 / extractInfo("SamplingInterval=", "num")
    chan <- chanInfo()
    #
    # import eeg data
    if (!is.null(import_options$eeg_ext)) {
        eeg <- readBin(file.path(file_path, 
                                 paste(file_name, 
                                       import_options$eeg_ext, sep = ".")),
                       what = "integer", n = eeg_length)
        eeg <- eeg / 1000
        if (eeg_orientation == "vectorized") {
            eeg <- t(matrix_(eeg, eeg_length/nr_chan, nr_chan))
        } else {
            matrix_(eeg, nr_chan, eeg_length/nr_chan)
        }
        setattr(eeg, "subject_id", as.character(id))
        procstep <- list(
            what = "import", call = match.call(),
            file_name = file_name, file_path = file_path,
            options = import_options)
    } else {
        eeg <- NULL
    }
    #
    # import markers
    markers.orig <- with(import_options, read.table(
        file.path(file_path, 
                  paste(file_name, marker_ext, sep = ".")),
        header = marker_header, 
        fill = marker_fill, 
        as.is = marker_asis, 
        sep = marker_sep, 
        skip = marker_skip))
    markers.orig[,3] <- suppressWarnings(as.numeric(markers.orig[, 3]))
    markers.orig[,5] <- suppressWarnings(as.integer(markers.orig[, 5]))
    markers.orig <- markers.orig[!is.na(markers.orig[, 3]), 1:5]
    markers.orig$segmind <- findInterval(
        markers.orig[,3],
        markers.orig[grepl("New Segment", 
                           markers.orig[, 1]), 3])
    if (eeg_length/nr_chan/nr_tpoints  != max(markers.orig$segmind))
        stop("Marker and data files do not match!")
    #
    # remove data in bad channel intervals
    if (!is.null(import_options$marker_badchan)) {
        chan_names <- rownames(chan)
        bch <- markers.orig[, 2]
        ind <- markers.orig[, 5] > 0
        bch[ind] <- paste(bch[ind],
                          chan_names[markers.orig[ind, 5]],
                          sep = "_")
        bch <- paste0(bch, "_")
        badchanFn <- function(i) {
            indstart <- grep(paste(import_options$marker_badchan[1], ".*_", 
                                   chan_names[i], "_", sep = ""), 
                             bch)
            indstart <- unique(markers.orig[indstart, 3])
            indstop <- grep(paste(import_options$marker_badchan[2], ".*_", 
                                  chan_names[i], "_", sep = ""), 
                            bch)
            indstop <- unique(markers.orig[indstop, 3])
            len1 <- length(indstart)
            len2 <- length(indstop)
            if (len1 == 0 & len2 == 0) {
                return(NULL)
            } else {
                lendiff <- len1 - len2
                if (lendiff == 1) {
                    message(
                        paste0("One badchan_stop marker was missing at channel ", chan_names[i], 
                               " and was automatically set to the last sampling point.")
                    )
                    indstop <- c(indstop, ncol(eeg))
                } else if (lendiff == -1) {
                    message(
                        paste0("One badchan_start marker was missing at channel ", chan_names[i],
                               " and was automatically set to the first sampling point.")
                    )
                    indstart <- c(1, indstart)    
                } else if (abs(lendiff) > 1) {
                    stop(
                        paste0("Badchan_start and badchan_stop markers do not match at channel ", chan_names[i], ".")
                    )
                }
                return(cbind(start = indstart, stop = indstop))
            }
            if (any((indstop - indstart) < 0)) {
                stop(
                    paste0("Badchan_start and badchan_stop markers do not match at channel ", chan_names[i], ".")
                )
            }
        }
        badchans <- lapply(seq_along(chan_names), badchanFn)
        names(badchans) <- chan_names
        badchans <- badchans[!sapply(badchans, is.null)]
        if (!is.null(import_options$eeg_ext)) {
            for (i in names(badchans)) {
                eeg[i, unlist(mapply(seq, 
                                     badchans[[i]][, "start"], 
                                     badchans[[i]][, "stop"], 
                                     SIMPLIFY = FALSE), 
                              use.names = FALSE)] <- NA
            }
        }
        procstep$bad_channels <- badchans
    } 
    # 
    # remove bad segments
    markers.orig$badsegm <- rep(
        sapply(
            split(markers.orig, markers.orig$segmind), 
            function(x) 
                any(mygrepl(import_options$marker_badcat, x = x[, 1]) | 
                        mygrepl(import_options$marker_badstim, x = x[, 2]))
        ),
        table(markers.orig$segmind))
    keepind <- with(import_options,
                    mygrepl(import_options$marker_segment, x = markers.orig[, 1]) & 
                        mygrepl(import_options$marker_keepstim, x = markers.orig[, 2]) &
                        rep(markers.orig[mygrepl("Time 0", x = markers.orig[, 1]), 3] %in% 
                                markers.orig[mygrepl(import_options$marker_time0, 
                                                     x = markers.orig[, 1]), 3],
                            table(markers.orig$segmind)) &
                        !markers.orig$badsegm)
    markers <- markers.orig[keepind, ]
    markers <- data.frame(
        segment = sapply(strsplit(markers[, 1], "="), "[[", 2),
        fullcode = markers[, 2],
        dpoint = markers[, 3],
        segmind = markers$segmind
    )
    #
    # format eeg data
    if (!is.null(import_options$eeg_ext)) {
        eeg <- eeg[, c(outer(import_options$segment_dpoints, 
                             markers$dpoint,"+"))]
        array_(eeg, 
               c(nr_chan, length(import_options$segment_dpoints), 
                 nrow(markers)),
               list(chan = rownames(chan),
                    time = import_options$segment_dpoints * 1000 / eeg_Hz,
                    trial = paste(markers$segment, 
                                  markers$fullcode, 
                                  sep = "_")))
    }
    # decorate
    setattr(eeg, "processing_steps", list(procstep))
    # return
    list(eeg = eeg, markers = markers, channels = chan)
}

#' Split concatenated strings. 
#' 
#' \code{splitMarker} is a conveniance function wrapping strsplit. It returns
#' a data.frame with columns of customizable classes. Useful for post-processing
#' markers.
#' @param marker a vector which contains concatenated elements
#' @param header a character vector which determines the column names of the
#' resulting data.frame
#' @param type a character vector which determines the class of each column in
#' the resulting vector. If its length is less then the number of columns, it
#' will be recycled. If set to NULL (default), character vectors are transformed to
#' factors.
#' @param splitchar a character vector of length one indicating the splitting
#' character (default: _). Splitting characters which are special characters in
#' R (e.g. "|", ".", etc.) should be given as in strsplit (e.g. "\\\|")
#' @export
#' @return A data.frame with header and appropriate classes, containing the 
#' spltted substrings of the original vector elements
#' @seealso \code{\link{strsplit}}
splitMarker <- function(marker, header, type = NULL, splitchar = "_") {
    out <- data.frame(
        matrix(unlist(strsplit(as.character(marker), splitchar), 
                      use.names = FALSE), 
               nrow = length(marker), ncol = length(header), byrow = TRUE))
    colnames(out) <- header
    if (!is.null(type)) {
        if (length(type) < length(header)) {
            type <- rep(type, length.out = length(header))
        }
        for (i in 1:length(type)) {
            if (is.factor(out[, i])) {
                out[, i] <- do.call(paste("as", type[i], sep = "."), 
                                    list(as.character(out[, i])))
            } else {
                out[, i] <- do.call(paste("as", type[i], sep = "."), list(out[, i]))
            }
        }
    }
    # return
    out
}

# base functions =========== 

#' Baseline correction 
#' 
#' \code{baselineCorr} performs baseline correction
#' @param dat numeric array containing the ERPs
#' @param along_dims numeric or character vector identifying the dimensions of 
#' dat along which separate baseline averaging should occur (default: "chan")
#' @param by_dims numeric or character vector identifying the dimensions of 
#' dat across which separate baseline averaging should occur (default: NULL).
#' See \code{\link{scaleArray}} for further explanation.
#' @param base_time numeric or character vector identifying time points which
#' form the baseline. If NULL (default), dat must have named dimension names,
#' and \code{base_time} is the vector of levels of the time dimension which 
#' are below 0.
#' @export
#' @return A numeric array with the same attributes as dat
#' @seealso \code{\link{scaleArray}} for the function behind the scenes
#' @examples
#' # example dataset
#' data(erps)
#' 
#' # remove baseline activity separately for each stimulus class, pairtype and 
#' # channel in each subject; now we use the by_dims argument because it needs 
#' # less typing
#' bc <- baselineCorr(erps, by_dims = "time")
#' 
#' # remove baseline activity separately for each channel in each subject;
#' # now we choose the along_dims argument 
#' bc <- baselineCorr(erps, along_dims = c("chan", "id"))
baselineCorr <- function(dat, along_dims = "chan", by_dims = NULL, 
                         base_time = NULL) {
    message("\n****\nPerform baseline correction ... ", appendLF = FALSE)
    origattr <- attributes(dat)
    if (is.null(along_dims)) {
        if (is.null(by_dims)) {
            stop("Both along_dims and by_dims are NULL")
        } else if (is.character(by_dims)) {
            by_dims <- match(by_dims, names(dimnames(dat)))
        }
        along_dims <- setdiff(seq_along(dim(dat)), by_dims)
    }
    if (is.null(base_time)) base_time <- as.numeric(dimnames(dat)$time)<0
    out <- scaleArray(dat, along_dims = along_dims, 
                      center_subset = list(time = base_time),
                      scale = FALSE)
    attributes(out) <- origattr
    setattr(out, "processing_steps",
            c(attr(out, "processing_steps"),
              list(list(what = "baseline correction", 
                        call = match.call(), along_dimensions = along_dims, 
                        base_time = base_time))))
    message("Done")
    # return
    out
}

#' Options for artifact rejection
#' 
#' \code{artrejOptions} allows to set the parameters of the artifact rejection
#' methods.
#' @param sampling_freq numeric value, the sampling frequency of the EEG-data
#' @param channels character vector containing the name or index of channels
#' which are subject to artifact rejection. If set to "all" (default), all 
#' channels are included.
#' @param apply_maxgrad logical value, if set to TRUE (default), the maximum 
#' gradient criterion is applied.
#' @param maxgrad_limit numeric value, the maximum gradient / millisecond 
#' (default: 50)
#' @param maxgrad_mark numeric vector of length 2; the placement of the Bad
#' Interval mark in milliseconds before and after the occurence of maxgrad_limit
#' violation (default: c(-200, 200))
#' @param apply_diffrange logical value, if set to TRUE (default), the difference 
#' range criterion is applied.
#' @param diffrange_limit numeric vector of length 2, the minimum and maximum
#' voltage difference in a given interval (default: 200)
#' @param diffrange_mark numeric vector of length 2; the placement of the Bad
#' Interval mark in milliseconds before and after the occurence of diffrange_limit
#' violation (default: c(-200, 200))
#' @param diffrange_interval numeric value, the length of interval for the 
#' difference range criterion in milliseconds (default: 200)
#' @param apply_amplrange logical value, if set to TRUE (default), the amplitude 
#' range criterion is applied.
#' @param amplrange_limit numeric vector of length 2, the minimum and maximum 
#' voltage in the whole segment (default: c(-200, 200))
#' @param amplrange_mark numeric vector of length 2; the placement of the Bad
#' Interval mark in milliseconds before and after the occurence of 
#' amplrange_limit violation (default: c(-200, 200))
#' @details The short definitions of the possible artifact rejection criteria 
#' are as follows:
#' \itemize{
#' \item{Maximum gradient:}{The absolute difference between the voltages 
#' measured at successive milliseconds.}
#' \item{Difference range:}{The minimum and maximum difference between the 
#' maximum and minimum voltages in a given sampling interval.}
#' \item{Amplitude range:}{The minimum and maximum voltages in the segments.}
#' }
#' @note The algorithm takes care of the sampling frequency for all parameters
#' which are provided in milliseconds (or /ms) and makes adjustments if needed.
#' However, *_mark parameters are not used since only segmented data can be
#' analyzed in the present version of artifactRejection().
#' @export
#' @return A list object with all parameters.
artrejOptions <- function(
    sampling_freq = 1000, channels = "all", 
    apply_maxgrad = TRUE, maxgrad_limit = 50, maxgrad_mark = c(-200, 200),
    apply_diffrange = TRUE, diffrange_limit = c(0.5, 100), 
    diffrange_mark = c(-200, 200), diffrange_interval = 200,
    apply_amplrange = TRUE, amplrange_limit = c(-200, 200), 
    amplrange_mark = c(-200, 200)) {
    #
    opt <- mget(setdiff(ls(), "sampling_freq"))
    freqmod <- sampling_freq/1000
    ind <- grep("mark|interval|limit", names(opt))
    lapply(opt[ind], function(x) is.numeric(x))
    opt[ind] <- lapply(opt[ind], sort)
    ind <- grep("mark|interval", names(opt))
    opt[ind] <- lapply(opt[ind], "*", freqmod)
    opt$maxgrad_limit <- opt$maxgrad_limit / freqmod
    # return
    opt
}

#' Artifact rejection
#' 
#' \code{artifactRejection} performs artifact rejection on segmented data.
#' @param dat numeric array (EEG-data) with the following named dimensions 
#' (dimension order does not matter): chan, time, trial
#' @param markers if not NULL (default), a matrix or data.frame containing the 
#' characteristics of the trials (markers)
#' @param artrej_options a named list containing the parameters for the 
#' artifact rejection criteria. See \code{\link{artrejOptions}} for details.
#' @param return_data logical value, if TRUE (default), dat and markers without
#' rejected trials are returned
#' @param return_details logical value, if TRUE (default), the full array of 
#' results (e.g., bad trials for each channel and for each criterion) is 
#' returned as an attribute of bad_trials (see Values section)
#' @param print_result logical value, if TRUE (default), a summary of the 
#' results is printed to the console
#' @export
#' @return A named list containing bad_trials (trials identified with artifacts)
#' and the modified input data (dat and markers without contaminated trials)
artifactRejection <- function(dat, markers = NULL, artrej_options = artrejOptions(),
                              return_data = TRUE, return_details = TRUE,
                              print_result = TRUE) {
    maxgradFn <- function(x, 
                          maxgrad_limit = artrej_options$maxgrad_limit,
                          maxgrad_mark = artrej_options$maxgrad_mark) {
        colAnys( abs(x[-1,]-x[-nrow(x),]) > maxgrad_limit )
    }
    diffrangeFn <- function(x, 
                            diffr_int = artrej_options$diffrange_interval,
                            diffr_lim = artrej_options$diffrange_limit) {
        ind <- rollFun(x, diffr_int, max, endrule = "NA") - 
            rollFun(x, diffr_int, min, endrule = "NA")
        ind <- ind[!rowAlls(is.na(ind)), ]
        #
        colAnys( abs(ind) < diffr_lim[1] ) |
            colAnys( abs(ind) > diffr_lim[2] )
    }
    amplrangeFn <- function(x, 
                            amplr_lim = artrej_options$amplrange_limit,
                            amplr_mark = artrej_options$amplrange_mark) {
        colMaxs(x) > amplr_lim[2] |
            colMins(x) < amplr_lim[1] 
    }
    # main function
    aRej <- function(x, details = return_details) {
        #assignList(artrej_options, verbose = FALSE)
        crits <- sub("apply_", "", 
                     names(artrej_options)[grep("apply", 
                                                names(artrej_options))])
        if (length(crits) == 0) {
            stop("No criterion to apply; check the apply_* parameters in artrej_options!")
        }
        attribs <- attr(x, "array_attributes")
        row_dim <- attribs$row_dim
        dims <- attribs$dim[-row_dim]
        dimn <- attribs$dimnames[-row_dim]
        names(dims) <- names(dimn)
        out <- matrix_(FALSE, ncol(x), length(crits))
        colnames(out) <- crits
        message("\n****\nStart artifact rejection / Criterion: ...\n")
        if (artrej_options$apply_maxgrad) {
            message("... Maximum gradient - ", appendLF = FALSE)
            out[, "maxgrad"] <- maxgradFn(x)
            message(" done\n")
        }
        if (artrej_options$apply_diffrange) {
            message("... Difference range - ", appendLF = FALSE)
            out[, "diffrange"] <- diffrangeFn(x)
            message(" done\n")
        }
        if (artrej_options$apply_amplrange) {
            message("... Amplitude range - ", appendLF = FALSE)
            out[, "amplrange"] <- amplrangeFn(x)
            message(" done\n")
        }
        array_(out, c(dims, ncol(out)), c(dimn, list(crit = crits)))
        artrej_summary <- matrix_(0, dims["chan"]+1, length(crits)+1,
                                  list(chan = c(dimn$chan, "all"),
                                       crit = c(crits, "all")))
        dimres <- dim(artrej_summary)
        artrej_summary[-dimres[1], -dimres[2]] <- avgDims(out, "trial")
        artrej_summary[dimres[1], -dimres[2]] <- 
            colMeans(apply(out, c("trial", "crit"), any))
        artrej_summary[-dimres[1], dimres[2]] <- 
            colMeans(apply(out, c("trial", "chan"), any))
        out.details <- out
        out <- apply(out, "trial", any)
        artrej_summary[dimres[1], dimres[2]] <- mean( out )
        setattr(out, "summary", artrej_summary)
        if (details) setattr(out, "details", out.details)
        return( out )
    }
    # input data check
    if (is.null(dimnames(dat)) || is.null(names(dimnames(dat))) ||
            !identical(sort(names(dimnames(dat))), 
                       c("chan", "time", "trial"))) {
        stop("Provide EEG data as an array with the following named 
             dimensions (dimension order does not matter): 
             chan, time, trial.")
    }
    if (is.null(markers)) {
        markers <- data.frame(fullcode = seq_along(dimnames(dat)$trial))
    }
    if (length(dimnames(dat)$trial) != nrow(markers)) {
        stop("EEG and marker data contain different number of trials!")
    }
    if (!is.list(artrej_options)) {
        stop("The artrej_options parameter must be a list; provide it through 
             artrejOptions() to avoid inconsistent results!")
    }
    keepchan <- 
        if (identical(artrej_options$channels, "all")) {
            dimnames(dat)$chan
        } else {
            artrej_options$channels
        }
    # run artifact rejection
    tempdat <- array2mat(subsetArray(dat, list(chan = keepchan), drop. = FALSE),
                         "time", keep_dimnames = FALSE)
    badtrials <- aRej(tempdat)
    if (print_result) {
        cat("\n----- Proportion of bad trials -----\n")
        print(attr(badtrials, "summary"))
        cat("\n------------------------------------\n")
    }
    if (return_data) {
        dat <- subsetArray(dat, list(trial = which(!badtrials)))
        setattr(dat, "processing_steps",
                c(attr(dat, "processing_steps"),
                  list(list(
                      what = "artifact rejection",
                      call = match.call(), results = badtrials, 
                      options = artrej_options)))
        )
        markers <- droplevels( markers[!badtrials, ] )
        # return
        list(bad_trials = badtrials, eeg = dat, markers = markers)
    } else {
        setattr(badtrials, "options", artrej_options)
        # return
        list(bad_trials = badtrials, eeg = NULL, markers = NULL)
    }
    }

#' Compute Global Field Power
#'
#' \code{compGfp} computes Global Field Power (the standard deviation of 
#' channel values for each sampling point)
#' @param dat numeric matrix or array, usually with named dimensions (one of 
#' which is "chan")
#' @param keep_channels logical value; if TRUE, the original channels are 
#' retained, and the GFP values are added with channel code "GFP"
#' (default = FALSE)
#' @param channel_dim a character value or numeric index indicating the
#' channel dimension of \code{dat} (default: "chan")
#' @export
#' @return The function returns a matrix or an array. Note that if
#' the original channels are not retained, the channel dimension is dropped.
compGfp <- function(dat, keep_channels = FALSE, channel_dim = "chan") {
    out <- fnDims(dat, channel_dim, colSds, vectorized = TRUE)
    if (keep_channels) {
        dimn <- names(dimnames(dat))
        if (!is.null(dimn) && is.numeric(channel_dim))
            channel_dim <- dimn[channel_dim]
        out <- 
            if (is.character(channel_dim)) {
                bindArrays(dat, GFP = out, along_name = channel_dim)
            } else {
                bindArrays(dat, GFP = out, along = channel_dim)
            }
        if (is.null(dimnames(dat)))
            setattr(out, "dimnames", NULL)
    }
    # return
    out
}

#' Scale channels
#' 
#' \code{scaleChan} normalizes the data across channels so that for
#' each sampling point, the mean of the channel amplitudes is zero and
#' the standard deviation is one.
#' @param dat numeric matrix or array with names dimensions, one of which
#' must be "chan"
#' @param keep_dimorder logical value; if TRUE (default), the order of the 
#' dimensions is kept intact, otherwise the channel dimension will be 
#' the first dimension in the resulting matrix or array
#' @export
#' @return A numeric matrix or array
scaleChan <- function(dat, keep_dimorder = TRUE) {
    scaleArray(dat, by_dims = "chan", keep_dimorder = keep_dimorder)
}

#' Compute centroids
#'
#' \code{centroid} computes the centroids (separately for the negative and
#' positive values)
#' @param dat numeric vector of amplitudes at a given sampling point
#' @param ch_pos matrix or data.frame containing channel positions (in the same
#' order as dat), see \code{\link{coordinates}}
#' @param proj2map logical value; if TRUE (default), the centroids are projected
#' onto a 2D plane
#' @param proj_unitsphere logical value; if TRUE, and proj2map is also TRUE, the 
#' projection assumes unit radius
#' @param ... additional parameters to \code{\link{project3dMap}}
#' @export
#' @return A data.frame containing the positions of the negative and positive 
#' centroids 
centroid <- function(dat, ch_pos, proj2map = TRUE, proj_unitsphere = FALSE, ...) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("x", "y", "z") %in% colnames(ch_pos))) {
        ch_pos <- sph2cart(ch_pos)
    }
    ind.neg <- dat < 0
    ind.pos <- dat >= 0
    centroids <- as.data.frame(
        rbind(colMeans(ch_pos[ind.neg, c("x", "y", "z")]),
              colMeans(ch_pos[ind.pos, c("x", "y", "z")]))
    )
    rownames(centroids) <- c("negative", "positive")
    if (proj2map) {
        pos <- cart2sph(centroids)
        if (proj_unitsphere) {
            pos$r <- 1
        }
        out <- project3dMap(pos, ...)
    }
    # return
    out
}

#' Average single-trials
#' 
#' \code{avgTrials} performs single-trial averaging
#' @param dat numeric array, the segmented ERPs
#' @param markers data.frame containing marker definitions
#' @param which_factors numeric or character vector indicating which columns
#' from markers should be used for averaging. If NULL (default), all columns
#' are used.
#' @export
#' @return numeric array
avgTrials <- function(dat, markers, which_factors = NULL) {
    stopifnot(is.numeric(dat))
    stopifnot(!is.null(markers))
    message("\n****\nAverage single-trials ... ", appendLF = FALSE)
    if (is.character(which_factors)) {
        if (!all(which_factors %in% colnames(markers))) {
            stop("markers and which_factors do not match")
        }
    } else if (is.numeric(which_factors)) {
        if (!all(which_factors %in% (1:ncol(markers)))) {
            stop("markers and which_factors do not match")
        }
    }
    if (!is.null(which_factors)) markers <- markers[, which_factors]
    groups <- interaction(markers, drop = TRUE, sep = "|")
    out <- fnDims(dat, "trial", 
                  function(x, g) sweep(rowsum(x, g), 1, table(g), "/"), 
                  list(g = groups), vectorized = TRUE, 
                  newdims = list(factor_level = levels(groups)))
    tempfac <- strsplit(dimnames(out)$factor_level, "\\|")
    setattr(out, "processing_steps",
            c(attr(out, "processing_steps"),
              list(list(
                  what = "averaging",
                  call = match.call(),
                  factors = as.data.frame(
                      matrix(unlist(tempfac, use.names = FALSE), 
                             nrow = length(tempfac), 
                             ncol = length(tempfac[[1]]), 
                             byrow = TRUE,
                             dimnames = list(1:length(tempfac), 
                                             colnames(markers))))
              ))
            )
    )
    message("Done")
    # return
    out
}
