#
# <<< plotting functions >>> -----------
#

#' Plot ERP curves
#' 
#' \code{plotERParray} is a generalization of \code{matplot} onto array inputs.
#' @param dat a numeric array with named dimensions
#' @param xdim character value; the name of the dimension of dat which defines 
#' the x axis (default: "time")
#' @param sepdim character value; the name of the dimension of dat which 
#' separates the lines (default: "chan")
#' @param title character value; the title of the plot
#' @param subtitle.col the colour of the subtitles
#' @param gfp_plot logical value; if TRUE (default), the GFP curves are also 
#' plotted
#' @param gfp_col the colour of the GFP curves (if plotted)
#' @param gfp_lwd the thickness of the GFP curves (if plotted)
#' @param minus_up logical value; if set to TRUE, minus values are plotted 
#' upwards. If set to FALSE, minus values are plotted downwards. If NULL 
#' (default), it is set to TRUE or FALSE depending on the gfp_plot parameter. 
#' NULL
#' @param grid_labels character vector; the name of the x and y axes
#' @param grid_dim integer vector giving the number of rows and columns
#' in which the plots are arranged; if set to NULL (default), the arrangement 
#' of the plot is set up automatically
#' @param ... additional parameters to be passed to \code{matlines}
#' @export
plotERParray <- function(dat, xdim = "time", sepdim = "chan", 
                         title = "", subtitle.col = "black",
                         gfp_plot = TRUE, gfp_col = "black", gfp_lwd = 1.3, 
                         minus_up = NULL, grid_labels = c("time", "ampl"), 
                         grid_dim = NULL, ...) {
    emptyplot <- function() {
        plot(0, 0, xlim = c(-1,1), type = "n", axes = FALSE, 
             frame.plot = FALSE, xlab = "", ylab = "")
    }
    if (length(dim(dat)) > 3) {
        dat <- mergeDims(dat, setdiff(names(dimnames(dat)), c(xdim, sepdim)))
    }
    wrapdim <- setdiff(names(dimnames(dat)), c(xdim, sepdim))
    dat <- aperm(dat, c(xdim, sepdim, wrapdim))
    subtitle.col <- rep(subtitle.col, length_out = dim(dat)[3])
    if (gfp_plot) gfpdat <- compGfp(dat)
    x <- as.numeric(as.character(dimnames(dat)[[1]]))
    xrange <- if (is.null(list(...)$xlim)) range(x) else list(...)$xlim 
    if (is.null(list(...)$ylim)) {
        yr <- range(dat)
        yrange <- mean(yr) + c(-1, 1)*(yr[2]-mean(yr))*1.02
        if (is.null(minus_up)) {
            minus_up <- if (!gfp_plot) TRUE else FALSE
        }
        if (minus_up) yrange <- -yrange
    } else {
        yrange <- list(...)$ylim
    } 
    if (is.null(grid_dim)) grid_dim = rep(ceiling(sqrt(dim(dat)[3])), 2)
    layoutmat <- cbind(
        c(0, rep(2, grid_dim[1]), 0, 0),
        rbind(rep(1, grid_dim[2]),
              matrix(1:prod(grid_dim), grid_dim[1], grid_dim[2], TRUE) + 3,
              rep(0, grid_dim[2]),
              rep(3, grid_dim[2])),
        c(rep(0, grid_dim[1] + 3)))
    layout(layoutmat, widths = c(0.3, rep(1, grid_dim[2]), 0.3),
           heights = c(0.3, rep(1, grid_dim[1]), 0.3, 0.3))
    par(mar = rep(0, 4))
    emptyplot(); text(0, 0, title, cex = 1.5)
    emptyplot(); text(0, 0, grid_labels[2], cex = 1.3, srt = 90)
    emptyplot(); text(0, 0, grid_labels[1], cex = 1.3)
    par(mar = c(0, 0, 2, 0))
    for (i in 1:dim(dat)[3]) {
        matplot(x, dat[,,i], xlim = xrange, ylim = yrange, yaxs = "i",  
                axes = FALSE, frame.plot = TRUE, type = "n")
        grid()
        mtext(dimnames(dat)[[3]][i], 3, cex = 0.7, col = subtitle.col[i])
        matlines(x, dat[,,i], axes = FALSE, ...)
        if (i == max(dim(dat)[3])) {
            axis(1)
            axis(4)
        }
        if (gfp_plot) lines(x, gfpdat[, i], col = gfp_col, lwd = gfp_lwd)
    }
}


#' Plot the topography of the signals
#' 
#' \code{plot2dview} creates 2D topoplots
#' @param dat a numeric vector or matrix. If a matrix is provided, it must 
#' contain the participants in rows and the channels in columns
#' @param ch_pos a data.frame with the channel coordinates
#' @param r a numeric value; the radius of the head (not used if ch_pos contains
#' a column \code{r})
#' @param timepoint an integer value
#' @param ampl_range numeric vector, the minimum and maximum value of the 
#' amplitudes (default: c(-5, 5))
#' @param resol integer value, the resolution of the projection (default: 100L)
#' @param resolcol integer value, the resolution of the colours (default: 1000L)
#' @param projection character value (default = "laea"). See 
#' \url{http://www.remotesensing.org/geotiff/proj_list/} for common projections
#' @param projref projection reference (pole [default] or equator)
#' @param origo a named character vector of lat and long coordinates of the 
#' origo
#' @param plot_centroid should centroids be plotted (default: TRUE)
#' @param centroid_circle numeric scalar (default=0.5), the ratio of individual
#' centroid coordinates which fall into the area of a shaded circle 
#' @param centroid_unitsphere should centroids be projected onto the 2d plane as
#' if they were located on the surface of a unit sphere (default: false)
#' @param plot_bar put color bar on the plot (default: TRUE)
#' @param plot_ch plot channels as dots (default: TRUE)
#' @param plot_chnames plot channel names (default: TRUE)
#' @param title title of the plot
#' @param gfp a numeric vector of the GFP values in the whole time window
#' @param gfp_max the upper limit if gfp is provided
#' @param color_scale function which creates a color scale with \code{n} values
#' (\code{n} is defined by \code{resolcol}). The default is to use the 
#' \code{\link[gplots]{colorpanel}} function, with blue > white > red 
#' color scale.
#' @param ... arguments to \code{\link{chanInterp}}
#' @export
#' @examples
#' # example data
#' data(erps)
#' chan_pos <- attr(erps, "chan")   # channel positions
#' 
#' # plot the grand average scalp distribution at 200 ms, for the stimulus 
#' # class "A" and pairtye "ident" condition; display also the GFP curve of the 
#' # whole epoch
#' plotdat <- transformArray(y ~ chan + time, erps,
#'                           subset = list(stimclass = "A", pairtype = "ident"),
#'                           datfr = FALSE)
#' gfpdat <- compGfp(plotdat)
#' plot2dview(subsetArray(plotdat, list(time = "200")),
#'            chan_pos, timepoint = 200, gfp = gfpdat)
plot2dview <- function(dat, ch_pos, r = 1, timepoint = NULL, 
                       ampl_range = NULL, 
                       resol = 100L, resolcol = 1000L,
                       projection = "laea", projref = c("pole", "equator"),
                       origo = c(lat = ifelse(projref == "pole", 90, 0),
                                 long = ifelse(projref == "pole", 270, 0)),
                       plot_centroid = TRUE, centroid_circle = 0.5, 
                       centroid_unitsphere = FALSE,
                       plot_bar = TRUE, plot_ch = TRUE, plot_chnames = TRUE, 
                       title = NULL,
                       gfp = NULL, gfp_max = NULL, 
                       color_scale = function(n) gplots::colorpanel(n, "blue", "white", "red"), 
                       ...) {
    #
    colors <- do.call(match.fun(color_scale), list(resolcol))
    ch_pos <- as.data.frame(ch_pos)
    if (all(c("x", "y", "z") %in% colnames(ch_pos))) {
        ch_pos <- ch_pos[, c("x", "y", "z")]
        posgeo <- cart2geo(ch_pos)
    } else if (all(c("theta", "phi") %in% colnames(ch_pos))) {
        if (!"r" %in% colnames(ch_pos)) ch_pos$r <- r
        ch_pos <- sph2cart(ch_pos)
        posgeo <- cart2geo(ch_pos)
    } else if (all(c("long", "lat") %in% colnames(ch_pos))) {
        if (!"r" %in% colnames(ch_pos)) ch_pos$r <- r
        posgeo <- ch_pos[, c("long", "lat", "r")]
        ch_pos <- geo2cart(ch_pos)
    } else {
        stop("Wrong coordinate object")
    }
    r <- posgeo$r
    if (is.data.frame(dat)) dat <- as.matrix(dat)
    if (is.matrix(dat)) {
        if (all(c("chan", "id") %in% names(dimnames(dat))))
            dat <- aperm(dat, c("id", "chan"))
        if (ncol(dat) != nrow(ch_pos)) 
            stop("Wrong data format! Should be: Rows = participants, Columns = channels")
        subjdat <- dat
        dat <- colMeans(dat, na.rm = TRUE)
    } else if (length(dim(dat)) > 2L) {
        stop("The input data may not have more than 2 dimensions")
    } else {
        centroid_circle <- NA
    }
    #
    projref <- match.arg(projref)
    #
    boundarypos <- 
        if (projref == "pole") {
            data.frame(
                long = seq(-180, 180, length.out = 180),
                lat = rep(max(-45/(1.11), min(posgeo$lat)), 180))
        } else {
            data.frame(
                long = c(rep(max(abs(posgeo$long)), 90),
                         rep(-max(abs(posgeo$long)), 90)),
                lat = c(seq(-90, 90, length.out = 90), 
                        seq(90, -90, length.out = 90)))
        }
    boundarypos <- unique(project3dMap(boundarypos, projection = projection, 
                                       projref = projref, origo = origo))
    #
    xlim <- range(boundarypos$x) * 1.1
    ylim <- range(boundarypos$y) * 1.1
    ylim.cm <- ylim.cp <- 0.1
    if (plot_bar) {
        ypos.bar <- c(0, ylim[1] * (1 + ylim.cm))
        ylim.cm <- ylim.cm + 0.1
        ypos.bar[1] <- 0.75 * ypos.bar[2] + 0.25 * ylim[1] * (1 + ylim.cm)    
    }
    if (!is.null(gfp)) {
        gfp <- drop(gfp)
        if (length(dim(gfp)) > 1L) {
            stop("The gfp argument can not have more than one dimension")
        }
        ypos.gfp <- c(ylim[2] * (1 + ylim.cp * 1.1), 0)
        ylim.cp <- ylim.cp + 0.3
        ypos.gfp[2] <- 0.1 * ypos.gfp[1] + 0.9 * ylim[2] * (1 + ylim.cp)
        if (is.null(gfp_max)) gfp_max <- max(gfp) * 1.05
        gfp <- ypos.gfp[1] + gfp * diff(ypos.gfp) / gfp_max
    } else {
        ypos.gfp <- c(ylim[2], ylim[2])
    }
    if (!is.null(title) | !is.null(timepoint)) {
        ylim.cp <- ylim.cp + 0.2
        ypos.title <- 0.3 * ypos.gfp[2] + 0.7 * ylim[2] * (1+ylim.cp)
    }
    ylim[1] <- ylim[1] * (1 + ylim.cm)
    ylim[2] <- ylim[2] * (1 + ylim.cp)
    #
    gridx <- seq(min(boundarypos$x), max(boundarypos$x), length.out = resol)
    gridy <- seq(min(boundarypos$y), max(boundarypos$y), length.out = resol)
    gridpos <- expand.grid(x = gridx, y = gridy)
    ind <- in.polygon(gridpos$x, gridpos$y, 
                      boundarypos$x * 1.1, boundarypos$y * 1.1)
    gridgeo <- project3dMap(gridpos[ind,], projection = projection, 
                            projref = projref, origo = origo, inverse = TRUE)
    gridcart <- geo2cart(gridgeo)
    z <- matrixIP(NA_real_, resol, resol)
    z[ind] <- chanInterp(dat, ch_pos, gridcart, ...)
    if (is.null(ampl_range)) ampl_range <- range(z, na.rm = TRUE)
    z[z > ampl_range[2]] <- ampl_range[2]
    z[z < ampl_range[1]] <- ampl_range[1]
    par(mar = rep(0, 4))
    image(gridx, gridy, z, useRaster = TRUE, col = colors,
          xlim = xlim, ylim = ylim,
          zlim = ampl_range, xlab = "", ylab = "", axes = FALSE)
    gridx <- seq(min(boundarypos$x), max(boundarypos$x), 
                 length.out = 1000) * 1.1
    gridy <- seq(min(boundarypos$y), max(boundarypos$y), 
                 length.out = 1000) * 1.1
    gridpos <- expand.grid(x = gridx, y = gridy)
    ind <- !in.polygon(gridpos$x, gridpos$y, 
                       boundarypos$x*1, boundarypos$y*1)
    z <- matrixIP(NA_integer_, 1000, 1000)
    z[ind] <- 1L
    image(gridx, gridy, z, useRaster = TRUE, col = "white", add = TRUE)
    lines(boundarypos, col = "white", lwd = 1)
    #
    if (plot_ch) {
        ch_pos_xy <- project3dMap(ch_pos)
        points(ch_pos_xy,, pch = 20)
        if (plot_chnames) {
            indleft <- ch_pos_xy$x <= 0
            indright <- ch_pos_xy$x > 0
            text(ch_pos_xy[indleft,], , rownames(ch_pos[indleft, ]), pos = 4)
            text(ch_pos_xy[indright,], , rownames(ch_pos[indright, ]), pos = 2)
        }
    }
    #
    if (plot_bar) {
        xpos.bar <- seq(min(boundarypos$x / 3), max(boundarypos$x / 3), 
                        length.out = resolcol)
        image(xpos.bar, ypos.bar,
              matrixIP(seq(ampl_range[1], ampl_range[2], length.out = resolcol), 
                       resolcol, 2),
              zlim = ampl_range, col = colors, add = TRUE)
        text(min(xpos.bar), ypos.bar[1], 
             substitute(paste(k, " ", mu, "V", sep = ""), 
                        list(k = formatC(ampl_range[1], 1, format = "f"))),
             pos = 2)
        text(max(xpos.bar), ypos.bar[1], 
             substitute(paste(k, " ", mu, "V", sep = ""), 
                        list(k = formatC(ampl_range[2], 1, format = "f"))), 
             pos = 4)
    }
    if (!is.null(title) | !is.null(timepoint)) {
        if (is.null(title)) {
            title <- paste("Time:", timepoint, "ms")
        }
        text(0, ypos.title, title, cex = 1.1)
    }
    if (!is.null(gfp)) {
        normFn <- function(x, xrange, lims) 
            (x - xrange[1])/diff(xrange) * diff(lims) + lims[1]
        lgfp <- length(gfp)
        gfp.x <- as.numeric(names(gfp))
        gfp.xr <- range(gfp.x)
        gfp.xlims <- xlim * 0.8
        gfp.xticks <- pretty(gfp.x)
        gfp.xticks <- gfp.xticks[gfp.xticks >= gfp.xr[1] & 
                                     gfp.xticks <= gfp.xr[2]]
        gfp.xticks.normed <- normFn(gfp.xticks, gfp.xr, gfp.xlims)
        gfp.x <- normFn(gfp.x, gfp.xr, gfp.xlims)
        gfp.xr <- range(gfp.x)
        rect(gfp.xr[1], ypos.gfp[1], gfp.xr[2], ypos.gfp[2], 
             border = NA, col = "grey60")
        axis(1, at = gfp.xticks.normed, 
             labels = gfp.xticks, 
             pos = ypos.gfp[1], cex.axis = 0.7, padj = -1)
        lines(gfp.x, gfp, col = "white")
        smpl <- 
            if (is.null(names(gfp))) {
                which(as.numeric(dimnames(gfp)[[1]]) == timepoint)
            } else {
                which(as.numeric(names(gfp)) == timepoint)
            }
        if (length(smpl) == 1) {
            points(gfp.x[smpl], gfp[smpl], pch = 16, col = "red", cex = 1.1)
        }
    }
    if (plot_centroid) {
        c1 <- centroid(dat, ch_pos, proj_unitsphere = centroid_unitsphere)
        if (!is.na(centroid_circle)) {
            if (reqFn("plotrix")) {
                cc <- apply(subjdat, 1, centroid, ch_pos, 
                            proj_unitsphere = centroid_unitsphere)
                temp <- matrixIP(unlist(lapply(cc, function(x) abs(x-c1))), 
                                 ncol = length(cc))
                rads <- apply(temp, 2, quantile, 
                              probs = centroid_circle, na.rm = TRUE)
                plotrix::draw.ellipse(x = c1$x, y = c1$y, 
                                      a = rads[c(1, 3)], b = rads[c(2, 4)], 
                                      col = c(rgb(0, 0, 0.5, 0.15), rgb(0.5, 0, 0, 0.15)), 
                                      border = NA)
            }
        }
        points(c1["neg", "x"], c1["neg", "y"], col = "white", 
               pch = 15, cex = 2.5)
        points(c1["neg", "x"], c1["neg", "y"], col = colors[1], 
               pch = 15, cex = 2)
        points(c1["pos", "x"], c1["pos", "y"], col = "white", 
               pch = 15, cex = 2.5)
        points(c1["pos", "x"], c1["pos", "y"], col = colors[length(colors)], 
               pch = 15, cex = 2)
    }
}


#' Plot topography in complex layout
#' 
#' \code{complexplot2dview} plots the scalp topographies of different groups or
#' conditions by calling \code{\link{plot2dview}} on a list of vectors or 
#' matrices
#' @param dat a list of numeric vectors or matrices, which can be plotted by
#' \code{\link{plot2dview}}
#' @param ch_pos channel position matrix
#' @param timepoint integer; the timepoint at which the scalp distributions are
#' compared. Can be a vector referring to the timepoints of each topomap.  
#' @param datgrid arrangement of the grid (number of rows X columns)
#' @param layout_matrix,heights,widths arguments to be passed to 
#' \code{\link{layout}}. If given, it overrides \code{datgrid}.
#' @param title_row,title_col character vectors of the titles of rows and 
#' columns, respectively. If NULL and dat has dimnames attribute, 
#' \code{rownames(dat)} and \code{colnames(dat)} are used.
#' @param title_main character; the title of the plot (default: the timepoint(s))
#' @param plot_bar logical; should a colour bar be plotted below the maps 
#' (default: TRUE)
#' @param ampl_range numeric vector of the amplitude range. If not provided,
#' the range is computed from the input data.
#' @param gfp a list of Global Field Power data along the whole epoch
#' @param ... additional arguments to be passed to \code{\link{plot2dview}}
#' @seealso \code{\link{plot2dview}}
#' @export
#' @examples
#' # example dataset
#' data(erps)
#' chan_pos <- attr(erps, "chan")   # channel positions
#' 
#' # plot the grand average scalp distribution at 180 ms, for each stimulus 
#' # class and pairtye condition; display also the GFP curves of the 
#' # whole epochs
#' timepoint <- 180
#' plotdat <- splitArray(avgDims(erps, "id"), 
#'                       c("stimclass", "pairtype"), 
#'                       drop = TRUE)
#' gfpdat <- lapply(plotdat, compGfp)
#' plotdat[] <- lapply(plotdat, subsetArray,
#'                     list(time = as.character(timepoint)))
#' complexplot2dview(plotdat,
#'                   chan_pos, timepoint, 
#'                   gfp = gfpdat, plot_ch = FALSE, 
#'                   resol = 50L, resolcol = 50L)
#
complexplot2dview <- function(dat, ch_pos, timepoint, 
                              datgrid = NULL, layout_matrix = NULL, 
                              heights = NULL, widths = NULL, 
                              title_row = NULL, title_col = NULL, 
                              title_main = paste("Time:", 
                                                 paste(timepoint, 
                                                       collapse = ", "), 
                                                 "ms"), 
                              plot_bar = TRUE, ampl_range = NULL, 
                              gfp = NULL, ...) {
    #
    emptyplot <- function() plot(0, 0, xlim = c(-1, 1), type = "n", 
                                 axes = FALSE, frame.plot = FALSE, 
                                 xlab = "", ylab = "")
    resolcol <- list(...)$resolcol
    if (is.null(resolcol)) resolcol <- formals(plot2dview)$resolcol
    if (is.null(layout_matrix)) {
        if (is.null(datgrid)) {
            datgrid <- c(floor(sqrt(length(dat))), 
                         ceiling(sqrt(length(dat))))
        }
        rownum <- datgrid[1]
        colnum <- datgrid[2]
        datnum <- prod(datgrid)
        layout_matrix <- rbind(
            matrix(c(
                rep(1, colnum), 0,
                2:(colnum + 1), 0), 2, colnum + 1, TRUE), 
            matrix(c(
                (colnum + rownum + 1) + 1:datnum,
                (colnum + 2):(colnum + rownum + 1)), rownum, colnum + 1))
        layout_matrix <- rbind(layout_matrix, 
                               c(rep(max(layout_matrix) + 1, colnum), 0))
        if (is.null(heights)) 
            heights <- c(0.3, 0.2, rep(1, rownum), 0.2)
        if (is.null(widths)) 
            widths <- c(rep(1, colnum), 0.2)
    } else {
        if (is.null(heights)) heights <- rep(1, nrow(layout_matrix))
        if (is.null(widths)) widths <- rep(1, ncol(layout_matrix))
    }
    if (is.null(title_row)) title_row <- rownames(dat)
    if (is.null(title_col)) title_col <- colnames(dat)
    layout(layout_matrix, widths = widths, heights = heights)
    par(mar = rep(0, 4))
    emptyplot(); text(0, 0, title_main, cex = 1.5)
    for (i in 1:length(unique(layout_matrix[2, -1]))) {
        emptyplot()
        text(0, 0, title_col[i], cex = 1.3)
    }
    for (i in 1:length(unique(layout_matrix[-c(1, 2, nrow(layout_matrix)),1]))) {
        emptyplot()
        text(0, 0, title_row[i], cex = 1.3, srt = -90L)
    }    
    if (is.null(ampl_range)) {
        ampl_range <- c(min(unlist(dat, use.names = FALSE)), 
                        max(unlist(dat, use.names = FALSE)))
        ampl_range <- round(c((1 - 0.2 * sign(ampl_range)[1]) * ampl_range[1],
                              (1 + 0.2 * sign(ampl_range)[2]) * ampl_range[2]), 
                            1)
    }
    if (length(unlist(gfp, use.names = FALSE)) == 1 && gfp == TRUE) {
        gfp <- lapply(dat, compGfp)    
    }
    gfp_max <- if (is.null(gfp)) NULL else max(unlist(gfp, use.names = FALSE))
    timepoint <- rep_len(timepoint, length(dat))
    for (i in 1:length(dat)) {
        plot2dview(dat[[i]], ch_pos = ch_pos, timepoint = timepoint[i], 
                   gfp = gfp[[i]], gfp_max = gfp_max,
                   plot_bar = FALSE, ampl_range = ampl_range, title = "", ...)
    }
    if (plot_bar) {
        xpos.bar <- seq(0.4, 0.6, length.out = resolcol)
        ypos.bar <- c(0.45, 0.55)
        image(xpos.bar, ypos.bar,
              matrixIP(seq(ampl_range[1], ampl_range[2], length.out = resolcol), 
                       resolcol, 2),
              xlim = c(0, 1), ylim = c(0, 1), zlim = ampl_range, 
              useRaster = TRUE, col = bluered(resolcol), 
              xlab = "", ylab = "", axes = FALSE 
        )
        text(min(xpos.bar), ypos.bar[1], 
             substitute(paste(k, " ", mu, "V", sep = ""), 
                        list(k = ampl_range[1])), pos = 2)
        text(max(xpos.bar), ypos.bar[1], 
             substitute(paste(k, " ", mu, "V", sep = ""), 
                        list(k = ampl_range[2])), pos = 4)
    }
}

# plot p-values with across-country similarity matrices 
reprPlot <- function(p, e, title = "Reproducibility", plim = c(0.05, 0.01), 
                     plim_shading = TRUE, plim_outline = 0.05,
                     corr_range = c(-0.99, 0.99), 
                     corr_method = c("pearson", "spearman")) {
    #
    emptyplot <- function() {
        plot(0, 0, xlim = c(-1, 1), type = "n", 
             axes = FALSE, frame.plot = FALSE, xlab = "", ylab = "")
    }
    pPolygon <- function(x, y, direction = "h") {
        plim <- sort(unique(c(0, plim, 1)))
        y_cat <- findInterval(y, plim)
        y_phase <- c(1, cumsum(abs(diff(y_cat))) + 1)
        yy <- 0.04 + -log(y) / (-log(p_min)) * 0.15
        cols <- grey(c(seq(0.4, 0.6, 
                           length.out = (length(plim)-2)), 0.95))
        if (direction == "h") {
            for (i in unique(y_phase)) {
                ind <- y_phase == i
                polygon(c(x[ind], rev(x[ind])), 
                        c(-yy[ind], rep(-0.04, sum(ind))),
                        col = cols[y_cat[ind][1]], border = "grey20")
                if (plim_shading && y_cat[ind]<(length(plim)-1)) {
                    acol <- col2rgb(cols[y_cat[ind][1]], TRUE) / 255
                    acol[4] <- 0.2
                    rect(x[min(which(ind))], 0,
                         x[max(which(ind))], 1,
                         col = rgb(acol[1], acol[2], acol[3], acol[4]), 
                         border = rgb(0, 1, 0, acol[4]), lty = "1F")
                }
            }
            if (!is.na(plim_outline)) {
                sig <- rep(NA_integer_, length(x))
                sig[y <= plim_outline] <- 1L
                lines(x, sig, col = "green", lwd = 3)
            }
        } else {
            for (i in unique(y_phase)) {
                ind <- y_phase == i
                polygon(c(-yy[ind], rep(-0.04, sum(ind))), 
                        c(x[ind], rev(x[ind])),
                        col = cols[y_cat[ind][1]], border = "grey20")
                if (plim_shading && y_cat[ind]<(length(plim)-1)) {
                    acol <- col2rgb(cols[y_cat[ind][1]], TRUE) / 255
                    acol[4] <- 0.2
                    rect(0, x[min(which(ind))],
                         1, x[max(which(ind))],
                         col = rgb(acol[1], acol[2], acol[3], acol[4]), 
                         border = rgb(0, 1, 0, acol[4]), lty = "1F")
                }
            }
            if (!is.na(plim_outline)) {
                sig <- rep(NA_integer_, length(x))
                sig[y <= plim_outline] <- 1L
                lines(sig, x, col = "green", lwd = 3)
            }
        }
    }                     
    implot <- function(corrmat, pp) {
        zlim <- atanh(corr_range)
        mat <- atanh(corrmat)
        mat[mat < zlim[1]] <- zlim[1]
        mat[mat > zlim[2]] <- zlim[2]
        image(mat, zlim = zlim, xlim = c(-0.2, 1), ylim = c(-0.2, 1),
              axes = FALSE, xlab = "", ylab = "", col = bluered(1000))
        lines(c(0, 1), c(0, 1), col = "grey70")
        tpoints <- as.numeric(rownames(corrmat))
        seqtpoints <- tpoints[seq(1, length(tpoints), 25)]
        at <- (seqtpoints - min(tpoints)) / (max(tpoints) - min(tpoints))
        x <- (tpoints - min(tpoints)) / (max(tpoints) - min(tpoints))
        for (ii in 1:2) {
            axis(ii + 2, at = at, labels = tpoints[tpoints %in% seqtpoints],
                 cex.axis = 0.8)
            y <- pp[, ii]
            pPolygon(x, y, c("h", "v")[ii])
        }
    }
    plim <- sort(unique(c(0, plim, 1)))
    corr_method <- match.arg(corr_method)
    e <- aperm(e, c("chan", "time", "nation"))
    p <- aperm(p, c("time", "nation"))
    p_min <- min(p)
    natnum <- ncol(p)
    nations <- colnames(p)
    imlayout <- matrixIP((2 * natnum + 1), natnum, natnum)
    diagpos <- diag(matrixIP(1:natnum^2, natnum, natnum))
    for (i in 1:length(imlayout)) {
        imlayout[i] <- 
            if (i %in% diagpos) 0 
        else max(imlayout)+1
    }
    layout_matrix <- 
        rbind(c(0, rep(1, natnum)),
              c(0, rep(2:(natnum + 1))),
              cbind((natnum + 2):(2 * natnum + 1),
                    imlayout))
    layout(layout_matrix, 
           widths = c(0.4, rep(1, natnum)),
           heights = c(0.2, 0.1, rep(1, natnum)))
    par(mar = rep(0, 4))
    emptyplot(); text(0, 0, title, cex = 1.5)
    for (i in 1:(2*natnum)) {
        emptyplot()
        if (i < 5) { 
            text(0, 0, nations[i], cex = 1.2)
        } else {
            text(-1, 0, nations[i-4], cex = 1.2, pos = 4)
        }
    }
    par(mar = c(1, 1, 3, 3))
    for (nat in 1:natnum) {
        for (i in 1:natnum) {
            if (nat != i) {
                corrs <- cor(e[,,nat], e[,,i], method = corr_method)
                implot(corrs, p[, c(nat, i)])
            }
        }
    }
}

#' Image plot of p-values
#' 
#' \code{imagePvalues} creates an image plot from a matrix or array of p-values
#' @param pvalues numeric matrix or array of p-values with named dimensions (at
#' least chan and time)
#' @param pcrit numeric vector of significancy limits 
#' (default: 0.001, 0.01, 0.05)
#' @param grid character vector or formula defining the layout of panels
#' @param wrap character vector or formula defining the dimension which 
#' separates panels (only considered if grid is NULL)
#' @param time_label character string; the label of the x (time) axis 
#' (default: "Time (ms)")
#' @param channel_label character string; the label of the y (channel) axis 
#' default: "Channels")
#' @param raster use raster image (TRUE, default) or not (FALSE)
#' @param cluster_order logical; if TRUE, channels are ordered with hierarchical
#' agglomerative clustering (default: FALSE)
#' @export
#' @return A ggplot object
#' @seealso \code{\link{imageValues}} for plotting effects or raw amplitudes
imagePvalues <- function(pvalues, pcrit = c(0.001, 0.01, 0.05),
                         grid = NULL, wrap = NULL,
                         time_label = "Time (ms)", channel_label = "Channels",
                         raster = TRUE, cluster_order = FALSE) {
    rp <- range(pvalues, na.rm = TRUE)
    if (rp[1] < 0 || rp[2] > 1) 
        stop("Data range should be between 0 and 1")
    chans <- dimnames(pvalues)$chan
    if (cluster_order) {
        temp <- -log(pvalues)
        temp <- avgDims(temp, setdiff(names(dimnames(temp)), c("chan", "time")))
        temp <- stats::hclust(dist(aperm(temp, c("chan", "time"))))
        chans <- chans[temp$order]
        pvalues <- subsetArray(pvalues, list(chan = chans))
    }
    timebreaks <- as.numeric( as.character( dimnames(pvalues)$time ) )
    timebreaks <- range( timebreaks %/% 100)
    timebreaks <- seq(timebreaks[1], timebreaks[2]) * 100
    pvalues_df <- transformArray(p ~ ., pvalues)
    pvalues_df$pcrit <- findInterval(pvalues_df$p, pcrit)
    pp <- ggplot(pvalues_df, aes_string(x = "time", y = "chan")) + 
        {if (raster) geom_raster(aes(fill = pcrit)) else 
            geom_tile(aes(fill = pcrit), size = 0)} + 
        scale_fill_gradient(guide = "legend", 
                            high = "white", 
                            low = "#a50f15",
                            name = "p-value",
                            breaks = seq_along(pcrit) - 1L,
                            labels = paste("< ", pcrit, sep = ""),
                            limits = c(0, length(pcrit))) + 
        scale_y_discrete(limits = rev(chans), expand = c(0.05, 0),
                         name = channel_label) + 
        scale_x_continuous(breaks = timebreaks, expand = c(0.05, 0.1),
                           name = time_label)
    if (!is.null(grid)) {
        pp <- pp + facet_grid( as.formula(grid) )
    } else if (!is.null(wrap)) {
        pp <- pp + facet_wrap( as.formula(wrap) )
    } else {
        griddims <- setdiff(names(dimnames(pvalues)), c("chan", "time"))
        if (length(griddims) > 0) {
            if (length(griddims) > 1) {
                grid <- paste(griddims, collapse = "~")  
                pp <- pp + facet_grid( as.formula(grid) )
            } else {
                grid <- paste0("~", griddims)
                pp <- pp + facet_wrap( as.formula(grid) )
            }
        }
    }
    # return
    pp
}


#' Image plot of channel values
#' 
#' \code{imageValues} creates an image plot from a matrix or array of channel x 
#' time values
#' @param dat numeric matrix or array of values with named dimensions (at least
#' chan and time)
#' @param grid character vector or formula defining the layout of panels
#' @param wrap character vector or formula defining the dimension which 
#' separates panels (only considered if grid is NULL)
#' @param bar_title character string; the title of the colour bar
#' @param time_label character string; the label of the x (time) axis 
#' (default: "Time (ms)")
#' @param channel_label character string; the label of the y (channel) axis 
#' default: "Channels")
#' @param raster use raster image (TRUE, default) or not (FALSE)
#' @param cluster_order logical; if TRUE, channels are ordered with hierarchical
#' agglomerative clustering (default: FALSE)
#' @param low,mid,high colour for low/mid/high end of gradient, respectively
#' @param midpoint the midpoint (in data value) of the diverging scale 
#' (default: 0)
#' @param ... other arguments passed to \code{\link{scale_fill_gradient2}}
#' @export
#' @return A ggplot object
#' @seealso \code{\link{imagePvalues}} for plotting p-values, and 
#' \code{\link{plotERParray}} for multiline (butterfly) plots
#' @examples
#' # example dataset
#' data(erps)
#' 
#' # plot grand averages
#' avgs <- avgDims(erps, "id")
#' imageValues(avgs)
#' 
#' # imageValues returns a ggplot object, which can be modified afterwards
#' implot <- imageValues(avgs)
#' 
#' # modify faceting and theme
#' library(ggplot2)
#' implot + facet_grid(pairtype ~ stimclass) + theme_bw()
imageValues <- function(dat, grid = NULL, wrap = NULL, bar_title = "effect",
                        time_label = "Time (ms)", channel_label = "Channels",
                        raster = TRUE, cluster_order = FALSE,
                        low = scales::muted("blue"), mid = "white", 
                        high = scales::muted("red"), midpoint = 0, ...) {
    chans <- dimnames(dat)$chan
    if (cluster_order) {
        temp <- avgDims(dat, setdiff(names(dimnames(dat)), c("chan", "time")))
        temp <- stats::hclust(dist(aperm(temp, c("chan", "time"))))
        chans <- chans[temp$order]
        dat <- subsetArray(dat, list(chan = chans))
    }
    timebreaks <- as.numeric( as.character( dimnames(dat)$time ) )
    timebreaks <- range( timebreaks %/% 100)
    timebreaks <- seq(timebreaks[1], timebreaks[2]) * 100
    dat_df <- transformArray(effect ~ ., dat)
    pp <- ggplot(dat_df, aes_string(x = "time", y = "chan")) + 
        {if (raster) geom_raster(aes_string(fill = "effect")) else 
            geom_tile(aes_string(fill = "effect"), size = 0)} + 
        scale_fill_gradient2(name = bar_title, 
                             low = low, mid = mid, high = high, ...) +               
        scale_y_discrete(limits = rev(chans), expand = c(0.05, 0),
                         name = channel_label) + 
        scale_x_continuous(breaks = timebreaks, expand = c(0.05, 0.1),
                           name = time_label)
    if (!is.null(grid)) {
        pp <- pp + facet_grid( as.formula(grid) )
    } else if (!is.null(wrap)) {
        pp <- pp + facet_wrap( as.formula(wrap) )
    } else {
        griddims <- setdiff(names(dimnames(dat)), c("chan", "time"))
        if (length(griddims) > 0) {
            if (length(griddims) > 1) {
                grid <- paste(griddims, collapse = "~")  
                pp <- pp + facet_grid( as.formula(grid) )
            } else {
                grid <- paste0("~", griddims)
                pp <- pp + facet_wrap( as.formula(grid) )
            }
        }
    }
    # return
    pp
}

#' Make image plots more colorful (only for ggplot objects)
#' 
#' \code{colorize} changes the color scale of \code{\link{imageValues}} plots
#' @param obj a ggplot object created by \code{\link{imageValues}}
#' @param low,mid,high colors to produce a low < mid < high color scale
#' @param ... additional arguments to be passed to \code{link[ggplot2]{scale_fill_gradient2}}
#' @export
#' @return a ggplot object
colorize <- function(obj, low=scales::muted("blue"), mid="white",
                     high=scales::muted("red"), ...) {
    if (!inherits(obj, "ggplot")) {
        stop("The object must be a ggplot")
    }
    suppressMessages(
        obj + scale_fill_gradient2(low = low, mid = mid, high = high, ...))
}


# simple function to plot TFCE effects
tfce.plot <- function(arraydat, breaks = c(0, 0.001, 0.01, 0.05), 
                      colors = rev(RColorBrewer::brewer.pal(length(breaks), "Reds")[-1]), title = "",
                      gridlines_step = 50) {
    extradim <- setdiff(names(dimnames(arraydat)), 
                        c("chan", "time", "modelterm", "nation"))
    if (!is.null(extradim)) arraydat <- avgDims(arraydat, extradim)
    arraydat <- aperm(arraydat, c("chan", "time", "modelterm", "nation"))
    dimnms <- dimnames(arraydat)
    tpoints <- as.numeric(dimnms$time)
    gridlines <- seq(
        min(tpoints) %/% gridlines_step * gridlines_step,
        max(tpoints) %/% gridlines_step * gridlines_step,
        gridlines_step)
    dims <- vapply(dimnms, length, 0L)
    emptyplot <- function() {
        plot(0, 0, xlim = c(-1, 1), 
             type = "n", axes = FALSE, frame.plot = FALSE, xlab = "", ylab = "")
    }
    layoutmat <- rbind(
        c(0, rep(1, dims[3])),
        c(0, (1:dims[3]) + dims[4] + 1),
        cbind(1:dims[4] + 1,
              matrix((1:prod(dims[3:4])) + sum(dims[3:4]) + 1, 
                     dims[4], dims[3], TRUE)))
    layout(layoutmat, 
           widths = c(0.3, rep(1, dims[3])),
           heights = c(0.5, 0.3, rep(1, dims[4])))
    par(mar = c(0, 0, 0, 0))
    # row and column labels
    emptyplot()
    text(0, 0, title, cex = 1.3)
    for (i in 1:dims[4]) {
        emptyplot()
        text(0, 0, dimnms[[4]][i], srt = 90)
    }
    for (i in 1:dims[3]) {
        emptyplot()
        text(0, 0, dimnms[[3]][i])
    }
    par(mar = c(2, 2, 0, 0.5))
    # plot images
    for (i in 1:dims[4]) {
        for (ii in 1:dims[3]) {
            temp <- aperm(arraydat[,,ii,i], c("time", "chan"))
            plot(0, 0,
                 xlim = range(tpoints), ylim = range(1:dims["chan"]),  
                 type = "n", 
                 xlab = "", ylab = "", yaxt = "n")
            abline(v = gridlines, lty = 1, col = "grey95")
            image(tpoints, 1:dims["chan"], temp, 
                  breaks = breaks, col = colors, add = TRUE,
                  xlab = "", ylab = "", yaxt = "n")
            axis(2, at = 1:dims["chan"], labels = FALSE, tick = FALSE)
            text(par("usr")[1] - 21, 1:dims["chan"], cex = 0.6, pos = 2,
                 labels = dimnms$chan, xpd = TRUE)
            abline(v = 0, col = "grey40")
        }
    }
}

###

#' Plot TANOVA results
#' 
#' \code{plotTanova} plots the result of the \code{\link{tanova}} function
#' @param results a list; the return value of the \code{\link{tanova}} function
#' @param grid character vector or formula defining the layout of panels
#' @param wrap character vector or formula defining the dimension which 
#' separates panels (only considered if grid is NULL)
#' @param plot_title character string; the title of the plot
#' @param time_label character string; the label of the x (time) axis 
#' (default: "Time (ms)")
#' @param only_p logical; if TRUE, p-values are plotted instead of combined 
#' (effect + p-value) plots (default: FALSE). See note.
#' @note Use only_p = TRUE if you want to check the uncorrected p-values. If 
#' only_p is set to FALSE (which is the default), \code{plotTanova} highlights
#' the significant phases of the effect curves based on the corrected p-values.
#' @export
#' @return A ggplot object
plotTanova <- function(results, grid = NULL, wrap = NULL, 
                       plot_title = "", time_label = "Time (ms)", 
                       only_p = FALSE) {
    reshapefn <- function(slot, headername) {
        x <- results[[slot]]
        x <- array2df(x, response_name = headername, 
                      dim_types = list(time = "numeric"))
        return(x)
    }
    #
    pcrit <- eval(as.list(results$call)$pcrit)
    if (is.null(pcrit)) pcrit <- formals(tanova)$pcrit
    pcrit <- union(sort(pcrit), 1)
    #
    dat <- transformArray(effect ~ ., results$effect)
    dat$pvalue <- -log(as.vector(results$perm_pvalues))
    dat$pcrit <- factor(results$perm_pvalues_consec, 
                        levels = as.character(pcrit), 
                        labels = c(as.character(pcrit[-length(pcrit)]), 
                                   "n.s."))
    final_pcrit <- levels(droplevels(dat)$pcrit)
    colour_pcrit <- 
        if (length(final_pcrit) > 2) {
            c(rev(brewer.pal(length(final_pcrit), "Reds")[-1]), "grey60")
        } else {
            c(brewer.pal(3, "Reds")[3], "grey60")
        }
    legendtitle <- "p-value"
    #
    if (only_p) {
        qp <- ggplot(dat[order(dat$time),], 
                     aes_string(x = "time", y = "pvalue", col = "pcrit", 
                                group = NA)) + 
            geom_hline(yintercept = -log(pcrit), lty = 3) + 
            ylab("-log(P-value)")
    } else {
        qp <- ggplot(dat[order(dat$time),], 
                     aes_string(x = "time", y = "effect", col = "pcrit", 
                                group = NA)) + 
            ylab("Effect")
    }
    if (!is.null(grid)) {
        qp <- qp + facet_grid( as.formula(grid) )
    } else if (!is.null(wrap)) {
        qp <- qp + facet_wrap( as.formula(wrap) )
    } else {
        griddims <- setdiff(names(dimnames(results$effect)), "time")
        if (length(griddims) > 0) {
            if (length(griddims) > 1) {
                grid <- paste(griddims, collapse = "~")  
                qp <- qp + facet_grid( as.formula(grid) )
            } else {
                grid <- paste0("~", griddims)
                qp <- qp + facet_wrap( as.formula(grid) )
            }
        }
    }
    qp <- qp + geom_path() + 
        xlab(time_label) + ggtitle(plot_title) + 
        scale_colour_manual(name = legendtitle,
                            values = colour_pcrit) + 
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(color = "grey95"),
              panel.grid.minor = element_line(color = "grey97"),
              panel.border = element_rect(color = "grey70", fill = NA))
    #
    qp
}


# plot neurodys tanova results
fastplot_tanova <- function(results, plot_title = "", pcrit = 0.05, 
                            only_p = FALSE) {
    reshapefn <- function(slot, headername) {
        x <- lapply(results, "[[", slot)
        #x <- rearrangeList(x, "nation")
        if (length(x) > 1) {
            temp <- matrix(unlist(strsplit(names(x), "-")), 
                           nrow = length(x), byrow = TRUE)
            names(x) <- temp[, 2]
            x <- rearrangeList(x, temp[1, 1])
        } else {
            x <- x[[1]]
        }
        x <- as.data.frame.table(x, responseName = headername)
        x$time <- as.numeric(as.character(x$time))
        return(x)
    }
    dimn <- names(results[[1]])
    dat <- reshapefn(dimn[grep("effect_", dimn)], "effect_size")
    dat$pvalue <- -log(reshapefn("perm_pvalues", "pvalue")$pvalue)
    dat$pvalue_consec <- reshapefn("perm_pvalues_consec", "pvalue")$pvalue
    dat$pcrit <- factor(dat$pvalue_consec<pcrit)
    legendtitle <- paste("pvalue <", substitute(pcrit), sep = " ")
    #
    if (only_p) {
        qp <- ggplot(dat[order(dat$time),], 
                     aes(x = time, y = pvalue, col = pcrit, group = NA)) + 
            geom_hline(yintercept = -log(pcrit), lty = 3) + 
            ylab("-log(P-value)")
    } else {
        qp <- ggplot(dat[order(dat$time),], 
                     aes(x = time, y = effect_size, col = pcrit, group = NA)) + 
            ylab("Effect size")
    }
    qp <- qp + geom_line() + facet_grid(nation~modelterm) +
        ggtitle(plot_title) + 
        scale_colour_manual(name = legendtitle,
                            values = c("FALSE"="grey70","TRUE"="red"))
    print(qp)
}


# plot peak results
fastplot_peak <- function(results, plot_title = "", 
                          pcrit = 0.05, only_p = FALSE) {
    reshapefn <- function(slot, headername, attr_slot = NULL) {
        if (is.null(attr_slot)) {
            x <- lapply(results, "[[", slot)
        } else {
            x <- lapply(results, 
                        function(x) attr(x[[slot]], attr_slot))
        }
        x <- rearrangeList(x, "nation")
        x <- as.data.frame.table(x, responseName = headername)
        x$peak <- factor(x$peak, levels = rev(levels(x$peak)))
        return(x)
    }
    dat_ampl <- reshapefn("F_obs", "pvalue", attr_slot = "pvalues")
    dat_ampl <- dat_ampl[with(dat_ampl, order(nation, peak, modelterm)),
                         c("nation", "peak", "modelterm", "pvalue")]
    dat_ampl$measure <- "amplitude"
    dat_lat <- reshapefn("lat_pvalues", "pvalue")
    dat_lat <- dat_lat[with(dat_lat, order(nation, peak, term)),
                       c("nation", "peak", "term", "pvalue")]
    colnames(dat_lat)[3] <- "modelterm"
    dat_lat$measure <- "latency"
    dat <- rbind(dat_ampl, dat_lat)
    dat$pcrit <- factor(dat$pvalue < pcrit)
    legendtitle <- paste("pvalue <", substitute(pcrit), sep = " ")
    qp <- ggplot(dat, aes(x = nation, y = peak, fill = pcrit)) + 
        geom_tile(col = "white") + facet_grid(modelterm ~ measure) + 
        ggtitle(plot_title) + 
        scale_fill_manual(name = legendtitle,
                          values = c("FALSE"="grey60","TRUE"="red3"))
    print(qp)
}

#
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow = 2, byrow = TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
# Author: Winston Chang
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
    require("grid")
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots == 1) {
        print(plots[[1]])   
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))        
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}
