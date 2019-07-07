#
# <<< plotting functions >>> -----------
#

#' Helper function to create customized ggplot2 theme
#' 
#' \code{theme_ndys} modifies the base ggplot2 theme for figures in Leppanen et 
#' al.
#' @param base_fontsize size of texts
#' @param legend_direction the direction of the various levels of a legend key
#' (default: "vertical")
#' @param panel_aspect_ratio the ratio of the y and x axes on the plot
#' @export
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(x = wt, y = mpg, colour = factor(cyl))) + 
#'     geom_point() + facet_grid(~ am) + 
#'     theme_ndys(12, "h", 1.5)
#'     
theme_ndys <- function(base_fontsize = 8, 
                       legend_direction = c("vertical", "horizontal"), 
                       panel_aspect_ratio = NULL) {
  # check argument
  assertNumber(base_fontsize, .var.name = "base_fontsize")
  legend_direction <- match.arg(legend_direction)
  if (!is.null(panel_aspect_ratio)) {
    assertNumber(panel_aspect_ratio, .var.name = "panel_aspect_ratio")
  }
  #
  theme_bw(base_size = base_fontsize) + 
    theme(axis.title.x = element_text(vjust = -0.3, size = rel(1)), 
          axis.title.y = element_text(vjust = -1, size = rel(1)), 
          strip.text = element_text(size = rel(1),
                                    lineheight = 0.9),
          axis.text.x = element_text(size = rel(0.9)),
          axis.text.y = element_text(size = rel(0.8)),
          plot.title = element_text(size = rel(1)),
          legend.title = element_text(size = rel(0.85)),
          legend.text = element_text(size = rel(0.85)),
          legend.key = element_blank(),
          legend.key.height = unit(11, "pt"),
          legend.key.width = unit(10, "pt"),
          legend.position = "top", 
          legend.box = "horizontal", 
          legend.direction = legend_direction,
          panel.border = element_rect(colour = "grey60", size = 0.3,
                                      linetype = 1),
          strip.background = element_rect(colour = NA, fill = "grey80"),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          aspect.ratio = panel_aspect_ratio)
}

#' Plot ERP curves
#' 
#' \code{plotERParray} is a generalization of \code{\link{matplot}} onto array 
#' inputs.
#' @param dat a numeric array with named dimensions
#' @param xdim character string; the name of the dimension of dat which defines 
#' the x axis (default: "time")
#' @param sepdim character string; the name of the dimension of dat which 
#' separates the lines (default: "chan")
#' @param gfp_plot logical value whether GFP curves should be also plotted
#' (only if 'sepdim' is "chan")
#' @param gfp_col the colour of the GFP curves (if plotted)
#' @param gfp_lwd the thickness of the GFP curves (if plotted)
#' @param minus_up logical value; if set to TRUE, minus values are plotted 
#' upwards. If FALSE, minus values are plotted downwards. By default, it 
#' is set to TRUE or FALSE depending on the 'gfp_plot' argument. 
#' @param title character string; the title of the plot
#' @param subtitle_col colour(s) of the subtitles (recycled if necessary)
#' @param xlab character string; the title of the x axis (defaults to the 
#' 'xdim' argument)
#' @param ylab character string; the title of the y axis (default: "amplitude")
#' @param grid_dim integer vector giving the number of rows and columns
#' in which the plots are arranged; if set to NULL (default), the arrangement 
#' of the plot is set up automatically
#' @param ... additional parameters to be passed to \code{\link{matlines}}
#' @seealso \code{\link{matlines}}
#' @examples
#' # example data
#' data(erps)
#' 
#' # compute grand averages
#' plotdat <- avgDims(erps, "id")
#' 
#' # plot the channels in light grey, draw only solid lines
#' plotERParray(plotdat, col = "grey70", lty = 1)
#' 
#' @export
plotERParray <- function(dat, xdim = "time", sepdim = "chan",  
                         gfp_plot = if (sepdim == "chan") TRUE else FALSE, 
                         gfp_col = "black", gfp_lwd = 1.3, 
                         minus_up = if (gfp_plot) FALSE else TRUE,
                         title = "", subtitle_col = "black",
                         xlab = xdim, ylab = "amplitude", 
                         grid_dim = NULL, ...) {
  # helper function
  emptyplot <- function() {
    plot(0, 0, xlim = c(-1,1), type = "n", axes = FALSE, 
         frame.plot = FALSE, xlab = "", ylab = "")
  }
  #
  assertArray(dat, mode = "numeric", min.d = 1L, .var.name = "dat")
  plot_args <- list(...)
  dat <- decorateDims(dat)
  if (!is.character(xdim)) xdim <- names(dimnames(dat))[xdim]
  if (!is.character(sepdim)) sepdim <- names(dimnames(dat))[sepdim]
  if (gfp_plot && !identical(sepdim, "chan")) {
    stop("If 'gfp_plot' is TRUE, 'sepdim' must be 'chan'")
  }
  if (length(dim(dat)) > 3L) {
    wrapdim <- setdiff(names(dimnames(dat)), c(xdim, sepdim))
    dat <- mergeDims(dat, list(xdim, sepdim, wrapdim))
  } else if (length(dim(dat)) < 3L) {
    dimn <- dimnames(dat)
    dim(dat) <- c(dim(dat), 1L)
    dimnames(dat) <- c(dimn, list(""))
    dat <- apermArray(dat, first = c(xdim, sepdim))
  } else {
    dat <- apermArray(dat, first = c(xdim, sepdim))
  }
  subtitle_col <- rep(subtitle_col, length_out = dim(dat)[3])
  if (gfp_plot) gfpdat <- compGfp(dat)
  suppressWarnings(x <- as.numeric(dimnames(dat)[[1]]))
  if (anyNA(x)) {
    x <- seq_along(x)
    setattr(x, "type", "char")
  } else {
    setattr(x, "type", "num")
  }
  xrange <- if (is.null(plot_args$xlim)) range(x) else plot_args$xlim 
  if (is.null(plot_args$ylim)) {
    yr <- range(dat)
    yrange <- mean(yr) + c(-1, 1)*(yr[2]-mean(yr))*1.1
    if (minus_up) {
      if (yrange[1] > -.Machine$double.eps ^ 0.5) {
        yrange[1] <- -.Machine$double.eps ^ 0.5
      }
      yrange <- rev(yrange)
    }
  } else {
    yrange <- plot_args$ylim
  } 
  if (is.null(grid_dim)) grid_dim = rep(ceiling(sqrt(dim(dat)[3L])), 2L)
  layoutmat <- matrix(0, grid_dim[1] + 2, grid_dim[2] + 2)
  layoutmat[-nrow(layoutmat), -c(1, ncol(layoutmat))] <- 
    rbind(rep(1, grid_dim[2]),
          matrix(1:prod(grid_dim), 
                 grid_dim[1], grid_dim[2], TRUE) + 1L)
  layout(layoutmat, 
         widths = c(0.1, rep(1, grid_dim[2]), 0.2),
         heights = c(0.2, rep(1, grid_dim[1]), 0.4))
  par(mar = rep(0, 4))
  emptyplot(); text(0, 0, title, cex = 1.5)
  par(mar = c(0, 0, 2, 0))
  for (i in 1:dim(dat)[3L]) {
    matplot(x, dat[,,i], xlim = xrange, ylim = yrange, yaxs = "i",  
            axes = FALSE, frame.plot = TRUE, type = "n")
    grid()
    mtext(dimnames(dat)[[3L]][i], 3, cex = 0.7, col = subtitle_col[i])
    matlines(x, dat[,,i], axes = FALSE, ...)
    if (i == max(dim(dat)[3L])) {
      mtext(xlab, 1, 3)
      mtext(ylab, 4, 3)
      if (attr(x, "type") == "char") {
        axis(1, at = x, labels = dimnames(dat)[[1]])
      } else {
        axis(1)
      }
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
#' @param ch_pos a data.frame of electrode positions, 
#' see \code{\link{coordinates}}
#' @param r a numeric value; the radius of the head (not used if ch_pos contains
#' a column \code{r})
#' @param timepoint an integer value representing the time point at which the 
#' amplitude values were recorded. This is used as the title of the topoplot if 
#' the 'title' argument is not provided, and for highlighting the time point on
#' the GFP curve if 'gfp' data are provided (with proper names or dimnames 
#' attribute).
#' @param ampl_range numeric vector, the minimum and maximum value of the 
#' amplitudes. If NULL (default), it is computed from the data. Data values
#' outside the \code{ampl_range} interval are winsorized (set to equal to the
#' minimum or maximum limit).
#' @param type the type of the topoplot, either "contour" (the default) or 
#' "raster". The latter is faster for very high resolution topoplots.
#' @param resol integer value, the resolution of the projection locations for
#' each axis (default: 30L) 
#' @param resolcol integer value, the number of levels for colouring (default: 
#' 101L)
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
#' @param colour_scale function which creates a colour scale with \code{n} values
#' (\code{n} is defined by \code{resolcol}). The default is to use the 
#' \code{\link[gplots]{colorpanel}} function, with blue > white > red 
#' colour scale.
#' @param colour_midpoint a number to which the middle colour of the colour scale
#' corresponds
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
                       ampl_range = NULL, type = c("contour", "raster"),
                       resol = 30L, resolcol = 101L, 
                       projection = "laea", projref = c("pole", "equator"),
                       origo = c(lat = ifelse(projref == "pole", 90, 0),
                                 long = ifelse(projref == "pole", 270, 0)),
                       plot_centroid = TRUE, centroid_circle = 0.5, 
                       centroid_unitsphere = FALSE,
                       plot_bar = TRUE, plot_ch = TRUE, plot_chnames = TRUE, 
                       title = NULL,
                       gfp = NULL, gfp_max = NULL, 
                       colour_scale = function(n) gplots::colorpanel(n, "blue", "white", "red"), 
                       colour_midpoint = 0, ...) {
  #
  type <- match.arg(type)
  projref <- match.arg(projref)
  assertFunction(colour_scale, .var.name = "color_scale")
  assertNumber(colour_midpoint, .var.name = "color_midpoint")
  #
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
  assertNumeric(dat, .var.name = "dat")
  if (is.matrix(dat)) {
    if (all(c("chan", "id") %in% names(dimnames(dat))))
      dat <- aperm(dat, c("id", "chan"))
    if (ncol(dat) != nrow(ch_pos)) 
      stop("Wrong data format! Should be: Rows = participants, Columns = channels")
    subjdat <- dat
    dat <- colMeans(dat, na.rm = TRUE)
  } else if (length(dim(dat)) > 2L) {
    stop("The input data may not have more than 2 dimensions")
  } else if (length(dim(dat)) == 1L) {
    dat <- as.vector(dat)
    centroid_circle <- NA
  } else {
    centroid_circle <- NA
  }
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
                    boundarypos$x*1.1, boundarypos$y*1.1)
  gridgeo <- project3dMap(gridpos[ind,], projection = projection, 
                          projref = projref, origo = origo, inverse = TRUE)
  gridcart <- geo2cart(gridgeo)
  z <- matrix_(NA_real_, resol, resol)
  z_ind <- suppressMessages(chanInterp(dat, ch_pos, gridcart, ...))
  if (is.null(ampl_range)) {
    ampl_range <- range(z_ind, na.rm = TRUE)
  } else {
    assert0(ampl_range, checkNumeric, len = 2L)
    sort(ampl_range)
  }
  z_ind[z_ind > ampl_range[2]] <- ampl_range[2]
  z_ind[z_ind < ampl_range[1]] <- ampl_range[1]
  z[ind] <- z_ind
  rm(z_ind)
  #
  colstep <- diff(ampl_range)/resolcol
  colrange <- sort(colour_midpoint + 
                     c(-1, 1) * max(abs(ampl_range - colour_midpoint)))
  coln <- ceiling(diff(colrange)/colstep/2)*2 + 1L
  colors <- match.fun(colour_scale)(coln)
  colindex <- floor((ampl_range[1] - colrange[1])/colstep) + 1L
  colors <- colors[colindex:(colindex + resolcol)]
  #
  par(mar = rep(0, 4))
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
  if (type == "contour") {
    .filled.contour(
      gridx, gridy, z, 
      levels = seq(ampl_range[1], ampl_range[2], length.out = resolcol + 1L), 
      col = colors)
  } else {
    image(gridx, gridy, z, useRaster = TRUE, col = colors,
          xlim = xlim, ylim = ylim,
          zlim = ampl_range, xlab = "", ylab = "", axes = FALSE, add = TRUE)
  }
  polypath(c(xlim, rev(xlim), xlim[1], NA_real_, boundarypos$x),
           c(rep(ylim, each = 2L), ylim[1], NA_real_, boundarypos$y),
           border = NA, col = "white", rule = "evenodd")
  #
  if (plot_ch) {
    ch_pos_xy <- project3dMap(ch_pos)
    points(ch_pos_xy, NULL, pch = 20)
    if (plot_chnames) {
      indleft <- ch_pos_xy$x <= 0
      indright <- ch_pos_xy$x > 0
      text(ch_pos_xy[indleft,], NULL, rownames(ch_pos[indleft, ]), 
           pos = 4)
      text(ch_pos_xy[indright,], NULL, rownames(ch_pos[indright, ]), 
           pos = 2)
    }
  }
  #
  if (plot_bar) {
    xpos.bar <- seq(min(boundarypos$x / 3), max(boundarypos$x / 3), 
                    length.out = resolcol)
    image(xpos.bar, ypos.bar,
          matrix_(seq(ampl_range[1], ampl_range[2], length.out = resolcol), 
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
        which(as.integer(dimnames(gfp)[[1]]) == timepoint)
      } else {
        which(as.integer(names(gfp)) == timepoint)
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
        temp <- matrix_(unlist(lapply(cc, function(x) abs(x-c1))), 
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
#' @param ch_pos a data.frame of electrode positions, 
#' see \code{\link{coordinates}}
#' @param timepoint integer; the timepoint at which the scalp distributions are
#' compared. Can be a vector referring to the timepoints of each topomap.  
#' @param datgrid arrangement of the grid (number of rows X columns)
#' @param layout_matrix,heights,widths arguments to be passed to 
#' \code{\link{layout}}. If given, it overrides \code{datgrid}. Note that the
#' layout matrix refers to the whole plot, including the main title, the 
#' column titles, the row titles, the scalp plots, and the colourbar plot 
#' (in this order). 
#' @param title_row,title_col character vectors of the titles of rows and 
#' columns, respectively. If NULL and dat has dimnames attribute, 
#' \code{rownames(dat)} and \code{colnames(dat)} are used.
#' @param title_map character vector; the titles of the individual scalp maps. 
#' if NULL and \code{dat} has \code{names} attribute, \code{names(dat)} is 
#' used.
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
#'                   gfp = gfpdat, plot_ch = FALSE)
#
complexplot2dview <- function(dat, ch_pos, timepoint, 
                              datgrid = NULL, layout_matrix = NULL, 
                              heights = NULL, widths = NULL, 
                              title_row = NULL, title_col = NULL, 
                              title_map = NULL,
                              title_main = paste("Time:", 
                                                 paste(timepoint, 
                                                       collapse = ", "), 
                                                 "ms"), 
                              plot_bar = TRUE, ampl_range = NULL, 
                              gfp = NULL, ...) {
  # helper function
  emptyplot <- function() plot(0, 0, xlim = c(-1, 1), type = "n", 
                               axes = FALSE, frame.plot = FALSE, 
                               xlab = "", ylab = "")
  # argument checks
  assertList(dat, .var.name = "dat")
  assertDataFrame(ch_pos, .var.name = "ch_pos")
  assertAtomicVector(timepoint, any.missing = FALSE, .var.name = "timepoint")
  assert0(title_main, checkString)
  assertFlag(plot_bar, .var.name = "plot_bar")
  assert0(ampl_range, checkNumeric, 
          finite = TRUE, any.missing = FALSE, len = 2L)
  assert0(gfp, checkList, types = "numeric", len = length(dat))
  # layout
  datnum <- length(dat)
  if (!is.null(layout_matrix)) {
    assert0(layout_matrix, checkMatrix, mode = "numeric")
    rownum <- abs(diff(range(layout_matrix[-(1:2), 1L]))) + 1L
    colnum <- abs(diff(range(layout_matrix[2L, -1L]))) + 1L
  } else if (!is.null(datgrid)) {
    assert0(datgrid, checkIntegerish, any.missing = FALSE, len = 2L)
    if (prod(datgrid) < datnum) {
      stop(paste0("in complexplot2dview: ",
                  "the product of 'datgrid' is less than the ",
                  "length of 'dat'"), call. = FALSE)
    }
    rownum <- datgrid[1]
    colnum <- datgrid[2]
  } else if (!is.null(dim(dat))) {
    dims <- dim(dat)
    rownum <- dims[1L]
    colnum <- datnum/rownum
  } else if (!is.null(title_row)) {
    assert0(title_row, checkCharacter, max.len = datnum)
    rownum <- length(title_row)
    colnum <- ceiling(datnum/rownum)
  } else if (!is.null(title_col)) {
    assert0(title_col, checkCharacter, max.len = datnum)
    colnum <- length(title_col)
    rownum <- ceiling(datnum/colnum)
  } else {
    rownum <- floor(sqrt(datnum))
    colnum <- ceiling(datnum/rownum)
  }
  if (is.null(layout_matrix)) {
    layout_matrix <- rbind(
      matrix(
        c(
          # main title
          rep_len(1L, colnum), 0L,
          # column titles
          2:(colnum + 1), 0L
        ), 
        2L, colnum + 1L, TRUE), 
      matrix(
        c(
          # scalp maps
          c(colnum + rownum + 1L + 1:datnum, 
            rep_len(0L, rownum*colnum - datnum)),
          # row titles
          (colnum + 2L):(colnum + rownum + 1L)
        ), 
        rownum, colnum + 1L))
    layout_matrix <- rbind(layout_matrix, 
                           c(rep(max(layout_matrix) + 1L, colnum), 0L))
  }
  # height and widths
  heights <- 
    if (is.null(heights)) {
      c(0.3, 0.2, rep(1, rownum), 0.2)
    } else {
      assert0(heights, checkNumeric)
      rep_len(heights, rownum + 3L)
    }
  widths <- 
    if (is.null(widths)) {
      widths <- c(rep(1, colnum), 0.2)
    } else {
      assert0(widths, checkNumeric)
      rep_len(widths, colnum + 1L) 
    }
  # row titles 
  if (is.null(title_row)) {
    rown <- rownames(dat)
    title_row <- if (is.null(rown)) rep_len("", rownum) else rown
  } else {
    assert0(title_row, checkCharacter, len = rownum)
  }
  # column titles
  if (is.null(title_col)) {
    dimn <- dimnames(dat)[-1L]
    title_col <- 
      if (length(dimn) > 0L) {
        dimn <- expand.grid(dimn, 
                            KEEP.OUT.ATTRS = FALSE, 
                            stringsAsFactors = FALSE)
        do.call(paste, c(dimn, list(sep = ".")))
      } else {
        rep_len("", colnum)
      }
  } else {
    assert0(title_col, checkCharacter, len = colnum)
  }
  # map titles
  if (is.null(title_map)) {
    n <- names(dat)
    title_map <- 
      if (!is.null(n) && is.null(title_row) && is.null(title_col)) {
        n
      } else {
        rep_len("", datnum)
      }
  } else {
    assert0(title_map, checkCharacter, len = length(dat))
  }
  #
  # start plotting
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
  for (i in 1:datnum) {
    plot2dview(dat[[i]], ch_pos = ch_pos, timepoint = timepoint[i], 
               gfp = gfp[[i]], gfp_max = gfp_max,
               plot_bar = FALSE, ampl_range = ampl_range, 
               title = title_map[i], ...)
  }
  if (plot_bar) {
    resolcol <- list(...)$resolcol
    if (is.null(resolcol)) resolcol <- formals(plot2dview)$resolcol
    xpos.bar <- seq(0.4, 0.6, length.out = resolcol)
    ypos.bar <- c(0.45, 0.55)
    image(xpos.bar, ypos.bar,
          matrix_(seq(ampl_range[1], ampl_range[2], length.out = resolcol), 
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
                         low = low, mid = mid, high = high, 
                         midpoint = midpoint, space = "Lab", 
                         ...) +               
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
colorize <- function(obj, low = scales::muted("blue"), mid = "white",
                     high = scales::muted("red"), ...) {
  if (!inherits(obj, "ggplot")) {
    stop("The object must be a ggplot")
  }
  suppressMessages(
    obj + 
      scale_fill_gradient2(low = low, mid = mid, high = high, 
                           space = "Lab", ...))
}

###

#' Plot results of t-test, ANOVA or TANOVA modeling functions
#' 
#' \code{modelplot} is a generic function for plotting
#' \code{\link{arrayTtest}}, \code{\link{arrayAnova}}, or \code{\link{tanova}}
#' objects, each produced by the corresponding function.
#' @inheritParams imageValues
#' @param results a list; the return value of the corresponding modeling 
#' functions
#' @param what a character string or vector indicating what should be plotted:
#' only the test statistics ("statistic"), only the p-values ("p-values"), or 
#' both (default). The names can be abbreviated.
#' @param type display "corrected" (default) or "uncorrected" statistics or
#' p-values. The names can be abbreviated.
#' @param pcrit numeric vector of significancy limits 
#' (default: 0.001, 0.01, 0.05) for highlighting significant areas
#' @param ... arguments passed to \code{\link{extract}}; e.g. use 
#' arguments 'time_window' and 'term' to display only a subset of the results
#' @export
#' @return A ggplot object or a list of such objects
#' @seealso See the examples for \code{\link{arrayAnova}} and 
#' \code{\link{tanova}}
modelplot <- function(...) UseMethod("modelplot")

#' Plot arrayTtest or arrayAnova results
#' 
#' \code{modelplot.default} plots the result of the \code{\link{arrayTtest}} or
#' \code{\link{arrayAnova}} function.
#' @export
#' @describeIn modelplot Default method
modelplot.default <- function(results, 
                              what = c("statistic", "p-value"), 
                              type = c("corrected", "uncorrected"),
                              pcrit = c(0.001, 0.01, 0.05),
                              grid = NULL, wrap = NULL, 
                              time_label = "Time (ms)", 
                              channel_label = "Channels", 
                              raster = TRUE, 
                              cluster_order = FALSE, 
                              low = scales::muted("blue"), mid = "white", 
                              high = scales::muted("red"), midpoint = 0, 
                              ...) {
  #
  if (!inherits(results, c("arrayTtest", "arrayAnova")))
    stop(paste0("'results' does not have class 'arrayTest' ",
                "or 'arrayAnova' or 'tanova'"))
  what <- match.arg(what, several.ok = TRUE)
  type <- match.arg(type)
  which_stat <- if (type == "corrected") "stat_corr" else "stat"
  which_p <- if (type == "corrected") "p_corr" else "p"
  dat <- out <- vector("list", 2)
  if ("statistic" %in% what) {
    dat <- extract(results, which_stat, ...)
    subplot_title <- attr(dat, "label")
    bar_title <- 
      if (type == "corrected") "TFCE statistic" else "Test statistic"
    out[[1L]] <- imageValues(
      dat, 
      grid = grid, wrap = wrap, bar_title = bar_title,
      time_label = time_label, channel_label = channel_label,
      raster = raster, cluster_order = cluster_order,
      low = low, mid = mid, high = high, midpoint = midpoint) +
      ggtitle(subplot_title)
    i <- 1L
  }
  if ("p-value" %in% what) {
    dat <- extract(results, which_p, ...)
    subplot_title <- attr(dat, "label")
    out[[2L]] <- imagePvalues(
      dat, 
      pcrit = pcrit, grid = grid, wrap = wrap,
      time_label = time_label, channel_label = channel_label,
      raster = raster, cluster_order = cluster_order) + 
      ggtitle(subplot_title)
    i <- 2L
  } 
  # return
  if (length(what) == 2L) {
    grid::grid.newpage()
    grid::grid.draw(
      gtable::gtable_row(
        "result",
        list(ggplotGrob(out[[1L]]), ggplotGrob(out[[2L]])),
        height = unit(1, "npc")))
    invisible(out)
  } else {
    out[[i]]
  }
}


#' Plot TANOVA results
#' 
#' \code{modelplot.tanova} plots the result of the \code{\link{tanova}} function.
#' @export
#' @keywords internal
#' @describeIn modelplot Method for \code{tanova} objects
modelplot.tanova <- function(results, 
                             what = c("statistic", "p-value"),
                             type = c("corrected", "uncorrected"), 
                             grid = NULL, wrap = NULL, 
                             time_label = "Time (ms)", 
                             ...) {
  #
  if (!inherits(results, "tanova"))
    stop("'results' does not have class 'tanova'")
  what <- match.arg(what, several.ok = TRUE)
  type <- match.arg(type)
  # if statistic is in 'what', the (corrected) test statistic will be 
  # displayed, potentially highlighted (if p-value is also in "what");
  # otherwise, the raw p-values will be displayed
  if ("statistic" %in% what) {
    extract_stat <- if (type == "corrected") "stat_corr" else "stat"
    dat <- extract(results, extract_stat, ...)
    ylabel <- attr(dat, "label")
    dat <- transformArray(y ~ ., dat)
    if ("p-value" %in% what) {
      highlight_p <- if (type == "corrected") "p_corr" else "p"
      hpvals <- extract(results, highlight_p, ...)
    }
  } else {
    hpvals <- extract(results, "p", ...)
    ylabel <- sprintf("-log(%s)", attr(hpvals, "label"))
    dat <- transformArray(y ~ ., -log(hpvals))
    # hpvals will be used for highlighting; overwrite with p_corr if needed
    if (type == "corrected") hpvals <- extract(results, "p_corr", ...)
  }
  # highlight with corrected or uncorrected p-values
  if ("p-value" %in% what) {
    # alpha-level
    pcrit <- try(eval(as.list(results$call)$pcrit), silent = TRUE)
    if (is.null(pcrit)) {
      pcrit <- formals(tanova)$pcrit
    }
    if (!is.numeric(pcrit) || is.na(pcrit)) {
      pcrit <- 
        if (type == "corrected") {
          unique(as.vector(hpvals))
        } else {
          0.05
        }
    }
    pcrit <- union(sort(pcrit), 1)
    #
    if (type == "uncorrected") {
      ind <- hpvals <= pcrit[1L]
      hpvals[ind] <- pcrit[1L] 
      for (i in 2:length(pcrit)) {
        ind <- !ind & hpvals <= pcrit[i]
        hpvals[ind] <- pcrit[i] 
      }
    }
    # shift edges
    hpvals <- fnDims(hpvals, "time", 
                     function(x) {
                       dx <- diff(x)
                       dx[dx < 0] <- 0
                       x[-nrow(x), ] <- x[-nrow(x), ] + dx
                       x
                     },
                     vectorized = TRUE,
                     keep_dimorder = TRUE)
    # transform to factor
    dat$pcrit <- factor(hpvals, 
                        levels = as.character(pcrit), 
                        labels = c(as.character(pcrit[-length(pcrit)]), 
                                   "n.s."))
    final_pcrit <- levels(droplevels(dat)$pcrit)
    colour_pcrit <- rev(brewer.pal(max(3L, length(final_pcrit)), 
                                   "Reds"))[seq_along(final_pcrit)]
    colour_pcrit[final_pcrit == "n.s."] <- "grey60"
    legendtitle <- "p-value"
    if (type == "corrected") {
      legendtitle <- paste(legendtitle, "(corrected)")
    }
  }
  #
  # basic plot
  scales <- "fixed"
  if (identical(what, "statistic")) {
    if (type == "uncorrected") scales <- "free_y"
    qp <- ggplot(dat[order(dat$time),], 
                 aes_string(x = "time", y = "y"))
  } else {
    qp <- ggplot(dat[order(dat$time),], 
                 aes_string(x = "time", y = "y", col = "pcrit", 
                            group = NA)) + 
      scale_colour_manual(name = legendtitle,
                          values = colour_pcrit)
    if (length(what) == 1L) {
      pcrit <- setdiff(pcrit, 1)
      qp <- qp + 
        geom_hline(yintercept = -log(pcrit), linetype = "dotted") + 
        annotate("text", x = max(dat$time), y = -log(pcrit),
                 size = 3, label = paste0("p < ", pcrit),
                 hjust = "right", vjust = -0.3,
                 colour = "grey50")
    }
  }
  #
  # facetting
  if (!is.null(grid)) {
    qp <- qp + facet_grid( as.formula(grid) )
  } else if (!is.null(wrap)) {
    qp <- qp + facet_wrap( as.formula(wrap) )
  } else {
    griddims <- names(fillMissingDimnames(results$dimnames, results$dim))
    griddims <- setdiff(griddims, "time")
    if ("modelterm" %in% griddims) {
      griddims <- c("modelterm", setdiff(griddims, "modelterm"))
    }
    if (length(griddims) > 0) {
      if (length(griddims) > 1) {
        grid <- paste(griddims, collapse = "~")  
        qp <- qp + facet_grid( as.formula(grid), scales = scales )
      } else {
        grid <- paste0("~", griddims)
        qp <- qp + facet_wrap( as.formula(grid), scales = scales )
      }
    }
  }
  #
  # finalize plot
  qp + 
    geom_path() + 
    ylab(ylabel) + 
    xlab(time_label) +  
    theme_bw()
}




#' Add caption to plot, and print it to the console
#' 
#' \code{addCaption} appends a string to a ggplot object as attribute 'caption',
#' and prints it to the console. No copy is made.
#' @param x a ggplot object
#' @param caption a character string
#' @keywords internal
addCaption <- function(x, caption) {
  assertClass(x, "ggplot", .var.name = "x")
  assertString(caption, .var.name = "caption")
  setattr(x, "caption", caption)
  cat(paste0("\nCaption: ", caption, "\n"), 
      file = stdout())
  invisible(x)
}

#' Create butterfly plot with extras: highlighting and peak topographies
#' 
#' \code{plotButterfly} creates a multi-channel (a.k.a butterfly) plot. 
#' Additionally, significant time windows can be highlighted, and peak 
#' topographies can be displayed above the ERP curves.
#' @inheritParams plotComplex
#' @param caption logical flag indicating if caption should be also returned
#' (default: TRUE)
#' @return \code{plotButterfly} returns a ggplot object.
#' @export
#' @examples
#' # load example data
#' data(erps)
#' 
#' # extract channel positions
#' chan_pos <- attr(erps, "chan")
#' 
#' # collapse pairtypes and participants
#' tempdat <- avgDims(erps, c("pairtype", "id"))
#' 
#' # plot butterfly with topo-maps at specified time points
#' plotButterfly(tempdat, topo_time = seq(24, 476, by = 50),
#'               chan_pos = chan_pos)
#' 
#' # plot butterfly with topo-maps at peaks which are selected 
#' # automatically; let's look for local maxima on the GFP curves between
#' # 0 and 480 ms
#' # 1) add GFP to the dataset
#' tempdat2 <- compGfp(tempdat, keep_channels = TRUE)
#' 
#' # 2) provide the peak definition
#' peak_def <- isLocalMaximum(
#'     subset. = list(time = isBetween(0, 480), chan = "GFP"),
#'     options. = list(along_dim = "time", n = 15))
#' 
#' # 3) find the peaks
#' peak_data <- selectValues(tempdat2, peak_def)
#' 
#' # 4) create plot
#' plotButterfly(tempdat2, topo_time = peak_data, chan_pos = chan_pos)
#' 
#' # highlight time windows where the effect of the 'stimclass' factor is
#' # statistically significant according to TANOVA
#' # 1) run TANOVA
#' result_tanova <- tanova(
#'     avgDims(erps, "pairtype"),
#'     list(within = "stimclass", w_id = "id"),
#'     parallel = .(ncores = 2),
#'     perm = .(n = 499))
#'     
#' # 2) extract p-values and bind them to a single array
#' pvalues <- extract(result_tanova, c("p", "p_corr"))
#' pvalues <- bindArrays(pvalues, along_name = "measure")
#' 
#' # 3) plot
#' plotButterfly(tempdat2, sig = pvalues, topo_time = peak_data, 
#'               chan_pos = chan_pos)
#'   
plotButterfly <- function(
  dat, sig = NULL, topo_time = NULL, 
  chan_pos = NULL, subset = list(),
  pcrit = 0.05, aspect_ratio = 0.5, scalp_ratio = 0.5, ampl_range = NULL,
  caption = TRUE, ...) {
  #
  # check arguments and prepare data
  #
  # subset (use also for creating the caption)
  assertList(subset, .var.name = "subset")
  #
  # caption
  assertFlag(caption, .var.name = "caption")
  if (caption) {
    cap <- "Grand averages and GFP curves"
    if (!is.null(topo_time)) {
      cap <- paste0(cap, 
                    ", and representative scalp topographies")
    }
    modterm <- if (!is.null(sig)) dimnames(sig)$modelterm else NULL
    cap <- paste0(cap,
                  " of the ",
                  if (is.null(modterm)) {
                    "effect"
                  } else if (grepl(":", modterm)) {
                    paste0(modterm, " interaction")
                  } else {
                    paste0(modterm, " main effect")
                  })
    if (!is.null(sig)) {
      cap <- paste0(
        cap,
        ". Statistically significant time windows are highlighted ",
        "(shaded rectangles depict time windows where the effect ",
        "was persistent according to the minimum duration criterion).")
    }
    if (length(subset)) {
      cap <- paste0(cap,
                    "Only a subset of the data is shown: ", 
                    deparse(substitute(subset)))
    }
  }
  #
  # create plot
  #
  qp <- plotComplex(dat, sig = sig, topo_time = topo_time, 
                    chan_pos = chan_pos, plot_curves = TRUE, 
                    subset = subset, pcrit = pcrit, 
                    aspect_ratio = aspect_ratio, scalp_ratio = scalp_ratio, 
                    ampl_range = ampl_range, ...)
  #
  # caption
  if (caption) addCaption(qp, cap)
  #
  # return
  qp
}



#' Plot scalp topographies 
#' 
#' \code{plotMap} takes an array (or a vector) of ERP amplitudes and plots 
#' scalp topographies at user-defined time points. 
#' @param dat either a numeric vector or a numeric array of ERP amplitudes. If
#' 'dat' is a named vector, its names are taken as channel names. If unnamed,
#' the order of the data must correspond to the order of channels in 'chan_pos'.
#' If 'dat' is an array, it must have at least 'time' and 'chan' dimensions.
#' @param aspect_ratio the ratio of \code{y} and \code{x} axes on the figure. 
#' If \code{NULL} (the default), it is set automatically. A user-defined 
#' 'aspect_ratio' might be adjusted to avoid overlapping scalp maps.
#' @inheritParams plotComplex
#' @return \code{plotMap} returns a ggplot object.
#' @export
#' @seealso \code{\link{plotButterfly}} for a more complex way of visualizaton
#' @examples
#' # load example data
#' data(erps)
#' 
#' # extract channel positions
#' chan_pos <- attr(erps, "chan")
#' 
#' # collapse pairtypes and participants
#' tempdat <- avgDims(erps, c("pairtype", "id"))
#' 
#' # plot topo-maps of stimclass differences from 0 to 500 in 50 ms steps
#' map <- plotMap(compareLevels(tempdat, "stimclass"),
#'                topo_time = seq(0, 500, by = 50),
#'                chan_pos = chan_pos)
#' 
#' # add title
#' library(ggplot2)
#' map + ggtitle("Pairwise differences between stimulus classes")
#'          
plotMap <- function(dat, chan_pos, topo_time = NULL, subset = list(),
                    aspect_ratio = NULL, ampl_range = NULL, 
                    map_marker_colour = "grey10", map_marker_shape = "|", 
                    map_marker_size = 1.2, ...) {
  assertNumeric(dat, finite = TRUE, any.missing = FALSE)
  assertDataFrame(chan_pos, min.rows = 10L, row.names = "unique")
  # if dat is a vector, it must be transformed to an array 
  if (is.vector(dat)) {
    if (testNamed(dat)) {
      common_channels <- intersect(rownames(chan_pos), names(dat))
      if (length(common_channels) <= 10L) {
        stop(paste0(
          "If 'dat' is a named vector, it must have at least ",
          "10 channel names as provided in 'chan_pos'"),
          call. = FALSE)
      }
      dat <- array(dat[common_channels], 
                   dim = length(common_channels),
                   dimnames = list(chan = common_channels))
    } else if (length(dat) == nrow(chan_pos)) {
      dat <- array(dat, 
                   dim = length(dat),
                   dimnames = list(chan = rownames(chan_pos)))
    } else {
      stop(paste0(
        "If 'dat' is an unnamed vector, its length must equal ",
        "the number of channels as provided in 'chan_pos'"), 
        call. = FALSE)
    }
  }
  # if dat has no "time" dimension, create it 
  # (without copy, dat might be large)
  if (!"time" %in% names(dimnames(dat))) {
    if (length(topo_time) == 0L) {
      topo_time <- "0"
    } else if (length(topo_time) > 1L) {
      stop(paste0("If 'dat' has no time dimension, ",
                  "'topo_time' must be a single value"),
           call. = FALSE)
    }
    origdims <- dim(dat)
    origdimnames <- dimnames(dat)
    on.exit({
      setattr(dat, "dim", origdims)
      setattr(dat, "dimnames", origdimnames)})
    setattr(dat, "dim", c(1L, origdims))
    setattr(dat, "dimnames", c(list(time = topo_time), origdimnames))
  }
  # call plotComplex
  out <- plotComplex(dat = dat, topo_time = topo_time, chan_pos = chan_pos,
                     subset = subset, ampl_range = ampl_range, 
                     plot_curves = FALSE, 
                     aspect_ratio = aspect_ratio,
                     scalp_ratio = 100, 
                     map_marker_colour = map_marker_colour, 
                     map_marker_shape = map_marker_shape, 
                     map_marker_size = map_marker_size,
                     ...)
  # return
  out
}


#' Compute projections for plotting scalp topographies
#' 
#' \code{topoCoord} computes the topographic projection of the amplitudes at 
#' the given peaks.
#' @param dat array of ERP data
#' @param peak_df data.frame of peaks (time and potentially other facetting
#' variables)
#' @param ch_pos data.frame of channel positions
#' @param size_x the width of the topoplot in time units
#' @param size_y the height of the topoplot in amplitude units
#' @param shift_y the vertical center of the topoplot in amplitude units
#' @param r radius of the channel positions if \code{ch_pos} does not contain 
#' it (default: 1)
#' @param ampl_range numeric vector of length 2 containing the minimum and 
#' maximum values. If NULL, it is computed from the data. Values outside the
#' \code{ampl_range} interval are set equal to the corresponding limit.
#' @param resol an integer value; the resolution of the grid topographic map in
#' horizontal and vertical direction (default: 30)
#' @param resolcol integer value, the number of levels (resolution) for 
#' colouring (default: 101L)
#' @param projection character value (default = "laea"). See 
#' \url{http://www.remotesensing.org/geotiff/proj_list/} for common projections.
#' @param projref projection reference (pole [default] or equator)
#' @param origo a named character vector of lat and long coordinates of the 
#' origo
#' @param ... arguments passed to \code{\link{chanInterp}}. You might consider
#' setting the argument \code{N} for dense electrode caps.
#' @note This function is called by \code{\link{plotButterfly}} and 
#' \code{\link{plotMap}}.
#' @return The function returns a named list of two data.table objects (topo
#' and boundary).
#' @keywords internal
topoCoord <- function(
  dat, peak_df, ch_pos, size_x, size_y, shift_y,
  r = 1,
  ampl_range = NULL, resol = 50L, resolcol = 101L, 
  projection = "laea", projref = c("pole", "equator"),
  origo = c(lat = ifelse(projref == "pole", 90, 0),
            long = ifelse(projref == "pole", 270, 0)),
  ...) {
  # helper function to get the interpolated values
  compInterp <- function(x, resol, grid_interp) {
    z <- matrix_(NA_real_, resol, resol)
    z_ind <- chanInterp(x, ch_pos, grid_interp, ...)
    if (is.null(ampl_range)) ampl_range <- range(z_ind, na.rm = TRUE)
    z_ind[z_ind > ampl_range[2]] <- ampl_range[2]
    z_ind[z_ind < ampl_range[1]] <- ampl_range[1]
    z[ind] <- z_ind
    rm(z_ind)
    # return
    as.vector(z)
  }
  #
  projref <- match.arg(projref)
  #
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
  corr_x <- size_x/diff(range(boundarypos$x)) / 1.1
  corr_y <- size_y/diff(range(boundarypos$y)) / 1.1
  xmin <- min(boundarypos$x)
  xmax <- max(boundarypos$x)
  ymin <- min(boundarypos$y)
  ymax <- max(boundarypos$y)
  boundary_polygon <- 
    data.table(xcoord = c(1.1*c(xmin, xmax, xmax, xmin, xmin), 
                          rev(boundarypos$x)),
               ycoord = c(1.1*c(ymin, ymin, ymax, ymax, ymin), 
                          rev(boundarypos$y)))
  len <- nrow(boundary_polygon)
  seq_subs <- 1:nrow(peak_df)
  boundary <- as.data.table(peak_df)[, peak := as.factor(seq_subs)]
  boundary <- boundary[rep(seq_subs, each = len), ]
  boundary[, time := as.numeric(as.character(time))]
  boundary[, xcoord := boundary_polygon$xcoord * corr_x + time]
  boundary[, ycoord := boundary_polygon$ycoord * corr_y + shift_y]
  #
  gridx <- seq(xmin, xmax, length.out = resol)
  gridy <- seq(ymin, ymax, length.out = resol)
  gridpos <- expand.grid(x = gridx, y = gridy)
  ind <- in.polygon(gridpos$x, gridpos$y, 
                    boundarypos$x*1.1, boundarypos$y*1.1)
  gridgeo <- project3dMap(gridpos[ind,], projection = projection, 
                          projref = projref, origo = origo, inverse = TRUE)
  gridcart <- geo2cart(gridgeo)
  out <- peak_df
  out <- setDT(out)[rep(1:nrow(out), each = nrow(gridpos))]
  temp <- lapply(1:nrow(peak_df), function(i) {
    x <- subsetArray(dat, subset. = peak_df[i, ])
    compInterp(as.vector(x), resol, gridcart)
  })
  out[, peak := as.factor(rep(seq_subs, each = resol^2))]
  out[, ampl := unlist(temp, use.names = FALSE)]
  out[, time := as.numeric(as.character(time))]
  out[, xcoord := gridpos$x * corr_x + time]
  out[, ycoord := gridpos$y * corr_y + shift_y]
  # return
  list(topo = out, boundary = boundary)
}




#' Visualize ERP curves and scalp topographies
#' 
#' \code{plotComplex} is a workhorse function called by 
#' \code{\link{plotButterfly}} and \code{\link{plotMap}}.
#' @param dat the array of ERP curves which must have at least 'time' and/or 
#' 'chan' dimensions
#' @param sig [optional] a corresponding array to \code{dat} which is used 
#' to highlight significant time windows. It must have at least 'time' and 
#' measure' dimensions. The 'measure' dimension must have 'p' and 'p_corr' 
#' levels indicating the uncorrected ('p') and corrected ('p_corr') p-values.
#' @param topo_time [optional] an object which describes the time points at 
#' which the scalp topographies should be plotted. Can be a simple atomic vector
#' or a data.frame if the time points are not identical across the faceting
#' dimension(s).
#' @param chan_pos a data.frame of channel positions. Obligatory if 
#' \code{topo_time} is not \code{NULL}.
#' @param subset a named list to subset the input arrays; see 
#' \code{\link[eegR]{subsetArray}} 
#' @param pcrit the level of alpha to highlight significant effects
#' @param aspect_ratio the ratio of \code{y} and \code{x} axes (default: 0.5) 
#' on the figure. If NULL, it is set automatically. If 'topo_time' is provided,
#' user-defined 'aspect_ratio' might be adjusted to avoid overlapping scalp
#' maps.
#' @param scalp_ratio the ratio of the diameter of the scalp and the vertical
#' range of the ERP curves on the figure
#' @param ampl_range the range of amplitudes to plot. If \code{NULL} (default), 
#' it is computed from the data.
#' @param map_marker_colour,map_marker_shape,map_marker_size the colour, shape,
#' and size of the marker used to mark the exact time points of the maps
#' @param ggplot_theme a function which produces a ggplot theme 
#' (default: theme_ndys)
#' @param ... additional arguments passed to \code{\link{topoCoord}}
#' @keywords internal
#' 
plotComplex <- function(
  dat, sig = NULL, topo_time = NULL, chan_pos = NULL, subset = list(),
  plot_curves = TRUE, pcrit = 0.05, aspect_ratio = 0.5, scalp_ratio = 0.5, 
  ampl_range = NULL, 
  map_marker_colour = "red", map_marker_shape = 16, map_marker_size = 1,
  ggplot_theme = theme_ndys,
  ...) {
  #
  # small helper functions
  #
  # subset and drop dimension(s)
  subsetAndDrop <- function(x, subset. = NULL, keep. = c("time", "chan")) {
    if (is.null(subset.)) {
      dropDims(dat, keep = keep.)
    } else {
      dropDims(
        subsetArray(dat, subset. = subset., drop. = FALSE),
        keep = keep.)
    }
  }
  # set xbreaks and ybreaks
  breakSetter <- function(limits) {
    breaks <- pretty(limits)
    if (max(breaks) > round(limits[2])) {
      breaks <- breaks[-length(breaks)]
    }
    if (min(breaks) < round(limits[1])) {
      breaks <- breaks[-1]
    }
    breaks_minor <- (breaks - diff(breaks)[1]/2)[-1]
    list(major = breaks, minor = breaks_minor)
  }
  #
  # check arguments and prepare data
  #
  # theme
  assertFunction(ggplot_theme)
  # aspect_ratio
  if (!is.null(aspect_ratio)) {
    assertNumber(aspect_ratio, lower = .Machine$double.eps,
                 .var.name = "aspect_ratio")
  }
  # subset
  assertList(subset, .var.name = "subset")
  if (!length(subset)) subset <- NULL
  # dat (curvedata)
  assertArray(dat, mode = "numeric", any.missing = FALSE, min.d = 2L,
              .var.name = "dat")
  if (!all(c("time", "chan") %in% names(dimnames(dat)))) {
    stop("no 'time' and 'chan' dimensions in 'dat'", call. = FALSE)
  }
  dat <- subsetAndDrop(dat, subset)
  singletime <- length(dimnames(dat)$time) == 1L
  if (plot_curves) {
    curvedata <- 
      if ("GFP" %in% toupper(dimnames(dat)$chan)) {
        transformArray(ampl ~ ., dat)
      } else {
        transformArray(compGfp(ampl, keep_channels = TRUE) ~ ., dat)
      }
    curvedata$channel <- factor(toupper(curvedata$chan) == "GFP", 
                                labels = c("electrodes", "GFP"))
  }
  #
  # faceting dimensions
  griddims <- setdiff(names(dimnames(dat)), c("time", "chan"))
  # concatenate
  facets <- 
    if (length(griddims) > 1L) {
      gr <- paste(griddims, collapse = "__+__")
      gsub("__+__", " + ", 
           sub("__+__", " ~ ", gr, fixed = TRUE), 
           fixed = TRUE)
    } else {
      paste0(griddims, "~.")
    }
  #
  # ampl_range
  if (is.null(ampl_range)) {
    ampl_range <- range(dat)
  } else {
    assertNumeric(ampl_range, any.missing = FALSE, len = 2L,
                  .var.name = "ampl_range")
    ampl_range <- sort(ampl_range)
  }
  ylim_orig <- ampl_range + 0.05 * c(-1, 1) * diff(ampl_range)
  ylim_ext <- ampl_range + 0.09 * c(-1, 1) * diff(ampl_range)
  #
  # time range
  xlim <- range(as.numeric(dimnames(dat)$time))
  #
  # sig (sigdata)
  if (!is.null(sig)) {
    assertArray(sig, mode = "numeric", any.missing = FALSE, min.d = 2L,
                .var.name = "sig")
    if (!all(c("time", "measure") %in% names(dimnames(sig)))) {
      stop("no 'time' and 'measure' dimensions in 'sig'",
           call. = FALSE)
    }
    if (!all(c("p", "p_corr") %in% dimnames(sig)$measure)) {
      stop(paste0("the 'measure' dimension of 'sig' must have ",
                  "both 'p' and 'p_corr' levels"), call. = FALSE)
    }
    assertNumber(pcrit, lower = 0, upper = 1, .var.name = "pcrit")
    if (length(subset)) {
      sig <- subsetAndDrop(sig, subset)
    }
    sigdata <- sigPhases(sig, pcrit, ylim_orig[1:2])
  }
  #
  # topo_time (topodata)
  if (!is.null(topo_time)) {
    if (is.null(chan_pos)) {
      stop("'chan_pos' must be provided", call. = FALSE)
    }
    assertDataFrame(chan_pos, row.names = "unique")
    coln <- colnames(chan_pos)
    if (!all(c("x", "y", "z") %in% coln) &&
        !all(c("theta", "phi") %in% coln) &&
        !all(c("long", "lat"))) {
      stop("'chan_pos' is not a valid channel position object", 
           call. = FALSE)
    }
    chan_pos <- chan_pos[rownames(chan_pos) %in% dimnames(dat)$chan, ]
    if (nrow(chan_pos) == 0) {
      stop("'chan_pos' contains no common channel with 'dat'")
    }
    topodata <- subsetAndDrop(dat,
                              list(chan = rownames(chan_pos)))
    topo_time <- 
      if (is.vector(topo_time) && is.atomic(topo_time)) {
        if (length(griddims)) {
          expand.grid(c(list(time = as.character(topo_time)),
                        autoConvert(dimnames(topodata)[griddims])),
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
        } else {
          data.frame(time = as.character(topo_time), 
                     stringsAsFactors = FALSE)
        }
      } else if (is.data.frame(topo_time)) {
        if (!"time" %in% colnames(topo_time)) {
          stop(paste0("if 'topo_time' is a data.frame, ",
                      "it must contain a 'time' variable"), 
               .call = FALSE)
        }
        topo_time$time <- as.character(topo_time$time)
        for (n in names(subset)) {
          topo_time <- topo_time[topo_time[[n]] %in% subset[[n]], ]
        }
        for (g in griddims) {
          if (!g %in% colnames(topo_time)) {
            lev <- autoConvert(dimnames(topodata)[[g]])
            topo_time <- topo_time[rep(1:nrow(topo_time), 
                                       each = length(lev)), ]
            topo_time[[g]] <- lev
          }
        }
        topo_time
      } else {
        stop(paste0(
          "'topo_time' must be an atomic object representing ",
          "time points or a data.frame"), call. = FALSE)
      }
    if (uniqueN(topo_time$time) == 1L) singletime <- TRUE
    if (singletime) {
      xlim <- c(xlim[1] - 10, xlim[2] + 10)
    } else if (!plot_curves && is.null(sig)) {
      xlim <- range(as.numeric(topo_time$time))
    }
    extra_y <- diff(ylim_orig)*scalp_ratio
    ylim_ext[2] <- ylim_orig[2] + extra_y * 1.05
    scalpcenter_y <- ylim_orig[2] + extra_y * 1.05/2
    scalpsize_y <- abs(extra_y)
    topotime_diff <- 
      if (singletime) {
        diff(xlim)
      } else if (length(griddims) == 0L) {
        min(diff(as.numeric(topo_time$time)))
      } else {
        min(tapply(as.numeric(topo_time$time), 
                   topo_time[, griddims], 
                   function(x) {
                     if (length(x) == 1L) {
                       Inf
                     } else {
                       min(diff(x))   
                     }
                   }))
      }
    ideal_aspect_ratio <- 
      topotime_diff/diff(xlim)*diff(ylim_ext)/scalpsize_y
    if (is.null(aspect_ratio) || (aspect_ratio > ideal_aspect_ratio)) {
      aspect_ratio <- ideal_aspect_ratio
    } 
    scalpsize_x <- diff(xlim)*aspect_ratio*scalpsize_y/diff(ylim_ext)
    topodata <- suppressMessages(
      topoCoord(topodata, 
                topo_time[, c("time", griddims), drop = FALSE], 
                chan_pos,
                scalpsize_x, scalpsize_y, scalpcenter_y,
                ampl_range = ampl_range, ...))
    if (!singletime) {
      topo_time$time <- as.numeric(topo_time$time)
      if (plot_curves) {
        topo_time$ampl <- NULL
        topo_time <- merge(topo_time, 
                           curvedata[curvedata$chan == "GFP",], 
                           by = c("time", griddims))
      } else {
        topo_time$ampl <- ylim_ext[1]
      }
    }
    xlim0 <- xlim
    xlim <- c(min(xlim[1], min(topodata$boundary$xcoord)),
              max(xlim[2], max(topodata$boundary$xcoord)))
    aspect_ratio <- aspect_ratio * diff(xlim0) / diff(xlim)
  }
  #
  # faceted plots
  #
  # set breaks on x and y axes
  xbreaks <- breakSetter(xlim)
  ybreaks <- breakSetter(ylim_orig)
  # start plot
  qp <- ggplot() 
  # highlight significant phases
  if (!is.null(sig)) {
    if (!is.null(sigdata)) {
      sig_fill <- alpha(brewer.pal(4, "Oranges")[2], 0.1)
      sig_col <- alpha(brewer.pal(4, "Oranges")[2], 0.9)
      ymin <- mean(c(ylim_orig[1], ylim_ext[1]))
      qp <- qp + 
        geom_rect(
          data = sigdata, 
          aes(xmin = time1, xmax = time2, ymax = ymax,
              colour = Significant), 
          fill = sig_fill,
          ymin = ymin) + 
        geom_segment(
          data = sigdata,
          aes(x = time1, xend = time2, linetype = Significant),
          colour = sig_col, size = 1,
          y = ymin, yend = ymin) + 
        scale_linetype_manual(
          name = paste("Significant at", "\u03B1", "=", 
                       sub("^0.", "\\.", 
                           formatC(pcrit, digits = 3))),
          values = c(1, 0), 
          drop = FALSE,
          guide = guide_legend(order = 2)) + 
        scale_colour_manual(
          name = paste("Significant at", "\u03B1", "=", 
                       sub("^0.", "\\.", 
                           formatC(pcrit, digits = 3))),
          values = c("white", sig_col),
          drop = FALSE,
          guide = guide_legend(order = 2))
    } else {
      message("no significant phase has been found")
    }
  }
  # butterfly plot
  if (plot_curves) {
    qp <- qp +
      geom_line(aes(x = time, y = ampl, group = chan, 
                    size = channel, alpha = channel),
                colour = "grey10",
                data = curvedata) + 
      scale_alpha_manual(
        name = "Channel", values = c(0.3, 1),
        guide = guide_legend(order = 1)) + 
      scale_size_manual(
        name = "Channel", values = c(0.2, 0.5),
        guide = guide_legend(order = 1)) + 
      scale_y_continuous(breaks = ybreaks$major,
                         minor_breaks = ybreaks$minor,
                         limits = ylim_ext) + 
      xlab("Time (ms)") + 
      ylab(paste0("Voltage (", "\u03BC", "V)")) + 
      ggplot_theme(panel_aspect_ratio = aspect_ratio)
  } else {
    qp <- qp + 
      ggplot_theme(panel_aspect_ratio = aspect_ratio) + 
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid = element_blank())
    if (singletime && is.null(sig)) {
      qp <- qp +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank())
    }
  }
  # facetting
  if (length(griddims) > 0L) {
    qp <- qp + facet_grid( as.formula(facets) )
  }
  #
  # plot maps
  if (!is.null(topo_time)) {
    if (!singletime) {
      qp <- qp + 
        geom_point(
          aes(x = time, y = ampl), data = topo_time,
          colour = map_marker_colour, 
          shape = map_marker_shape, 
          size = map_marker_size)
    }
    qp <- qp + 
      geom_raster(
        aes(xcoord, ycoord, fill = ampl), 
        data = topodata$topo) + 
      geom_polygon(
        aes(xcoord, ycoord, group = peak), fill = "white", 
        data = topodata$boundary) + 
      scale_fill_gradient2(
        name = paste0("Scalp voltage (", "\u03BC", "V)"),
        guide = guide_colourbar(direction = "horizontal",
                                title.position = "top",
                                barheight = unit(5, "pt"),
                                order = 3),
        low = "blue", mid = "grey97", high = "red", 
        midpoint = 0, 
        space = "Lab")
  }
  #
  # add x scale
  qp <- qp + 
    scale_x_continuous(breaks = xbreaks$major, limits = xlim) + 
    coord_cartesian(xlim = xlim, ylim = ylim_ext)
  #
  # return
  qp
}






#' Multiple plot function
#'
#' \code{multiplot} plots multiple ggplot objects on one page
#' @param ... ggplot objects
#' @param plot_list a list of ggplot objects
#' @param cols the number of columns in layout
#' @param a matrix specifying the layout. If present, 'cols' is ignored.
#' @author Winston Chang
#' @export
#' @keywords internal
multiplot <- function(..., plot_list = NULL, cols = 1L, layout = NULL) {
  requireNamespace("grid")
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plot_list)
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
    grid::grid.newpage()
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))        
      print(plots[[i]], 
            vp = grid::viewport(layout.pos.row = matchidx$row,
                                layout.pos.col = matchidx$col))
    }
  }
}
