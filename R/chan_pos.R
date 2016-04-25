#
# <<< channel positions >>> --------
#

#' Transformation between spherical, cartesian and geographical coordinates
#'
#' This package contains several functions which require a data.frame of the 
#' exact channel positions (e.g., \code{\link{plot2dview}}, 
#' \code{\link{chanNb}}). To facilitate the transformation between various
#' coordinate systems, \code{sph2cart}, \code{sph2geo}, \code{cart2sph}, 
#' \code{cart2geo}, \code{geo2sph}, \code{geo2cart} transform spherical (sph), 
#' cartesian (cart), or geographical (geo) coordinates into the respective 
#' coordinate system.
#' @name coordinates
#' @param ch_pos a data frame or matrix containing the spherical, cartesian,
#' or geographical coordinates of the electrode positions. It should contain
#' at least the following (named) columns, respectively: theta and phi;
#' x, y, and z; or lat and long. The rownames attribute of \code{ch_pos} is
#' interpreted as the names of the electrodes.
#' @param r radius (default = 1)
#' @param deg logical variable indicating whether spherical or geographical
#' coordinates are or should be given in degrees (TRUE, default)
#' @param long360 logical variable; if TRUE (default), longitudes range from
#' 0 to 360 (otherwise from -180 to 180)
#' @param orient a character value specifying the orientation in the
#' geographical coordinate system; "northpole" (default) or "equatorial"
#' @return A data.frame with converted coordinates
NULL

#' @rdname coordinates
#' @export
# Transform cartesian to spherical coordinates
sph2cart <- function(ch_pos, r = 1, deg = TRUE) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("theta", "phi") %in% colnames(ch_pos)))
        stop("Either theta or phi angles are missing!")
    if (!is.null(ch_pos$r)) r <- ch_pos$r
    if (deg)
        ch_pos[, c("theta", "phi")] <- ch_pos[, c("theta", "phi")]/180*pi
    theta <- ch_pos$theta
    phi <- ch_pos$phi
    out <- data.frame(
        x = r * sin(theta) * cos(phi),
        y = r * sin(theta) * sin(phi),
        z = r * cos(theta)
    )
    rownames(out) <- rownames(ch_pos)
    # return
    out
}

#' @rdname coordinates
#' @export
# Transform cartesian to spherical coordinates
cart2sph <- function(ch_pos, deg = TRUE) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("x", "y", "z") %in% colnames(ch_pos)))
        stop("Either x, y, or z coordinates are missing!")
    x <- ch_pos$x
    y <- ch_pos$y
    z <- ch_pos$z
    r <- sqrt(x^2 + y^2 + z^2)
    theta <- acos(z/r)
    phi <- atan2(y, x)
    out <- data.frame(theta, phi, r)
    rownames(out) <- rownames(ch_pos)
    if (deg) out[, 1:2] <- out[, 1:2] * 180/pi
    # return
    out
}

#' @rdname coordinates
#' @export
# Transform spherical to geographical coordinates
sph2geo <- function(ch_pos, r = 1, deg = TRUE, long360 = TRUE,
                    orient = c("northpole", "equatorial")) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("theta", "phi") %in% colnames(ch_pos)))
        stop("Either theta or phi angles are missing!")
    if (!is.null(ch_pos$r)) r <- ch_pos$r
    if (!deg) ch_pos <- ch_pos * 180/pi
    orient <- match.arg(orient)
    if (orient == "equatorial") {
        temp <- sph2cart(ch_pos)[, c("z", "x", "y")]
        colnames(temp) <- c("x", "y", "z")
        temp <- cart2sph(temp)
        ch_pos$theta <- temp$theta
        ch_pos$phi <- temp$phi
    }
    theta <- ch_pos$theta
    phi <- ch_pos$phi
    long <- ifelse(theta < 0, 180, 0) + phi
    long[long < 0] <- 360 + long[long < 0]
    if (!long360) long[long > 180] <- long[long > 180] - 360
    lat <- 90 - abs(theta)
    out <- data.frame(long, lat, r)
    rownames(out) <- rownames(ch_pos)
    # return
    out
}

#' @rdname coordinates
#' @export
# Transform geographical to spherical coordinates
geo2sph <- function(ch_pos, r = 1, deg = TRUE) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("lat", "long") %in% colnames(ch_pos)))
        stop("Either latitude [lat] or longitude [long] coordinates are missing!")
    if (!is.null(ch_pos$r)) r <- ch_pos$r
    neglong <- ch_pos$long < 0
    long360 <- if (any(neglong)) F else T
    theta <- (90 - ch_pos$lat) * sign(90 - ch_pos$long) * sign(270 - ch_pos$long)
    theta <- ifelse(ch_pos$lat == 45 & theta == 0, 45, theta)
    phi <- -(ifelse(theta < 0, 180, 0) - ch_pos$long)
    phi <- ifelse(phi > 180, phi - 360, phi)
    out <- data.frame(r = r, theta = theta, phi = phi)
    backcheck <- sum(abs(as.matrix(
        sph2geo(out, long360 = long360)[, c("long", "lat")] -
            ch_pos[, c("long", "lat")]))) < 1e-8
    if (backcheck) {
        if (!deg) out <- out/180 * pi
        return(out)
    } else {
        stop("Oops, something went wrong with the transformation and I don't know why.")
    }
}

#' @rdname coordinates
#' @export
# Transform geographical to cartesian coordinates
geo2cart <- function(ch_pos, r = 1, deg = TRUE) {
    #sph2cart( geo2sph(ch_pos, r, deg) )
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("long", "lat") %in% colnames(ch_pos)))
        stop("Either longitudes (long) or latitudes (lat) are missing!")
    if (!is.null(ch_pos$r)) r <- ch_pos$r
    if (deg) {
        ch_pos[, c("long", "lat")] <- ch_pos[, c("long", "lat")] * pi/180
    }
    long <- ch_pos$long
    lat <- ch_pos$lat
    x <- r * cos(long) * cos(lat)
    y <- r * sin(long) * cos(lat)
    z <- r * sin(lat)
    out <- data.frame(x = x, y = y, z = z)
    rownames(out) <- rownames(ch_pos)
    # return
    out
}

#' @rdname coordinates
#' @export
# Transform cartesian to geographical coordinates
cart2geo <- function(ch_pos, deg = TRUE) {
    ch_pos <- as.data.frame(ch_pos)
    if (!all(c("x", "y", "z") %in% colnames(ch_pos)))
        stop("Either x, y, or z coordinates are missing!")
    x <- ch_pos$x
    y <- ch_pos$y
    z <- ch_pos$z
    r <- sqrt(x^2 + y^2 + z^2)
    long <- atan2(y, x)
    lat <- asin(z/r)
    lat[r == 0] <- 0
    out <- data.frame(long, lat, r)
    rownames(out) <- rownames(ch_pos)
    if (deg) out[, 1:2] <- out[, 1:2] * 180/pi
    # return
    out
}


#' Find channel neighbours
#'
#' \code{chanNb} finds neighbouring channels
#' @param ch_pos a data.frame of electrode positions, 
#' see \code{\link{coordinates}}
#' @param check_alpha a two-element numeric vector defining the range which is
#' supposed to contain the optimal value of alpha
#' @param alpha a numeric value which influences the allowed distance between
#' neighbouring electrodes; if other than NULL (the default), check_alpha is
#' ignored
#' @param ... parameters to \code{\link{sph2cart}}
#' @export
#' @return An electrode neighbourhood matrix
chanNb <- function(ch_pos, check_alpha = c(0.1, 10), alpha = NULL, ...) {
    options(rgl.useNULL = TRUE)
    reqFn(c("alphashape3d", "geometry"))
    if (!all(c("x", "y", "z") %in% colnames(ch_pos))) {
        if (all(c("theta", "phi") %in% colnames(ch_pos))) {
            ch_pos <- sph2cart(ch_pos, ...)
        } else {
            stop(paste0("Channel coordinates should be spherical (polar)",
                        "or cartesian coordinates!"))
        }
    }
    ch_pos <- ch_pos[, c("x", "y", "z")]
    channames <- paste(1:nrow(ch_pos), rownames(ch_pos), sep = ". ")
    if (is.null(alpha)) {
        if (!(requireNamespace("shiny") & requireNamespace("shinyRGL"))) {
            stop(paste0(
                "If you wish to determine alpha by an interactive plot, ",
                "install the 'shinyRGL' package first."), call. = FALSE)
        }
        alpha <- shiny::runApp(list(
            ui = shiny::pageWithSidebar(
                # Application title
                shiny::headerPanel("Find channel neighbours"),
                # Sidebar with a slider input for number of points
                shiny::sidebarPanel(
                    shiny::sliderInput("alpha",
                                "Alpha value: ",
                                min = min(check_alpha),
                                max = max(check_alpha),
                                value = 1, step = 0.1),
                    shiny::actionButton("submit", "Use selected alpha")
                ),
                # Show the generated 3d scatterplot
                shiny::mainPanel(
                    shinyRGL::webGLOutput("chanPlot", height = "700px")
                )
            ),
            server = function(input, output) {
                output$chanPlot <- shinyRGL::renderWebGL({
                    #bg3d("grey40")
                    a <- suppressWarnings(alphashape3d::ashape3d(as.matrix(ch_pos),
                                                                 input$alpha, pert = TRUE))
                    plot(a, walpha = TRUE, transparency = 0.95,
                         col = c("red", "red2", "red"), shininess = 100)
                    rgl::text3d(ch_pos$x, ch_pos$y, ch_pos$z, cex = 0.6,
                                channames, col = "black")
                })
                shiny::observe({
                    if (input$submit == 0)
                        return()
                    shiny::stopApp(input$alpha)
                })
            }
        ))
    }
    a <- suppressWarnings(alphashape3d::ashape3d(as.matrix(ch_pos),
                                                 alpha, pert = TRUE))$edge
    a <- a[,c(1, 2, ncol(a))]
    out <- matrix_(0, nrow(ch_pos), nrow(ch_pos))
    rownames(out) <- rownames(ch_pos)
    mirror <- out > 0
    for (i in 1:nrow(ch_pos)) {
        nb <- sort(unique(c(
            which(mirror[i,]),
            i,
            a[a[,1]==i & a[,3]>0, 2])))
        out[i, 1:length(nb)] <- nb
        mirror[nb, i] <- T
    }
    out <- out[, apply(out > 0, 2, any)]
    # return
    out
}

#' Cosine of the angles between electrodes
#'
#' \code{cosAngle} computes the cosine of angles between electrodes or
#' between any N-dimensional vectors, or between two sets of electrodes or
#' any N-dimensional vectors.
#' @param x a matrix or data.frame with at least two columns and rows
#' @param y If not NULL (default), a matrix or data.frame with at least two
#' columns and rows
#' @param coords a logical variable; if TRUE (default), data in x are
#' coordinates. Coordinates should be named as "x", "y", "z" or "theta", "phi",
#' and if names are present, only coordinate columns are used in the
#' computations.
#' @param units_in_rows a logical variable indicating whether the units (vectors)
#' are in the rows of x (TRUE, default) or in the columns.
#' @param check_params logical; if TRUE (default), the appropriateness of
#' input data are checked before computing the cosine of angles. Set it to
#' FALSE if you really know what you are doing.
#' @export
#' @return A symmetric matrix
cosAngle <- function(x, y = NULL, coords = TRUE, units_in_rows = TRUE,
                     check_params = TRUE) {
    checkFn <- function(xx) {
        if (!units_in_rows) xx <- t(xx)
        if (length(dim(xx)) <= 1 || ncol(xx)<2) {
            stop("Provide a matrix or data.frame with at least 2 columns.")
        }
        if (coords) {
            if (!is.null(colnames(xx))) {
                if (!all(c("x", "y", "z") %in% colnames(xx)) &&
                        all(c("theta", "phi") %in% colnames(xx))) {
                    xx <- sph2cart(xx)
                }
                ind <- na.omit(match(c("x", "y", "z"), colnames(xx)))
                if (length(ind) < 2) {
                    stop("Provide a matrix or data.frame with x and y coordinates.")
                } else {
                    xx <- xx[, ind]
                }
            }
        }
        if (is.data.frame(xx)) xx <- as.matrix(xx)
        return( xx )
    }
    if (is.null(y)) {
        if (check_params) {
            x <- checkFn(x)
        }
        len <- sqrt(rowSums(x^2))
        out <- (x %*% t(x)) / outer(len, len)
        rownames(out) <- colnames(out) <- rownames(x)
    } else {
        if (check_params) {
            x <- checkFn(x)
            y <- checkFn(y)
        }
        lenx <- sqrt(rowSums(x^2))
        leny <- sqrt(rowSums(y^2))
        out <- (x %*% t(y)) / outer(lenx, leny)
        rownames(out) <- rownames(x)
        colnames(out) <- rownames(y)
    }
    # return
    out
}


#' Spherical spline interpolation
#'
#' \code{chanInterp} performs spherical spline interpolation. It can impute
#' bad channels (i.e. channels with missing values) or interpolate to
#' artbitrary positions on the scalp.
#' @param dat a numeric vector, matrix, data.frame or array with named dimnames
#' @param ch_pos channel (electrode) positions; should be a matrix or data.frame
#' with the following column names (order is not important): "theta", "phi"
#' (spherical coordinates), or "x", "y", "z" (cartesian coordinates)
#' @param interp_pos interpolation positions, if they differ from ch_pos;
#' should be in the same format as ch_pos (default: NULL)
#' @param maxNA numeric value given as ratio (0 < maxNA < 1) or integer
#' (maxNA >= 1 ); it determines the maximum ratio or number of non-missing
#' channels in a time sample. If exceeded, no interpolation occurs in the given
#' time sample.
#' @param m integer value (default: 4); the flexibility of the spline
#' @param N integer value (default: 7); number of terms in the Legendre
#' polynomial. You should probably increase it (even up to 100L) for
#' realistic (i.e. non-spherical) or high-density electrode arrangements.
#' @param lambda numeric value (default: 1e-10); smoothing factor
#' @param type character value, do not use yet
#' @param alarm_tolerance numeric value (default: 1e-2); if the maximal absolute
#' interpolation error at any time sample exceeds this limit, a message is
#' shown or an error is thrown depending on \code{error_on_alarm}. If set to
#' NULL, no check is performed.
#' @param error_on_alarm defaults to TRUE, see \code{alarm_tolerance}
#' @export
#' @return An object having the same attributes as dat
chanInterp <- function(dat, ch_pos, interp_pos = NULL, maxNA = 0.3,
                       m = 4L, N = 7L, lambda = 1e-10,
                       type = c("voltage", "laplacian", "scd"),
                       alarm_tolerance = 1e-2, error_on_alarm = TRUE) {
    message("\n****\nStart interpolation...\n")
    # Function to compute G matrices:
    # G = g(cos(ch_pos)) / interpG = g(cos(ch_pos, interp_pos))
    compGmat <- function(cos_angles, m, N) {
        n <- seq.int(N)
        series <- (2*n + 1) / (n^m * (n+1)^m)
        k <- 1/4 / pi
        P <- polynomial.values(
            legendre.polynomials(N, normalized = FALSE)[-1],
            cos_angles)
        G <- k * Reduce("+", mapply("*", P, series, SIMPLIFY = FALSE))
        return( G )
    }
    # Function to compute C coefficients
    compCmat <- function(x, G, tol) {
        Gx <- array_(1, dim(G) + 1)
        Gx[-1, -1] <- G
        Gx[1, 1] <- 0
        diag(Gx) <- diag(Gx) + tol
        invG <- corpcor::pseudoinverse(Gx, tol)[, -1]
        return( invG %*% x )
    }
    # Function to perform spherical spline interpolation
    interpFn <- function(y, m, N, lambda, type) {
        na_ind <- is.na(y)
        message("...Looking for missing value patterns - ", appendLF = FALSE)
        if (is.null(interp_pos)) {
            temp <- colMeans(na_ind)
            keep_columns <- temp <= maxNA & temp > 0
            if (!any(keep_columns)) return(y)
            yorig <- y
            y <- y[, keep_columns, drop = FALSE]
            na_ind <- na_ind[, keep_columns, drop = FALSE]
            missing_patterns <- fastUnique(na_ind, margin = 2L)
            mischan <- which(rowAnys(missing_patterns))
            mischan <-
                if (!is.null(rownames(ch_pos))) {
                    paste(rownames(ch_pos)[mischan], collapse = ", ")
                } else {
                    paste(seq_along(mischan), collapse = ", ")
                }
            message("Done")
            message("...There are missing values in the following channels: ",
                    mischan)
        } else {
            missing_patterns <- fastUnique(na_ind, margin = 2L)
            message("Done")
            yy <- array_(0, c(nrow(interp_pos), ncol(y)),
                         list(rownames(interp_pos), colnames(y)))
            names(dimnames(yy)) <- names(dimnames(y))
        }
        if (!is.null(alarm_tolerance)) {
            dev <- rep(FALSE, ncol(y))
            names(dev) <- colnames(y)
            maxdev <- 0
        }
        message("...Perform interpolation - ", appendLF = FALSE)
        G0 <- compGmat(cosAngle(ch_pos, check_params = FALSE), m, N)
        for (i in 1:ncol(missing_patterns)) {
            na_vec <- missing_patterns[, i]
            y_ind <- colSums(na_ind==na_vec) == length(na_vec)
            ch_good <- ch_pos[!na_vec, , drop = F]
            ch_interp <-
                if (is.null(interp_pos)) ch_pos[na_vec,,drop = F] else interp_pos
            #
            G <- G0[!na_vec, !na_vec]
            Coef <- compCmat(y[!na_vec, y_ind, drop = F], G, lambda)
            interpG <- compGmat(cosAngle(ch_good, ch_interp, check_params = FALSE),
                                m, N)
            if (is.null(interp_pos)) {
                res <- crossprod(Coef[-1, , drop = F], interpG) + Coef[1, ]
                y[na_vec, y_ind] <- t(res)
            } else {
                yy[, y_ind] <- t(crossprod(Coef[-1, , drop = F], interpG) +
                                     Coef[1, ])
            }
            if (!is.null(alarm_tolerance)) {
                ch_interp <- ch_pos
                interpG <- compGmat(cosAngle(ch_good,
                                             ch_interp, check_params = FALSE),
                                    m, N)
                res <- crossprod(Coef[-1, , drop = F], interpG) + Coef[1, ]
                mdevs <- colMaxs(abs(y[, y_ind] - t(res)), na.rm = TRUE)
                maxdev <- max(c(maxdev, mdevs))
                dev[y_ind] <- temp <- ( maxdev > alarm_tolerance )
                if (error_on_alarm && temp)
                    stop(sprintf("Deviance %f exceeds threshold %f",
                                 maxdev, alarm_tolerance))
            }
        }
        if (!is.null(alarm_tolerance)) {
            if (any(dev)) {
                devind <- which(dev)
                devind <- devind[1:min(c(20, length(devind)))]
                message(
                    "****\nThe maximum absolute interpolation error exceeds the limit
                    at the following time samples (only the first 20 are shown):\n-----\n",
                    paste(names(dev)[devind], collapse = " ")
                )
                message(paste("Maximal deviation = ", maxdev, sep = ""))
            }
        }
        message("Done")
        if (is.null(interp_pos)) {
            yorig[, keep_columns] <- y
            return(yorig)
        } else {
            return(yy)
        }
    }
    #
    type <- match.arg(type)
    if (is.null(lambda)) {
        lambda <- if (type == "voltage") 1e-7 else 1e-5
    }
    #
    # check if ch_pos has an appropriate format
    if (!all(c("x", "y", "z") %in% colnames(ch_pos))) {
        if (all(c("theta", "phi") %in% colnames(ch_pos))) {
            ch_pos <- sph2cart(ch_pos)
        } else {
            stop("Provide channel locations in spherical or cartesian coordinates!")
        }
    }
    ch_pos <- as.matrix(ch_pos[, c("x", "y", "z")])
    #
    # do the same for the interpolation locations
    if (!is.null(interp_pos)) {
        if (!all(c("x", "y", "z") %in% colnames(interp_pos))) {
            if (all(c("theta", "phi") %in% colnames(interp_pos))) {
                interp_pos <- sph2cart(interp_pos)
            } else {
                stop("Provide interpolation locations in spherical or cartesian coordinates!")
            }
        }
        interp_pos <- as.matrix(interp_pos[, c("x", "y", "z")])
    }
    #
    # set up interpolation infos
    argnames <- setdiff(names(as.list(args(chanInterp))), "dat")
    argnames <- argnames[argnames  != ""]
    procstep <- list(
        what = "interpolation", call = match.call(),
        options = mget(argnames), nr_of_missings = sum(is.na(dat)))
    #
    # check input data and run interpolation
    dim_names <- dimnames(dat)
    if (is.list(dat) || !is.numeric(dat)) {
        stop("Provide a numeric vector, matrix, data.frame or array as input!")
    } else if (is.null(interp_pos) &&
                   procstep$nr_of_missings == 0) {
        setattr(dat, "processing_steps",
                c(attr(dat, "processing_steps"), list(procstep)))
        message("...No missing data to interpolate.\nInterpolation finished.")
        return(dat)
    } else if (is.vector(dat)) {
        if (length(dat) != nrow(ch_pos)) {
            stop("Number of channels and number of data points do not match!")
        }
        out <- interpFn(as.matrix(dat), m, N, lambda, type)
    } else {
        if (!any(dim(dat) == nrow(ch_pos))) {
            stop("Number of channels and data size do not match!")
        }
        target_dim <- if (!is.null(dim_names$chan)) "chan" else 1
        arg_list <- list(m = m, N = N, lambda = lambda, type = type)
        if (length(arg_list) == 0) arg_list <- NULL
        newdims <-
            if (is.null(interp_pos)) list() else list(chan = rownames(interp_pos))
        out <- fnDims(dat, "chan", interpFn, arg_list = arg_list,
                      newdims = newdims, vectorized = TRUE)
        out <- apermArray(out, names(dim_names),
                          keep_attributes. = TRUE)
    }
    #
    if (is.null(interp_pos)) {
        NAs_left <- sum(is.na(out))
        procstep$nr_of_interp <-
            procstep$nr_of_missings - NAs_left
        message("Number of NAs after bad channel interpolation: ",
                NAs_left)
    }
    message("Interpolation finished.")
    setattr(out, "processing_steps",
            c(attr(dat, "processing_steps"), list(procstep)))
    # return
    out
    }

#' Project channel positions onto 2D plane
#'
#' \code{project3dMap} projects 3D coordinates onto a 2D plane or vice versa
#' @param pos electrode positions (matrix or data.frame)
#' @param r radius
#' @param projection character string (default = "laea"). See
#' \url{http://www.remotesensing.org/geotiff/proj_list/} for common projections
#' @param projref projection reference (pole [ = default] or equator)
#' @param origo a named character vector of lat ( = latitude) and long (longitude)
#' @param inverse if set to TRUE, back-projection is performed (default = FALSE)
#' @export
#' @return A data.frame containing the projected coordinates
project3dMap <- function(pos, r = 1,
                         projection = "laea", projref = c("pole", "equator"),
                         origo = c(lat = ifelse(projref == "pole", 90, 0),
                                   long = ifelse(projref == "pole", 270, 0)),
                         inverse = FALSE) {
    #
    xyfn <- function(dat) data.frame(x = dat[,1], y = dat[,2])
    projfn <- function(p, projcall, ...) {
        res <- lapply(unique(projcall), function(x) {
            ind <- projcall == x
            out <- project(as.data.frame(p[ind, , drop = F]), proj = x, ...)
            out <- as.data.frame(out)
            rownames(out) <- which(ind)
            return( out )
        })
        res <- Reduce("rbind", res)
        res <- res[order(as.numeric(rownames(res))), , drop = F]
        return( xyfn(res) )
    }
    #
    if (!is.data.frame(pos) & is.list(pos)) {
        stop("Pos can be a vector, matrix, or data.frame.")
    }
    if (is.atomic(pos) && is.vector(pos)) {
        if (length(pos)  != 2) {
            stop("If pos is a vector, it must contain exactly two elements!")
        }
        pos <- cbind(x = pos[1], y = pos[2])
    }
    pos <- as.data.frame(pos)
    if (is.null(pos$r)) pos$r <- 1
    projref <- match.arg(projref)
    proj4call <- paste0("+proj=", projection,
                        " +lat_0=", origo[1],
                        " +lon_0=", origo[2],
                        paste0(" +R=", pos$r))
    #
    if (inverse) {
        if (is.null(colnames(pos))) colnames(pos)[1:2] <- c("x", "y")
        if (!all(c("x", "y") %in% colnames(pos))) {
            stop("Provide a matrix or data.frame with x and y coordinates!")
        }
        mapxy <- projfn(as.matrix(pos[, c("x", "y")]), proj4call, inverse = TRUE)
        colnames(mapxy) <- c("long", "lat")
    } else {
        if (!all(c("long", "lat") %in% colnames(pos))) {
            if (!all(c("theta", "phi") %in% colnames(pos)))
                pos <- cart2sph(pos)
            if (projref == "equator") {
                temp <- sph2cart(pos)
                colnames(temp) <- c("y", "z", "x")
                pos <- cart2sph(temp)
            }
            pos <- sph2geo(pos, long360 = FALSE)
        }
        mapxy <- projfn(pos[, c("long", "lat")], proj4call)
    }
    rownames(mapxy) <- rownames(pos)
    # return
    mapxy
}
