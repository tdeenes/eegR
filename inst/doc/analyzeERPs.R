## ----echo = FALSE--------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 8, fig.height = 8)

## ------------------------------------------------------------------------
library(eegR)
data(erps)
?erps

## ------------------------------------------------------------------------
id_dat <- attributes(erps)$id
chan_pos <- attributes(erps)$chan

## ------------------------------------------------------------------------
# subset erps 
dat03 <- subsetArray(erps, list(id = "03", stimclass = c("A", "B")))

# plot curves
plotERParray(dat03, lty = 1, col = "grey70", grid_dim = c(3, 2))

## ------------------------------------------------------------------------
# prepare data
tempdat <- transformArray(
    compGfp(amplitude) ~ . | id, 
    erps, list(readgroup = id_dat$group))

# create plot
library(ggplot2)
ggplot(tempdat, aes(x = time, y = amplitude, col = readgroup, lty = pairtype)) + 
    geom_line() + 
    facet_grid(stimclass ~ .)

## ------------------------------------------------------------------------
# subset the data (pairtype = "ident") and average across id based on 
# reading group membership
tempdat <- transformArray(y ~ . | id, erps, 
                          group = list(readgroup = id_dat$group), 
                          subset = list(pairtype = "ident"), 
                          datfr = FALSE)

# split the array along the readgroup and stimclass dimensions
tempdat <- splitArray(tempdat, c("readgroup", "stimclass"), drop = TRUE)

# compute GFP
gfpdat <- lapply(tempdat, compGfp)

# subset tempdat to include only the time samples at 300 ms
# use [] to keep the attributes of tempdat
plotdat <- tempdat
plotdat[] <- lapply(tempdat, subsetArray, list(time = "300"))

# plot the topographies and the centroids at 300 ms, and add the GFP curves 
# of the whole segments; do not display the channel positions and labels
complexplot2dview(
    plotdat, chan_pos, 300, 
    gfp = gfpdat, plot_ch = FALSE)

## ------------------------------------------------------------------------
# first call chanNb without the alpha parameter
#chanNb(chan_pos)

# now that you know which alpha parameter gives good results, call the 
# function again with the given value
ChN <- chanNb(chan_pos, alpha = 0.7)

## ---- fig.height=14------------------------------------------------------
# do bin averaging
tempdat <- avgBin(erps, "time", bin_length = 6)

# perform the analysis
results <- arrayAnova(
    tempdat, 
    factordef = list(between = "group",
                     w_id = "id", 
                     within = c("stimclass", "pairtype"),
                     observed = "group"),
    bwdat = id_dat, 
    perm = .(n = 99), 
    tfce = .(ChN = ChN),
    parallel = .(ncores = 2))

# plot the results (for now, only the p-values of the TFCE method)
modelplot(results, "p-value")

# print the summary
summary(results)

