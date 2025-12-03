#!/usr/bin/env Rscript

# Import libraries.
library(HiCExperiment)
library(HiContacts)
library(InteractionSet)
library(ggplot2)
library(tibble)
library(RColorBrewer)

# Import variables
args <- commandArgs(trailingOnly = TRUE)

#cool file
cf <- args[1]

#resolution to plot the contact map
res <- as.numeric(args[2])

#area to zoom in --> could be B.caecimuris or B.caecimuris:190000-210000 
zoom <- args[3]

#outfile (with extension such as .png or .pdf)
outfile <- args[4]

percentile = .99

# Import HiC matrices
hic <- HiCExperiment::import(CoolFile(cf),resolution=res)

hic3 <- hic[c(zoom1)]
hic4 <- hic3[pairdist(interactions(hic3)) > 2 * res]

max_lim2 <- quantile(scores(hic4, "balanced"), percentile, na.rm = TRUE)


p1 <- plotMatrix(
    hic3,
    scale = "linear",
    use.scores = "balanced",
    cmap = brewer.pal(n = 9, name = "Reds"),
    limits = c(0, max_lim2),
    caption = FALSE,
)
ggsave(outfile1, p1)