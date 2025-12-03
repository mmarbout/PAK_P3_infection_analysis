#!/usr/bin/env Rscript

# Import libraries.
library(HiCExperiment)
library(HiContacts)
library(InteractionSet)
library(ggplot2)
library(tibble)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)

cf1 <- args[1]
cf2 <- args[2]
outfile1 <- args[3]
zoom1 <- args[4]
res <- as.numeric(args[5])

cond1 <- cf1
cond2 <- cf2

percentile = .99

# Import HiC matrices
hic0 <- HiCExperiment::import(CoolFile(cond1), resolution = res)
hic1 <- HiCExperiment::import(CoolFile(cond2), resolution = res)

hic3 <- hic0[c(zoom1)]
hic4 <- hic1[c(zoom1)]


# LOG RATIO

div_contacts <- divide(hic4, by = hic3)


p5 <- plotMatrix(
    div_contacts,
    use.scores = "balanced.fc",
    scale = 'log2', 
    limits = c(-2, 2),
    cmap = bwrColors(),
    caption = FALSE,
)    

ggsave(outfile1, p5)
