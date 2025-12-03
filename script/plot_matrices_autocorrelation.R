#!/usr/bin/env Rscript

# Import libraries.
library(HiCExperiment)
library(rtracklayer)
library(ggplot2)
library(GenomicRanges)
library(tidyr)
library(dplyr)
library(purrr)
library(GenomeInfoDb)
library(Biostrings)
library(stringr)
library(tibble)
library(HiContacts)
library(plyranges)
library(InteractionSet)
library(RColorBrewer)
library(cowplot)
library(patchwork)

#### parameters ####

args <- commandArgs(trailingOnly = TRUE)

#cool file
cf <- args[1]

#resolution to plot the contact map
res <- as.numeric(args[2])

#area to zoom in --> could be B.caecimuris or B.caecimuris:190000-210000 
zoom <- args[3]

#outfile (with extension such as .png or .pdf)
outfile <- args[4]



#### handle hic matrix ####

# Import HiC matrices
hic <- HiCExperiment::import(CoolFile(cf), resolution = res, focus=zoom)
hic <- autocorrelate(hic)

p1 <- plotMatrix(hic, use.scores = 'autocorrelated', scale = 'linear', limits = c(-1, 1), caption = FALSE)


ggsave(outfile,p1)