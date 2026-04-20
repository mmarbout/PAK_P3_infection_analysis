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
library(tidyr)
library(dplyr)

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

hic_compts <- getCompartments(hic,neigens = 5)

eigen <- coverage(metadata(hic_compts)$eigens, weight = 'E1')[[1]]
eigen_df <- tibble(pos = cumsum(runLength(eigen)), eigen = runValue(eigen))

eigen2 <- coverage(metadata(hic_compts)$eigens, weight = 'E2')[[1]]
eigen2_df <- tibble(pos = cumsum(runLength(eigen2)), eigen = runValue(eigen2))

eigen3 <- coverage(metadata(hic_compts)$eigens, weight = 'E3')[[1]]
eigen3_df <- tibble(pos = cumsum(runLength(eigen3)), eigen = runValue(eigen3))

eigen4 <- coverage(metadata(hic_compts)$eigens, weight = 'E4')[[1]]
eigen4_df <- tibble(pos = cumsum(runLength(eigen4)), eigen = runValue(eigen4))

eigen5 <- coverage(metadata(hic_compts)$eigens, weight = 'E5')[[1]]
eigen5_df <- tibble(pos = cumsum(runLength(eigen5)), eigen = runValue(eigen5))


p1 <- plotMatrix(hic, use.scores = 'autocorrelated', scale = 'linear', limits = c(-1, 1), caption = FALSE)

p2 <- ggplot(eigen_df, aes(x = pos, y = eigen)) + 
    geom_area() + 
   # theme_void() + 
    coord_cartesian(expand = FALSE) + 
    labs(x = "Genomic position", y = "Eigenvector value n1")

p3 <- ggplot(eigen2_df, aes(x = pos, y = eigen)) + 
    geom_area() + 
 #   theme_void() + 
    coord_cartesian(expand = FALSE) + 
    labs(x = "Genomic position", y = "Eigenvector value n2")

p4 <- ggplot(eigen3_df, aes(x = pos, y = eigen)) + 
    geom_area() + 
  #  theme_void() + 
    coord_cartesian(expand = FALSE) + 
    labs(x = "Genomic position", y = "Eigenvector value n3")

p5 <- ggplot(eigen4_df, aes(x = pos, y = eigen)) + 
    geom_area() + 
   # theme_void() + 
    coord_cartesian(expand = FALSE) + 
    labs(x = "Genomic position", y = "Eigenvector value n4")

p6 <- ggplot(eigen5_df, aes(x = pos, y = eigen)) + 
    geom_area() + 
 #   theme_void() + 
    coord_cartesian(expand = FALSE) + 
    labs(x = "Genomic position", y = "Eigenvector value n5")

#final_plot <- cowplot::plot_grid(p1, p2,p3,p4,p5,p6,  ncol = 1, align = "v", axis = "tb")
final_plot <- wrap_plots(p1,p2,p3,p4,p5,p6, ncol = 1, heights = c(20,1,1,1,1,1))


ggsave(outfile,final_plot)