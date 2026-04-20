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
library(data.table)

#### parameters ####

args <- commandArgs(trailingOnly = TRUE)

#cool file
cf <- args[1]

#resolution to plot the contact map
res <- as.numeric(args[2])

#area to zoom in --> could be B.caecimuris or B.caecimuris:190000-210000 
zoom <- args[3]

#outfile (with extension .tsv)
outfile <- args[4]

#number of pixels to compute the HiC signal at short scale
window_size <- as.numeric(args[5])

#### handle hic matrix ####

# Import HiC matrices
hic <- HiCExperiment::import(CoolFile(cf), resolution = res, focus=zoom)
contact_map <- as.matrix(hic, use.scores = "balanced", sparse = TRUE)

#extract diagonals signals
# Assuming 'contact_map' is your square matrix
n <- nrow(contact_map)  # Number of rows in the matrix
num_diagonals <-  window_size # Number of diagonals above and below the main diagonal to consider

# Initialize the sum vector with zeros, length equal to the main diagonal
diagonal_sum <- numeric(n)

# Extract and sum the diagonals manually
for (offset in -num_diagonals:num_diagonals) {
  if (offset == 0) {
    diag_values <- contact_map[cbind(1:n, 1:n)]
  } else if (offset > 0) {
    row_indices <- 1:(n - offset)
    col_indices <- (1 + offset):n
    diag_values <- contact_map[cbind(row_indices, col_indices)]
  } else {  # offset < 0
    row_indices <- (1 - offset):n
    col_indices <- 1:(n + offset)
    diag_values <- contact_map[cbind(row_indices, col_indices)]
  }

  # Add to the diagonal sum vector, adjusting the length for shorter diagonals
  if (offset >= 0) {
    diagonal_sum[(1 + abs(offset)):n] <- diagonal_sum[(1 + abs(offset)):n] + diag_values
  } else {
    diagonal_sum[1:(n + offset)] <- diagonal_sum[1:(n + offset)] + diag_values
  }
}


# Assuming smoothListGaussian function is defined to handle numeric vectors
smoothListGaussian <- function(x, degree = 5) {
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector")
  }
  window <- degree*2-1
  weights <- rep(1/window, window)
  stats::filter(x, filter = weights, circular = TRUE, sides = 2)
}

# Smoothing the extracted diagonal
smoothed_contact <- smoothListGaussian(diagonal_sum)

hic_signal <- list(smoothed_contact)

fwrite(hic_signal, file = outfile, sep = "\n", col.names = TRUE)
