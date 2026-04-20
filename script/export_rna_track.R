#!/usr/bin/env Rscript

# Import libraries.
library(purrr)
library(GenomicRanges)
library(GenomeInfoDb)
library(Biostrings)
library(stringr)
library(HiCExperiment)
library(tibble)
library(HiContacts)
library(plyranges)
library(InteractionSet)
library(ggplot2)
library(RColorBrewer)
library(rtracklayer)
library(cowplot)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

#track file RNA in bigwig format
rna <- args[1]

#outfile (with extension .tsv)
outfile <- args[2]

#area to zoom in --> could be B.caecimuris or B.caecimuris:190000-210000 
zoom <- args[3]

# Define bin size (e.g., 100 bp) for RNA track
bin_size <- as.numeric(args[4])

# Function to bin data

bin_signal_rle <- function(track, zooming, bin_size) {
  # Ensure the Rle object has data for the full zooming range
  track_zoomed <- track[zooming][[1]]

  # Debug: Check the zoomed data range
  zoom_start <- start(zooming)
  zoom_end <- end(zooming)
  cat("Zoomed track length:", length(track_zoomed), "\n")
  cat("Zoomed range:", zoom_start, "-", zoom_end, "\n")

  # Define bin coordinates (absolute coordinates)
  binned_coords <- seq(zoom_start, zoom_end, by = bin_size)

  # Convert absolute coordinates to relative indices
  binned_scores <- sapply(seq_along(binned_coords), function(i) {
    range_start <- binned_coords[i]
    range_end <- min(range_start + bin_size - 1, zoom_end)

    # Convert to relative indices within track_zoomed
    relative_start <- range_start - zoom_start + 1
    relative_end <- range_end - zoom_start + 1

    # Handle out-of-bounds indices
    if (relative_start > length(track_zoomed) || relative_end < 1) {
      cat("Skipping out-of-bound bin:", range_start, "-", range_end, "\n")
      return(NA)
    }
    relative_start <- max(relative_start, 1)
    relative_end <- min(relative_end, length(track_zoomed))

    # Extract data for the bin
    bin_data <- track_zoomed[relative_start:relative_end]
    cat("Processing bin (relative indices):", relative_start, "-", relative_end, "\n")

    if (length(bin_data) > 0) {
      mean(bin_data, na.rm = TRUE)
    } else {
      NA
    }
  })

  # Debug: Check final scores
  cat("Binned scores:", binned_scores, "\n")

  tibble(coord = binned_coords, score = binned_scores)
}

## -- track (bigwig format)

zooming <-GRanges(zoom)

# Apply binning
track_rna <- rtracklayer::import(rna, format = "bw", as = "Rle")


library(tidyr)
library(dplyr)

# Binning and renaming for each track
a <- bin_signal_rle(track_rna, zooming, bin_size) %>%
  rename(score_fw = score)



fwrite(a, file = outfile, sep = "\t", col.names = TRUE)