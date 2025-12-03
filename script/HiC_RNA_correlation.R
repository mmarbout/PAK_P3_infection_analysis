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

#### parameters ####

args <- commandArgs(trailingOnly = TRUE)

#cool file
cf <- args[1]

#resolution to plot the contact map
res <- as.numeric(args[2])

#area to zoom in --> could be B.caecimuris or B.caecimuris:190000-210000 
zoom <- args[3]

#bigwig file for rna
bw <- args[4]

# Define bin size (e.g., 100 bp) for RNA track
bin_size <- as.numeric(args[5])

#outfile (with extension such as .png or .pdf)
outfile <- args[6]

#number of pixels to compute the HiC signal at short scale
window_size <- as.numeric(args[7])

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

#### handle RNA data ####

# Import the RNA-seq data as Rle
zooming <-GRanges(zoom)
track_rna <- rtracklayer::import(bw, format = "bw", as = "Rle")

# Function to bin data

bin_signal_rle <- function(track, zooming, bin_size) {
  # Ensure the Rle object has data for the full zooming range
  track_zoomed <- track[zooming][[1]]
  
  # Define bin coordinates
  binned_coords <- seq(start(zooming), end(zooming), by = bin_size)
  binned_scores <- sapply(seq_along(binned_coords), function(i) {
    # Define bin range
    range_start <- binned_coords[i]
    range_end <- min(range_start + bin_size - 1, end(zooming))
    
    # Debug: Print bin range
    cat("Processing bin:", range_start, "-", range_end, "\n")
    
    # Extract data for the bin
    bin_data <- track_zoomed[range_start:range_end]
    
    # Debug: Print bin data summary
    cat("Bin data length:", length(bin_data), "Mean:", mean(bin_data, na.rm = TRUE), "\n")
    
    # Compute the mean, handle cases with no data
    if (length(bin_data) > 0) {
      mean(bin_data, na.rm = TRUE)
    } else {
      NA
    }
  })
  
  tibble(coord = binned_coords, score = binned_scores)
}

# extract rna scores

rna_track <- bin_signal_rle(track_rna, zooming, bin_size) %>%
  rename(score_fw = score)

rna_scores <- rna_track$score_fw


smoothed_rna <- smoothListGaussian(rna_scores)


# Function to process and normalize signals
normalize_signal <- function(signal) {
  # Handle NA values by replacing them with 0
  signal[is.na(signal)] <- 0
  # Normalize the signal to have a mean of 0 and SD of 1
  return((signal - mean(signal)) / sd(signal))
}

# Normalize RNA and HiC signals
normalized_rna <- normalize_signal(smoothed_rna)
normalized_hic <- normalize_signal(smoothed_contact)

# Ensure RNA and HiC signals are of the same length
common_length <- min(length(normalized_rna), length(normalized_hic))
normalized_rna <- normalized_rna[1:common_length]
normalized_hic <- normalized_hic[1:common_length]

# Remove NA values and align the data
aligned_data <- na.omit(data.frame(
  RNA = normalized_rna,
  Contact = normalized_hic,
  Index = 1:common_length
))

# Calculate correlation metrics
correlation_pearson <- cor(aligned_data$RNA, aligned_data$Contact, method = "pearson")
correlation_spearman <- cor(aligned_data$RNA, aligned_data$Contact, method = "spearman")

# Create the plot
p <- ggplot(aligned_data, aes(x = Index)) +
  geom_line(aes(y = RNA, colour = "RNA Signal"), na.rm = TRUE) +
  geom_line(aes(y = Contact, colour = "Contact Signal"), na.rm = TRUE) +
  scale_colour_manual(values = c("RNA Signal" = "red", "Contact Signal" = "blue")) +
  labs(
    title = "Normalized Contact and RNA Signals",
    x = "Index",
    y = "Normalized Signal",
    colour = "Signal"
  ) +
  theme_minimal() +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("Pearson: %.2f", correlation_pearson),
    hjust = 1.1,
    vjust = 1.5,
    size = 5,
    color = "black"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = sprintf("Spearman: %.2f", correlation_spearman),
    hjust = 1.1,
    vjust = 3,
    size = 5,
    color = "black"
  )

# Save the plot to a file
ggsave(outfile, plot = p, width = 10, height = 8)

# Print correlation metrics to the console
cat("Pearson correlation:", correlation_pearson, "\n")
cat("Spearman correlation:", correlation_spearman, "\n")
