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

######################################

args <- commandArgs(trailingOnly = TRUE)

#cool file
cf <- args[1]

#track file RNA in bigwig format
rna1 <- args[2]
rna2 <- args[3]

#genome (in gff)
genome <- args[4]

#outfile (with extension such as .png or .pdf)
outfile1 <- args[5]

#resolution to plot the contact map
res <- as.numeric(args[6])

#area to zoom in --> could be B.caecimuris or B.caecimuris:190000-210000 
zoom <- args[7]

#area to zoom in --> specify both coordinates (zoom1=lower limit , zoom2 =upper limit , same as zoom area above)
zoom1 <- as.numeric(args[8])
zoom2 <- as.numeric(args[9])

# Define bin size (e.g., 100 bp) for RNA track
bin_size <- as.numeric(args[10])

percentile = .99


######################################

# Import HiC matrices
hic <- HiCExperiment::import(CoolFile(cf), resolution = res, focus=zoom)
hic2 <- hic[pairdist(interactions(hic)) > 2 * res]
max_lim <- quantile(scores(hic2, "balanced"), percentile, na.rm = TRUE)


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
track_rna1 <- rtracklayer::import(rna1, format = "bw", as = "Rle")
track_rna2 <- rtracklayer::import(rna2, format = "bw", as = "Rle")
#track_rna3 <- rtracklayer::import(rna3, format = "bw", as = "Rle")


library(tidyr)
library(dplyr)

# Binning and renaming for each track
a <- bin_signal_rle(track_rna1, zooming, bin_size) %>%
  rename(score_fw = score)

b <- bin_signal_rle(track_rna2, zooming, bin_size) %>%
  rename(score_rv = score)

#c <- bin_signal_rle(track_rna3, zooming, bin_size) %>%
 # rename(score_fw = score)


#import genome track
gtf <- rtracklayer::import(genome, format="gff")
start_pos <- zoom1
end_pos <- zoom2
gr <- as.data.frame(gtf) %>%
  #filter(seqnames == "IV_SynEC") %>%
  filter(start > start_pos) %>%
  filter(start < end_pos) %>%
  filter(end > start_pos) %>%
  filter(end < end_pos) 


#create plot and do the figure

p0 <- plotMatrix(
  hic,
  maxDistance = 300000,
  scale="linear",
  use.scores = "balanced",
  cmap = brewer.pal(n = 9, name = "Reds"),
  limits = c(0, max_lim),
  caption = FALSE
) + coord_cartesian(xlim = c(zoom1, zoom2))  # Set x-axis limit

track_plot_rna1 <- ggplot(data = a, aes(x = coord, y = score_fw)) + 
  geom_area(col = NA, fill = '#ee766f') +
  #coord_cartesian(xlim = c(zoom1, zoom2), ylim = c(0, 250), expand = FALSE) +  # Align x-axis
  coord_cartesian(xlim = c(zoom1, zoom2), expand = FALSE) +  # Align x-axis
  labs(x = NULL)

track_plot_rna2 <- ggplot(data = b, aes(x = coord, y = score_rv)) + 
  geom_area(col = NA, fill = '#5ab5e6') +
  coord_cartesian(xlim = c(zoom1, zoom2), expand = FALSE) +  # Align x-axis
  labs(x = NULL)

genes <- ggplot(data = gr, aes(xmin = start, xmax = end, ymin = 0, ymax = 1)) + 
  geom_rect(aes(color = strand, fill = strand), alpha = 0.4) +
  coord_cartesian(xlim = c(zoom1, zoom2), expand = FALSE)  # Align x-axis

# Combine all plots into one figure
final_plot <- cowplot::plot_grid(p0, track_plot_rna1, track_plot_rna2, genes, ncol = 1, align = "v", axis = "tb")

# Save the output
ggplot2::ggsave(outfile1, plot = final_plot, width = 10, height = 8)
