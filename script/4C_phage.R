#!/usr/bin/env Rscript

# Import libraries
library(data.table)
library(Matrix)
library(HiCExperiment)
library(ggplot2)

# Define the 4C-like function
fourc_like <- function(M, frags, reg1, reg2, bin_size, window = 10000, stride = 1000, out_file = NULL) {
    # Convert window and stride to bin units
    window_bins <- ceiling((window / bin_size) / 2)
    stride_bins <- ceiling(stride / bin_size)
    
    # Helper function to parse UCSC regions
    parse_ucsc <- function(region, bins) {
        chrom <- strsplit(region, ":")[[1]][1]
        start_end <- strsplit(strsplit(region, ":")[[1]][2], "-")[[1]]
        start <- as.numeric(start_end[1])
        end <- as.numeric(start_end[2])
        region_start <- which(bins$chromosome == chrom & bins$start >= start)[1]
        region_end <- which(bins$chromosome == chrom & bins$end <= end)[length(which(bins$chromosome == chrom & bins$end <= end))]
        c(region_start, region_end)
    }
    
    # Parse the regions
    reg1_bins <- parse_ucsc(reg1, frags)
    reg2_bins <- parse_ucsc(reg2, frags)
    
    # Debug parsed regions
    cat("Parsed reg1 bins:", reg1_bins, "\n")
    cat("Parsed reg2 bins:", reg2_bins, "\n")
    
    # Check if regions are valid
    if (is.na(reg1_bins[1]) || is.na(reg2_bins[1])) {
        stop("Region parsing failed. Check reg1 and reg2 definitions.")
    }
    
    # Create the 4C-like table
    fourc <- data.table(bin = reg1_bins[1]:reg1_bins[2], val = 0)
    
    # Ensure the matrix is symmetric
    if (!isSymmetric(M)) {
        cat("Matrix is not symmetric. Forcing symmetry...\n")
        M <- M + t(M)
    }
    
    # Compute normalization proportion
    proportion <- sum(M[reg2_bins[1]:reg2_bins[2], reg2_bins[1]:reg2_bins[2]]) /
                  sum(M[reg1_bins[1]:reg1_bins[2], reg1_bins[1]:reg1_bins[2]])
    cat("Normalization proportion:", proportion, "\n")
    
    # Sliding window calculation
    for (i in seq(reg1_bins[1], reg1_bins[2], by = stride_bins)) {
        reg <- c(max(reg1_bins[1], i - window_bins), min(reg1_bins[2], i + window_bins))
        count <- sum(M[reg[1]:reg[2], reg2_bins[1]:reg2_bins[2]])
        count_ref <- sum(M[reg[1]:reg[2], reg1_bins[1]:reg1_bins[2]])
        if (count_ref > 0) {  # Avoid division by zero
            normalized_val <- ((count / length(reg[1]:reg[2])) / (count_ref / length(reg[1]:reg[2]))) / proportion
        } else {
            normalized_val <- 0
        }
        fourc[bin == i, val := normalized_val]
    }
    
    # Write output to file if specified
    if (!is.null(out_file)) {
        fwrite(fourc, file = out_file, sep = "\t", col.names = TRUE)
    }
    
    return(fourc)
}

# Import variables
args <- commandArgs(trailingOnly = TRUE)
mcool_files <- strsplit(args[1], ",")[[1]]  # Comma-separated list of mcool files
resolution <- as.numeric(args[2])  # Resolution
reg1 <- args[3]  # Region 1 in UCSC format
reg2 <- args[4]  # Region 2 in UCSC format
time_points <- strsplit(args[5], ",")[[1]]  # Comma-separated list of time points
out_file <- args[6]  # Output file path
out_graph <- args[7]  # Output graph path

# Initialize an empty list to store results for each time point
results_list <- list()

# Process each mcool file
for (i in seq_along(mcool_files)) {
    # Import Hi-C data
    hic <- HiCExperiment::import(CoolFile(mcool_files[i]), resolution = resolution)
    
    # Extract fragments (bins) and contact matrix
    frags <- as.data.table(bins(hic))
    colnames(frags)[1:3] <- c("chromosome", "start", "end")  # Adjust to match fragment data format
    interactions_data <- as.data.table(interactions(hic))

    print(interactions_data)
    contact_matrix <- sparseMatrix(
        i = interactions_data$bin_id1 + 1,
        j = interactions_data$bin_id2 + 1,
        x = interactions_data$count,
        dims = c(nrow(frags), nrow(frags))
    )
    
    # Run 4C-like analysis
    fourc_result <- fourc_like(contact_matrix, frags, reg1, reg2, resolution, window = 10000, stride = 1000)
    fourc_result[, Time := time_points[i]]  # Add time point information
    
    # Store the result
    results_list[[i]] <- fourc_result
}

# Combine all results into a single data table
combined_results <- rbindlist(results_list)

# Write the combined results to a file
fwrite(combined_results, file = out_file, sep = "\t", col.names = TRUE)

# Generate overlay plot
p <- ggplot(combined_results, aes(x = bin, y = val, color = Time, group = Time)) +
    geom_line(linewidth = 0.3) +  # Line plot for each time point
    theme_minimal() +
    labs(
        x = "Bin Index",
        y = "Normalized 4C Signal",
        color = "Time Point",
        title = "4C-Like Signal Across Time Points"
    ) +
    scale_color_brewer(palette = "Set1")  # Optional: Use a nice color palette

# Save the plot
ggsave(out_graph, plot = p, width = 10, height = 6)