#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Function to read and prepare data
prepare_data <- function(file, time_point) {
  dat <- read.table(file, header = FALSE, sep = "\t")
  colnames(dat) <- c("X", "Y", "Category")
  dat$TimePoint <- time_point  # Adds a column for the time point
  return(dat)
}

# Function to normalize data for each specific category and time point
normalize_by_category_and_time <- function(data) {
  data <- data %>% 
    group_by(Category, TimePoint) %>% 
    mutate(Y = Y / sum(Y[X >= 3000], na.rm = TRUE)) %>% 
    ungroup()
  return(data)
}

# Get file paths, output file name, specific category, and time points from command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if enough arguments were provided
if (length(args) < 4 || length(args) %% 2 != 0) {
  stop("Usage: Rscript script_name.R output_file.pdf specific_category time1 file1.txt time2 file2.txt ...")
}

# Separate output file name and specific category
output_file <- args[1]
specific_category <- args[2]

# Process time points and file paths
input_args <- args[-c(1, 2)]
time_points <- input_args[seq(1, length(input_args), by = 2)]
file_paths <- input_args[seq(2, length(input_args), by = 2)]

# Ensure the number of time points matches the number of file paths
if (length(time_points) != length(file_paths)) {
  stop("Number of time points must match the number of file paths.")
}

# Read, combine, and normalize files with associated time points
data_list <- mapply(prepare_data, file_paths, time_points, SIMPLIFY = FALSE)
combined_data <- bind_rows(data_list)

# Normalize data for each specific category and time point
normalized_data <- normalize_by_category_and_time(combined_data)

# Filter for the specific category
filtered_data <- normalized_data %>% filter(Category == specific_category)

# Convert factors to avoid issues with ggplot coloring
filtered_data$TimePoint <- as.factor(filtered_data$TimePoint)

# Determine appropriate x-axis upper limit
x_max <- max(filtered_data$X, na.rm = TRUE)

# Plotting and saving to PDF
pdf(output_file, width = 8, height = 6)
ggplot(filtered_data, aes(x = X, y = Y, group = TimePoint)) +
  geom_line(aes(color = TimePoint), linewidth = 0.5) +
  scale_x_log10(limits = c(3000, x_max)) +  # Set x-axis limits dynamically
  scale_y_log10() +
  labs(title = paste("Log-Log Line Plot for Category:", specific_category),
       x = "distance (log scale)",
       y = "Normalized P(s) (log scale)",
       color = "Time Point") +
  theme_minimal()
dev.off()
