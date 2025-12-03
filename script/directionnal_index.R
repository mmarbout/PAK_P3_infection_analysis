#!/usr/bin/env Rscript

# Import libraries.
#library(HiContactsData)
#library(HiCool)
library(GenomicRanges)
library(InteractionSet)
library(HiCExperiment)
library(ggplot2)
library(HiContacts)


# Specify the path to your .mcool file and the desired resolution
args <- commandArgs(trailingOnly = TRUE)
mcool <- args[1] #mcool file
zoom <- args[2] #i.e. LR657304.1:1-6395872
res <- as.numeric(args[3]) # Change this to your desired resolution
nw <- as.numeric(args[4]) #Define the window size
outfile <- args[5] #define the path of your output plotfile

# Load the contact matrix at specified resolution
hic <- import(mcool, resolution= res)
hic2 <- hic[c(zoom)]

# Convert to a regular matrix (if needed and feasible, depending on size)

contact_map <- as.matrix (hic2, use.scores = "balanced", sparse = TRUE)



# Define the directional function
directional <- function(A, nw) {
  n1 <- nrow(A)
  print(paste("Size of the matrix entered for the directional index:", n1))
  signal1 <- matrix(0, nrow = n1, ncol = 1)
  
  for (i in 1:n1) {
    vect_left <- numeric()
    vect_right <- numeric()
    
    for (k in (i-1):(i-nw-1)) {
      kp <- k
      if (k < 1) {
        kp <- n1 + k
      }
      # Check for NA or non-positive values
      if (!is.na(A[i, kp]) && A[i, kp] > 0) {
        vect_left <- c(vect_left, log(A[i, kp]))
      } else {
        vect_left <- c(vect_left, 0)
      }
    }
    
    for (k in (i+1):(i+nw+1)) {
      kp <- k
      if (k > n1) {
        kp <- k - n1
      }
      # Check for NA or non-positive values
      if (!is.na(A[i, kp]) && A[i, kp] > 0) {
        vect_right <- c(vect_right, log(A[i, kp]))
      } else {
        vect_right <- c(vect_right, 0)
      }
    }
    
    if (sum(vect_left) != 0 && sum(vect_right) != 0) {
      signal1[i, 1] <- t.test(vect_right, vect_left, paired = TRUE)$statistic
    } else {
      signal1[i, 1] <- 0
    }
  }
  
  return(signal1)
}


#compute the directionnal index
directional_indices <- directional(contact_map, nw)

# Convert the matrix to a data frame for plotting
directional_data <- data.frame(Index = 1:nrow(directional_indices), Value = directional_indices[,1])

# Create a new variable in the data frame for color mapping
directional_data$Color <- ifelse(directional_data$Value >= 0, "red", "green")

# Plotting using ggplot2
p <- ggplot(directional_data, aes(x = Index, y = Value, fill = Color)) +
  geom_bar(stat = "identity", position = "identity") +  # Using bars to represent indices
  scale_fill_manual(values = c("red" = "red", "green" = "green")) +
  theme_minimal() +
  labs(title = "Directional Indices",
       x = "Index",
       y = "Directional Value") +
  theme(legend.position = "none",  # Hide the legend
        panel.border = element_rect(colour = "black", fill=NA, size=1), # Add black border
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points")) + # Adjust plot margins if needed
  coord_cartesian(ylim = c(-2, 2))  # Strictly enforce y-limits without clipping data

  
ggsave(outfile, p)
