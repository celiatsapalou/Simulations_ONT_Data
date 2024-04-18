# Set the CRAN mirror
options(repos = "http://cran.us.r-project.org")

#install.packages("rlang")
#install.packages("readxl")
#install.packages("ggplot2")
library("Biostrings")
library("readxl")
library("ggplot2")
library("ggbeeswarm")
library("ggbeeswarm")


# Read the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract the prefix from the command-line arguments
prefix <- args[1]

# Construct the file path with the prefix
file_path <- paste0(prefix, "_for_R.bed")

####After running bedtools
inversions_sniffles_sim <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#inversions_sniffles_sim <- as.data.frame(read.table("/g/korbel2/tsapalou/SURVIVOR-master/Debug/for_R.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
colnames(inversions_sniffles_sim) <- c("chrom_ONT", "start_ONT", "end_ONT", "chrom_SNIFFLES", "start_SNIFFLES", "end_SNIFFLES","overlap")

# Create new column based on start values
inversions_sniffles_sim$presence <- ifelse(inversions_sniffles_sim$start_SNIFFLES == -1, "no", "yes")
inversions_sniffles_sim$Length_ONT <- inversions_sniffles_sim$end_ONT - inversions_sniffles_sim$start_ONT

# Define the colors for "Yes" and "No"
colors <- c("yes" = "white", "no" = "pink")

# Create the boxplot with colors using ggplot2
plot <- ggplot(inversions_sniffles_sim, aes(x = presence, y = Length_ONT, fill = presence)) +
  geom_boxplot() + geom_beeswarm() +
  scale_fill_manual(values = colors) +
  xlab("Sniffles") +
  ylab("Length") + scale_y_log10()

# Save the plot as a PNG file with the prefix in the filename
ggsave(paste0(prefix, "_plot.png"), plot, width = 8.5, height = 8, units = "in", dpi = 300)

