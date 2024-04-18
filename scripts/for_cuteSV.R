# Set the CRAN mirror
options(repos = "http://cran.us.r-project.org")

# install.packages("rlang")
# install.packages("readxl")
# install.packages("ggplot2")
#library("readxl")
library("Biostrings")
library("ggplot2")
library("ggbeeswarm")


# Read the command-line arguments
args <- commandArgs(trailingOnly = TRUE)


# Construct the file path with the prefix
infile_path  <- args[1]
outfile_path  <- args[2]

####After running bedtools
inversions_cuteSV_sim <- read.table(infile_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#inversions_sniffles_sim <- as.data.frame(read.table("/g/korbel2/tsapalou/SURVIVOR-master/Debug/for_R.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
colnames(inversions_cuteSV_sim) <- c("chrom_ONT", "start_ONT", "end_ONT", "chrom_cuteSV", "start_cuteSV", "end_cuteSV","overlap")

# Create new column based on start values
inversions_cuteSV_sim$presence <- ifelse(inversions_cuteSV_sim$start_cuteSV== -1, "no", "yes")
inversions_cuteSV_sim$Length_ONT <- inversions_cuteSV_sim$end_ONT - inversions_cuteSV_sim$start_ONT

# Define the colors for "Yes" and "No"
colors <- c("yes" = "white", "no" = "pink")

# Create the boxplot with colors using ggplot2
plot <- ggplot(inversions_cuteSV_sim, aes(x = presence, y = Length_ONT, fill = presence)) +
  geom_boxplot() + geom_beeswarm() +
  scale_fill_manual(values = colors) +
  xlab("cuteSV") +
  ylab("Length") + scale_y_log10()

# Get the directory path
dir_path <- dirname(outfile_path)

print(outfile_path)
# Check if the directory exists
if (!dir.exists(dir_path)) {
  # Create the directory if it doesn't exist
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# Save the plot as a PNG file with the prefix in the filename
ggsave(outfile_path, plot, width = 8.5, height = 8, units = "in", dpi = 300)
