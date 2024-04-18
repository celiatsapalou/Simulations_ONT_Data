library(Biostrings)

introduce_mutations <- function(bed_file, fasta_file, output_file, mutation_length = 100) {
  # Load the sequences
  sequences <- readDNAStringSet(fasta_file)
  
  # Read the bed file
  bed_data <- read.table(bed_file, header=FALSE, stringsAsFactors=FALSE)
  
  for (i in 1:nrow(bed_data)) {
    chrom <- bed_data[i, 1]
    sv_start <- bed_data[i, 2]
    sv_end <- bed_data[i, 3]
    
    sequence <- as.character(sequences[chrom])
    
    # Decide whether to introduce insertion or deletion
    mutation_type <- sample(c("insert", "delete"), 1)
    
    if (mutation_type == "insert") {
      # Generate random DNA sequence for insertion
      insertion_seq <- paste(sample(c("A", "C", "G", "T"), mutation_length, replace=TRUE), collapse="")
      position <- sample(c(sv_start, sv_end + 1), 1)  # Decide whether to insert before or after the SV
      
      # Insert the sequence
      sequence <- paste0(substr(sequence, 1, position - 1), insertion_seq, substr(sequence, position, nchar(sequence)))
    } else {
      position <- sample(c(sv_start - mutation_length, sv_end + 1), 1)  # Decide whether to delete before or after the SV
      
      # Delete the sequence
      sequence <- paste0(substr(sequence, 1, position - 1), substr(sequence, position + mutation_length, nchar(sequence)))
    }
    
    sequences[chrom] <- DNAStringSet(sequence)
  }
  
  # Write modified sequences to the output fasta file
  writeXStringSet(sequences, output_file)
}

# Example usage:
bed_file <- "path_to_bed_file.bed"
fasta_file <- "path_to_input_fasta.fa"
output_file <- "path_to_output_fasta.fa"

introduce_mutations(bed_file, fasta_file, output_file)
