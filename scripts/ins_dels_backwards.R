library(Biostrings)

#fasta_file <- "/Users/tsapalou/Downloads/for_ins_dels_invs_sim.fasta"
#bed_file <- "/Users/tsapalou/Downloads/for_ins_dels_invs_sim.fasta.bed"
#output_file_modifications <- "/Users/tsapalou/Downloads/output_modifications_4.fa"

# Getting the command line arguments passed to the script
args <- commandArgs(trailingOnly = TRUE)

# Construct the file path with the prefix
infile_fasta <- args[1]
infile_bed <- args[2]
outfile_fasta <- args[3]
outfile_modifications_bed <- args[4]

introduce_insertions <- function(sequence, sv_start, sv_end, mutation_length) {
  # Generate random DNA sequence for insertion
  insertion_seq <- paste(sample(c("A", "C", "G", "T"), mutation_length, replace=TRUE), collapse="")
  
  # Decide whether to insert before or after the inversion
  if (sample(c(TRUE, FALSE), 1)) {
    # Insert before the inversion
    position <- sv_start - mutation_length
    sequence <- paste0(substr(sequence, 1, position), insertion_seq, substr(sequence, position + 1, nchar(sequence)))
    mod_start <- position + 1
    mod_end <- sv_start - 1
  } else {
    # Insert after the inversion
    position <- sv_end
    sequence <- paste0(substr(sequence, 1, position), insertion_seq, substr(sequence, position + 1, nchar(sequence)))
    mod_start <- position + 1
    mod_end <- position + mutation_length
  }
  return(list(sequence=sequence, mod_start=mod_start, mod_end=mod_end, type="insertion"))
}

introduce_deletions <- function(sequence, sv_start, sv_end, mutation_length) {
  # Decide randomly whether to delete before or after the structural variation (SV).
  if (sample(c(TRUE, FALSE), 1)) {

    # Delete right before the SV Start if TRUE
    position <- sv_start - mutation_length - 1 # Here -1 ensures we capture right before the SV start.
    sequence <- paste0(substr(sequence, 1, position), substr(sequence, sv_start, nchar(sequence)))
    mod_start <- position + 1
    mod_end <- sv_start - 1
    return(list(sequence=sequence, mod_start=mod_start, mod_end=mod_end, type="deletion"))
  } else {

    # Delete right after the SV end if FALSE
    position <- sv_end
    sequence <- paste0(substr(sequence, 1, position), substr(sequence, sv_end + mutation_length + 1, nchar(sequence)))
    mod_start <- sv_end + 1
    mod_end <- sv_end + mutation_length
    return(list(sequence=sequence, mod_start=mod_start, mod_end=mod_end, type="deletion"))
  }
}
  

process_bed_file <- function(bed_file, fasta_file, output_file, mutation_length=8000) {
  sequences <- readDNAStringSet(fasta_file)
  
  bed_data <- read.table(bed_file, header=FALSE, stringsAsFactors=FALSE, colClasses=c("character", "numeric", "character", "numeric", "character"))
  
  # Initialize the modifications data frame outside of any loops or conditions
  modifications <- data.frame(chrom=character(), start=integer(), end=integer(), type=character(), stringsAsFactors=FALSE)
  
  # Process from the end of the BED data toward the beginning
  for (i in nrow(bed_data):1) {  # <-- Note the change here
    chrom <- bed_data[i, 1]
    sv_start <- as.numeric(bed_data[i, 2])
    sv_end <- as.numeric(bed_data[i, 4])
    sequence <- as.character(sequences[chrom])
    
 
    # Flip a coin to decide between insertion or deletion
    if (sample(c(TRUE, FALSE), 1)) {
      result <- introduce_insertions(sequence, sv_start, sv_end, mutation_length)
    } else {
      result <- introduce_deletions(sequence, sv_start, sv_end, mutation_length)
    }
    
    sequence <- result$sequence
    sequences[chrom] <- DNAStringSet(sequence)
    
    modifications <- rbind(modifications, data.frame(chrom=chrom, start=result$mod_start, end=result$mod_end, type=result$type))
  }


  # Check the sequences content before writing
  cat("Checking sequences content before writing:\n")
  print(sequences)
  writeXStringSet(sequences, output_file)

  
  # Write the modifications to a BED file
  #write.table(modifications, file=paste0(output_file, "_modifications.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(modifications, file=outfile_modifications_bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  
}

process_bed_file(infile_bed, infile_fasta, outfile_fasta, outfile_modifications_bed)
