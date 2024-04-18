library(Biostrings)

fasta_file <- "/Users/tsapalou/Downloads/for_ins_dels_invs_sim.fasta"
bed_file <- "/Users/tsapalou/Downloads/for_ins_dels_invs_sim.fasta.bed"
output_file <- "/Users/tsapalou/Downloads/output_test5.fa"


introduce_mutations <- function(bed_file, fasta_file, output_file, mutation_length = 10000) {
  # Load the sequences
  sequences <- readDNAStringSet(fasta_file)
  
  # bed file
  bed_data <- read.table(bed_file, header=FALSE, stringsAsFactors=FALSE)
  
  #IDENTIFY THE RIGHT COLUMNS IN THE BED FILE
  for (i in 1:nrow(bed_data)) {
    chrom <- bed_data[i, 1]
    sv_start <- as.numeric(bed_data[i, 2])
    sv_end <- as.numeric(bed_data[i, 4])
    
    sequence <- as.character(sequences[chrom])
    
    # Decide whether to introduce insertion or deletion
    mutation_type <- sample(c("insert", "delete"), 1)
    
    if (mutation_type == "insert") {
      
      # Generate random DNA sequence for insertion
      insertion_seq <- paste(sample(c("A", "C", "G", "T"), mutation_length, replace=TRUE), collapse="")
      
      #randomly choose either the SV start or thee SV end position for the insertion
      position <- sample(c(sv_start, sv_end + 1), 1) 
      
      # Insertion - Insert the sequence in the position decided above
      #The position - 1 ensures that you're taking all characters up to, but not including the position.
     
       sequence <- paste0(substr(sequence, 1, position - 1), insertion_seq, substr(sequence, position, nchar(sequence)))
    } else {
      if (sample(c(TRUE, FALSE), 1)) {
        
        ### Delete Before the SV Start if condition is TRUE
        sequence <- paste0(substr(sequence, 1, sv_start - 1 - mutation_length), substr(sequence, sv_start, nchar(sequence)))
      } else {
        #### Delete right after the SV end if condition is FALSE
        sequence <- paste0(substr(sequence, 1, sv_end), substr(sequence, sv_end + mutation_length + 1, nchar(sequence)))
      }
    }
    
    sequences[chrom] <- DNAStringSet(sequence)
  }
  
  # Write modified sequences to the output fasta
  writeXStringSet(sequences, output_file)
}

#### function
introduce_mutations(bed_file, fasta_file, output_file)


