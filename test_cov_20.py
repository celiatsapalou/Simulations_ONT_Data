import argparse
import random
from Bio import SeqIO

def split_fasta(fasta_file, length_start, length_end, output_file):
    coverage = 20  # Fixed coverage value
    records = SeqIO.parse(fasta_file, "fasta")
    sequences = []

    for record in records:
        seq_len = len(record.seq)
        if length_start <= seq_len < length_end:
            sequences.append(record)

    # Calculate the target read count based on coverage
    target_read_count = int(coverage * len(sequences))

    # Shuffle the sequences to ensure random selection
    random.shuffle(sequences)

    # Select the required number of sequences
    selected_sequences = sequences[:target_read_count]

    with open(output_file, "w") as output:
        SeqIO.write(selected_sequences, output, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split a FASTA file based on sequence lengths while maintaining coverage.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file.')
    parser.add_argument('length_start', type=int, help='Start of length range.')
    parser.add_argument('length_end', type=int, help='End of length range.')
    parser.add_argument('output_file', type=str, help='Output FASTA file.')

    args = parser.parse_args()

    split_fasta(args.fasta_file, args.length_start, args.length_end, args.output_file)
