import argparse
import random
from Bio import SeqIO
import pdb

def split_fasta(fasta_file, length_start, length_end, output_file, GENOME_SIZE):
    records = SeqIO.parse(fasta_file, "fasta")
    sequences = []

    for record in records:
        seq_len = len(record.seq)
        if length_start <= seq_len < length_end:
            sequences.append(record)

    # Shuffle the sequences to ensure random selection
    random.shuffle(sequences)

    # Calculate the target read count based on the average read length and genome size
    average_read_length = sum(len(record.seq) for record in sequences) / len(sequences)
    target_read_count = int((GENOME_SIZE * 20) / average_read_length)

    # Select the required number of sequences
    selected_sequences = sequences[:target_read_count]

    with open(output_file, "w") as output:
        SeqIO.write(selected_sequences, output, "fasta")

    # Calculate the approximate genome size based on the number of reads and read length
    approximate_genome_size = (target_read_count * average_read_length) / 20

    desired_total_length = int((GENOME_SIZE * 20) / 1000)
    real_total_length = int(sum(len(record.seq) for record in selected_sequences) / 1000)

    print(desired_total_length)
    print(real_total_length)
    print("Desired: {}, Real: {}".format(desired_total_length, real_total_length))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Split a FASTA file based on sequence lengths while maintaining coverage.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file.')
    parser.add_argument('length_start', type=int, help='Start of length range.')
    parser.add_argument('length_end', type=int, help='End of length range.')
    parser.add_argument('output_file', type=str, help='Output FASTA file.')

    args = parser.parse_args()

    GENOME_SIZE = 223680
    split_fasta(args.fasta_file, args.length_start, args.length_end, args.output_file, GENOME_SIZE)


