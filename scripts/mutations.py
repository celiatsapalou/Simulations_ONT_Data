import sys
import random
from Bio import SeqIO

def mutate_sequence(sequence, start, end):
    mutation_length = 100
    # Decide randomly if insertion or deletion
    mutation_type = random.choice(["insert", "delete"])

    # Random sequence generator for insertion
    def random_sequence(length):
        return ''.join(random.choice('ACGT') for _ in range(length))
    
    if mutation_type == "insert":
        insertion = random_sequence(mutation_length)
        # Decide to insert before or after
        if random.choice([True, False]):
            sequence = sequence[:start] + insertion + sequence[start:]
        else:
            sequence = sequence[:end+1] + insertion + sequence[end+1:]
    
    elif mutation_type == "delete":
        # Decide to delete before or after
        if random.choice([True, False]):
            sequence = sequence[:start-mutation_length] + sequence[start:]
        else:
            sequence = sequence[:end+1] + sequence[end+1+mutation_length:]
    
    return sequence

fasta_file = sys.argv[1]
bed_file = sys.argv[2]
output_file = sys.argv[3]

# Read the fasta file
sequences = list(SeqIO.parse(fasta_file, "fasta"))

# Go through the bed file and perform mutations
with open(bed_file, 'r') as bed:
    for line in bed:
        chrom, start, end = line.strip().split()[:3]
        start, end = int(start), int(end)
        # For simplicity, we assume the fasta has only one sequence. For multiple sequences, loop through them.
        sequences[0].seq = mutate_sequence(sequences[0].seq, start, end)

# Write the mutated sequences back to a new fasta file
SeqIO.write(sequences, output_file, "fasta")
