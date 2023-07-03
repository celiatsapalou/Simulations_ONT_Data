from Bio import SeqIO

fasta_file = "sample1_chr1_cov_20.fasta"
output_files = [
    ("reads_0_3000.fasta", (0, 3000)),
    ("reads_3000_6000.fasta", (3000, 6000)),
    ("reads_6000_10000.fasta", (6000, 10000)),
    ("reads_10000_15000.fasta", (10000, 15000)),
    ("reads_above_15000.fasta", (15000, float("inf")))
]

records = SeqIO.parse(fasta_file, "fasta")

for record in records:
    seq_len = len(record.seq)
    for file, length_range in output_files:
        if length_range[0] <= seq_len < length_range[1]:
            with open(file, "a") as output:
                SeqIO.write(record, output, "fasta")
