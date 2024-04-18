# Simulations for Benchmarking Small Inversion Detection

## Overview
This project utilizes a Snakefile to simulate Oxford Nanopore Technologies (ONT) data on a genome mutated with small inversions. The objective is to compare various genome mappers and inversion callers to identify the optimal approach for detecting genomic inversions under 1kb.


### Prerequisites
- Snakemake
- Python 3.10.8
- R 4.2.0 or higher
- bedtools
- samtools
- sniffles
- delly
