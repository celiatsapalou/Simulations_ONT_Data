# Simulations for Benchmarking Small Inversion Detection

## Overview
This project utilizes a Snakefile to simulate Oxford Nanopore Technologies (ONT) data on a genome mutated with small inversions. The objective is to compare various genome mappers (NGMLR vs. minimap2) and inversion callers (Sniffles vs. Delly) to identify the optimal approach for detecting genomic inversions under 1kb.

### System requirement 
- Snakemake 7.32.4 (or higher)
- Python 3.10.8
- SURVIVOR (https://github.com/fritzsedlazeck/SURVIVOR)
- R 4.2.0 or higher
- bedtools
- samtools
- Sniffles2
- Delly

## Reproduction instructions - Input/Output of the pipeline

### 1) Example input data
The `data_simulation` folder contains example data that can be used as an input to run SURVIVOR:
1. The `parameter_file` contains the list of parameters to simulate inversions of selected sizes (and other SVs).
2. The `ref_chr.fa` file is an example reference genome that can be used to simulate the genomic locations of inversions.
3. The `error_profile` file effectively prepares a position-based error profile that will be utilized to simulate reads with realistic error characteristics in subsequent steps of the workflow. This is crucial for generating meaningful and practical data.

### 2) Output
The output of the pipeline will produce 4 plots, each one comparing how many inversions of the selected size have been found using a different combination each time of the mappers and callers that are being tested.
