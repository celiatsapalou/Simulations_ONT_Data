rule all:
    input:
        expand("results/{sample}/{sample}.{cell}.bam1", sample=["sample1", "sample2"], cell=["cell1", "cell2"])

rule fastq_to_bam:
    input:
        "cell{number}.fastq1"
    output:
        "cell{number}.bam1"
    shell:
        "echo 'Converting {input} to BAM' && touch {output}"

rule demultiplexing:
    input:
        expand("cell{number}.bam1", number=[1, 2])
    output:
        "demultiplex_table.txt"
    shell:
        """
        echo -e 'cell\\tPool\\t1KG_identified_sample\\tz-score_value\\tSNP_nb\\tTrustable\\n\
        cell1\\tpool1\\tsample1\\t5.95255829\\t19\\tTRUE\\n\
        cell2\\tpool1\\tsample2\\t5.986672013\\t15\\tTRUE' > {output}
        """

checkpoint parse_demultiplex_table:
    input:
        "demultiplex_table.txt"
    output:
        touch("parsed_table.txt")
    script:
        "scripts/parse_table.py"

def gather_samples(wildcards):
    checkpoint_output = checkpoints.parse_demultiplex_table.get(**wildcards).output[0]
    # This script should parse the demultiplex table and return a list of dicts with 'sample' and 'cell' keys
    import pandas as pd
    table = pd.read_csv(checkpoint_output, sep='\t')
    return expand("results/{{sample}}/{{sample}}.{{cell}}.bam1", zip, sample=table['1KG_identified_sample'], cell=table['cell'])

rule organize_samples:
    input:
        gather_samples
    output:
        "results/{sample}/{sample}.{cell}.bam1"
    shell:
        """
        mkdir -p results/{wildcards.sample}
        mv {input} {output}
        """

