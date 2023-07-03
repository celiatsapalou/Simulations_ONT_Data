rule all:
    input:
        expand("{prefix}_{length_start}-{length_end}_plot.png", prefix=config["prefix"], length_start=config["length_start"], length_end=config["length_end"])

rule sim_SVs:
    output:
        fasta="{prefix}_invs_sim.fasta",
        bed="{prefix}_invs_sim.fasta.bed"
    shell:
        """
        ./SURVIVOR simSV ref_chr1.fa parameter_file 0.1 1 {output.fasta}
        mv {output.fasta}.fasta {output.fasta}
        """

rule sim_reads:
    input:
        fasta = rules.sim_SVs.output.fasta
    output:
        fasta_out="{prefix}_chr1_cov_20.fasta"
    shell:
        "./SURVIVOR simreads {input.fasta} /g/korbel2/tsapalou/SURVIVOR-master/NA12878_nano_error_profile_bwa.txt 500 {output.fasta_out}"


rule split_fasta_by_length:
    input:
        fasta_file="{prefix}_chr1_cov_20.fasta"
    output:
        processed_fasta="{prefix}_reads_{length_start}-{length_end}.fasta"
    shell:
        """
        ~/mambaforge/bin/python test_cov_20.py {input.fasta_file} {wildcards.length_start} {wildcards.length_end} {output.processed_fasta}
        """

rule align_ref:
    input:
        fasta="{prefix}_reads_{length_start}-{length_end}.fasta"
    output:
        sam="{prefix}_{length_start}-{length_end}_minimap_aligned.sam"
    shell:
        "minimap2 -ax map-ont ref_chr1.fa {input.fasta} > {output.sam}"

rule sam_to_bam:
    input:
        "{prefix}_{length_start}-{length_end}_minimap_aligned.sam"
    output:
        "{prefix}_{length_start}-{length_end}_minimap.bam"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort_bam:
    input:
        "{prefix}_{length_start}-{length_end}_minimap.bam"
    output:
        "{prefix}_{length_start}-{length_end}_minimap_sorted.bam"
    shell:
        "samtools sort {input} > {output}"

rule index_bam:
    input:
        "{prefix}_{length_start}-{length_end}_minimap_sorted.bam"
    output:
        "{prefix}_{length_start}-{length_end}_minimap_sorted.bam.bai"
    shell:
        "samtools index {input}"


rule run_sniffles:
    input:
        bam="{prefix}_{length_start}-{length_end}_minimap_sorted.bam",
        index="{prefix}_{length_start}-{length_end}_minimap_sorted.bam.bai"
    output:
        "{prefix}_{length_start}-{length_end}_variants_invs.vcf"
    shell:
        "sniffles --input {input.bam} -v {output}"

rule run_delly:
    input:
        bamfile="{prefix}_{length_start}-{length_end}_minimap_sorted.bam",
        index="{prefix}_{length_start}-{length_end}_minimap_sorted.bam.bai"
    output:
        "{prefix}_{length_start}-{length_end}_inv_delly_lr.vcf"
    shell:
        "/g/korbel/shared/software/delly/bin/delly lr -g ref_chr1.fa {input.bamfile} > {output}"

rule sniffles_to_bed:
    input:
        "{prefix}_{length_start}-{length_end}_variants_invs.vcf"
    output:
        "{prefix}_{length_start}-{length_end}_sniffles.bed"
    shell:
        r"""
        awk -F '\t' '!/^#/ {{split($8,a,";"); for(i in a){{if(a[i] ~ /END=/){{split(a[i],b,"="); print $1 "\t" $2 "\t" b[2]}}}}}}' "{input}" > "{output}"
        """

rule sim_bed_format:
    input:
        "{prefix}_invs_sim.fasta.bed"
    output:
        "{prefix}_invs_ONT.bed"
    shell:
        """
        cut -f1,2,4 {input} > {output}
        """

rule bedtools_intersect:
    input:
        "{prefix}_invs_ONT.bed",
        "{prefix}_{length_start}-{length_end}_sniffles.bed"
    output:
        "{prefix}_{length_start}-{length_end}_for_R.bed"
    envmodules: "BEDTools/2.30.0-GCC-11.2.0"
    shell:
        "bedtools intersect -a {input[0]} -b {input[1]} -f 0.7 -r -wao > {output}"

rule run_R_script:
    input:
        "{prefix}_{length_start}-{length_end}_for_R.bed"
    output:
        "{prefix}_{length_start}-{length_end}_plot.png"
    envmodules: "R/4.2.2-foss-2022b"
    shell:
        r"""
        Rscript for_snakemake.R {wildcards.prefix}_{wildcards.length_start}-{wildcards.length_end} > Rscript.log 2>&1 || (echo "R script failed, see Rscript.log for details" && exit 1)
        """
