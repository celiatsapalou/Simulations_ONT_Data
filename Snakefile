rule all:
    input:
        "{prefix}_plot.png"

rule sim_SVs:
    output:
       fasta="{prefix}_invs_sim.fasta"
    shell:
       """ 
       ./SURVIVOR simSV ref_chr1.fa parameter_file 0.1 1 {output.fasta} 
       mv {output.fasta}.fasta {output.fasta}
       """ 

rule sim_reads:
    input:
        fasta = rules.sim_SVs.output.fasta,
    output:
        fasta_out="{prefix}_chr1_cov_20.fasta"
    shell:
        "./SURVIVOR simreads {input.fasta} /g/korbel2/tsapalou/SURVIVOR-master/NA12878_nano_error_profile_bwa.txt 20 {output.fasta_out}"


rule align_ref:
    input:
        fasta="{prefix}_chr1_cov_20.fasta"
    output:
        sam="{prefix}_minimap_aligned.sam"
    shell:
        "minimap2 -ax map-ont ref_chr1.fa {input.fasta} > {output.sam}"

rule sam_to_bam:
    input:
        "{prefix}_minimap_aligned.sam"
    output:
        "{prefix}_minimap.bam"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort_bam:
    input:
        "{prefix}_minimap.bam"
    output:
        "{prefix}_minimap_sorted.bam"
    shell:
        "samtools sort {input} > {output}"

rule index_bam:
    input:
        "{prefix}_minimap_sorted.bam"
    output:
        "{prefix}_minimap_sorted.bam.bai"
    shell:
        "samtools index {input}"

rule run_sniffles:
    input:
        bam="{prefix}_minimap_sorted.bam",
        index="{prefix}_minimap_sorted.bam.bai"
    output:
        "{prefix}_variants_invs.vcf"
    shell:
        "sniffles --input {input.bam} -v {output}"

rule run_delly:
    input:
        bamfile="{prefix}_minimap_sorted.bam",
        index="{prefix}_minimap_sorted.bam.bai"
    output:
        "{prefix}_inv_delly_lr.vcf"
    shell:
        "/g/korbel/shared/software/delly/bin/delly lr -g ref_chr1.fa {input.bamfile} > {output}"

rule sniffles_to_bed:
    input:
        "{prefix}_variants_invs.vcf"
    output:
        "{prefix}_sniffles.bed"
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
        "{prefix}_sniffles.bed"
    output:
        "{prefix}_for_R.bed"
    envmodules: "BEDTools/2.30.0-GCC-11.2.0"
    shell:
        "bedtools intersect -a {input[0]} -b {input[1]} -f 0.7 -r -wao > {output}"

rule run_R_script:
    input:
        "{prefix}_for_R.bed"
    output:
        "{prefix}_plot.png"
    envmodules: "R/4.2.2-foss-2022b"
    shell:
        r"""
        Rscript for_snakemake.R {wildcards.prefix}  > Rscript.log 2>&1 || (echo "R script failed, see Rscript.log for details" && exit 1)
        """
