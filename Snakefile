rule all:
    input:
        expand("res/R_script_output/{prefix}_{length_start}-{length_end}_Delly_plot.png", prefix=config["prefix"], length_start=config["length_start"], length_end=config["length_end"]),
        expand("res/R_script_output/{prefix}_{length_start}-{length_end}_sniffles_plot.png", prefix=config["prefix"], length_start=config["length_start"], length_end=config["length_end"])

rule sim_SVs:
    input:
        parameter_file="data/parameter_file",
        ref_fa="data/ref_chr1.fa"
    output:
        fasta="res/sim_SVs_output/{prefix}_invs_sim.fasta",
        bed="res/sim_SVs_output/{prefix}_invs_sim.fasta.bed"
    shell:
        """
        scripts/SURVIVOR simSV {input.ref_fa} {input.parameter_file} 0.1 1 {output.fasta}
        mv {output.fasta}.fasta {output.fasta}
        """

rule error_profile:
    input:
        survivor_config="data/NA12878_nano_error_profile_bwa.txt"
    output:
        "data/{prefix}_reads_{length_start}-{length_end}.txt"
    shell:
         "awk -v min={wildcards.length_start} -v max={wildcards.length_end} 'BEGIN{{OFS=FS=\"\t\"}} NR==1{{print $0; next}} NR<min{{$2=0; print $0}} NR>=min && NR<=max{{$2=(NR-min)/(max-min); print $0}} NR>max{{$2=1; print $0}}' {input} > {output}"

rule sim_reads:
    input:
        fasta = "res/sim_SVs_output/{prefix}_invs_sim.fasta",
        survivor_config = "data/{prefix}_reads_{length_start}-{length_end}.txt" 
    output:
        fasta_out= "res/sim_reads_output/{prefix}_{length_start}-{length_end}_chr1_cov_20.fasta"
    params:
        simread_param = 20
    shell:
        "scripts/SURVIVOR simreads {input.fasta} {input.survivor_config} {params.simread_param} {output.fasta_out}"
        

rule split_fasta_by_length:
    input:
        all_sim_reads= "res/sim_reads_output/{prefix}_{length_start}-{length_end}_chr1_cov_20.fasta"
    output:
        processed_fasta= "res/split_fasta_output/{prefix}_reads_{length_start}-{length_end}.fasta"
    shell:
        """
        cp {input} {output}
        #~/mambaforge/bin/python scripts/split_reads.py {input.all_sim_reads} {wildcards.length_start} {wildcards.length_end} {output.processed_fasta}
        """

rule align_ref:
    input:
        fasta= "res/split_fasta_output/{prefix}_reads_{length_start}-{length_end}.fasta",
        ref_fa="data/ref_chr1.fa"
    output:
        sam= "res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap_aligned.sam"
    shell:
        "minimap2 -ax map-ont {input.ref_fa} {input.fasta} > {output.sam}"

rule sam_to_bam:
    input:
         "res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap_aligned.sam"
    output:
         "res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap.bam"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort_bam:
    input:
         "res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap.bam"
    output:
         "res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap_sorted.bam"
    shell:
        "samtools sort {input} > {output}"

rule index_bam:
    input:
         "res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap_sorted.bam"
    output:
         "res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap_sorted.bam.bai"
    shell:
        "samtools index {input}"

rule run_sniffles:
    input:
        bam="res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap_sorted.bam",
        index="res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap_sorted.bam.bai"
    output:
        "res/sniffles_output/{prefix}_{length_start}-{length_end}_variants_invs.vcf"
    shell:
        "sniffles --input {input.bam} -v {output}"

rule run_delly:
    input:
        bamfile="res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap_sorted.bam",
        index="res/align_ref_output/{prefix}_{length_start}-{length_end}_minimap_sorted.bam.bai",
        ref_fa = "data/ref_chr1.fa"
    output:
        "res/delly_output/{prefix}_{length_start}-{length_end}_inv_delly_lr.vcf"
    shell:
        "/g/korbel/shared/software/delly/bin/delly lr -g {input.ref_fa} {input.bamfile} > {output}"

rule sniffles_to_bed:
    input:
        "res/sniffles_output/{prefix}_{length_start}-{length_end}_variants_invs.vcf"
    output:
        "res/sniffles_to_bed_output/{prefix}_{length_start}-{length_end}_sniffles.bed"
    shell:
        r"""
        awk -F '\t' '!/^#/ {{split($8,a,";"); for(i in a){{if(a[i] ~ /END=/){{split(a[i],b,"="); print $1 "\t" $2 "\t" b[2]}}}}}}' "{input}" > "{output}"
        """

rule sim_bed_format:
    input:
        "res/sim_SVs_output/{prefix}_invs_sim.fasta.bed"
    output:
        "res/simulated_invs_bed/{prefix}_invs_ONT.bed"
    shell:
        """
        cut -f1,2,4 {input} > {output}
        """

rule delly_to_bed:
    input:
        "res/delly_output/{prefix}_{length_start}-{length_end}_inv_delly_lr.vcf"
    output:
        "res/delly_to_bed_output/{prefix}_{length_start}-{length_end}_delly.bed"
    shell:
        r"""
        awk -F'\t' 'BEGIN {{OFS=FS}} $1 !~ /^#/ {{split($8,a,";"); for(i in a) {{split(a[i],b,"="); if(b[1]=="END") print $1, $2, b[2]}}}}' "{input}" > "{output}"
        """


rule bedtools_intersect_delly:
     input:
         "res/simulated_invs_bed/{prefix}_invs_ONT.bed",
         "res/delly_to_bed_output/{prefix}_{length_start}-{length_end}_delly.bed"
     output:
        "res/bedtools_intersect_delly/{prefix}_{length_start}-{length_end}_for_delly.bed"
     envmodules: "BEDTools/2.30.0-GCC-11.2.0"
     shell:
         "bedtools intersect -a {input[0]} -b {input[1]} -f 0.7 -r -wao > {output}"

rule bedtools_intersect_sniffles:
    input:
        "res/simulated_invs_bed/{prefix}_invs_ONT.bed",
        "res/sniffles_to_bed_output/{prefix}_{length_start}-{length_end}_sniffles.bed"
    output:
        "res/bedtools_intersect_output/{prefix}_{length_start}-{length_end}_for_R.bed"
    envmodules: "BEDTools/2.30.0-GCC-11.2.0"
    shell:
        "bedtools intersect -a {input[0]} -b {input[1]} -f 0.7 -r -wao > {output}"


rule run_R_script_sniffles:
    input:
        "res/bedtools_intersect_output/{prefix}_{length_start}-{length_end}_for_R.bed"
    output:
        "res/R_script_output/{prefix}_{length_start}-{length_end}_sniffles_plot.png"
    envmodules: "R/4.2.2-foss-2022b"
    shell:
        r"""
        /g/korbel/hoeps/anaconda_2022/miniconda3/envs/nahrwhals-big/bin/Rscript scripts/plot_sniffles_results.R {input} {output}
        """

rule run_R_script_delly:
    input:
        "res/bedtools_intersect_delly/{prefix}_{length_start}-{length_end}_for_delly.bed"
    output:
        "res/R_script_output/{prefix}_{length_start}-{length_end}_Delly_plot.png"
    envmodules: "R/4.2.2-foss-2022b"
    shell:
        r"""
        /g/korbel/hoeps/anaconda_2022/miniconda3/envs/nahrwhals-big/bin/Rscript scripts/for_Delly.R {input} {output}
        """