
#rule copy_ref:
#    input:
#        config["reference"],
#    output:
#        OUT + "/reference/reference.fasta",
#    message:
#        "Copying reference genome to output directory"
#    shell:
#        """
#cp {input} {output}
#        """


# TODO: exact copy from apollo-mapping/workflow/rules/map_clean_reads.smk
rule bwa_index_ref:
    input:
        OUT + "/reference/reference.fasta",
    output:
        OUT + "/reference/reference.fasta.sa",
    message:
        "Indexing reference genome using bwa"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/bwa:0.7.17"
    log:
        OUT + "/log/bwa_index_ref.log",
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
bwa index {input} 2>&1>{log}
       """
# TODO: exact copy from apollo-mapping/workflow/rules/map_clean_reads.smk
rule bwa_mem:
    input:
        ref_index=OUT + "/mmidx/{reference_reference.fasta.sa",
        ref=OUT + "/reference/reference.fasta",
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
    output:
        sam=temp(OUT + "/mapped_reads/raw/{sample}.sam"),
    message:
        "Mapping reads for {wildcards.sample}"
    params:
        bases_per_batch="100000000",
        verbosity="3",
        softclip_supp_aln="-Y",
        rgid="'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'",
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/bwa:0.7.17"
    log:
        OUT + "/log/bwa_mem/{sample}.log",
    threads: config["threads"]["bwa"]
    resources:
        mem_gb=config["mem_gb"]["bwa"],
    shell:
        """
bwa mem \
-K {params.bases_per_batch} \
-v {params.verbosity} \
-t {threads} \
{params.softclip_supp_aln} \
-R {params.rgid} \
{input.ref} \
{input.r1} {input.r2} 2>{log} 1>{output}
        """

# TODO: exact copy from apollo-mapping/workflow/rules/map_clean_reads.smk
rule sam_to_sorted_bam:
    input:
        sam=OUT + "/mapped_reads/raw/{sample}.sam",
    output:
        bam=OUT + "/mapped_reads/sorted/{sample}.bam",
    message:
        "Convert sam to sorted bam for {wildcards.sample}"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    log:
        OUT + "/log/sam_to_sorted_bam/{sample}.log",
    threads: config["threads"]["samtools"]
    resources:
        mem_gb=config["mem_gb"]["samtools"],
    shell:
        """
samtools view -b -@ {threads} {input.sam} 2>{log} | \
samtools sort -@ {threads} - 1> {output.bam} 2>>{log}
        """


rule sam_to_sorted_bam:
    input:
        sam=OUT + "/mapped_reads/sorted/{sample}.bam",
    output:
        bam=OUT + "/mapped_reads/sorted/{sample}.bam",
    message:
        "Convert sam to sorted bam for {wildcards.sample}"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    log:
        OUT + "/log/sam_to_sorted_bam/{sample}.log",
    threads: config["threads"]["samtools"]
    resources:
        mem_gb=config["mem_gb"]["samtools"],
    shell:
        """
samtools view -b -@ {threads} {input.sam} 2>{log} | \
samtools sort -@ {threads} - 1> {output.bam} 2>>{log}
        """


rule samtools_stat:
    input:
        OUT + "/mapped_reads/sorted/{sample}.bam"
    output:
        OUT + "/mapped_reads/sorted/{sample}.bam.stats"
    message:
        "samtools stats on mappings of {sample}"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    log:
        OUT + "/log/samtools_stat.log",
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
samtools stats {input} > {output} 2>&1>{log}
        """



rule delimit_to_proper_paired_bam:
    input:
        bam=OUT + "/mapped_reads/sorted/{sample}.bam"
        stats=OUT + "/mapped_reads/sorted/{sample}.bam.stats"
    output:
        OUT + "/mapped_reads/filtered/{sample}.bam"
    message:
        "samtools stats on mappings of {sample}"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    log:
        OUT + "/log/delimit_to_proper_paired_bam.log",
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
echo "todo";
        """

rule bam_to_softclipped_bed:
    input:
        OUT + "/mapped_reads/filtered/{sample}.bam"
    output:
        OUT + "/varia/{sample}.softclipped.bed"
    message:
        "convert SAM/BAM of {sample} to softclipped bed"
    conda:
        "../envs/sam2sofclippedbed.yaml"
    ##container:
    ##    "docker://staphb/samtools:1.17"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
../scripts/convert-sam-to-softclipped-bed.py {input} > {output}
        """

rule bed_to_coverage_bw:
    input:
        bed=OUT + "/varia/{sample}.softclipped.bed"
        fal=OUT + "/reference/reference.fasta.fal"
    output:
        #bg=OUT + "/varia/{sample}.softclipped.bg"
        bw=OUT + "/varia/{sample}.softclipped.bw"
    message:
        "convert sorted bed features into a BigWig file with genomic coverage"
    #conda:
    #    "../envs/TODO-genomecov-TODO-containing-samtools-bedtools-bedGraphToBigWig.yaml"
    ##container:
    ##    "TODO"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
bedtools genomecov -i {input.bed} -g {input.fal} -bg > {output.bw}.bg;
~/software/ucsc/bedGraphToBigWig {output.bw}.bg {input.fal} {output.bw}
rm -f {output.bw}.bg;
        """

rule bam_to_coverage_bw:
    input:
        bed=OUT + "/varia/{sample}.softclipped.bed"
        fal=OUT + "/reference/reference.fasta.fal"
    output:
        #bg=OUT + "/varia/{sample}.softclipped.bg"
        bw=OUT + "/varia/{sample}.softclipped.bw"
    message:
        "convert sorted BAM of a BigWig file with genomic coverage"
    #conda:
    #    "../envs/TODO-genomecov-TODO-containing-samtools-bedtools-bedGraphToBigWig.yaml"
    ##container:
    ##    "TODO"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
bedtools bamtobed -i {input.bam}  | bedtools genomecov -i stdin -g {input.fal} -bg > {output.bw}.bg
~/software/ucsc/bedGraphToBigWig {output.bw}.bg {input.fal} {output.bw}
rm -f {output.bw}.bg;
        """




#rule bam_to_coverage_bw:
#rule softclipped_bed_to_coverage_bw:
#rule subtract_bigwig:
#rule bigwig_to_bed:

# 2. Aggregation rule at the end of the pipeline
rule aggregate_all_softclipped_beds:
    input:
        all_samples = expand(OUT + "/varia/{sample}.softclipped.bed", sample=SAMPLES)
    #output:
    #    final_summary = "results/final_calculation_summary.txt"
    message:
        "aggregate all per-clade softclip-excess regions into a single file"
    conda:
        "../envs/bedtools.yaml"
    #container:
    #    "docker://staphb/samtools:1.17"
    log:
        OUT + "/log/aggregate_all_softclipped_beds.log",
    threads:
        config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        # {input} will be replaced by a space-separated list of all files
        "my_calculation_tool {input.all_samples} > {output.final_summary}"