configfile: "config.yaml"

SAMPLES = config["samples"]
REF = config["reference"]

rule all:
    input:
        expand("results/qc/fastqc/{sample}_fastqc.zip", sample=SAMPLES),
        "results/qc/multiqc/multiqc_report.html",
        expand("results/{sample}.vcf.gz", sample=SAMPLES)

rule fastqc:
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/qc/fastqc/{sample}_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    shell:
        "fastqc --outdir results/qc/fastqc {input}"

rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_fastqc.zip", sample=SAMPLES)
    output:
        "results/qc/multiqc/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc results/qc/fastqc -o results/qc/multiqc"

rule trim:
    input:
        "data/{sample}.fastq.gz"
    output:
        "results/trim/{sample}.trim.fastq.gz"
    conda:
        "envs/fastp.yaml"
    shell:
        "fastp -i {input} -o {output} --json results/trim/{wildcards.sample}.fastp.json --html results/trim/{wildcards.sample}.fastp.html"

rule map:
    input:
        reads="results/trim/{sample}.trim.fastq.gz",
        ref=REF
    output:
        bam="results/{sample}.sorted.bam"
    threads: 4
    conda:
        "envs/bwa_samtools.yaml"
    shell:
        "bwa mem {input.ref} {input.reads} | samtools sort -o {output.bam}"

rule index_bam:
    input:
        "results/{sample}.sorted.bam"
    output:
        "results/{sample}.sorted.bam.bai"
    conda:
        "envs/bwa_samtools.yaml"
    shell:
        "samtools index {input}"

rule coverage:
    input:
        "results/{sample}.sorted.bam"
    output:
        "results/{sample}.coverage.txt"
    conda:
        "envs/bwa_samtools.yaml"
    shell:
        "samtools depth -a {input} > {output}"

rule mpileup_call:
    input:
        bam="results/{sample}.sorted.bam",
        ref=REF
    output:
        vcf="results/{sample}.vcf.gz"
    threads: 4
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools mpileup -f {input.ref} {input.bam} | bcftools call -mv -Oz -o {output.vcf} && bcftools index {output.vcf}"
