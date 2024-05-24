include:"setup.smk"

def get_r1_fastq(wildcards):
    out = {} # Dictionary that will hold two reads with r1 in index 0 and r2 in index 1
    reads = get_subsample_attributes(wildcards.subsample, "reads", pep)
    r1=[x for x in reads if ("_R1" in x or ".R1" in x or ".r1" in x or "_r1" in x or "_1.fq" in x or "_1.fastq" in x)]
    return r1[0]

def get_r2_fastq(wildcards):
    out = {} # Dictionary that will hold two reads with r1 in index 0 and r2 in index 1
    reads = get_subsample_attributes(wildcards.subsample, "reads", pep)
    r2=[x for x in reads if ("_R2" in x or ".R2" in x or ".r2" in x or "_r2" in x or "_2.fq" in x or "_2.fastq" in x)]
    return r2[0]

rule sra_human_scrubber:
    input:
        r1=get_r1_fastq,
        r2=get_r2_fastq
    output:
        r1="results/{project}/clean/sra_human_scrubber/{subsample}_r1.fastq",
        r2="results/{project}/clean/sra_human_scrubber/{subsample}_r2.fastq"
    log:
    conda: "../envs/clean.yml"
    shell:
        """
        """

rule fastp:
    input:
        r1=rules.sra_human_scrubber.output.r1,
        r2=rules.sra_human_scrubber.output.r2
    output:
        r1="results/{project}/clean/fastp/{subsample}_r1.fastq",
        r2="results/{project}/clean/fastp/{subsample}_r2.fastq",
        html="results/{project}/clean/fastp/{subsample}.html"
    log: "logs/{project}/clean/fastp/{subsample}.log"
    conda: "../envs/clean.yml"
    threads:config["fastp"]["threads"]
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -w {threads} 2> {log}
        """

rule bowtie2:
    input:
        index=rules.bowtie2_index.output,
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2
    output: "results/{project}/clean/bowtie2/{subsample}.sam"
    log: "logs/{project}/clean/bowtie2/{subsample}.log"
    conda: "../envs/clean.yml"
    threads:config["bowtie2"]["threads"]
    shell:
        """
        bowtie2 -x {input.index} -1 {input.r1} -2 {input.r2} -S {output} -p {threads} 2> {log}
        """

rule retain_unmapped:
    input: rules.bowtie2.output
    output: "results/{project}/clean/retain_unmapped/{subsample}.sam"
    log: "logs/{project}/clean/retain_unmapped/{subsample}.log"
    conda: "../envs/clean.yml"
    threads:config["bowtie2"]["threads"]
    shell:
        """
        samtools view -b -f 4 {input} > {output} 2> {log}
        """

