""" First remove the contamination using the official SRA-human-scrubber tool. Then align 
using the short read based reference to a combined gencode and T2T. Retain unmapped reads only
and classify them with kraken against human. Visualize all of the samples in krona, so that 
they can be viewed in one place to check for remaining human contamination.

@Author Kathryn Kananen
"""
include:"setup.smk"

# Remove adaptors and low quality reads.
rule fastp:
    input:
        r1=rules.gunzip_r1.output,
        r2=rules.gunzip_r2.output
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

# Run the sra human scrubber for each read to mask human regions.
rule sra_human_scrubber_r1:
    input:
        r1=rules.fastp.output.r1,
        db=rules.download_ncbi_human_db.output
    output:"results/{project}/clean/sra_human_scrubber/{subsample}_r1.fastq"
    log:"logs/{project}/clean/sra_human_scrubber/{subsample}.log"
    conda: "../envs/clean.yml"
    threads:config["sra_human_scubber"]["threads"]
    shell:
        """
        scrub.sh -i {input.r1} -d {input.db} -o - -p {threads} > {output} 2> {log}
        """
rule sra_human_scrubber_r2:
    input:
        r2=rules.fastp.output.r2,
        db=rules.download_ncbi_human_db.output
    output:"results/{project}/clean/sra_human_scrubber/{subsample}_r2.fastq"
    log:"logs/{project}/clean/sra_human_scrubber/{subsample}.log"
    conda: "../envs/clean.yml"
    threads:config["sra_human_scubber"]["threads"]
    shell:
        """
        scrub.sh -i {input.r2} -d {input.db} -o - -p {threads} > {output} 2> {log}
        """

# Align against the cmobined human reference genome.
rule bowtie2:
    input:
        index=rules.bowtie2_index.output,
        ref=rules.gunzip_fasta.output,
        r1=rules.sra_human_scrubber_r1.output,
        r2=rules.sra_human_scrubber_r2.output
    output: "results/{project}/clean/bowtie2/{subsample}.sam"
    params:
        prefix=rules.bowtie2_index.params.prefix
    log: "logs/{project}/clean/bowtie2/{subsample}.log"
    conda: "../envs/clean.yml"
    threads:config["bowtie2"]["threads"]
    shell:
        """
        bowtie2 -x {params.prefix} -1 {input.r1} -2 {input.r2} -S {output} -p {threads} 2> {log}
        """

# Only keep reads that are unmapped. If reads are partially mapped they are removed. Also if they are 
# multi-aligned reads due to using a merged (potentially duplicate sequence) reference they will still
# not be retained in this use-case.
rule retain_unmapped:
    input: rules.bowtie2.output
    output: "results/{project}/clean/retain_unmapped/{subsample}.sam"
    log: "logs/{project}/clean/retain_unmapped/{subsample}.log"
    conda: "../envs/clean.yml"
    shell:
        """
        samtools view -b -f 4 {input} > {output} 2> {log}
        """

rule sam_to_fastq:
    input:rules.retain_unmapped.output
    output:
        r1="results/{project}/clean/sam_to_fastq/{subsample}_r1.fastq",
        r2="results/{project}/clean/sam_to_fastq/{subsample}_r2.fastq",
        singleton="results/{project}/clean/sam_to_fastq/{subsample}_singleton.fastq"
    log: "logs/{project}/clean/sam_to_fastq/{subsample}.log"
    conda: "../envs/clean.yml"
    shell:
        """
        samtools fastq -1 {output.r1} -2 {output.r2} -s {output.singleton} {input} 2> {log}
        """

# Classify against the human kraken2 database.
rule kraken2:
    input:
        db=rules.kraken_build_db.output,
        r1=rules.sam_to_fastq.output.r1,
        r2=rules.sam_to_fastq.output.r2
    output: 
        out="results/{project}/clean/kraken2/{subsample}_kraken2_out.txt",
        report="results/{project}/clean/kraken2/{subsample}_kraken2_report.txt",
        unclassified="results/{project}/clean/kraken2/{subsample}_filtered-reads.fq"
    params:
        db=rules.kraken_build_db.params.dbdir
    log: "logs/{project}/clean/retain_unmapped/{subsample}.log"
    conda: "../envs/clean.yml"
    threads:config["bowtie2"]["threads"]
    shell:
        """
        kraken2 --db {params.db} --threads {threads} --output {output.out} --report {output.report} --unclassified-out {output.unclassified} {input.r1} {input.r2} 2> {log}
        """

# Output the kraken report needed for Krona visualizations.
rule kraken2_percentage_report:
    input:
        db=rules.kraken_build_db.output,
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2
    output: "results/{project}/clean/kraken2_percentage_report/samples_kraken2_out.txt"
    log: "logs/{project}/clean/kraken2_percentage_report/samples_kraken2.log"
    conda: "../envs/clean.yml"
    threads:config["bowtie2"]["threads"]
    shell:
        """
        printf "sample-1\nsample-2\n" > {output}
        scripts/kraken_reads_removed.sh >> {output}
        """
