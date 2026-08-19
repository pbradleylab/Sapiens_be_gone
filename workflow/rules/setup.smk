""" Setup the databases to be the latest. This will automatciall download any of the resources needed
for removal of human contaminants.

@Author Kathryn Kananen; Nia Tran
"""

READ1_PATTERNS = ("_R1", ".R1", ".r1", "_r1", "_1.fq", "_1.fastq")
READ2_PATTERNS = ("_R2", ".R2", ".r2", "_r2", "_2.fq", "_2.fastq")
SEQ_METHODS = ("paired_end", "single_end")


def get_seq_method(wildcards):
    seq_method = get_subsample_attributes(wildcards.subsample, "seq_method", pep)
    if isinstance(seq_method, list):
        if len(seq_method) != 1:
            raise ValueError(
                f"Subsample '{wildcards.subsample}' should have exactly one seq_method, "
                f"found {len(seq_method)}: {seq_method}"
            )
        seq_method = seq_method[0]
    if seq_method not in SEQ_METHODS:
        raise ValueError(
            f"Unsupported seq_method '{seq_method}' for subsample '{wildcards.subsample}'. "
            f"Expected one of: {', '.join(SEQ_METHODS)}."
        )
    return seq_method


def get_reads(wildcards):
    reads = get_subsample_attributes(wildcards.subsample, "reads", pep)
    if isinstance(reads, str):
        reads = [reads]
    reads = list(reads)
    if not reads:
        raise ValueError(f"No FASTQ files found for subsample '{wildcards.subsample}'.")
    return reads


# Pattern match for reads assuming conventional naming conventions. Paired-end samples must
# include recognizable R1 and R2 files. Single-end samples must resolve to exactly one FASTQ.
def get_r1_fastq(wildcards):
    reads = get_reads(wildcards)
    if get_seq_method(wildcards) == "single_end":
        if len(reads) != 1:
            raise ValueError(
                f"Single-end subsample '{wildcards.subsample}' should resolve to exactly "
                f"one FASTQ file, found {len(reads)}: {reads}"
            )
        return reads[0]
    r1 = [x for x in reads if any(pattern in x for pattern in READ1_PATTERNS)]
    if not r1:
        raise ValueError(
            f"Could not find an R1 FASTQ for paired-end subsample '{wildcards.subsample}'. "
            f"Recognized patterns are: {', '.join(READ1_PATTERNS)}."
        )
    return r1[0]


def get_r2_fastq(wildcards):
    if get_seq_method(wildcards) == "single_end":
        raise ValueError(
            f"Single-end subsample '{wildcards.subsample}' does not have an R2 FASTQ."
        )
    reads = get_reads(wildcards)
    r2 = [x for x in reads if any(pattern in x for pattern in READ2_PATTERNS)]
    if not r2:
        raise ValueError(
            f"Could not find an R2 FASTQ for paired-end subsample '{wildcards.subsample}'. "
            f"Recognized patterns are: {', '.join(READ2_PATTERNS)}."
        )
    return r2[0]

  
# Download the latest gencode reference that contains all contigs including unmapped. 
# This can be changed to a static version in the config file by changing the config.
rule download_gencode_hg38:
    output: "resources/gencode/GRCh38.p14.genome.fa.gz"
    params:
        url=config["hg38"]
    conda: "../envs/setup.yml"
    shell:
        """
        wget -c --no-http-keep-alive {params.url} -O {output} 
        """
        
rule download_t2t_chm13:
    output: "resources/t2t/T2T-CHM13v2.0_genomic.fna.gz"
    params:
        url=config["t2t"]
    conda: "../envs/setup.yml"
    shell: 
        """
        wget -c --no-http-keep-alive {params.url} -O {output}
        """

# Used for the ncbi-scrubber
rule download_ncbi_human_db:
    output: "resources/ncbi/human_filter.db"
    params:
        url=config["human_filter"]
    shell:
        """
        wget -c --no-http-keep-alive {params.url} -O {output}
        """

# Needed to decompress for rules in clean.smk (unfortunately).
rule gunzip_grch:
    input: rules.download_gencode_hg38.output
    output:"resources/gencode/GRCh38.p14.genome.fa"
    conda: "../envs/setup.yml"
    shell:
        """
        gunzip {input} > {output}
        """

rule gunzip_t2t:
    input:rules.download_t2t_chm13.output
    output:"resources/t2t/T2T-CHM13v2.0_genomic.fa"
    conda: "../envs/setup.yml"
    shell:
        """
        gunzip {input} > {output}
        """

# Merge the fastas together
rule merge_fastas:
    input:
        in1=rules.gunzip_t2t.output,
        in2=rules.gunzip_grch.output
    output: "resources/zcat/merged_ref.fa"
    shell:
        """
        cat {input.in1} {input.in2} > {output}
        """
        
rule gunzip_r1:
    input: get_r1_fastq
    output: temporary("resources/{project}/gunzip/{subsample}_r1.fastq")
    conda: "../envs/setup.yml"
    shell:
        """
        gunzip -c {input} > {output}
        """
rule gunzip_r2:
    input: get_r2_fastq
    output: temporary("resources/{project}/gunzip/{subsample}_r2.fastq")
    conda: "../envs/setup.yml"
    shell:
        """
        gunzip -c {input} > {output}
        """

# Generate an index for bowtie2
rule bowtie2_index:
    input: rules.merge_fastas.output
    output: "resources/zcat/merged_ref.3.bt2"
    params:
        prefix="resources/zcat/merged_ref"
    conda: "../envs/setup.yml"
    threads:config["bowtie2"]["threads"]
    shell:
        """
        bowtie2-build --threads {threads} -f -q {input} {params.prefix}
        """

# Download all of the necessary kraken resources. These are human libraries and NOT
# prokaryotic. They are used soley for confirmation of the removal of human contamination.
rule kraken_download_library:
    output: "resources/kraken/kraken2_human_db/library/human/assembly_summary.txt"
    conda: "../envs/setup.yml"
    params:
        dbdir="resources/kraken/kraken2_human_db/"
    threads:config["kraken2"]["threads"]
    shell:
        """
        kraken2-build --download-library human --db {params.dbdir} --threads {threads}
        """
        
rule kraken_download_tax:
    output: "resources/kraken/kraken2_human_db/taxonomy/gc.prt"
    conda: "../envs/setup.yml"
    params:
        dbdir="resources/kraken/kraken2_human_db/"
    threads:config["kraken2"]["threads"]
    shell:
        """
        kraken2-build --download-taxonomy --db {params.dbdir}  --threads {threads}
        """
        
rule kraken_build_db:
    input:
        tax=rules.kraken_download_tax.output,
        library=rules.kraken_download_library.output
    output: "resources/kraken/kraken2_human_db/taxo.k2d"
    params:
        dbdir="resources/kraken/kraken2_human_db"
    conda: "../envs/setup.yml"
    threads:config["kraken2"]["threads"]
    shell:
        """
        kraken2-build --build --db {params.dbdir} --threads {threads}
        """
# Setup krona with the standard microbe and non-prokaryotic databases.
rule krona_setup:
    output: "results/temporary/done.krona"
    conda: "../envs/visualize.yml"
    shell:
        """
        ktUpdateTaxonomy.sh
        touch {output}
        """
