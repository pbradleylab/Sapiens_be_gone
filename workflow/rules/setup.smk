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


rule download_gencode_hg38:
    output: "resources/gencode/GRCh38.p14.genome.fa.gz"
    params:
        url=config["hg38"]
    conda: "../envs/setup.yml"
    shell:
        """
        wget -c --no-http-keep-alive {params.url} -O {output} 
        """

rule gunzip_fasta:
    input: rules.download_gencode_hg38.output
    output: "resources/gencode/GRCh38.p14.genome.fa"
    conda: "../envs/setup.yml"
    shell:
        """
        gunzip {input}
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

rule bowtie2_index:
    input: rules.gunzip_fasta.output
    output: "resources/gencode/GRCh38.p14.genome.fa.bt2"
    params:
        prefix="GRCh38.p14.genome"
    conda: "../envs/setup.yml"
    shell:
        """
        bowtie2-build -f -q {input} {params.prefix}
        """

rule kraken_download_library:
    output: "resource/kraken/kraken2_human_db/tmp.db"
    conda: "../envs/setup.yml"
    params:
        dbdir="resource/kraken/kraken2_human_db/"
    threads:config["kraken2"]["threads"]
    shell:
        """
        kraken2-build --download-library human --db {params.dbdir} --threads {threads}
        """

rule kraken_download_tax:
    output: "resource/kraken/kraken2_human_db/tax.db"
    conda: "../envs/setup.yml"
    params:
        dbdir="resource/kraken/kraken2_human_db/"
    threads:config["kraken2"]["threads"]
    shell:
        """
        kraken2-build --download-taxonomy --db {params.dbdir}  --threads {threads}
        """

rule kraken_build_db:
    input:
        tax=rules.kraken_download_tax.output,
        library=rules.kraken_download_library.output
    output: "resource/kraken/kraken2_human_db/see.db"
    params:
        dbdir="resource/kraken/kraken2_human_db/"
    conda: "../envs/setup.yml"
    threads:config["kraken2"]["threads"]
    shell:
        """
        kraken2-build --build --db {params.dbdir} --threads {threads}
        """

rule krona_setup:
    output: directory("resource/krona/")
    conda: "../envs/setup.yml"
    shell:
        """
        ktUpdateTaxonomy.sh
        """