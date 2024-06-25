""" Setup the databases to be the latest. This will automatciall download any of the resources needed
for removal of human contaminants. We assume paired-end reads are used here as infrustructre is not
setup for single-end at this moment. 

@Author Kathryn Kananen; Nia Tran
"""

#pattern match for any of the reads assuming conventional naming conventions. If a pattern is not
#in these function pairs, then it will not be found. Remember to add a new pattern to both functions
#below for paired reads or you will get an error of a sample not existing and/or not found.
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
rule gunzip_fasta:
    input: 
        in1=rules.download_gencode_hg38.output,
        in2=rules.download_t2t_chm13.output
    output: 
        out1= "resources/gencode/GRCh38.p14.genome.fa",
        out2= "resources/t2t/T2T-CHM13v2.0_genomic.fa"
    conda: "../envs/setup.yml"
    shell:
        """
        gunzip {input.in1} > {output.out1}
        gunzip {input.in2} > {output.out2}
        """

rule zcat_all:
    input:
        in1=rules.gunzip_fasta.output.out1,
        in2=rules.gunzip_fasta.output.out2 
    output:"resources/zcat/merged_ref.fa"
    shell:
        """
        zcat {input.in1} {input.in2} > {output}
        gunzip -c {input} > {output}
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
    input: rules.gunzip_fasta.output
    output: "resources/gencode/merged_ref.3.bt2"
    params:
        prefix="resources/gencode/merged_ref"
    conda: "../envs/setup.yml"
    threads:config["bowtie2"]["threads"]
    shell:
        """
        bowtie2-build --threads {threads} -f -q {input} {params.prefix}
        mv {params.prefix}* $(dirname {output})
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
    output: "/tmp/done.krona"
    conda: "../envs/visualize.yml"
    shell:
        """
        ktUpdateTaxonomy.sh
        touch {output}
        """
