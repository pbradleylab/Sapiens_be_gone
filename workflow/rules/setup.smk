rule download_gencode_hg38:
    output: "resources/gencode/GRCh38.p14.genome.fa.gz"
    params:
        url=config["hg38"]
    conda: "../envs/setup.yml"
    shell:
        """
        wget -c --no-http-keep-alive {params.url} -O {output} 
        """

rule gunzip:
    input: rules.download_gencode_hg38.output
    output: "resources/gencode/GRCh38.p14.genome.fa"
    conda: "../envs/setup.yml"
    shell:
        """
        gunzip {input}
        """

rule bowtie2_index:
    input: rules.gunzip.output
    output: "resources/gencode/GRCh38.p14.genome.fa.bt2"
    params:
        prefix="GRCh38.p14.genome"
    conda: "../envs/setup.yml"
    shell:
        """
        bowtie2-build -f -q {input} {params.prefix}
        """