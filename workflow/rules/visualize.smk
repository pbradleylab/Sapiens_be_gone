from scripts.utils import *
include:"clean.smk"

rule sam_to_krona:
    input:rules.retain_unmapped.output
    output:"results/{project}/clean/sam_to_krona/{subsample}.xml"
    conda: "../envs/visualize.yml"
    log: "logs/{project}/visualize/sam_to_krona/{subsample}.log"
    shell:
        """
        ktImportSAM {input} -o {output} 2> {log}
        """

rule krona:
    input:rules.sam_to_krona.output
    output:"results/{project}/clean/sam_to_krona/{subsample}.html"
    conda: "../envs/visualize.yml"
    log: "logs/{project}/visualize/sam_to_krona/{subsample}.log"
    shell:
        """
        ktImportText {input} -o {output} 2> {log}
        """