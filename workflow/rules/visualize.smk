from scripts.utils import *
include:"clean.smk"

def get_reports(wildcards):
    samlst = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        samlst.append(rules.kraken2.output[1].format(project=project, subsample=subsample))
    return(samlst)

rule krona:
    input:
        reports=get_reports,
        db=rules.krona_setup.output
    output:"results/{project}/visualize/sam_to_krona/multi-krona.html"
    conda: "../envs/visualize.yml"
    log: "logs/{project}/visualize/sam_to_krona/multi-krona.log"
    shell:
        """
        ktImportTaxonomy -t 5 -m 3 -o {output} *.report 2> {log}
        """