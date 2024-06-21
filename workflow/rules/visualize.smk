""" For visualization of contamination.

@Author Kathryn Kananen
"""
from scripts.utils import *
include:"clean.smk"

# This function is important. This allows for a single rule to be called in the rule all
# IF other visualizations are needed, you can add them here to the list. To keep the 
# code clean and easy to maintain, only update the visualizations or other reports here
# if adding additional steps.
def get_reports(wildcards):
    samlst = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        samlst.append(rules.kraken2.output[1].format(project=project, subsample=subsample))
    return(samlst)

# All samples will be visualized with krona here. Krona will generate a nice html report
# that should show all classification as unclassifed. If something is classified, then 
# that means there is still existing and detectable human sequencing information in the
# sample that man need to be removed.
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
