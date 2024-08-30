import os
import json
import pandas as pd

# create log directory - needed for cluster slurm execution
os.makedirs("log", exist_ok=True)


configfile: "config.yml"


def load_acc_nums(dset_name):
    fname = "datasets/" + dset_name + "/acc_nums.txt"
    with open(fname) as f:
        accnums = f.read().splitlines()
    # remove empty lines and strip whitespaces
    accnums = [acc.strip() for acc in accnums if acc.strip()]
    return accnums


# load config file entries
dsets_config = config["datasets"]
dset_names = list(dsets_config.keys())

# read accession numbers from dataset files
dset_chrom_accnums = {}
for dset_name, dset_info in dsets_config.items():
    dset_chrom_accnums[dset_name] = load_acc_nums(dset_name)

print("Datasets:", dset_names)


wildcard_constraints:
    dset="(" + "|".join(dset_names) + ")",
    acc=r"[^/]+",


rule gbk_to_fa:
    input:
        gbk="data/gbk/{acc}.gbk",
    output:
        fa="data/fa/{acc}.fa",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/utils/gbk_to_fa.py --gbk {input.gbk} --fa {output.fa}
        """


include: "rules/pangraph.smk"
include: "rules/distances.smk"
include: "rules/backbone_joints.smk"
include: "rules/annotations.smk"
include: "rules/rates.smk"
include: "rules/hotspots.smk"
include: "rules/figs.smk"


rule all:
    input:
        expand(rules.FG_homoplasies.output, dset=dset_names),
        expand(rules.FG_recombination_filter.output, dset=dset_names),
        expand(rules.FG_block_distr_fig.output, dset=dset_names),
        expand(rules.FG_large_tree.output, dset=dset_names),
        expand(rules.FG_distances.output, dset=dset_names),
        expand(rules.FG_coresynt.output, dset=dset_names),
        expand(rules.FG_circle_synteny.output, dset=dset_names),
        expand(rules.FG_junctions_survey.output, dset=dset_names),
        expand(rules.FG_junctions_overview.output, dset=dset_names),
        expand(rules.FG_junctions_ann.output, dset=dset_names),
        expand(rules.FG_rates.output, dset=dset_names),


localrules:
    Dfinder_models_download,
    GM_download_db,
