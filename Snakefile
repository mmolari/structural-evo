import os
import json
import pandas as pd

# create log directory - needed for cluster slurm execution
os.makedirs("log", exist_ok=True)


configfile: "config.yml"


# load config file entries
dsets_config = config["datasets"]
dset_names = list(dsets_config.keys())

# read accession numbers from dataset files
dset_chrom_accnums = {}
for dset_name, dset_info in dsets_config.items():
    fname = "datasets/" + dset_name + "/assembly_to_chrom.tsv"
    df = pd.read_csv(fname, sep="\t")
    acc_nums = df["chromosome_acc"].tolist()
    dset_chrom_accnums[dset_name] = acc_nums

# load accession numbers of excluded isolates
excluded = {k: [] for k in dset_names}
for dset_name, dset_info in dsets_config.items():
    fname = f"datasets/{dset_name.strip()}/excluded.txt"
    acc_nums = []
    try:
        with open(fname, "r") as f:
            acc_nums = f.readlines()
        acc_nums = [an.strip() for an in acc_nums]
        acc_nums = [an for an in acc_nums if len(an) > 0]
    except:
        pass
    excluded[dset_name] = acc_nums


wildcard_constraints:
    dset="(" + "|".join(dset_names) + ")",
    acc=r"[^/]+",


include: "rules/downloads.smk"
include: "rules/pangraph.smk"
include: "rules/distances.smk"
include: "rules/backbone_joints.smk"
include: "rules/annotations.smk"
include: "rules/rates.smk"
include: "rules/figs.smk"

rule all:
    input:
        expand(rules.FG_homoplasies.output, dset=dset_names),
        expand(rules.FG_recombination_filter.output, dset=dset_names),
        expand(rules.FG_block_distr_fig.output, dset=dset_names),
        expand(rules.FG_distances.output, dset=dset_names),
        expand(rules.FG_coresynt.output, dset=dset_names),
        expand(rules.FG_circle_synteny.output, dset=dset_names),
        expand(rules.FG_junctions_survey.output, dset=dset_names),
        expand(rules.FG_junctions_overview.output, dset=dset_names),
        expand(rules.FG_junctions_ann.output, dset=dset_names),
        expand(rules.FG_rates.output, dset=dset_names),

localrules:
    download_gbk,
    Dfinder_models_download,
    GM_download_db,
