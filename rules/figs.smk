rule FG_recombination_filter:
    input:
        info_idxs=rules.PG_filtered_corealignment.output.info_idxs,
        info_size=rules.PG_filtered_corealignment.output.info_size,
    output:
        full="figs/{dset}/corealn_remove_recombination/full.pdf",
        reduced="figs/{dset}/corealn_remove_recombination/reduced.pdf",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/recombination_filter.py \
            --idxs {input.info_idxs} \
            --size {input.info_size} \
            --fig_full {output.full} \
            --fig_reduced {output.reduced}
        """


rule FG_homoplasies:
    input:
        tree=rules.PG_coregenome_tree.output.nwk,
        aln=rules.PG_corealignment.output.fa,
        aln_info=rules.PG_corealignment.output.json,
        filt_tree=rules.PG_filtered_coregenome_tree.output.nwk,
        filt_aln=rules.PG_filtered_corealignment.output.fa,
        filt_aln_info=rules.PG_filtered_corealignment.output.info_size,
    output:
        out_fld=directory("figs/{dset}/homoplasies"),
    conda:
        "../conda_env/tree_inference.yml"
    shell:
        """
        python3 scripts/figs/homoplasies.py \
            --tree {input.tree} \
            --aln {input.aln} \
            --aln_info {input.aln_info} \
            --filt_tree {input.filt_tree} \
            --filt_aln {input.filt_aln} \
            --filt_aln_info {input.filt_aln_info} \
            --out_fld {output.out_fld}
        """


rule FG_block_distr_fig:
    input:
        rules.PG_polish.output,
    output:
        "figs/{dset}/pangraph/block_distr.pdf",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/pangraph/plot_block_distr.py \
            --pangraph {input} --fig {output}
        """


rule FG_distances:
    input:
        csv=rules.DST_merge.output,
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
    output:
        directory("figs/{dset}/distances"),
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/distances.py \
            --dist_df {input.csv} \
            --tree {input.tree} \
            --fig_fld {output}
        """


rule FG_coresynt:
    input:
        pg=rules.PG_polish.output,
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
    output:
        directory("figs/{dset}/coresynt"),
    params:
        len_thr=config["backbone-joints"]["len-thr"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/fig_core_synteny.py \
            --pangraph {input.pg} \
            --tree {input.tree} \
            --len_thr {params.len_thr} \
            --fig_fld {output}
        """


rule FG_junctions_survey:
    input:
        df=rules.BJ_extract_joints_df.output.dfl,
    output:
        directory("figs/{dset}/junctions_survey"),
    params:
        len_thr=config["backbone-joints"]["len-thr"],
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/backbone_joints/fig_junctions_survey.py \
            --joints_df {input.df} \
            --fig_fld {output}
        """


rule FG_junctions_overview:
    input:
        jdf=rules.BJ_junct_stats.output.stats,
        epg=rules.BJ_extract_pangenome_info.output.info,
    output:
        ff=directory("figs/{dset}/junctions_overview/"),
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/figs/junctions_simple.py \
            --junctions_stats {input.jdf} \
            --edge_pangenome {input.epg} \
            --out_fld {output.ff}
        """


rule FG_all:
    input:
        expand(rules.FG_homoplasies.output, dset=dset_names),
        expand(rules.FG_recombination_filter.output, dset=dset_names),
        expand(rules.FG_block_distr_fig.output, dset=dset_names),
        expand(rules.FG_distances.output, dset=dset_names),
        expand(rules.FG_coresynt.output, dset=dset_names),
        expand(rules.FG_junctions_survey.output, dset=dset_names),
        expand(rules.FG_junctions_overview.output, dset=dset_names),
