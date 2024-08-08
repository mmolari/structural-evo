rule HS_extract_annotations:
    input:
        dfinder=expand(
            rules.AN_assign_positions.output,
            tool="defensefinder",
            K="real",
            allow_missing=True,
        ),
        genomad=expand(
            rules.AN_assign_positions.output,
            tool="genomad",
            K="real",
            allow_missing=True,
        ),
        isescan=expand(
            rules.AN_assign_positions.output,
            tool="ISEScan",
            K="real",
            allow_missing=True,
        ),
        integrons=expand(
            rules.AN_assign_positions.output,
            tool="integronfinder",
            K="real",
            allow_missing=True,
        ),
        gbk_fld="data/gbk",
    output:
        gbk="results/{dset}/hotspots/{edge}/annotations.csv",
        tools="results/{dset}/hotspots/{edge}/tools.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/hotspots/extract_joint_annotations.py \
            --junction {wildcards.edge} \
            --dfinder {input.dfinder} \
            --genomad {input.genomad} \
            --isescan {input.isescan} \
            --integrons {input.integrons} \
            --out_gbk {output.gbk} \
            --out_tools {output.tools} \
            --gbk_fld {input.gbk_fld}
        """


rule HS_plot:
    input:
        gbk=rules.HS_extract_annotations.output.gbk,
        tools=rules.HS_extract_annotations.output.tools,
        pan=rules.BJ_pangraph.output.pan,
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
    output:
        fig_bl="figs/{dset}/hotspots/{edge}/blocks.pdf",
        fig_ann="figs/{dset}/hotspots/{edge}/annotations.pdf",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/hotspots/draw_junction.py \
            --junction {wildcards.edge} \
            --tree {input.tree} \
            --gbk_ann {input.gbk} \
            --tool_ann {input.tools} \
            --pan {input.pan} \
            --fig_blocks {output.fig_bl} \
            --fig_ann {output.fig_ann}
        """


def HS_all_plots(wildcards):
    files = []
    for dset, opt in itt.product(dset_names, kernel_opts):
        wc = {"dset": dset}
        # define list of edges
        edge_count_file = checkpoints.BJ_extract_joints_df.get(**wc).output["dfc"]
        edges = read_edge_count(edge_count_file)
        # add desired output files
        files += expand(rules.HS_plot.output, edge=edges, **wc)
    return files


rule HS_all:
    input:
        HS_all_plots,
