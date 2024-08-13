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
    for dset in dset_names:
        wc = {"dset": dset}
        # define list of edges
        df_file = checkpoints.HS_dataframe.get(**wc).output["df"]
        df = pd.read_csv(df_file, index_col=0)
        # select edges
        mask = df["df"] > 0
        edges = df[mask].index.to_list()
        edges += df.sort_values("n_categories", ascending=False).index[:5].to_list()
        edges += df.sort_values("pangenome_len", ascending=False).index[:3].to_list()
        edges += df.sort_values("gm", ascending=False).index[:3].to_list()
        edges += df.sort_values("is", ascending=False).index[:3].to_list()
        edges = list(set(edges))
        # add desired output files
        files += expand(rules.HS_plot.output, edge=edges, **wc)
    return files


rule HS_all:
    input:
        HS_all_plots,
