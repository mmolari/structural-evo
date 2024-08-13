import pandas as pd


checkpoint HS_dataframe:
    input:
        jdf=rules.BJ_junct_stats.output.stats,
        epg=rules.BJ_extract_pangenome_info.output.info,
        gm=expand(
            rules.AN_assign_positions.output,
            tool="genomad",
            K="real",
            allow_missing=True,
        ),
        ig=expand(
            rules.AN_assign_positions.output,
            tool="integronfinder",
            K="real",
            allow_missing=True,
        ),
        IS=expand(
            rules.AN_assign_positions.output,
            tool="ISEScan",
            K="real",
            allow_missing=True,
        ),
        df=expand(
            rules.AN_assign_positions.output,
            tool="defensefinder",
            K="real",
            allow_missing=True,
        ),
        df_hmm=expand(
            rules.AN_assign_positions.output,
            tool="dfinder_hmm",
            K="real",
            allow_missing=True,
        ),
    output:
        df="results/{dset}/hotspots/hs_df.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/hotspots/junctions_df.py \
            --junctions_stats {input.jdf} \
            --edge_pangenome {input.epg} \
            --ann_gm {input.gm} \
            --ann_if {input.ig} \
            --ann_is {input.IS} \
            --ann_df {input.df} \
            --ann_df_hmm {input.df_hmm} \
            --out_df {output.df}
        """


rule HS_extract_annotations:
    input:
        dfinder=expand(
            rules.AN_assign_positions.output,
            tool="defensefinder",
            K="real",
            allow_missing=True,
        ),
        dfinder_hmm=expand(
            rules.AN_assign_positions.output,
            tool="dfinder_hmm",
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
        gen_len=rules.PG_genome_lengths.output,
        j_pos=rules.BJ_extract_joints_pos.output.pos,
    output:
        gbk="results/{dset}/hotspots/hs/{edge}/annotations.csv",
        tools="results/{dset}/hotspots/hs/{edge}/tools.csv",
    conda:
        "../conda_env/bioinfo.yml"
    shell:
        """
        python3 scripts/hotspots/extract_joint_annotations.py \
            --junction {wildcards.edge} \
            --dfinder {input.dfinder} \
            --dfinder_hmm {input.dfinder_hmm} \
            --genomad {input.genomad} \
            --isescan {input.isescan} \
            --integrons {input.integrons} \
            --out_gbk {output.gbk} \
            --out_tools {output.tools} \
            --gbk_fld {input.gbk_fld} \
            --genome_len {input.gen_len} \
            --joint_pos {input.j_pos}
        """


rule HS_plot:
    input:
        gbk=rules.HS_extract_annotations.output.gbk,
        tools=rules.HS_extract_annotations.output.tools,
        pan=rules.BJ_pangraph.output.pan,
        tree=rules.PG_filtered_coregenome_tree.output.nwk,
    output:
        fig_bl="results/{dset}/hotspots/hs/{edge}/blocks.pdf",
        fig_ann="results/{dset}/hotspots/hs/{edge}/annotations.pdf",
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
        # add desired output files
        files += expand(rules.HS_plot.output, edge=edges, **wc)
    return files


rule HS_all:
    input:
        HS_all_plots,