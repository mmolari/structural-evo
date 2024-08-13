#  %%
# draw paths with gene annotations

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
import pypangraph as pp
import argparse
import pathlib
from collections import defaultdict

left_clr = "#2db7a7"
right_clr = "#cf4b36"


def parse_args():
    parser = argparse.ArgumentParser(description="Draw linear junctions")
    parser.add_argument("--junction", type=str, help="junction name")
    parser.add_argument("--tree", type=str, help="coregenome tree")
    parser.add_argument("--gbk_ann", type=str, help="genbank annotations")
    parser.add_argument("--tool_ann", type=str, help="tool annotations")
    parser.add_argument("--pan", type=str, help="junction graph")
    parser.add_argument("--fig_blocks", type=str, help="output figure blocks")
    parser.add_argument("--fig_ann", type=str, help="output figure annotations")
    return parser.parse_args()


def despine(ax, y=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if not y:
        ax.spines["left"].set_visible(False)
        ax.set_yticks([])
        ax.set_ylabel("")

def load_data(args):
    # load tree
    tree = Phylo.read(args.tree, "newick")

    #  load gene df
    gdf = pd.read_csv(args.gbk_ann)

    # load graph
    pan = pp.Pangraph.load_json(args.pan)

    # load tool annotations
    tdf = pd.read_csv(args.tool_ann)

    return tree, gdf, pan, tdf


def color_generator():
    cm = plt.get_cmap("tab20")
    for i in range(20):
        yield cm(i)
    cm = plt.get_cmap("tab20b")
    for i in range(20):
        yield cm(i)
    cm = plt.get_cmap("tab20c")
    for i in range(20):
        yield cm(i)
    cm = plt.get_cmap("hsv")
    while True:
        yield cm(np.random.rand())


def draw_gene(ax, b, e, y, strand, w=1, c="black", arr_l=100):
    if strand:
        a = max(b, e - arr_l)
        X = [b, a, e, a, b, b]
        Y = [w / 2, w / 2, 0, -w / 2, -w / 2, w / 2]
    else:
        a = min(b + arr_l, e)
        X = [e, a, b, a, e, e]
        Y = [-w / 2, -w / 2, 0, w / 2, w / 2, -w / 2]
    Y = np.array(Y) + y
    ax.plot(X, Y, color=c, lw=1)


def draw_genes(ax, gdf, strain_y):

    colors = {"CDS": "black", "tRNA": "orange", "ncRNA": "green", "tmRNA": "orange", "rRNA": "cyan"}

    for i, row in gdf.iterrows():
        b, e, strand = row.start, row.end, row.strand
        kind, iso, gene = row.kind, row.iso, row.gene
        if kind == "gene":
            continue
        y = strain_y[iso]
        c = colors[kind]
        draw_gene(ax, b, e, y, strand, c=c, w=0.7)
        # add gene name if not NaN
        if not pd.isna(gene):
            gx = (b + e) / 2
            ax.text(gx, y, gene, fontsize=3, va="center", ha="left")


def mark_tool_ann(ax, tdf, strain_y):
    colors = {
        "ISEScan": "C0",
        "defensefinder": "C2",
        "dfinder_hmm": "palegreen",
        "integronfinder": "C3",
        "genomad": "C4",
    }

    # draw annotations as barh in the background
    for i, row in tdf.iterrows():
        y = strain_y[row.iso]

        ax.barh(
            y,
            row.end - row.start,
            left=row.start,
            color=colors[row.kind],
            alpha=0.5,
            height=0.5,
            zorder=-1,
        )

        if row.kind == "defensefinder":
            x = 0.5*(row.end + row.start)
            ax.text(x, y, row.idx, fontsize=4, va="center", ha="center", color="darkgreen")
        


def mark_blocks(ax, strain_y, pan, block_colors):

    ax = axs[1]
    for iso, y in strain_y.items():
        if iso not in pan.strains():
            continue
        path = pan.paths[iso]
        Bs = path.block_ids
        Ss = path.block_strands
        Ps = path.block_positions
        for n_block in range(len(Bs)):
            bid = Bs[n_block]
            strand = Ss[n_block]
            x0, x1 = Ps[n_block], Ps[n_block + 1]
            color = block_colors[bid]
            edgecolor = "none"
            if not strand:
                edgecolor = "k"
            ax.barh(
                y,
                x1 - x0,
                left=x0,
                height=0.8,
                color=color,
                edgecolor=edgecolor,
                linewidth=0.5,
                zorder=-1,
            )

def create_fig(n, L):
    X = max(24, L/3000)
    Y = n / 8
    fig, axs = plt.subplots(
        1, 2, figsize=(X, Y), sharey=True, gridspec_kw={"width_ratios": [1, 10]}
    )
    return fig, axs

def draw_tree(tree, ax):
    Phylo.draw(
        tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: ""
    )
    for n, l in enumerate(tree.get_terminals()):
        y = n+1
        d = tree.distance(tree.root, l)
        ax.text(d, y, l.name, fontsize=3, va="center", ha="left") 

if __name__ == "__main__":
    args = parse_args()

    tree, gdf, pan, tdf = load_data(args)
    strain_y = {l.name: i + 1 for i, l in enumerate(tree.get_terminals())}

    gen = color_generator()
    block_colors = defaultdict(lambda: next(gen))

    # assign colors to first/last blocks
    strains = pan.strains()
    Bs = pan.paths[strains[0]].block_ids
    block_colors[Bs[0]] = left_clr
    block_colors[Bs[-1]] = right_clr

    L = max(max(p.block_positions) for p in pan.paths)

    fig, axs = create_fig(len(strain_y), L)

    ax = axs[0]
    draw_tree(tree, ax)
    despine(ax, y=False)
    ax.set_xlabel("branch length")
    ax.set_ylabel("isolates")

    ax = axs[1]
    draw_genes(ax, gdf, strain_y)
    mark_tool_ann(ax, tdf, strain_y)
    ax.set_xlabel("junction position (bp)")
    ax.set_xlim(0, max([pan.paths[iso].block_positions[-1] for iso in pan.strains()]))

    despine(ax)
    plt.tight_layout()
    plt.savefig(args.fig_ann)
    plt.close()

    fig, axs = create_fig(len(strain_y), L)

    ax = axs[0]
    draw_tree(tree, ax)
    despine(ax, y=False)

    ax = axs[1]
    draw_genes(ax, gdf, strain_y)
    mark_blocks(ax, strain_y, pan, block_colors)
    ax.set_xlabel("junction position (bp)")
    ax.set_xlim(0, max([pan.paths[iso].block_positions[-1] for iso in pan.strains()]))

    despine(ax)
    plt.tight_layout()
    plt.savefig(args.fig_blocks)
    plt.close()