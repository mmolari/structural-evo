import pathlib as pthl
import matplotlib.pyplot as plt
from Bio import Phylo
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="plot large core-genome tree")
    parser.add_argument("--coregenome_tree", type=str)
    parser.add_argument("--fig", type=str)
    return parser.parse_args()

def large_tree(tree, svname):
    leaves = [l.name for l in tree.get_terminals()]
    N = len(leaves)
    fig, ax = plt.subplots(1, 1, figsize=(10, N * 0.14 + 1))
    Phylo.draw(
        tree,
        label_func=lambda x: x.name if x.name in leaves else "",
        do_show=False,
        axes=ax,
    )
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tight_layout()
    plt.savefig(svname)
    plt.close(fig)


if __name__ == "__main__":
    # load input data
    args = parse_args()
    tree = Phylo.read(args.coregenome_tree, "newick")
    tree.root_at_midpoint()
    tree.ladderize()
    svname = pthl.Path(args.fig)
    large_tree(tree, svname)
