import numpy as np
import pandas as pd
import pathlib
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--junctions_stats", type=str)
    parser.add_argument("--edge_pangenome", type=str)
    parser.add_argument("--ann_df", type=str)
    parser.add_argument("--ann_gm", type=str)
    parser.add_argument("--ann_if", type=str)
    parser.add_argument("--ann_is", type=str)
    parser.add_argument("--out_df", type=str)
    return parser.parse_args()


def add_ann_info(args, df):
    fnames = {
        "df": args.ann_df,
        "gm": args.ann_gm,
        "if": args.ann_if,
        "is": args.ann_is,
    }

    for k, fname in fnames.items():
        df2 = pd.read_csv(fname, index_col=0)
        df2 = df2["junction"].value_counts()
        df[f"{k}"] = df.index.map(df2)
        df[f"{k}"] = df[f"{k}"].fillna(0)
        df[f"{k}"] = df[f"{k}"].astype(int)

    return df


def load_df(args):
    df = pd.read_csv(args.junctions_stats, index_col=0)

    df2 = pd.read_csv(args.edge_pangenome, index_col=0)

    df = pd.merge(df, df2, on="edge", validate="one_to_one")
    df["delta_L"] = df["max_length"] - df["min_length"]

    df = add_ann_info(args, df)

    df = df[
        [
            "n_iso",
            "n_categories",
            "min_length",
            "max_length",
            "n_all_cores",
            "transitive",
            "nonempty_freq",
            "pangenome_len",
            "df",
            "gm",
            "if",
            "is",
        ]
    ]

    assert np.all(
        df["transitive"] == (df["n_categories"] == 1)
    ), "transitive means only one category"
    assert np.all(df[df["transitive"]]["nonempty_freq"] == 0), "transitive means empty"
    df = df[~df["transitive"]]

    return df


if __name__ == "__main__":

    args = parse_args()

    df = load_df(args)

    df.to_csv(args.out_df)
    