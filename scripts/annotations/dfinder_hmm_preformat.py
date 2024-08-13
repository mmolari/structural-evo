import pandas as pd
import numpy as np
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_hmm", required=True)
    parser.add_argument("--output_df", required=True)
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    hdf = pd.read_csv(args.input_hmm, sep="\t")

    # for each entry, combine hit_id and gene_name in one string
    hdf["hmm_id"] = hdf["hit_id"] + "|" + hdf["gene_name"]

    cols = {"hmm_id": "id", "replicon": "iso", "gene_beg": "beg", "gene_end": "end"}
    df = hdf.reset_index()[list(cols.keys())].rename(columns=cols)
    df["type"] = "dfinder_hmm"
    df.to_csv(args.output_df, index=False)
