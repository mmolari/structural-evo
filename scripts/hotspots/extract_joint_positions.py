# %%
from Bio import SeqIO
import pandas as pd
from dataclasses import dataclass
import pathlib
import json
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--jpos", type=str)
    parser.add_argument("--edge", type=str)
    parser.add_argument("--out", type=str)
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()
    
    with open(args.jpos, "r") as f:
        jpos = json.load(f)[args.edge]

    df = []
    for iso, p in jpos.items():
        cb, ab, ae, ce, strand = p
        df.append({
            "iso": iso,
            "start": cb,
            "end": ce,
            "strand": strand
        })
    df = pd.DataFrame(df)
    df.to_csv(args.out, index=False)