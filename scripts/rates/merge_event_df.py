# %%
import pandas as pd
import numpy as np
import argparse
import utils as ut


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--terminal_df", type=str)
    parser.add_argument("--internal_df", type=str)
    parser.add_argument("--branch_df", type=str)
    parser.add_argument("--out_df", type=str)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    tdf = pd.read_csv(args.terminal_df, index_col=0)
    idf = pd.read_csv(args.internal_df, index_col=0)
    bdf = pd.read_csv(args.branch_df, index_col=0)

    tdf = ut.assign_mge_category(tdf)
    idf = ut.assign_mge_category(idf)

    terminals = bdf[bdf["terminal"]].index.to_list()

    events = []

    # process terminal events
    for j, row in tdf.iterrows():
        info = {
            "junction": j,
            "type": row["event_type"],
            "mge_cat": row["cat"],
            "branch": row["event_iso"],
            "terminal": True,
            "singleton": True,
        }
        events.append(info)

    # process internal events
    for j, row in idf.iterrows():
        for ev in ["gain", "loss", "other"]:
            N = row[ev]
            if N == 0:
                continue
            branches = row[f"{ev}_branches"].split("|")
            n = 0
            for b in branches:
                info = {
                    "junction": j,
                    "type": ev,
                    "mge_cat": row["cat"],
                    "branch": b,
                    "terminal": b in terminals,
                    "singleton": False,
                }
                events.append(info)
                n += 1
            assert n == N, f"n. branches ({n}) does not match n. events ({N})"
    df = pd.DataFrame(events)
    df.to_csv(args.out_df, index=False)

# %%
