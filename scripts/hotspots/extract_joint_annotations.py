# %%
from Bio import SeqIO
import pandas as pd
from dataclasses import dataclass
import pathlib
import json
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--junction", type=str)
    parser.add_argument("--dfinder", type=str)
    parser.add_argument("--dfinder_hmm", type=str)
    parser.add_argument("--genomad", type=str)
    parser.add_argument("--isescan", type=str)
    parser.add_argument("--integrons", type=str)
    parser.add_argument("--out_gbk", type=str)
    parser.add_argument("--out_tools", type=str)
    parser.add_argument("--gbk_fld", type=str)
    parser.add_argument("--genome_len", type=str)
    parser.add_argument("--joint_pos", type=str)
    return parser.parse_args()


@dataclass
class Ann:
    b: int
    e: int
    s: bool
    iso: str
    idx: str
    kind: str
    product: str = None
    gene: str = None
    locus_tag: str = None

    def to_dict(self):
        return {
            "start": self.b,
            "end": self.e,
            "strand": self.s,
            "iso": self.iso,
            "idx": self.idx,
            "kind": self.kind,
            "product": self.product,
            "gene": self.gene,
            "locus_tag": self.locus_tag,
        }


def load_tools_ann(args):
    dfs = []
    for tool, fname in [
        ("defensefinder", args.dfinder),
        ("dfinder_hmm", args.dfinder_hmm),
        ("genomad", args.genomad),
        ("ISEScan", args.isescan),
        ("integronfinder", args.integrons),
    ]:
        df = pd.read_csv(fname)
        df["tool"] = tool
        dfs.append(df)
    dfs = pd.concat(dfs)
    mask = dfs["junction"] == args.junction
    dfs = dfs[mask]
    dfs.reset_index(inplace=True, drop=True)
    return dfs


def dist_fwd(s, e, L):
    return (e - s) % L


def to_junction_coords(Jb, Je, Js, Ab, Ae, As, L):
    assert (Ab < Ae) or (Ab - Ae > 2e6)
    assert (Jb < Je) or (Jb - Je > 2e6)

    if Js:
        B = dist_fwd(Jb, Ab, L)
        E = dist_fwd(Jb, Ae, L)
    else:
        B = dist_fwd(Ae, Je, L)
        E = dist_fwd(Ab, Je, L)
    if As is None:
        S = None
    else:
        S = not (As ^ Js)

    if B > 1e6:
        B -= L
    if E > 1e6:
        E -= L
    assert -B < 1e6, f"{B=} {Ab=} {Ae=} {Je=} {L=}"
    assert -E < 1e6, f"{E=} {Ab=} {Ae=} {Je=} {L=}"

    return B, E, S


def extract_tool_annotations(dfs, Ls):

    Adf = []
    for idx, row in dfs.iterrows():
        Jb, Je, Js = row["jcb"], row["jce"], row["j_strand"]
        Ab, Ae, As = row["ib"], row["ie"], None
        iso = row["iso"]
        L = Ls[iso]
        B, E, S = to_junction_coords(Jb, Je, Js, Ab, Ae, None, L)
        idx = row["id"]
        kind = row["tool"]
        ann = Ann(B, E, S, iso, idx, kind)
        Adf.append(ann.to_dict())
    Adf = pd.DataFrame(Adf)
    return Adf


def parse_genbank(fname):
    with open(fname) as f:
        records = list(SeqIO.parse(f, "genbank"))
    assert len(records) == 1
    record = records[0]
    return record


def parse_genbank(fname):
    with open(fname) as f:
        records = list(SeqIO.parse(f, "genbank"))
    assert len(records) == 1
    record = records[0]
    return record


def isin(Js, Je, pos):
    if Js < Je:
        return Js <= pos < Je
    else:
        return (Js <= pos) or (pos < Je)


def extract_features(record, iso, Jb, Je, Js, L):
    A = []
    for feat in record.features:
        Ab = feat.location.start
        Ae = feat.location.end

        if isin(Jb, Je, Ab) or isin(Jb, Je, Ae):
            kind = feat.type
            if kind == "source" or kind == "gene" or kind == "misc_feature":
                continue
            As = True if feat.location.strand == 1 else False
            B, E, S = to_junction_coords(Jb, Je, Js, Ab, Ae, As, L)
            if B == E:
                continue
            assert (
                B < E
            ), f"{B=} {E=} {Ab=} {Ae=} {Jb=} {Je=} {L=} {iso=} {kind=} {feat.qualifiers}"
            labels = ["product", "gene", "locus_tag"]
            if not "locus_tag" in feat.qualifiers:
                print(f"locus_tag not in {feat.qualifiers}")
                print(feat)
                continue
            lt = feat.qualifiers["locus_tag"][0]
            kwargs = {l: feat.qualifiers[l][0] for l in labels if l in feat.qualifiers}
            ann = Ann(B, E, S, iso, lt, kind, **kwargs)
            A.append(ann.to_dict())
    A = pd.DataFrame(A)
    return A


def print_loading_bar(k, N):
    L = 50
    f = k / N
    done = int(f * L)
    missing = L - done
    print(f"|{'#'*done}{'.'*(missing)}|   {f*100:.0f}%", end="\r")


def parse_gbk_annotations(Ls, jpos, gbk_fld):
    # parse genbank annotations
    As = []
    strains = list(Ls.keys())
    for k, iso in enumerate(strains):

        if not iso in jpos:
            continue

        # loading bar
        print_loading_bar(k, len(strains))

        fname = gbk_fld / f"{iso}.gbk"
        record = parse_genbank(fname)
        Jb, _, _, Je, Js = jpos[iso]
        L = Ls[iso]
        # assert record.id == iso, f"{record.id=}, {iso=}"
        A = extract_features(record, iso, Jb, Je, Js, L)
        As.append(A)
    As = pd.concat(As)
    As.reset_index(inplace=True, drop=True)
    return As


def set_cols_if_empty(df):
    if len(df) == 0:
        df = pd.DataFrame(
            columns=[
                "start",
                "end",
                "strand",
                "iso",
                "idx",
                "kind",
                "product",
                "gene",
                "locus_tag",
            ]
        )
    return df


if __name__ == "__main__":

    args = parse_args()
    hs = args.junction
    gbk_fld = pathlib.Path(args.gbk_fld)

    # load data
    dfs = load_tools_ann(args)
    Ls = pd.read_csv(args.genome_len, index_col=0)["length"].to_dict()
    with open(args.joint_pos, "r") as f:
        jpos = json.load(f)[hs]

    # parse tool annotations
    ann_df = extract_tool_annotations(dfs, Ls)
    ann_df = set_cols_if_empty(ann_df)
    ann_df.to_csv(args.out_tools, index=False)

    # parse genbank annotations
    ann_gbk_df = parse_gbk_annotations(Ls, jpos, gbk_fld)
    ann_gbk_df = set_cols_if_empty(ann_gbk_df)
    ann_gbk_df.to_csv(args.out_gbk, index=False)

# %%