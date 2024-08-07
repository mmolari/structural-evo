# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# dset = "ecoli_ST69"
dset = "paeruginosa_ST235"

fname = f"../../results/{dset}/backbone_joints/junctions_stats.csv"
df = pd.read_csv(fname)
fname = f"../../results/{dset}/backbone_joints/edge_pangenome.csv"
df2 = pd.read_csv(fname)
df = df.merge(df2, on="edge", validate="one_to_one")
# %%
df["n. paths > 20"] = df["n_categories"] > 20
sns.scatterplot(
    data=df,
    y="n_categories",
    x="pangenome_len",
    hue="n. paths > 20",
    alpha=0.5,
)
plt.xlabel("Pangenome length")
plt.ylabel("Number of categories")
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()
plt.savefig("hotspot_selection.png")
plt.show()

# %%
sdf = df[df["n. paths > 20"]]
sdf

# %%
x = (
    sdf[["edge", "n_categories", "pangenome_len"]]
    .sort_values("n_categories", ascending=False)
    .to_markdown()
)
print(x)
# %%
import json

fname = "../../results/ST131_ABC/backbone_joints/joints_pos.json"
with open(fname) as f:
    pos = json.load(f)

new_pos = {}
for k in sdf["edge"].values:
    new_pos[k] = {
        iso: (p[0], p[3], "fwd" if p[4] else "rev") for iso, p in pos[k].items()
    }

with open("hotspot_positions.json", "w") as f:
    json.dump(new_pos, f)

# %%
