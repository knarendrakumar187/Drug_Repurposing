"""Paper 2 Task 1: master progression-gene table.

Builds Paper2_Results/Task1/Paper2_progression_genes.csv from the corrected
6,222-gene union (Paper1_HighConfidence/00_all_progression_genes.csv) and
annotates it with Paper 1 trajectory, transition, tier and shortlist evidence.
No gene is dropped. Run from the repository root:

    python Paper2_Results/code/task1_master_gene_table.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Paper2_Results" / "Task1"
OUT.mkdir(parents=True, exist_ok=True)

SRC_ALL = ROOT / "Paper1_HighConfidence/00_all_progression_genes.csv"
SRC_CLUSTERS = ROOT / "Paper1_Results/Task2/task2_gene_clusters.csv"
SRC_MASTER = ROOT / "Paper1_HighConfidence/01_master_table.csv"
SRC_SHORTLIST = ROOT / "Paper1_HighConfidence/11_FINAL_shortlist_25_genes.csv"
SRC_P1_SUMMARY = ROOT / "Paper1_Results/Task1/task1_summary.csv"

EXPECTED_GENES = 6222
DIRECTION = {"increasing": "positive", "decreasing": "negative"}
CLASS = {"fibrosis_only": "fibrosis_only", "inflammation_only": "inflammation_only", "both": "shared"}


def paper1_group(row):
    """Paper 1 Task 1 group definition (class x direction)."""
    cls = row.progression_class
    if cls == "fibrosis_only":
        return f"fibrosis_only_{row.fibrosis_direction}"
    if cls == "inflammation_only":
        return f"inflammation_only_{row.inflammation_direction}"
    if row.fibrosis_direction == row.inflammation_direction:
        return f"shared_{row.fibrosis_direction}"
    return f"shared_discordant_fib_{row.fibrosis_direction}"


genes = pd.read_csv(SRC_ALL)
assert len(genes) == EXPECTED_GENES and genes.ensembl_id.is_unique, "unexpected source gene set"

df = pd.DataFrame({
    "Gene": genes.symbol,
    "ensembl_id": genes.ensembl_id,
    "fibrosis_rho": genes.fibrosis_rho,
    "fibrosis_q": genes.fibrosis_q,
    "inflammation_rho": genes.inflammation_rho,
    "inflammation_q": genes.inflammation_q,
    "progression_class": genes.ladder.map(CLASS),
    "fibrosis_direction": genes.fibrosis_direction.map(DIRECTION),
    "inflammation_direction": genes.inflammation_direction.map(DIRECTION),
    "F0": genes.F0, "F1": genes.F1, "F2": genes.F2, "F3": genes.F3, "F4": genes.F4,
})
assert df.progression_class.notna().all()
df["paper1_group"] = df.apply(paper1_group, axis=1)

# Trajectory clusters (fibrosis genes only)
clusters = pd.read_csv(SRC_CLUSTERS, encoding="utf-8-sig")
clusters = clusters[["ensembl_id", "cluster", "cluster_name", "membership_confidence"]].copy()
clusters["cluster"] = "C" + clusters.cluster.astype(int).astype(str)
clusters = clusters.rename(columns={"cluster": "trajectory_cluster"})
df = df.merge(clusters, on="ensembl_id", how="left")

# Stage transitions and tier (fibrosis genes only)
master = pd.read_csv(SRC_MASTER)[["ensembl_id", "largest_transition", "transition_beta", "transition_q", "tier"]]
master = master.rename(columns={"tier": "Tier"})
df = df.merge(master, on="ensembl_id", how="left")

# High-confidence shortlist
shortlist = set(pd.read_csv(SRC_SHORTLIST).Ensembl_ID)
df["in_shortlist"] = df.ensembl_id.isin(shortlist)

assert len(df) == EXPECTED_GENES and df.ensembl_id.is_unique, "merge changed the gene set"
df.to_csv(OUT / "Paper2_progression_genes.csv", index=False)

# ---------------- summary ----------------
fib = df[df.progression_class != "inflammation_only"]
counts_class = df.progression_class.value_counts().reindex(["fibrosis_only", "inflammation_only", "shared"])
counts_group = df.paper1_group.value_counts().reindex([
    "fibrosis_only_positive", "fibrosis_only_negative",
    "shared_positive", "shared_negative", "shared_discordant_fib_negative", "shared_discordant_fib_positive",
    "inflammation_only_positive", "inflammation_only_negative",
]).fillna(0).astype(int)
counts_shared_fibdir = df[df.progression_class == "shared"].fibrosis_direction.value_counts()

# Check against Paper 1 Task 1 summary
p1 = pd.read_csv(SRC_P1_SUMMARY, encoding="utf-8-sig")
p1_key = {
    ("fibrosis_only", "all"): counts_class["fibrosis_only"],
    ("fibrosis_only", "positive"): counts_group["fibrosis_only_positive"],
    ("fibrosis_only", "negative"): counts_group["fibrosis_only_negative"],
    ("inflammation_only", "all"): counts_class["inflammation_only"],
    ("inflammation_only", "positive"): counts_group["inflammation_only_positive"],
    ("inflammation_only", "negative"): counts_group["inflammation_only_negative"],
    ("shared", "all"): counts_class["shared"],
    ("shared", "positive (both up)"): counts_group["shared_positive"],
    ("shared", "negative (both down)"): counts_group["shared_negative"],
    ("shared", "discordant (fib up, inf down)"): counts_group["shared_discordant_fib_positive"],
    ("shared", "discordant (fib down, inf up)"): counts_group["shared_discordant_fib_negative"],
}
check = []
for _, r in p1.iterrows():
    ours = int(p1_key[(r.group, r.direction)])
    check.append((r.group, r.direction, int(r.n_genes), ours, ours == int(r.n_genes)))
check = pd.DataFrame(check, columns=["group", "direction", "paper1_n", "paper2_n", "match"])
check.to_csv(OUT / "count_check_vs_paper1.csv", index=False)
all_match = bool(check.match.all())

annot = {
    "trajectory_cluster": df.trajectory_cluster.notna().sum(),
    "largest_transition": df.largest_transition.notna().sum(),
    "Tier": df.Tier.notna().sum(),
    "in_shortlist": int(df.in_shortlist.sum()),
}
tiers = df.Tier.value_counts()
clust = df.trajectory_cluster.value_counts().sort_index()


def md_table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(lines)


summary = f"""# Paper 2 — Task 1: master progression-gene table

This table uses the corrected Paper 1 union of **{len(df):,} genes** (`Paper1_HighConfidence/00_all_progression_genes.csv`, final combined version 10 Oct 2026). No gene was removed.

**Naming note:** columns and classes called "inflammation" hold Paper 1's **ordered disease-group** axis (control/NAFL/NASH labels; not a direct inflammation measurement). The plan's column names are kept so later tasks can read the file; the manuscript should say "disease-group".

Output: `Paper2_progression_genes.csv` ({len(df):,} rows). Script: `Paper2_Results/code/task1_master_gene_table.py`.


## Genes per progression class

{md_table(["Progression class", "Genes", "% of 6,222"], [(k, f"{v:,}", f"{100*v/len(df):.1f}") for k, v in counts_class.items()])}

Fibrosis-associated genes (fibrosis_only + shared): **{len(fib):,}**.

## Genes per class x direction (Paper 1 group definitions)

Direction = fibrosis direction for fibrosis_only, inflammation direction for inflammation_only, and agreement of both for shared (as in Paper 1 Task 1).

{md_table(["Group", "Genes"], [(k, f"{v:,}") for k, v in counts_group.items()])}

Shared genes split by **fibrosis** direction only: positive {counts_shared_fibdir.get('positive', 0):,}, negative {counts_shared_fibdir.get('negative', 0):,}.

## Paper 1 annotations available

{md_table(["Annotation", "Genes annotated", "Source"], [
    ("Trajectory cluster (C1–C5) + membership confidence", f"{annot['trajectory_cluster']:,}", "Paper1_Results/Task2/task2_gene_clusters.csv"),
    ("Largest stage transition + beta + q", f"{annot['largest_transition']:,}", "Paper1_HighConfidence/01_master_table.csv"),
    ("Tier", f"{annot['Tier']:,}", "Paper1_HighConfidence/01_master_table.csv"),
    ("In 23-gene shortlist", f"{annot['in_shortlist']:,}", "Paper1_HighConfidence/11_FINAL_shortlist_25_genes.csv"),
])}

These annotations exist only for the {len(fib):,} fibrosis-associated genes; inflammation-only genes are blank by design (Paper 1 clustered and scored fibrosis genes only).

Trajectory clusters: {", ".join(f"{k} {v:,}" for k, v in clust.items())}.
Tiers: {", ".join(f"{k} {v:,}" for k, v in tiers.items())}.

## Check against Paper 1 Task 1

{md_table(["Group", "Direction", "Paper 1", "Paper 2", "Match"], check.itertuples(index=False))}

**All counts match Paper 1: {"YES" if all_match else "NO — investigate"}.**
"""
(OUT / "summary.md").write_text(summary, encoding="utf-8")
print(summary)
