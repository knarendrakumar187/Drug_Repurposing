"""Paper 2: repurposing candidates whose drug action opposes the Paper 1 progression direction.

Input: the Task 4 main list (mechanism-annotated repurposing candidates) and their Task 2
drug-gene links. A link is
  * direction-consistent   : the drug switches the gene product DOWN and the gene RISES with
                             progression, or switches it UP and the gene FALLS;
  * direction-inconsistent : the opposite;
  * not decided            : the drug's direction is mixed or unknown.
Gene direction = Paper 1 fibrosis direction for fibrosis-associated genes, and the ordered
disease-group direction for disease-group-only genes (Task 2 column gene_direction).

A drug is kept for the metapath analysis when it has at least one direction-consistent link.
Drugs with both consistent and inconsistent links are kept and marked "mixed".
This is a mechanism-consistency filter, not evidence that a drug helps.

No Neo4j needed. Run from the repository root:

    python Paper2_Results/code/direction_consistent_drugs.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Paper2_Results" / "Direction_consistent_drugs"
OUT.mkdir(parents=True, exist_ok=True)

EDGES = ROOT / "Paper2_Results/Task2/drug_gene_edges.csv"
MAIN = ROOT / "Paper2_Results/Task4/repurposing_candidates_mechanism_annotated.csv"

edges = pd.read_csv(EDGES, low_memory=False)
main = pd.read_csv(MAIN, low_memory=False)

e = edges[edges.Drug_ID.isin(main.Drug_ID)].copy()
known = e.Action_direction.isin(["down", "up"]) & e.gene_direction.isin(["positive", "negative"])
e["direction_check"] = "not decided"
e.loc[known & (((e.Action_direction == "down") & (e.gene_direction == "positive"))
               | ((e.Action_direction == "up") & (e.gene_direction == "negative"))), "direction_check"] = "direction-consistent"
e.loc[known & (((e.Action_direction == "down") & (e.gene_direction == "negative"))
               | ((e.Action_direction == "up") & (e.gene_direction == "positive"))), "direction_check"] = "direction-inconsistent"
e.to_csv(OUT / "candidate_links_direction_check.csv", index=False)

g = e.groupby("Drug_ID")
drugs = pd.DataFrame({
    "n_consistent": g.direction_check.apply(lambda s: int((s == "direction-consistent").sum())),
    "n_inconsistent": g.direction_check.apply(lambda s: int((s == "direction-inconsistent").sum())),
    "n_not_decided": g.direction_check.apply(lambda s: int((s == "not decided").sum())),
    "consistent_targets": e[e.direction_check == "direction-consistent"].groupby("Drug_ID").Gene.agg(lambda s: ";".join(sorted(set(s)))),
    "inconsistent_targets": e[e.direction_check == "direction-inconsistent"].groupby("Drug_ID").Gene.agg(lambda s: ";".join(sorted(set(s)))),
}).reset_index()
drugs["has_known_direction"] = (drugs.n_consistent + drugs.n_inconsistent) > 0
drugs["direction_group"] = "no decided link"
drugs.loc[(drugs.n_consistent > 0) & (drugs.n_inconsistent == 0), "direction_group"] = "only consistent"
drugs.loc[(drugs.n_consistent > 0) & (drugs.n_inconsistent > 0), "direction_group"] = "mixed"
drugs.loc[(drugs.n_consistent == 0) & (drugs.n_inconsistent > 0), "direction_group"] = "only inconsistent"
cols = ["Drug_ID", "Drug", "max_phase", "is_approved", "added_drug_node", "liver_contraindicated",
        "has_liver_indication", "liver_indications", "n_progression_genes", "n_fibrosis_genes"]
drugs = drugs.merge(main[[c for c in cols if c in main.columns]], on="Drug_ID", how="left")
drugs.to_csv(OUT / "candidate_drugs_direction_summary.csv", index=False)

keep = drugs[drugs.n_consistent > 0].sort_values(["n_consistent", "n_fibrosis_genes"], ascending=False)
keep_links = e[e.Drug_ID.isin(keep.Drug_ID) & (e.direction_check == "direction-consistent")]
keep.to_csv(OUT / "direction_consistent_drugs.csv", index=False)
keep.to_excel(OUT / "direction_consistent_drugs.xlsx", index=False)
keep_links.to_csv(OUT / "direction_consistent_links.csv", index=False)

counts = drugs.direction_group.value_counts()
lc = e.direction_check.value_counts()
summary = f"""# Direction-consistent repurposing candidates

Input: the Task 4 main list ({len(main):,} mechanism-annotated repurposing candidates) and their {len(e):,} drug–gene links (Task 2, present graph).
Script: `Paper2_Results/code/direction_consistent_drugs.py`.

**Rule:** a link is direction-consistent when the drug switches the gene product **down** and the gene **rises** with progression, or switches it **up** and the gene **falls**. Gene direction = Paper 1 fibrosis direction; disease-group direction for disease-group-only genes.

## Drugs

| Group | Drugs |
|---|---|
| Main list | {len(main):,} |
| …with ≥1 link of known direction (down/up) | {int(drugs.has_known_direction.sum()):,} |
| …only direction-consistent links | {counts.get("only consistent", 0):,} |
| …both consistent and inconsistent links (mixed) | {counts.get("mixed", 0):,} |
| …only direction-inconsistent links | {counts.get("only inconsistent", 0):,} |
| …no decided link | {counts.get("no decided link", 0):,} |
| **Kept for the metapath analysis (≥1 direction-consistent link)** | **{len(keep):,}** |

## Links of the main-list drugs

| Direction check | Links |
|---|---|
""" + "\n".join(f"| {k} | {v:,} |" for k, v in lc.items()) + f"""

Kept drugs' direction-consistent links: {len(keep_links):,} on {keep_links.Gene_ID.nunique():,} genes.

## Outputs
- `direction_consistent_drugs.csv` / `.xlsx`: the kept drugs (input for the metapath runs).
- `direction_consistent_links.csv`: their direction-consistent drug–gene links.
- `candidate_drugs_direction_summary.csv`, `candidate_links_direction_check.csv`: every main-list drug and link with its check.

A direction-consistent link is a mechanism-consistency check, not evidence that a drug would help (repurposing hypotheses only).
"""
(OUT / "summary.md").write_text(summary, encoding="utf-8")
print(summary)
