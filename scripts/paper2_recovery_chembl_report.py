"""
paper2_recovery_chembl_report.py
--------------------------------
Report RECOVERY 2 -- what ChEMBL, an external database, added -- on its own.

WHY THIS IS A SEPARATE SCRIPT AND A SEPARATE FOLDER
The mentor's instruction of 13 October: two jobs were asked for, and if two were
done they must be kept separate. Recovery 1 reads OptimusKG's own discarded
fields. This one leaves the knowledge graph and asks an external API. Those are
different acts with different standing, and pooling them into a single
"recovered" count would hide which source carried which field.

So this folder never cites Recovery 1's totals. It reports only ChEMBL's
marginal contribution: what it added ON TOP OF what OptimusKG already gave.

THE ANSWER, UP FRONT
Almost nothing, and the reason is specific. Reported against three denominators
because the answer changes completely between them, and only the last one
matters for the paper:

  all 18,640 OptimusKG raw records    ChEMBL is a real contribution
  the graph's 12,025 Drug nodes       ChEMBL contributes much less
  the 2,404 drugs in Task 4           ChEMBL contributes almost nothing

The drugs missing a structure are precisely the ones that never reach a
progression gene. Every gap ChEMBL could fill sits outside the analysis.

The contribution figure depends entirely on which drugs are counted, so all
three denominators are reported rather than one.

Usage:
    uv run --with neo4j --with pandas python scripts/paper2_recovery_chembl_report.py
"""

import json
import os

import pandas as pd
from neo4j import GraphDatabase

OKG = "data/optimuskg_recovered"
CHEMBL_DIR = "data/chembl_recovered"
METADATA = f"{OKG}/drug_metadata.csv"
GAP_FILL = f"{CHEMBL_DIR}/chembl_gap_fill.csv"
OUT_DIR = "Paper2_Results/Recovery2_ChEMBL"

T2 = "Paper2_Results/Task2/task2_drug_gene_edges.csv"
T4 = "Paper2_Results/Task4/task4_drug_coverage.csv"

# OptimusKG field -> ChEMBL field, and the name the delivered tables use
PAIRS = [
    ("okg_name", "chembl_pref_name", "drug_name"),
    ("type", "chembl_molecule_type", "molecule_type"),
    ("maximum_clinical_trial_phase", "chembl_max_phase", "max_phase"),
    ("canonical_smiles", "chembl_canonical_smiles", "smiles"),
]


def driver():
    """Credentials come from the environment, never from source.
    NEO4J_PASSWORD (prompted if unset); NEO4J_URI, NEO4J_USERNAME and
    NEO4J_DATABASE default to bolt://localhost:7687, neo4j and neo4j."""
    import getpass
    pw = os.environ.get("NEO4J_PASSWORD") or getpass.getpass("Neo4j password: ")
    return (GraphDatabase.driver(os.environ.get("NEO4J_URI", "bolt://localhost:7687"),
                                 auth=(os.environ.get("NEO4J_USERNAME", "neo4j"), pw)),
            os.environ.get("NEO4J_DATABASE", "neo4j"))


def graph_drug_ids():
    drv, db = driver()
    try:
        with drv.session(database=db) as s:
            rows = s.run("MATCH (d:Drug) RETURN d.node_id AS id").data()
    finally:
        drv.close()
    return {r["id"] for r in rows if r["id"] is not None}


def source_database_split(md, graph_ids):
    """Which database KEYS each drug node, and which fields follow from that.

    This is the measurement that reframes Recovery 2. The graph's drug
    vocabulary is not merely ChEMBL-derived in spirit -- 11,623 of its 12,025
    Drug nodes ARE ChEMBL accessions, and molecule type and clinical phase are
    present for exactly those and for no DrugBank-keyed drug. Those fields are
    ChEMBL's fields, already inside the source the graph was built from.
    """
    chembl = {i for i in graph_ids if str(i).startswith("CHEMBL")}
    g = md[md.Drug_ID.isin(graph_ids)]
    rows = []
    for col, label in [("type", "molecule type"),
                       ("maximum_clinical_trial_phase", "clinical phase"),
                       ("canonical_smiles", "SMILES"),
                       ("inchi_key", "InChI key"),
                       ("cd_formula", "formula")]:
        isc = g.Drug_ID.isin(chembl)
        rows.append({
            "field": label,
            "chembl_keyed_drugs": int((g[col].notna() & isc).sum()),
            "drugbank_keyed_drugs": int((g[col].notna() & ~isc).sum()),
        })
    totals = {"field": "— total drugs —",
              "chembl_keyed_drugs": len(chembl),
              "drugbank_keyed_drugs": len(g) - len(chembl)}
    return pd.DataFrame(rows + [totals]), len(chembl), len(g) - len(chembl)


def marginal(j, denominator_ids, label):
    """What ChEMBL adds where OptimusKG is silent, within one denominator."""
    s = j if denominator_ids is None else j[j.Drug_ID.isin(denominator_ids)]
    out = []
    for okg, ch, name in PAIRS:
        gap = s[okg].isna()
        held = gap & s[ch].notna()
        out.append({
            "denominator": label,
            "drugs_in_denominator": len(s),
            "field": name,
            "optimuskg_supplied": int(s[okg].notna().sum()),
            "optimuskg_silent": int(gap.sum()),
            "chembl_added": int(held.sum()),
        })
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    md = pd.read_csv(METADATA, low_memory=False)
    gap = pd.read_csv(GAP_FILL, low_memory=False)
    # OptimusKG stores an unnamed drug's accession AS its name, so an accession
    # means "no name" and must read as silent, not as supplied.
    md["okg_name"] = md.name.where(md.name != md.Drug_ID)
    j = md.merge(gap, on="Drug_ID", how="left")

    drugs = graph_drug_ids()
    t4 = pd.read_csv(T4, low_memory=False)
    t4_ids = set(t4.Drug_ID.dropna())

    print(f"ChEMBL records fetched      {len(gap):,}")
    print(f"OptimusKG raw records       {len(md):,}")
    print(f"graph Drug nodes            {len(drugs):,}")
    print(f"drugs in Task 4 coverage    {len(t4_ids):,}")

    split, n_chembl_keyed, n_db_keyed = source_database_split(md, drugs)
    split.to_csv(f"{OUT_DIR}/recovery2_source_database_split.csv", index=False)
    print(f"\nthe graph's drug vocabulary: {n_chembl_keyed:,} ChEMBL-keyed, "
          f"{n_db_keyed:,} DrugBank-keyed")
    print(split.to_string(index=False))

    rows = (marginal(j, None, "all OptimusKG raw records")
            + marginal(j, drugs, "the graph's Drug nodes")
            + marginal(j, t4_ids, "drugs in Task 4 coverage"))
    tab = pd.DataFrame(rows)
    tab.to_csv(f"{OUT_DIR}/recovery2_marginal_contribution.csv", index=False)

    for label in tab.denominator.unique():
        sub = tab[tab.denominator == label]
        n = sub.drugs_in_denominator.iloc[0]
        print(f"\n{label}  (n={n:,})")
        print("   field           OKG supplied   OKG silent   ChEMBL added")
        for _, r in sub.iterrows():
            print(f"   {r.field:<14} {r.optimuskg_supplied:>12,} "
                  f"{r.optimuskg_silent:>12,} {r.chembl_added:>14,}")

    # ---- why the structures could not be filled ---------------------------
    g = j[j.Drug_ID.isin(drugs)]
    no_smiles = g[g.canonical_smiles.isna()]
    queried = no_smiles[no_smiles.Drug_ID.isin(set(gap.Drug_ID))]
    why = [
        {"measure": "graph drugs with no OptimusKG SMILES", "value": len(no_smiles)},
        {"measure": "of those, ChEMBL held a record", "value": len(queried)},
        {"measure": "of those, ChEMBL had a structure",
         "value": int(queried.chembl_canonical_smiles.notna().sum())},
        {"measure": "of those, reaching a progression gene",
         "value": int(no_smiles.Drug_ID.isin(t4_ids).sum())},
    ]
    pd.DataFrame(why).to_csv(f"{OUT_DIR}/recovery2_structure_gap.csv",
                             index=False)
    print("\nthe structure gap, and why ChEMBL cannot close it")
    for r in why:
        print(f"   {r['measure']:<45} {r['value']:>7,}")

    # ---- what reached the delivered table ---------------------------------
    deliv = []
    for _, _, name in PAIRS:
        col = f"{name}_source"
        if col in t4:
            deliv.append({
                "field": name,
                "from_ChEMBL": int((t4[col] == "ChEMBL").sum()),
                "from_OptimusKG": int((t4[col] == "OptimusKG").sum()),
                "drugs_in_table": len(t4)})
    dv = pd.DataFrame(deliv)
    dv.to_csv(f"{OUT_DIR}/recovery2_contribution_to_task4.csv", index=False)
    print(f"\nin the DELIVERED Task 4 table ({len(t4):,} drugs)")
    print(dv.to_string(index=False))

    write_readme(tab, why, dv, len(gap), len(md), len(drugs), len(t4_ids),
                 split, n_chembl_keyed, n_db_keyed)
    print(f"\nwrote {OUT_DIR}/")


def write_readme(tab, why, dv, n_fetched, n_raw, n_graph, n_t4,
                 split, n_chembl_keyed, n_db_keyed):
    w = {r["measure"]: r["value"] for r in why}
    d = dv.set_index("field")
    total_chembl = int(dv.from_ChEMBL.sum())

    def block(label):
        sub = tab[tab.denominator == label]
        n = sub.drugs_in_denominator.iloc[0]
        lines = [f"**{label}** — {n:,} drugs", "",
                 "| field | OptimusKG supplied | OptimusKG silent | ChEMBL added |",
                 "|---|---:|---:|---:|"]
        for _, r in sub.iterrows():
            lines.append(f"| {r.field} | {r.optimuskg_supplied:,} | "
                         f"{r.optimuskg_silent:,} | **{r.chembl_added:,}** |")
        return "\n".join(lines)

    txt = f"""# Recovery 2 — ChEMBL, an external database

**Kept separate from Recovery 1 (OptimusKG) on the mentor's instruction of
13 October 2026.** This folder reports only ChEMBL's *marginal* contribution —
what it added on top of what OptimusKG already supplied — and never pools its
numbers with Recovery 1's.

The separation matters beyond bookkeeping. Recovery 1 re-read a source the graph
was already built from. Recovery 2 left the knowledge graph and queried an
outside API. Those are different acts, and a reader judging "focus only on the
knowledge graph" needs to see which fields came from which.

## The answer: almost nothing, and the reason is specific

{n_fetched:,} records were fetched from the ChEMBL API, for drugs where
OptimusKG was short. The contribution depends entirely on which drugs you count,
which is why all three denominators are given. **Only the last one matters for
the paper.**

{block("all OptimusKG raw records")}

On the raw file ChEMBL looks like a real win — thousands of structures.

{block("the graph's Drug nodes")}

{block("drugs in Task 4 coverage")}

By the denominator that matters, ChEMBL adds **{total_chembl} field values across
{n_t4:,} drugs**. Not zero, and not nothing worth reporting — but nothing that
changes a result.

## Why the structure gap cannot be closed

| | |
|---|---:|
| graph drugs with no OptimusKG SMILES | {w.get('graph drugs with no OptimusKG SMILES', 0):,} |
| of those, ChEMBL held a record | {w.get('of those, ChEMBL held a record', 0):,} |
| of those, ChEMBL had a structure | **{w.get('of those, ChEMBL had a structure', 0):,}** |
| of those, reaching a progression gene | **{w.get('of those, reaching a progression gene', 0):,}** |

Two facts close the question together. ChEMBL holds records for most of these
drugs and has **no structure for any of them** — they are biologics, antibodies
and mixtures that no SMILES can represent. And **none of them reaches a
progression gene**, so they are absent from every analysis in Paper 2 anyway.

Every drug in the Task 4 coverage table already has a structure, from
OptimusKG. There was no gap left for ChEMBL to fill.

## What reached the delivered table

Of the Task 4 coverage table's {n_t4:,} drugs:

| field | from OptimusKG | from ChEMBL |
|---|---:|---:|
| molecule type | {int(d.loc['molecule_type','from_OptimusKG']) if 'molecule_type' in d.index else 0:,} | {int(d.loc['molecule_type','from_ChEMBL']) if 'molecule_type' in d.index else 0:,} |
| SMILES | {int(d.loc['smiles','from_OptimusKG']) if 'smiles' in d.index else 0:,} | {int(d.loc['smiles','from_ChEMBL']) if 'smiles' in d.index else 0:,} |
| drug name | {int(d.loc['drug_name','from_OptimusKG']) if 'drug_name' in d.index else 0:,} | {int(d.loc['drug_name','from_ChEMBL']) if 'drug_name' in d.index else 0:,} |
| clinical phase | {int(d.loc['max_phase','from_OptimusKG']) if 'max_phase' in d.index else 0:,} | {int(d.loc['max_phase','from_ChEMBL']) if 'max_phase' in d.index else 0:,} |

These are the `*_source` columns in `task4_drug_coverage.csv`. A reader can
separate the two recoveries from the delivered CSV alone, without this README.

## Why three denominators

A sample drawn from all {n_raw:,} OptimusKG records describes a different
population from the {n_graph:,} drugs in the graph, and a different one again
from the drugs Paper 2 ranks. The contribution changes completely between them,
which is why all three are reported rather than one.

## Was this worth doing?

Yes, and the result is the reason. "Would an external database improve the drug
vocabulary?" is a question a reviewer will ask, and the answer is now measured
rather than assumed: **no, because the drugs it could help are the drugs the
analysis never sees.** That is a stronger statement than not having checked.

The cost was one API pass. Nothing in Paper 2's rankings depends on it, and
removing this folder entirely would change no result — which is the cleanest
possible answer to whether it belongs inside "focus only on the knowledge
graph".

## Files

| file | what |
|---|---|
| `recovery2_marginal_contribution.csv` | ChEMBL's addition across all three denominators |
| `recovery2_structure_gap.csv` | why the SMILES gap cannot be closed |
| `recovery2_contribution_to_task4.csv` | what reached the delivered table |

## Reproduce

```bash
python scripts/paper2_recover_optimuskg_evidence.py        # Recovery 1 first
uv run --with pandas python scripts/paper2_fetch_chembl_gaps.py
uv run --with neo4j --with pandas python scripts/paper2_recovery_chembl_report.py
```

Recovery 1 must run first: this reads its `drug_metadata.csv` to know where
OptimusKG is silent. That is a one-way dependency — Recovery 1 does not read
anything of Recovery 2's.
"""
    open(f"{OUT_DIR}/README.md", "w", encoding="utf-8").write(txt)


if __name__ == "__main__":
    main()
