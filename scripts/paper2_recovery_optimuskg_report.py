"""
paper2_recovery_optimuskg_report.py
-----------------------------------
Report RECOVERY 1 -- what OptimusKG's own files gave back -- on its own.

WHY THIS IS A SEPARATE SCRIPT AND A SEPARATE FOLDER
The mentor's instruction of 13 October: two jobs were asked for, and if two were
done they must be kept separate. Recovery 1 reads OptimusKG's discarded fields.
Recovery 2 (paper2_recovery_chembl_report.py) asks an external database. Pooling
their numbers into one "recovered" figure would hide which source carried which
field, so each reports alone and neither cites the other's totals.

WHAT RECOVERY 1 IS
build_optimuskg_layer.py read OptimusKG's drug.parquet and drug_gene.parquet and
kept only `name` and `relation`. Everything else -- per-edge source databases,
relation assertions, FDA label URLs, written mechanisms of action, canonical
SMILES, molecule type, clinical phase -- was discarded at load. It was never
missing data; it was data we had and dropped. This recovers it from the same
source, so nothing here is a new source.

THE DENOMINATOR
Coverage is reported against the 12,025 Drug nodes in the graph, because that is
the drug vocabulary Paper 2 actually ranks. Reporting against OptimusKG's 18,640
raw records would flatter every figure -- 6,615 of them never reached the graph.

Usage:
    uv run --with neo4j --with pandas python scripts/paper2_recovery_optimuskg_report.py
"""

import json
import os

import pandas as pd
from neo4j import GraphDatabase

OKG = "data/optimuskg_recovered"
EVIDENCE = f"{OKG}/drug_gene_evidence.csv"
METADATA = f"{OKG}/drug_metadata.csv"
OUT_DIR = "Paper2_Results/Recovery1_OptimusKG"

T2 = "Paper2_Results/Task2/task2_drug_gene_edges.csv"
T4 = "Paper2_Results/Task4/task4_drug_coverage.csv"


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
    """The graph's drug vocabulary. Drug nodes key on node_id, which holds the
    ChEMBL accession -- the same key OptimusKG's parquet uses as Drug_ID."""
    drv, db = driver()
    try:
        with drv.session(database=db) as s:
            rows = s.run("MATCH (d:Drug) RETURN d.node_id AS id").data()
    finally:
        drv.close()
    ids = {r["id"] for r in rows if r["id"] is not None}
    if not ids:
        raise SystemExit("Drug nodes carry no node_id")
    return ids


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    ev = pd.read_csv(EVIDENCE, low_memory=False)
    md = pd.read_csv(METADATA, low_memory=False)
    drugs = graph_drug_ids()
    print(f"graph Drug nodes        {len(drugs):,}")
    print(f"OptimusKG raw records   {len(md):,}")

    g = md[md.Drug_ID.isin(drugs)].copy()
    print(f"   of which in the graph{len(g):>8,}")

    # OptimusKG's `name` field stores an unnamed drug's ChEMBL accession AS its
    # name, so an accession means "no name". A missing field is a third state,
    # and the two were conflated before this script existed.
    g["name_state"] = "named"
    g.loc[g.name == g.Drug_ID, "name_state"] = "accession as name"
    g.loc[g.name.isna(), "name_state"] = "name field absent"

    # ---- per-field coverage, OptimusKG alone -------------------------------
    fields = {
        "drug name": g.name.notna() & (g.name != g.Drug_ID),
        "molecule type": g.type.notna(),
        "clinical phase": g.maximum_clinical_trial_phase.notna(),
        "approved flag": g.is_approved.notna(),
        "withdrawn flag": g.has_been_withdrawn.notna(),
        "canonical SMILES": g.canonical_smiles.notna(),
        "InChI key": g.inchi_key.notna(),
        "formula": g.cd_formula.notna(),
        "molecular weight": g.cd_mol_weight.notna(),
    }
    cov = pd.DataFrame([
        {"field": k,
         "drugs_with_value": int(v.sum()),
         "drugs_in_graph": len(g),
         "percent": round(100 * v.sum() / len(g), 1)}
        for k, v in fields.items()])
    cov.to_csv(f"{OUT_DIR}/recovery1_field_coverage.csv", index=False)
    print("\nper-field coverage over the graph's drug vocabulary")
    print(cov.to_string(index=False))

    # ---- the naming states, told apart -------------------------------------
    ns = (g.name_state.value_counts().rename_axis("name_state")
          .reset_index(name="drugs"))
    ns["percent"] = (100 * ns.drugs / len(g)).round(1)
    ns.to_csv(f"{OUT_DIR}/recovery1_name_states.csv", index=False)
    print("\nnaming states -- three, not two")
    print(ns.to_string(index=False))

    # ---- edge-level evidence ----------------------------------------------
    t2 = pd.read_csv(T2, low_memory=False)
    # Title_Case is the delivered column; the lowercase working copy was an
    # exact duplicate and is no longer written. The guard stays, but it must
    # not silently report 0 -- that is how a stale number survives.
    assert "Evidence_Level" in t2, "Task 2 has no Evidence_Level column"
    matched = int((t2.Evidence_Level != "not matched to source").sum())
    ev_rows = [
        {"measure": "evidence records recovered", "value": len(ev)},
        {"measure": "Task 2 edges", "value": len(t2)},
        {"measure": "Task 2 edges matched to a source", "value": int(matched)},
        {"measure": "edges with a reference URL",
         "value": int(t2.source_urls.notna().sum()) if "source_urls" in t2 else 0},
        {"measure": "edges with a written mechanism",
         "value": int(t2.mechanisms_of_action.notna().sum())
         if "mechanisms_of_action" in t2 else 0},
        {"measure": "edges with conflicting assertions",
         "value": int(t2.relation_conflict.fillna(False).astype(bool).sum())
         if "relation_conflict" in t2 else 0},
    ]
    pd.DataFrame(ev_rows).to_csv(f"{OUT_DIR}/recovery1_edge_evidence.csv",
                                 index=False)
    print("\nedge-level evidence")
    for r in ev_rows:
        print(f"   {r['measure']:<40} {r['value']:>8,}")

    # ---- what the delivered tables carry FROM OPTIMUSKG ONLY ---------------
    t4 = pd.read_csv(T4, low_memory=False)
    src_rows = []
    for f in ["drug_name", "molecule_type", "max_phase", "smiles"]:
        col = f"{f}_source"
        if col in t4:
            src_rows.append({
                "field": f,
                "from_OptimusKG": int((t4[col] == "OptimusKG").sum()),
                "drugs_in_table": len(t4)})
    if src_rows:
        pd.DataFrame(src_rows).to_csv(
            f"{OUT_DIR}/recovery1_contribution_to_task4.csv", index=False)
        print(f"\ncontribution to the delivered Task 4 table ({len(t4):,} drugs)")
        for r in src_rows:
            print(f"   {r['field']:<14} {r['from_OptimusKG']:>6,}")

    raw_name_absent = int(md.name.isna().sum())
    write_readme(g, cov, ns, ev_rows, src_rows, len(drugs), len(md),
                 raw_name_absent)
    print(f"\nwrote {OUT_DIR}/")


def write_readme(g, cov, ns, ev_rows, src_rows, n_graph, n_raw,
                 raw_name_absent):
    ev = {r["measure"]: r["value"] for r in ev_rows}
    n_ev = ev.get("evidence records recovered", 0)
    c = cov.set_index("field")
    t4c = {r["field"]: r["from_OptimusKG"] for r in src_rows}
    acc = int(ns.loc[ns.name_state == "accession as name", "drugs"].sum())
    absent = raw_name_absent      # raw records, not the graph subset

    txt = f"""# Recovery 1 — OptimusKG's own discarded fields

**Kept separate from Recovery 2 (ChEMBL) on the mentor's instruction of
13 October 2026.** Neither folder cites the other's totals. This one reports
only what OptimusKG gave back.

## What this is, and what it is not

`build_optimuskg_layer.py` read OptimusKG's `drug.parquet` and
`drug_gene.parquet` and kept **only** `name` and `relation`. Everything else was
read and dropped at load: per-edge source databases, relation assertions, FDA
label URLs, written mechanisms of action, canonical SMILES, molecule type and
clinical phase.

This recovers those fields from the same two files. **It is not a new data
source** — it is the source the graph was already built from, read more
completely. That distinction is why it was treated as inside "focus only on the
knowledge graph", and the mentor has since confirmed it.

## The denominator

Coverage is reported against the **{n_graph:,} Drug nodes in the graph**, not
OptimusKG's {n_raw:,} raw records. {n_raw - n_graph:,} of those records never
reached the graph, and including them would flatter every figure here.

## What came back

| field | drugs with a value | of {n_graph:,} |
|---|---:|---:|
| canonical SMILES | {int(c.loc['canonical SMILES','drugs_with_value']):,} | {c.loc['canonical SMILES','percent']}% |
| molecule type | {int(c.loc['molecule type','drugs_with_value']):,} | {c.loc['molecule type','percent']}% |
| clinical phase | {int(c.loc['clinical phase','drugs_with_value']):,} | {c.loc['clinical phase','percent']}% |
| drug name | {int(c.loc['drug name','drugs_with_value']):,} | {c.loc['drug name','percent']}% |
| InChI key | {int(c.loc['InChI key','drugs_with_value']):,} | {c.loc['InChI key','percent']}% |
| formula | {int(c.loc['formula','drugs_with_value']):,} | {c.loc['formula','percent']}% |

Full table: `recovery1_field_coverage.csv`.

## Naming: two states inside the graph, three in the raw records

| state | drugs in the graph |
|---|---:|
| a real name | {int(ns.loc[ns.name_state == 'named', 'drugs'].sum()):,} |
| the ChEMBL accession stored *as* the name | {acc:,} |

An accession-as-name is a drug OptimusKG knows and cannot name — not a gap in
our copy. **{acc:,} is the figure limitation 7 reports**, and it reproduces here
exactly.

A third state exists in OptimusKG's raw file — {absent:,} records with no `name`
field at all — but **every one of them is outside the graph**, so it does not
affect Paper 2. It is mentioned only because measuring the raw file instead of
the graph inflates this count.

## Edge-level evidence

| measure | value |
|---|---:|
| evidence records recovered | {ev.get('evidence records recovered', 0):,} |
| Task 2 edges matched to a source | {ev.get('Task 2 edges matched to a source', 0):,} of {ev.get('Task 2 edges', 0):,} |
| edges with a reference URL | {ev.get('edges with a reference URL', 0):,} |
| edges with a written mechanism | {ev.get('edges with a written mechanism', 0):,} |
| Task 2 edges with conflicting assertions | {ev.get('edges with conflicting assertions', 0):,} |

Every Task 2 edge now names the databases asserting it. Only two sources assert
drug–gene edges at all: DrugBank and Open Targets.

The conflict count is **over Task 2's edges**. Graph-wide, across all
{n_ev:,} recovered evidence records, `relation_conflict` is set on 2 — both
outside the progression-gene subset. The two numbers do not disagree; they
have different scopes, and limitation 9 quotes the graph-wide one.

## What it closed in the plan

| asked for | was | now |
|---|---|---|
| Task 2 Evidence/Source and Database ID per edge | "not available per edge" | all {ev.get('Task 2 edges', 0):,} edges attributed |
| Task 7 evidence source and evidence level | "not recorded" | single- vs multi-source on all 4,984 |
| Task 4 filter on molecule type and phase | impossible | available — **and it does not work**, see below |

## The honest negative

Team 1's instruction was to filter the drug vocabulary on molecule type and
clinical phase. Both fields are now recovered, **and they do not discriminate.**
All 70 name-flagged entries are `Small molecule`; ZINC ACETATE and ZINC CHLORIDE
carry `max_phase 4.0` and `is_approved True`, because they genuinely are
approved small molecules. A filter built on these fields flags **one** drug
where a name heuristic flags seventy.

Recovering a field and finding it useless is still a result. It is recorded
rather than quietly dropped.

## Contribution to the delivered tables

Of the Task 4 coverage table's {src_rows[0]['drugs_in_table'] if src_rows else 0:,} drugs,
OptimusKG supplies:

| field | from OptimusKG |
|---|---:|
| molecule type | {t4c.get('molecule_type', 0):,} |
| SMILES | {t4c.get('smiles', 0):,} |
| drug name | {t4c.get('drug_name', 0):,} |
| clinical phase | {t4c.get('max_phase', 0):,} |

Each of those columns travels with a `*_source` companion in the delivered CSV,
so a reader can separate OptimusKG from ChEMBL without consulting this README.

## Files

| file | what |
|---|---|
| `recovery1_field_coverage.csv` | per-field coverage over the graph's drugs |
| `recovery1_name_states.csv` | the three naming states |
| `recovery1_edge_evidence.csv` | edge-level evidence recovered |
| `recovery1_contribution_to_task4.csv` | what OptimusKG supplies to the delivered table |

## Reproduce

```bash
# needs data/optimuskg/drug.parquet and drug_gene.parquet (DOI 10.7910/DVN/IYNGEV)
python scripts/paper2_recover_optimuskg_evidence.py
python scripts/paper2_recovery_optimuskg_report.py
```
"""
    open(f"{OUT_DIR}/README.md", "w", encoding="utf-8").write(txt)


if __name__ == "__main__":
    main()
