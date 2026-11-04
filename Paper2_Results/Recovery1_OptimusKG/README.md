# Recovery 1 — OptimusKG's own discarded fields

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

Coverage is reported against the **12,025 Drug nodes in the graph**, not
OptimusKG's 18,640 raw records. 6,615 of those records never
reached the graph, and including them would flatter every figure here.

## What came back

| field | drugs with a value | of 12,025 |
|---|---:|---:|
| canonical SMILES | 9,036 | 75.1% |
| molecule type | 11,623 | 96.7% |
| clinical phase | 9,902 | 82.3% |
| drug name | 10,862 | 90.3% |
| InChI key | 9,069 | 75.4% |
| formula | 1,965 | 16.3% |

Full table: `recovery1_field_coverage.csv`.

## Naming: two states inside the graph, three in the raw records

| state | drugs in the graph |
|---|---:|
| a real name | 10,862 |
| the ChEMBL accession stored *as* the name | 1,163 |

An accession-as-name is a drug OptimusKG knows and cannot name — not a gap in
our copy. **1,163 is the figure limitation 7 reports**, and it reproduces here
exactly.

A third state exists in OptimusKG's raw file — 2,699 records with no `name`
field at all — but **every one of them is outside the graph**, so it does not
affect Paper 2. It is mentioned only because measuring the raw file instead of
the graph inflates this count.

## Edge-level evidence

| measure | value |
|---|---:|
| evidence records recovered | 36,611 |
| Task 2 edges matched to a source | 8,948 of 8,948 |
| edges with a reference URL | 1,028 |
| edges with a written mechanism | 1,030 |
| Task 2 edges with conflicting assertions | 0 |

Every Task 2 edge now names the databases asserting it. Only two sources assert
drug–gene edges at all: DrugBank and Open Targets.

The conflict count is **over Task 2's edges**. Graph-wide, across all
36,611 recovered evidence records, `relation_conflict` is set on 2 — both
outside the progression-gene subset. The two numbers do not disagree; they
have different scopes, and limitation 9 quotes the graph-wide one.

## What it closed in the plan

| asked for | was | now |
|---|---|---|
| Task 2 Evidence/Source and Database ID per edge | "not available per edge" | all 8,948 edges attributed |
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

Of the Task 4 coverage table's 2,404 drugs,
OptimusKG supplies:

| field | from OptimusKG |
|---|---:|
| molecule type | 2,404 |
| SMILES | 2,404 |
| drug name | 1,756 |
| clinical phase | 1,432 |

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
