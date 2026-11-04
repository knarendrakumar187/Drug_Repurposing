# Recovery 2 — ChEMBL, an external database

**Kept separate from Recovery 1 (OptimusKG) on the mentor's instruction of
13 October 2026.** This folder reports only ChEMBL's *marginal* contribution —
what it added on top of what OptimusKG already supplied — and never pools its
numbers with Recovery 1's.

The separation matters beyond bookkeeping. Recovery 1 re-read a source the graph
was already built from. Recovery 2 left the knowledge graph and queried an
outside API. Those are different acts, and a reader judging "focus only on the
knowledge graph" needs to see which fields came from which.

## The answer: almost nothing, and the reason is specific

7,665 records were fetched from the ChEMBL API, for drugs where
OptimusKG was short. The contribution depends entirely on which drugs you count,
which is why all three denominators are given. **Only the last one matters for
the paper.**

**all OptimusKG raw records** — 18,640 drugs

| field | OptimusKG supplied | OptimusKG silent | ChEMBL added |
|---|---:|---:|---:|
| drug_name | 14,518 | 4,122 | **938** |
| molecule_type | 13,183 | 5,457 | **2,500** |
| max_phase | 11,007 | 7,633 | **688** |
| smiles | 10,262 | 8,378 | **2,638** |

On the raw file ChEMBL looks like a real win — thousands of structures.

**the graph's Drug nodes** — 12,025 drugs

| field | OptimusKG supplied | OptimusKG silent | ChEMBL added |
|---|---:|---:|---:|
| drug_name | 10,862 | 1,163 | **16** |
| molecule_type | 11,623 | 402 | **0** |
| max_phase | 9,902 | 2,123 | **28** |
| smiles | 9,036 | 2,989 | **0** |

**drugs in Task 4 coverage** — 2,404 drugs

| field | OptimusKG supplied | OptimusKG silent | ChEMBL added |
|---|---:|---:|---:|
| drug_name | 1,756 | 648 | **8** |
| molecule_type | 2,404 | 0 | **0** |
| max_phase | 1,432 | 972 | **17** |
| smiles | 2,404 | 0 | **0** |

By the denominator that matters, ChEMBL adds **25 field values across
2,404 drugs**. Not zero, and not nothing worth reporting — but nothing that
changes a result.

## Why the structure gap cannot be closed

| | |
|---|---:|
| graph drugs with no OptimusKG SMILES | 2,989 |
| of those, ChEMBL held a record | 2,654 |
| of those, ChEMBL had a structure | **0** |
| of those, reaching a progression gene | **0** |

Two facts close the question together. ChEMBL holds records for most of these
drugs and has **no structure for any of them** — they are biologics, antibodies
and mixtures that no SMILES can represent. And **none of them reaches a
progression gene**, so they are absent from every analysis in Paper 2 anyway.

Every drug in the Task 4 coverage table already has a structure, from
OptimusKG. There was no gap left for ChEMBL to fill.

## What reached the delivered table

Of the Task 4 coverage table's 2,404 drugs:

| field | from OptimusKG | from ChEMBL |
|---|---:|---:|
| molecule type | 2,404 | 0 |
| SMILES | 2,404 | 0 |
| drug name | 1,756 | 8 |
| clinical phase | 1,432 | 17 |

These are the `*_source` columns in `task4_drug_coverage.csv`. A reader can
separate the two recoveries from the delivered CSV alone, without this README.

## Why three denominators

A sample drawn from all 18,640 OptimusKG records describes a different
population from the 12,025 drugs in the graph, and a different one again
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
