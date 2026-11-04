# Paper 2 — Task 2: drug–gene relationships

Gene set: the corrected Paper 1 union of **6,222 progression genes** (Task 1). "Inflammation" columns hold Paper 1's ordered disease-group axis.

**Source:** the present Neo4j graph only (live, read-only; no CSV fallback). Query: `MATCH (d:Drug)-[r:ACTS_ON]->(g:Gene) WHERE g.ensembl_id IN $ids`.
Script: `Paper2_Results/code/task2_drug_gene_edges.py`. Output: `drug_gene_edges.csv` (14,164 rows).

## Graph check

| Check | Result |
|---|---|
| Drug nodes in graph | 12,544 (write log after 21 Oct: 12,544) |
| ACTS_ON relationships in graph | 25,890 (write log: 25,890; of which added 21 Oct: 5,216) |
| Graph matches the 21 Oct write log | YES |
| Duplicate drug–gene pairs among progression-gene edges | 0 |

## Headline numbers

| Measure | Value |
|---|---|
| Progression genes with ≥1 drug | 1,391 / 6,222 (22.4%) |
| Unique drugs | 5,118 (1,776 flagged approved) |
| Total drug–gene relationships | 14,164 |
| Relationships with known direction (down/up) | 5,644 (39.8%) |
| …comparable with Paper 1 gene direction (positive/negative) | 5,644 (39.8%) |
| …comparable with an F0–F4 trajectory cluster | 3,908 (27.6%) |
| 23-gene shortlist genes with ≥1 drug | 17 / 23 |

"Comparable" means the drug pushes the gene product down or up **and** Paper 1 gives the gene a direction (fibrosis direction; disease-group direction for disease-group-only genes) or an F0–F4 trajectory cluster. Whether each comparison is direction-consistent is Task 7. This is an early mechanism-consistency check, not therapeutic proof.

## By link origin

| Origin | Edges | Genes | Drugs | Known direction | % |
|---|---|---|---|---|---|
| ChEMBL/DGIdb/GtoPdb/OpenTargets (added 21 Oct 2026) | 5,216 | 767 | 2,623 | 4,346 | 83.3% |
| OptimusKG | 8,948 | 1,072 | 3,027 | 1,298 | 14.5% |

## Drug action direction

| Direction | Edges | % |
|---|---|---|
| down | 4,604 | 32.5% |
| up | 1,040 | 7.3% |
| mixed | 43 | 0.3% |
| unknown | 8,477 | 59.8% |

- OptimusKG edges: direction from the consolidated `action` (down: ALLOSTERIC_ANTAGONIST, ANTAGONIST, BLOCKER, DEGRADER, INHIBITOR, INVERSE_AGONIST, NEGATIVE_ALLOSTERIC_MODULATOR, NEGATIVE_MODULATOR; up: ACTIVATOR, AGONIST, OPENER, PARTIAL_AGONIST, POSITIVE_ALLOSTERIC_MODULATOR, POSITIVE_MODULATOR, STABILISER; mixed: BINDING_AGENT, MODULATOR, RELEASING_AGENT, SUBSTRATE or sources disagree; unknown: CARRIER, ENZYME, TARGET, TRANSPORTER or blank). OptimusKG `relation_conflict` edges: 0.
- Added edges: direction from the edge's `push` field (down / up; "conflicting" counted as mixed; unknown when no action label).
- `Action_raw` keeps the action text exactly as stored on the edge; blank actions are written as "unknown".
- OptimusKG edges with no direction take the GtoPdb/DGIdb/ChEMBL direction stored on the same edge (`external_action_direction`): 297 edges. It never overrides a direction OptimusKG gives.

| Direction source | Edges |
|---|---|
| none | 8,477 |
| added edge (push) | 4,361 |
| OptimusKG action | 1,029 |
| GtoPdb/DGIdb/ChEMBL on OptimusKG edge | 297 |

## By progression class

| Class | Genes | Genes with drug | % | Unique drugs | Edges |
|---|---|---|---|---|---|
| fibrosis_only | 2,023 | 471 | 23.3% | 2,904 | 5,540 |
| shared | 1,822 | 476 | 26.1% | 2,815 | 5,124 |
| inflammation_only | 2,377 | 444 | 18.7% | 2,057 | 3,500 |

(Statistics comparing these groups are Task 3.)

## Notes

- `Database_ID` is the stable pair `drug ID | Ensembl gene ID`; the graph holds no per-assertion accession IDs.
- No binding strength, dose or tissue context is stored on the edges.
- Action labels such as TARGET, ENZYME, TRANSPORTER and CARRIER say what the gene product is for the drug, not whether the drug switches it down or up; they count as unknown direction.
