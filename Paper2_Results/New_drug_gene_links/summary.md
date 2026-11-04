# Drug–gene links added to the graph from ChEMBL, DGIdb, GtoPdb and Open Targets (graph preparation, step 0.4)

Gene set: the 6,222 Paper 1 progression genes. The headline links below were written to the graph on 21 Oct 2026 (`write_new_links_to_graph.py`; log in `graph_write/`).

## In one sentence
Outside databases list **5,216 drug–gene links for approved or clinical-trial medicines** (curated drug-target evidence) on our genes that the graph does not have; they reach **319 progression genes that have no drug at all in the graph** (genes with a drug: 1,072 → 1,391).

## Sources searched (all free)
| Source | Version | How genes were matched |
|---|---|---|
| IUPHAR/BPS Guide to Pharmacology | 2026.3.1 | Ensembl ID (human, endogenous ligands removed) |
| DGIdb | 5 (interactions.tsv, downloaded Oct 2026) | HGNC symbol (whole database, kept for our genes) |
| ChEMBL drug mechanisms | REST API (Oct 2026) | gene symbol of human targets |
| Open Targets Platform | Open Targets Platform 26.09 (latest) | Ensembl ID |

Each drug is reduced to its ChEMBL parent compound (salts count once). A link is **new** if the graph has no Drug–[ACTS_ON]→Gene link for that drug and gene (all 8,948 graph links on our genes were compared).

## Does the matching work? (positive control)
- The graph's 1,030 links that came from Open Targets: Open Targets finds **1,019 (99%)**, ChEMBL finds **1,028 (100%)**.
- Neo4j spot check (live graph, 21 Oct 2026): 120 random "new" links → **0** have any Drug–Gene relationship in Neo4j; 150 random new-link drugs → 121 already exist as Drug nodes but only 7 have any ACTS_ON link. So the graph knows most of these drugs but is missing their targets.

| source | graph links tested | n | found by source | % found |
|---|---|---|---|---|
| GtoPdb | all graph links | 8,948 | 736 | 8.2 |
| GtoPdb | pharmacological graph links | 4,987 | 696 | 14.0 |
| DGIdb | all graph links | 8,948 | 1,539 | 17.2 |
| DGIdb | pharmacological graph links | 4,987 | 1,167 | 23.4 |
| ChEMBL | all graph links | 8,948 | 1,035 | 11.6 |
| ChEMBL | pharmacological graph links | 4,987 | 1,033 | 20.7 |
| OpenTargets | all graph links | 8,948 | 1,025 | 11.5 |
| OpenTargets | pharmacological graph links | 4,987 | 1,023 | 20.5 |

## Two errors found in DGIdb and removed
1. **DGIdb's copy of Guide to Pharmacology is on the wrong genes.** Guide to Pharmacology target numbers were read as NCBI gene numbers (target 215 = dopamine D2 receptor, but NCBI gene 215 = ABCD1). All 6,095 such rows on our genes were dropped; Guide to Pharmacology is used directly instead.
2. **DGIdb's copy of ChEMBL keeps virus and bacteria targets** matched by short names (hepatitis B polymerase "P" → human OCA2). These 4,748 rows were dropped; ChEMBL is used directly with human targets only.

## How good is the evidence? (all new links)
| evidence | approved | clinical trial | research compound | all |
|---|---|---|---|---|
| curated drug-target | 2,218 | 2,998 | 4,478 | 9,694 |
| curated, TTD only (gene matched by protein name) | 216 | 143 | 837 | 1,196 |
| lab binding screen | 2,002 | 887 | 5,914 | 8,803 |
| literature / biomarker only | 3,616 | 283 | 1,739 | 5,638 |
| all | 8,052 | 4,311 | 12,968 | 25,331 |

- **curated drug-target**: Guide to Pharmacology, ChEMBL, Open Targets, or DGIdb sources that give the gene by symbol/UniProt (TEND, TALC, FDA, clinical-trial target lists). **Used for the headline numbers.**
- **TTD only**: curated, but TTD matches genes from protein names, which sometimes goes wrong (e.g. "KIR" → GEM). Kept, not in the headline.
- **lab binding screen** (Drug Target Commons): measured binding in an assay; mostly research compounds.
- **literature / biomarker only** (NCI, PharmGKB, CIViC, cancer knowledge bases): "drug is related to gene" (response, resistance, expression), **not** "drug acts on gene". Not counted as a target.

## By source
| source | drug–gene pairs on our genes | already in graph | NEW pairs | NEW pairs with a medicine | found ONLY by this source |
|---|---|---|---|---|---|
| GtoPdb | 6,471 | 734 | 5,737 | 1,703 | 4,760 |
| DGIdb | 19,443 | 1,528 | 17,915 | 8,885 | 16,579 |
| ChEMBL | 4,755 | 1,029 | 3,726 | 3,720 | 134 |
| OpenTargets | 4,629 | 1,019 | 3,610 | 3,607 | 21 |

## Headline: new links for medicines with curated evidence
| | links | drugs | genes | genes with NO drug in the graph before |
|---|---|---|---|---|
| approved or clinical trial | 5,216 | 2,623 | 767 | **319** |
| …approved only | 2,218 | 1,168 | 508 | 187 |

- 3,732 of 5,216 links (72%) are reported by two or more sources.
- 230 of the 319 newly reached genes are fibrosis genes (fibrosis-only or shared); **49 are Tier 1** (high confidence).
- Push direction is known for 4,346 links: **2,829 direction-consistent**, 1,517 inconsistent (direction rule: the drug pushes down a gene that rises, or up a gene that falls).

| progression_class | pairs | genes | genes_with_no_drug_before |
|---|---|---|---|
| fibrosis_only | 1,813 | 268 | 113 |
| inflammation_only | 1,623 | 219 | 89 |
| shared | 1,780 | 280 | 117 |

### Drugs already tested against fibrosis that the graph was missing (sanity check)
| Gene | Drug | found | push | direction rule |
|---|---|---|---|---|
| CCN2 (CTGF) | pamrevlumab | yes | down | direction-consistent |
| LGALS3 (galectin-3) | belapectin | yes | down | direction-consistent |
| ITGB5 (αv integrin) | abituzumab | yes | down | direction-consistent |
| ITGB5 (αv integrin) | cilengitide | yes | down | direction-consistent |
| ITGB8 (αv integrin) | abituzumab | yes | down | direction-consistent |
| ITGA4 | natalizumab | yes | down | direction-consistent |
| ITGA4 | vedolizumab | yes | down | direction-consistent |
| TGFB1 (TGF-β1) | luspatercept | yes | down | direction-consistent |

These are drugs that have been tested against fibrosis (pamrevlumab in lung fibrosis, belapectin in NASH cirrhosis, αv-integrin blockers), which suggests the new layer adds real biology, not noise.

**Limit of the simple direction rule:** soluble guanylate cyclase stimulators (riociguat, vericiguat, praliciguat) are found on GUCY1A1 / GUCY1B1 but are labelled *direction-inconsistent*: these genes rise with fibrosis and the drugs push them up, yet sGC stimulation reduces fibrosis in animal models. A rising gene is not always a harmful gene; direction labels should be read with biology in mind.

## Caveats
- Protein-complex targets are expanded to every subunit gene (e.g. proteasome, tubulin, CD3), as the graph itself does; antibody–drug conjugates (…vedotin) also hit tubulin through their payload.
- A few drugs appear twice (salt or prodrug with its own ChEMBL parent, or a DGIdb RxNorm ID that could not be matched to ChEMBL).
- Literature-only and TTD-only links are listed but not counted in the headline.

## Files
- `new_links_medicines_curated.csv`: headline links (one row per drug–gene, with sources, action, push direction, gene tier).
- `new_links_medicines.csv`: all medicine links, any evidence. `new_links_all.csv`: everything incl. research compounds.
- `genes_reached_for_first_time.csv`, `medicines_by_drug.csv`, `by_source.csv`, `by_evidence_and_stage.csv`, `medicines_direction_consistency.csv`, `positive_control.csv`, `numbers.json`.
- `raw/`: cached downloads; `raw_rows_all_sources.csv.gz`: every source row used.
- Scripts: `code/new_links_search.py` (search), `code/new_links_report.py` (this page + diagram).
