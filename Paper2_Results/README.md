# Paper 2: knowledge-graph drug repurposing for liver fibrosis progression

This folder holds the code and results of the Paper 2 analysis: preparing the knowledge graph (including drug information from ChEMBL, DGIdb, GtoPdb and Open Targets), then Tasks 1–4, the direction filter, the metapath B/BB analysis and the final ranking of liver-treating repurposing candidates. Every analysis step uses the same, final graph.

Input from Paper 1 (final combined version, 10 October 2026): 6,222 progression genes (3,845 fibrosis-associated, 4,199 ordered disease-group-associated, 1,822 shared; fibrosis-only 2,023, disease-group-only 2,377). Tier 1/2/3 = 447 / 1,403 / 1,995; 23-gene shortlist. The second axis is the **ordered disease-group** axis (control/NAFL/NASH); some file and column names call it "inflammation".

## Step 0: graph preparation

The graph is the Paper 1 Neo4j graph (OptimusKG knowledge layer + expression layer + `TRACKS_FIBROSIS` / `TRACKS_INFLAMMATION` progression links), prepared for Paper 2 in these steps, in this order. Every write step has a dry run, a backup or undo file, and a write log.

| # | What was added to the graph | Code | Outputs / logs |
|---|---|---|---|
| 0.1 | **Drug metadata and drug–gene evidence recovered from OptimusKG** (fields not loaded when the graph was built: SMILES, molecule type, phase, sources per edge, mechanism text, reference URLs) | `scripts/paper2_recover_optimuskg_evidence.py`, `scripts/paper2_recovery_optimuskg_report.py` | `Recovery1_OptimusKG/` |
| 0.2 | **Missing drug metadata filled from ChEMBL** (kept separate from 0.1; OptimusKG values take precedence) | `scripts/paper2_fetch_chembl_gaps.py`, `scripts/paper2_recovery_chembl_report.py` | `Recovery2_ChEMBL/` |
| 0.3 | 0.1 + 0.2 written onto Drug nodes and `ACTS_ON` edges (new properties only, never overwrites) | `scripts/paper2_integrate_recovered_evidence.py`, `scripts/paper2_write_recovered_to_graph.py` (`--apply`, `--rollback`) | |
| 0.4 | **New drug–gene links from ChEMBL, DGIdb, GtoPdb and Open Targets** (approved/clinical drugs with a curated drug–target source): **5,216 `ACTS_ON` links and 519 new Drug nodes**, written 21 Oct 2026. The links carry `layer = 'external (Idea 4)'` (the label stored in the graph), `sources`, `action_labels`, `push` (direction) and `mechanism`. Graph went from 126,244 nodes / 35,132,831 relationships to 126,763 / 35,138,047. | `code/new_links_search.py` → `code/new_links_report.py` → `code/write_new_links_to_graph.py` | `New_drug_gene_links/` (`summary.md`, `numbers.json` with source releases: GtoPdb 2026.3.1, DGIdb 5, ChEMBL REST API Oct 2026, Open Targets 26.09; `graph_write/write_log.json`, `undo_new_links.cypher`) |
| 0.5 | **Drug action directions from GtoPdb/DGIdb/ChEMBL** for OptimusKG links that had none (TARGET/ENZYME/TRANSPORTER/CARRIER): added as `external_action_*` properties on 297 existing `ACTS_ON` links (the `action` field is unchanged), written 21 Oct 2026 | `code/external_action_check.py` → `code/write_external_directions_to_graph.py` | `External_action_check/` (lookup, positive control, `graph_write/`), `Task7_external_layer/external_layer_links.csv` (the 297 links written) |
| 0.6 | **KEGG (legacy) and Hallmark pathways** from MSigDB 2024.1, the same files Paper 1 used for enrichment (the graph had only Reactome) | [`scripts/add_msigdb_pathways.py`](https://github.com/nissankararaoThanvitha/liver_graph/blob/main/scripts/add_msigdb_pathways.py) (Paper 1 repository) | [`data/graph_msigdb/`](https://github.com/nissankararaoThanvitha/liver_graph/tree/main/data/graph_msigdb) in the Paper 1 repository (Hallmark rows; KEGG rows are not redistributed) |
| 0.7 | **Pathway → NAFLD links** (`IMPLICATED_IN`), kept only when Paper 1 enrichment (BH-FDR < 0.05, fibrosis groups) **and** CTD (NAFLD, MESH:D065626) agree: 89 links | [`scripts/build_pathway_nafld_links.py`](https://github.com/nissankararaoThanvitha/liver_graph/blob/main/scripts/build_pathway_nafld_links.py) (Paper 1 repository; needs `CTD_diseases_pathways.csv.gz` in `data/ctd/`) | [`data/graph_pathway_disease/`](https://github.com/nissankararaoThanvitha/liver_graph/tree/main/data/graph_pathway_disease) in the Paper 1 repository |
| 0.8 | **Drug → disease indications for the drugs added in 0.4** (they came with drug–gene links only): ChEMBL_37 indications, matched by ChEMBL ID or exact name, written as **2,087 `TREATS` links** (138 drugs, 866 diseases already in the graph; 172 indications whose disease ID is not in the graph skipped). Same meaning as OptimusKG's TREATS (an indication at any clinical phase); each link carries `source = 'ChEMBL indication'`, `layer = 'ChEMBL indication (added drugs)'` and its highest phase. Written 3 Nov 2026. | `code/task4_added_drug_indications.py` → `code/write_added_drug_indications_to_graph.py` (`--apply`, `--undo`) | `Task4/added_drug_*`, `Task4/graph_write/` |

**Final graph (checked 3 Nov 2026):** 126,999 nodes; 35,160,254 relationships; 15 relationship types (the 14 of Paper 1 plus `IMPLICATED_IN`). Paper 1's figure (126,763 / 35,138,047) is the graph after step 0.4; steps 0.6–0.8 add 236 Pathway nodes, 20,031 `IN_PATHWAY` links, 89 `IMPLICATED_IN` links and 2,087 `TREATS` links.

**Sources (for Methods):** none of the 5,216 links of step 0.4 is TTD-only; all have curated drug-target evidence from ChEMBL, DGIdb, GtoPdb and/or Open Targets. In 701 links DGIdb names TTD as its own upstream source, so Methods may say "DGIdb (which aggregates sources including TTD)".

## Analysis workflow and results — use this for the manuscript

All steps use the final graph (live Neo4j, read-only). Run from the repository root with `NEO4J_PASSWORD` set; `code/run_tasks.ps1 <script> ...` runs scripts with one password prompt.

| Step | What was done | Result | Code | Output |
|---|---|---|---|---|
| Task 1 | Master table of all Paper 1 progression genes (no gene removed) | 6,222 genes: fibrosis-only 2,023 · disease-group-only 2,377 · shared 1,822; all counts match Paper 1 | `code/task1_master_gene_table.py` | `Task1/` |
| Task 2 | Every Drug -ACTS_ON-> Gene link on these genes, with action and source. Direction from the link's own action; where OptimusKG gives only TARGET/ENZYME/TRANSPORTER/CARRIER, the GtoPdb/DGIdb/ChEMBL direction stored on the same link (step 0.5) | 14,164 links (8,948 OptimusKG + 5,216 from step 0.4); 1,391 genes (22.4%); 5,118 drugs; 5,644 links (39.8%) with known direction | `code/task2_drug_gene_edges.py` | `Task2/` |
| Task 3 | Druggability by Paper 1 group vs the 14,794-gene background: Fisher + BH-FDR, pairwise tests, 1,000 degree-matched permutations (seed 42) | Fibrosis-only (21.9%, 1.47×) and shared (24.4%, 1.59×) genes more druggable than degree-matched background, both rising and falling; disease-group-only not (17.3%, 1.07×) | `code/task3_druggability.py` | `Task3/` (bar plot `druggability_barplot.png`) |
| Task 4 | Per-drug coverage (by class, direction, cluster), progression-target fraction, top 50, annotated drug–gene network; candidate filter | 4,623 drugs with a pharmacological link → phase ≥ 2: 3,360 → not withdrawn: 3,279 → no fibrosis-type liver use (graph + ChEMBL indications for the added drugs): **3,163 candidates** → mechanism annotated: **2,587** | `code/task4_added_drug_indications.py`, `code/task4_drug_coverage.py` | `Task4/` |
| Direction | Keep drugs whose action opposes the Paper 1 direction (switch down a rising gene / up a falling gene) on every link with a known direction | 2,402 with a known direction → **1,131 only direction-consistent** (used); 501 mixed and 770 only inconsistent not used | `code/direction_consistent_drugs.py` | `Direction_consistent_drugs/` |
| Metapaths | For each Gene* set: metapath **B** (Drug → Gene → Gene* → Pathway → NAFLD), then **BB** (Drug → Gene → Gene → Gene* → Pathway → NAFLD) for the drugs with B routes. Rules: routes start only from correct-direction targets; in-between genes ≤ 500 interaction partners; pathway 5–200 genes; drug has no TREATS/OFF_LABEL_FOR to NAFLD, NASH or cirrhosis | see table below | `scripts/run_six_metapaths.py --metapaths B`, `scripts/run_two_gene_metapath.py`, `code/metapath_pipeline.py` | `Metapath_runs/` |
| Ranking | Direct diseases of the BB drugs; keep drugs with ≥ 1 TREATS link to a liver disease, remove drugs contraindicated in a liver disease, rank by BB routes; stage table for ranks 1–25 | see table below | `code/metapath_followup.py`, `code/compare_top25.py` | `Metapath_runs/` |

| Run (Gene* set) | Drugs with B routes | BB routes | Drugs with BB routes | Treat a liver disease | Removed (liver-contraindicated) | Final ranked drugs |
|---|---|---|---|---|---|---|
| Run1 top-23 shortlist | 97 | 10,177 | 96 | 13 | 0 | 13 |
| Run2 Tier 1 (447 genes) | 801 | 963,625 | 801 | 101 | 7 | 94 |
| Run3 Tier 2 (1,403 genes) | 827 | 1,762,572 | 826 | 111 | 6 | 105 |
| **Run4 Tier 3 (1,995 genes) — main result** | **728** | **1,251,496** | **717** | **95** | **6** | **89** |

**Main result:** `Metapath_runs/Run4_onlyconsistent_tier3_BB_liver_treating_drugs_final.xlsx` (89 drugs ranked by BB routes; top: Quizartinib, Alvocidib, Milciclib, Zotiraciclib, Trastuzumab, Zanidatamab). For all 25 top drugs, most BB routes reach Gene* whose largest change is at F3 → F4 (`Run4_…_top_drugs_stages.xlsx`).
**Robustness across tiers:** 24 drugs are in the top 25 of Tier 1, Tier 2 and Tier 3; the Tier 2 and Tier 3 top 25 are identical (`Metapath_runs/top25_comparison.md`).

Files per run (`<run>` = e.g. `Run4_onlyconsistent_tier3`): `<run>.xlsx` (kept B routes), `<run>_all_routes.csv` (all B routes with drop reasons), `<run>_BB.xlsx` (Summary + every BB route; long lists split over BB, BB (2)), `<run>_BB_direct_diseases.xlsx`, `<run>_BB_liver_treating_drugs_final.xlsx`, `<run>_BB_top_drugs_stages.xlsx`, `<run>_BB_followup.json`. The full BB route CSVs (up to 320 MB) are not in git; rebuild them with `python Paper2_Results/code/metapath_pipeline.py Run1 Run2 Run3 Run4`.

### Points to state in the manuscript
- Direction is checked only at the drug → target step; gene–gene (INTERACTS_WITH) links have no sign. Routes show a biological connection, not a demonstrated effect.
- Only drugs whose every link with a known direction opposes the progression direction are used; drugs that also push some progression gene the wrong way ("mixed", e.g. kinase inhibitors that also block KDR or EGFR, which fall with fibrosis) are excluded.
- For some top drugs the liver disease they treat is the general label "liver disease" (e.g. Quizartinib, Siremadlin, Bardoxolone methyl, Pacritinib); check the specific indication before describing them.
- All results are repurposing hypotheses, not treatments.

## Scope (from the mentor, 13 October 2026)

- "Repurposing, not reversal": no LINCS/CMap transcriptomic reversal analysis.
- Evidence comes from the knowledge graph (including the drug information added in step 0) and Paper 1 results.
- OptimusKG and ChEMBL recoveries are kept separate (separate scripts, folders and `*_source` columns).

## Standing rules

- Do not trim the graph; filter at read time.
- Keep the progression classes separate (fibrosis-only, disease-group-only, shared; up/down; clusters C1–C5).
- The graph records association direction, not benefit. Outputs are repurposing hypotheses, not treatments.
- Paper 1's tier, trajectory and transition evidence covers only the 3,845 fibrosis-associated genes.

## Graph instance

Instance **LiverDiseaseDB**, database `neo4j`, Neo4j 2026.09.0 Enterprise.
