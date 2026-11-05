# Knowledge-graph drug repurposing for liver fibrosis progression

Code and results for Paper 2: knowledge-graph-based drug repurposing for liver fibrosis progression. The study ranks approved and clinical-stage drugs that reach stage-associated liver fibrosis genes through the knowledge graph, keeping only drugs whose action opposes the direction of the disease-associated change.

**To rebuild from the public Paper 1 graph, start with [docs/REBUILD.md](docs/REBUILD.md).** All required graph-restoration code is in this repository; KEGG/CTD source inputs remain local.

**For the analysis steps, see [Paper2_Results/README.md](Paper2_Results/README.md).** It lists every step (graph preparation, Tasks 1–4, the direction filter, the metapath B/BB analysis and the final ranking) with its numbers, code and output files.

## Built on Paper 1

This work uses the progression genes and the knowledge graph from the companion study (Paper 1):

- **Paper 1 repository:** https://github.com/nissankararaoThanvitha/liver_graph
- **Knowledge graph:** GitHub release [graph-export-2026-11-05](https://github.com/nissankararaoThanvitha/liver_graph/releases/tag/graph-export-2026-11-05) of that repository (126,480 nodes and 35,152,951 relationships in the starting full graph; the public release leaves out the KEGG pathways, as explained there).

The eight Paper 1 result files that the Paper 2 code reads are copied here at the same paths, unchanged, so the code runs as it is:

| File | Used by |
|---|---|
| `Paper1_HighConfidence/00_all_progression_genes.csv` | Task 1 |
| `Paper1_HighConfidence/01_master_table.csv` | Task 1 |
| `Paper1_HighConfidence/11_final_all_3845_genes.csv` | metapath runs (Tier 1–3 gene sets) |
| `Paper1_HighConfidence/11_FINAL_shortlist_25_genes.csv` | Task 1, metapath runs (23-gene shortlist; renamed `11_FINAL_shortlist_23_genes.csv` in the Paper 1 repository, same content) |
| `Paper1_HighConfidence/10_loso_all_eligible_runs.csv` | Task 3 (14,794-gene background) |
| `Paper1_Results/Task1/task1_summary.csv` | Task 1 (count check) |
| `Paper1_Results/Task2/task2_gene_clusters.csv` | Task 1 (clusters C1–C5) |
| `Paper1_Results/Task4/task4_all_significant.csv` | Task 4 |

The current reconstruction restores the omitted KEGG content using `scripts/restore_kegg.py`, called by `scripts/rebuild_graph.py`. No Paper 1 checkout or OptimusKG parquet download is needed for this route.

## Layout

| Folder | Contents |
|---|---|
| `Paper2_Results/` | Results of every step, and the analysis code in `Paper2_Results/code/` |
| `scripts/` | Graph-preparation scripts for drug metadata recovered from OptimusKG and ChEMBL (steps 0.1–0.3), and the metapath runners used by `Paper2_Results/code/metapath_pipeline.py` |
| `Paper1_HighConfidence/`, `Paper1_Results/` | The Paper 1 input files listed above |

Use the cleaned Paper 1 release as the starting graph; Paper 2 adds its drug evidence to a separate database. The schema-correct loader and legacy-conversion commands are documented under “Final drug schema” in `Paper2_Results/README.md`.

## Running

Python 3.12 with the packages in `requirements.txt`. Run scripts from the repository root. Most steps query the knowledge graph in Neo4j (load it from the Paper 1 release); set `NEO4J_PASSWORD` (and `NEO4J_URI` / `NEO4J_USERNAME` if they differ from `bolt://localhost:7687` / `neo4j`), or let the script prompt for it. On Windows, `Paper2_Results/code/run_tasks.ps1 <script> ...` runs several steps with one password prompt.

Raw downloads (OptimusKG parquet files, CTD) are not included; the Paper 1 repository's `docs/provenance/DATA_SOURCES.md` says where each comes from. The full BB route CSVs (up to 320 MB) are not in git; `python Paper2_Results/code/metapath_pipeline.py Run1 Run2 Run3 Run4` rebuilds them, and the same routes are in the `*_BB.xlsx` workbooks.

## Data terms

Results derived from the knowledge graph keep the terms of its sources. The largest part comes from OptimusKG (CC BY-NC-SA 4.0); drug data come from ChEMBL (CC BY-SA 3.0), DGIdb, GtoPdb (CC BY-SA 4.0) and Open Targets (CC0). Use them non-commercially, credit the sources, and share derived data under the same terms. See the Paper 1 repository's `docs/graph_export/DATA_TERMS.md`.
