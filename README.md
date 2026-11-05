# Knowledge-graph drug repurposing for liver fibrosis progression

Code and results for Paper 2: knowledge-graph-based drug repurposing for liver fibrosis progression. The study ranks approved and clinical-stage drugs that reach stage-associated liver fibrosis genes through the knowledge graph, keeping only drugs whose action opposes the direction of the disease-associated change.

**To rebuild from the public Paper 1 graph, start with [docs/REBUILD.md](docs/REBUILD.md).** All required graph-restoration code is in this repository; KEGG/CTD source inputs remain local.

**For the analysis steps, see [Paper2_Results/README.md](Paper2_Results/README.md).** It maps the saved results to Tasks 1–4, direction filtering, B/BB metapaths and ranking. Its graph-preparation table is historical provenance; use `docs/REBUILD.md` for the current execution order.

## Built on Paper 1

This work uses the progression genes and the knowledge graph from the companion study (Paper 1):

- **Paper 1 repository:** https://github.com/nissankararaoThanvitha/liver_graph
- **Knowledge graph:** GitHub release [graph-export-2026-11-05](https://github.com/nissankararaoThanvitha/liver_graph/releases/tag/graph-export-2026-11-05) of that repository. The downloaded public graph contains **126,294 nodes and 35,140,166 relationships**. Restoring its omitted KEGG portion gives the full Paper 1 graph (126,480 nodes / 35,152,951 relationships); Paper 2 then adds its drug evidence.

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
| `scripts/` | Current release verification, CTD download, KEGG restoration, complete graph rebuild and analysis runners; historical recovery tools and optional legacy-schema conversion |
| `docs/REBUILD.md` | Current import → restore → enrich → verify → analyse instructions |
| `Paper1_HighConfidence/`, `Paper1_Results/` | The Paper 1 input files listed above |

Use the cleaned Paper 1 release as the starting graph and `scripts/rebuild_graph.py` for the complete Paper 2 additions in a separate database. The older component loaders and optional legacy conversion are explained in `Paper2_Results/README.md`; they are not extra steps after the complete rebuild.

## Running

Python 3.12 with the packages in `requirements.txt`. Run commands from the repository root. The current rebuild and analysis runners require an explicit Bolt URI and database name; the analysis runner also requires the matching HTTP address. Set `NEO4J_PASSWORD` and, if needed, `NEO4J_USERNAME` for that separate Paper 2 instance. Follow [docs/REBUILD.md](docs/REBUILD.md) for exact commands and the two required local pathway inputs.

The current saved-evidence rebuild needs no fresh drug API queries or OptimusKG parquet download. Older recovery/search scripts are retained for provenance or deliberate data updates; their historical raw files are not all bundled. Do not run the old metadata `--rollback` on a graph imported from Paper 1: it removes imported metadata globally.

Use `scripts/run_analysis.py --stage core`, then `--stage metapaths`, with the connection arguments in the guide. The older `run_tasks.ps1` is a component-script launcher, not the complete reconstruction workflow. Analysis commands intentionally regenerate outputs; use a separate checkout or retain existing results before rerunning. The full BB route CSVs are not tracked; the metapath runner regenerates them, and the saved `*_BB.xlsx` workbooks also contain the routes.

## Data terms

Results derived from the knowledge graph keep the terms of its sources. The largest part comes from OptimusKG (CC BY-NC-SA 4.0); drug data come from ChEMBL (CC BY-SA 3.0), DGIdb, GtoPdb (CC BY-SA 4.0) and Open Targets (CC0). Use them non-commercially, credit the sources, and share derived data under the same terms. See the Paper 1 repository's `docs/graph_export/DATA_TERMS.md`.
