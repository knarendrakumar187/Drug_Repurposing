# Rebuild Paper 2 from the public Paper 1 release

This is the current workflow. It uses this repository's saved drug-evidence tables, not changing database APIs. It does not require a Paper 1 checkout or the original OptimusKG parquet files. Work in a **separate Neo4j instance/database for Paper 2**. Do not modify the published Paper 1 graph.

## 1. Install and obtain the starting graph

Use Python 3.12 and install `requirements.txt`. Neo4j 5.26 or newer is required. Download the CSV ZIP or database dump from:

https://github.com/nissankararaoThanvitha/liver_graph/releases/tag/graph-export-2026-11-05

The public graph has **126,294 nodes and 35,140,166 relationships**. It includes original OptimusKG drugs, metadata and evidence, Hallmark, Reactome, expression data and the Paper 1 gene associations. It excludes KEGG and Paper 2's additional drug records.

The quickest route is the release's `liver_kg_2026-11-05.dump`: follow its `IMPORT.md` to load it into a new, separate instance. For CSV import, verify and extract the ZIP with the tool in this repository:

```bash
python scripts/prepare_graph_release.py --archive /downloads/liver_kg_graph_export_2026-11-05.zip --out /private/paper1_csv
```

The tool pins the published ZIP's SHA-256, verifies all 25 data-file hashes, rejects unsafe ZIP paths, and requires a new output directory. It does not alter a database. With the separate Paper 2 server stopped, import using that instance's `neo4j-admin` from the extracted folder:

```bash
args=()
for f in nodes/*.csv.gz; do args+=("--nodes=$f"); done
for f in relationships/*.csv.gz; do args+=("--relationships=$f"); done
HEAP_SIZE=2G /path/to/paper2/bin/neo4j-admin database import full neo4j \
  --id-type=string --array-delimiter=";" --multiline-fields=true \
  --ignore-empty-strings=false --bad-tolerance=0 --max-off-heap-memory=1G "${args[@]}"
```

No overwrite option is supplied: choose an empty database. On Community edition use the default database `neo4j` in the separate instance. The enclosed release `IMPORT.md` also explains the PowerShell route. Start the separate server and run the enclosed `constraints.cypher` in the imported database; these indexes are needed for efficient matching. The release dump already includes the tested indexes/constraints; running the file again is harmless.

## 2. Obtain the two local pathway inputs

1. Your authorised **MSigDB 2024.1 human-symbol KEGG legacy GMT**, named `c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt`. Obtain it from the MSigDB 2024.1 release directory, subject to the provider's access and use terms:
   https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/
2. **CTD disease–pathway CSV**, `CTD_diseases_pathways.csv.gz`, from CTD's official reports. The helper downloads it to a new local file:

```bash
python scripts/fetch_rebuild_inputs.py --out /private/inputs/CTD_diseases_pathways.csv.gz
```

CTD source: https://ctdbase.org/reports/CTD_diseases_pathways.csv.gz
CTD terms: https://ctdbase.org/about/legal.jsp

The tested input checksums are:

- KEGG GMT: `06a5d7f82d90a0a6a2e464579d049fc42d279fc23fa9f79a956e56c4a442447f`
- CTD CSV gzip: `c48cf9bbc70cf601f08573233b7c5de92af8303483a4c0fdc69fd917bb44a067`

These identify retained files, not official release identifiers. Retain matching inputs for exact reproduction; future source downloads may differ.

Keep both inputs outside version control (or under the ignored `data/` directory). **Do not upload KEGG data or generated KEGG graph exports.** The rebuild reads the GMT directly; it does not need the Hallmark GMT because Hallmark is already in the release. It matches CTD against the committed Paper 1 enrichment output. Downloaded CTD data may change; retain the file used. The preflight refuses inputs that do not reproduce 186 KEGG pathways, 12,724 gene memberships and 61 KEGG–NAFLD links. Input SHA-256 values are recorded in the checkpoint, but are not currently enforced against the reference hashes above. Matching counts alone do not prove identical inputs; compare both hashes before applying an exact-reproduction run.

## 3. Plan and rebuild all additions

Run from this repository's root. Set `NEO4J_PASSWORD` for the separate Paper 2 database and optionally `NEO4J_USERNAME` (default `neo4j`). Never put the password in source or a committed configuration. The URI/database are required explicitly, so no default can silently select Paper 1.

First run the read-only plan:

```bash
python scripts/rebuild_graph.py --uri bolt://localhost:7687 --database neo4j \
  --kegg-gmt /private/inputs/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt \
  --ctd /private/inputs/CTD_diseases_pathways.csv.gz
```

Replace the example address with the separate Paper 2 instance's actual address. Stop other writers, then apply with a **new** checkpoint path:

```bash
python scripts/rebuild_graph.py --uri bolt://localhost:7687 --database neo4j \
  --kegg-gmt /private/inputs/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt \
  --ctd /private/inputs/CTD_diseases_pathways.csv.gz \
  --apply --backup /private/backups/paper2-before-additions.json
```

The checkpoint contains pre-existing drug evidence, graph counts and input hashes, not credentials. It is **not a full database backup**: keep the original imported dump/ZIP for full restoration. All additions and checks run in one transaction. A failure rolls back the entire rebuild. A repeat against a complete graph verifies it without writing.

| Reconstructed content | Input / implementation |
|---|---|
| 186 KEGG pathways; 12,724 gene memberships | Local GMT; `scripts/restore_kegg.py` |
| 61 KEGG–NAFLD links; existing 28 Reactome links retained | Local CTD + `Paper1_Results/Task4/task4_all_significant.csv` |
| 519 Drug nodes; 5,216 ACTS_ON evidence links | Saved curated table, alias rows and source versions in `Paper2_Results/New_drug_gene_links/`; schema-correct loader |
| 297 action annotations on original links | `Paper2_Results/Task7_external_layer/external_layer_links.csv` |
| 2,087 TREATS indication links | Saved accepted pairs in `Paper2_Results/Task4/graph_write/indication_pairs_plan.csv` |
| ChEMBL gap fill on 28 original drugs: 28 phase values and 16 names | `Paper2_Results/Recovery2_ChEMBL/graph_gap_properties.csv`; OptimusKG values retained |

The gap-fill CSV preserves values from the verified earlier graph, with field provenance; it does not claim a newly fetched or inferred ChEMBL version. See its directory README. Drug attribution remains in properties. No `ExternalDrug` labels or `EXTERNAL_ACTS_ON` relationships are created.

Expected completed graph: **126,999 nodes; 35,160,254 relationships; 12,544 Drug nodes; 25,890 ACTS_ON; 59,688 TREATS; 89 IMPLICATED_IN**. Original expression, gene associations, Reactome/Hallmark and non-drug knowledge remain in place. The tool checks counts, restored KEGG contents, drug evidence, metadata and directions before committing.

For an older Paper 2 graph, use `scripts/unify_drug_schema.py` to migrate its legacy labels/types; the new full rebuild accepts only the clean public import or a complete graph, rather than guessing how to repair partial prior writes. Do not replay the historical steps 0.1–0.8 over the new release or invoke the old global metadata rollback.

## 4. Run downstream analyses

Only after the complete rebuild verifies successfully, run Tasks 1–4 and direction filtering:

```bash
python scripts/run_analysis.py --uri bolt://localhost:7687 --http http://localhost:7474 --database neo4j --stage core
```

Then run the B/BB metapaths and ranking:

```bash
python scripts/run_analysis.py --uri bolt://localhost:7687 --http http://localhost:7474 --database neo4j --stage metapaths
```

Use the actual Bolt **and HTTP addresses of the same Paper 2 instance**. `--stage all` runs both. These analysis commands intentionally regenerate outputs under `Paper2_Results/`; preserve existing results or use a separate checkout before recomputing. No scientific filters, gene sets or ranking definitions are changed by the rebuild tools. Large metapath runs require time, memory and disk space.

## 5. Checks and limitations

Ordinary tests: `python -m unittest discover -s tests -v`. Database integration tests require the disposable public-release import at `TEST_NEO4J_URI=bolt://127.0.0.1:17689`, plus `TEST_KEGG_GMT` and `TEST_CTD` pointing to local source files; never point them at a production graph. The full test intentionally leaves that disposable database rebuilt after checking rollback and repeat verification. Historical raw API caches are not required for the current saved-evidence rebuild. Re-querying drug databases is a separate data update and can change outputs; this guide deliberately replays retained study evidence.
