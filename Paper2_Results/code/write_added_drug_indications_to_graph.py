"""Write the ChEMBL indications of the drugs added on 21 Oct 2026 into Neo4j as TREATS links.

The 519 Drug nodes added from ChEMBL, DGIdb, GtoPdb and Open Targets have drug-gene links
but no drug-disease links. task4_added_drug_indications.py fetched their indications from
ChEMBL (Paper2_Results/Task4/added_drug_indications_chembl.csv). This script adds them to the
graph so later steps (metapaths, rule "drug has no TREATS link to NAFLD/NASH/cirrhosis") see them.

Meaning: OptimusKG's TREATS links are drug indications at any clinical phase (Open Targets
clinical-trial indications, phase 0.5-4, and DrugCentral approved uses). ChEMBL indications are
the same kind of record, so they are written as TREATS, one link per drug-disease pair, with:
    source = 'ChEMBL indication', layer = 'ChEMBL indication (added drugs)',
    max_phase_for_ind (highest phase for the pair), efo_terms, mesh_headings, drugind_ids,
    chembl_release, chembl_release_date, retrieved, added
Only diseases already in the graph are used (matched on Disease.node_id, e.g. EFO_0000274); no
node is created. Pairs that already have a TREATS link are skipped. Nothing existing is changed.

Default is a DRY RUN that changes nothing.
    python Paper2_Results/code/write_added_drug_indications_to_graph.py           # dry run
    python Paper2_Results/code/write_added_drug_indications_to_graph.py --apply   # write
    python Paper2_Results/code/write_added_drug_indications_to_graph.py --undo    # remove what this script added
Password from NEO4J_PASSWORD (prompted if unset).
"""
import base64
import getpass
import json
import os
import sys
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "Paper2_Results/Task4/added_drug_indications_chembl.csv"
OUT = ROOT / "Paper2_Results/Task4/graph_write"
OUT.mkdir(parents=True, exist_ok=True)
LAYER = "ChEMBL indication (added drugs)"
PW = os.environ.get("NEO4J_PASSWORD") or getpass.getpass("Neo4j password: ")
MODE = "apply" if "--apply" in sys.argv else "undo" if "--undo" in sys.argv else "dry run"


def cy(stmt, params=None, write=False):
    h = {"Authorization": "Basic " + base64.b64encode(f"neo4j:{PW}".encode()).decode(),
         "Content-Type": "application/json", "Accept": "application/json"}
    body = json.dumps({"statement": stmt, "parameters": params or {},
                       "accessMode": "Write" if write else "Read"}).encode()
    r = json.load(urllib.request.urlopen(urllib.request.Request(
        "http://127.0.0.1:7474/db/neo4j/query/v2", data=body, headers=h, method="POST"), timeout=900))
    if r.get("errors"):
        raise RuntimeError(r["errors"])
    return pd.DataFrame(r["data"]["values"], columns=r["data"]["fields"])


def snapshot():
    s = cy("""CALL () { MATCH (n) RETURN count(n) AS nodes }
              CALL () { MATCH ()-[r]->() RETURN count(r) AS rels }
              CALL () { MATCH ()-[r:TREATS]->() RETURN count(r) AS treats }
              CALL () { MATCH ()-[r:TREATS {layer: $layer}]->() RETURN count(r) AS treats_added }
              RETURN nodes, rels, treats, treats_added""", {"layer": LAYER})
    return {k: int(v) for k, v in s.iloc[0].items()}


if MODE == "undo":
    before = snapshot()
    n = cy("MATCH ()-[r:TREATS {layer: $layer}]->() DELETE r RETURN count(r) AS n", {"layer": LAYER}, write=True)
    log = {"date": date.today().isoformat(), "mode": "undo", "removed": int(n.n.iloc[0]),
           "before": before, "after": snapshot()}
    (OUT / "undo_log.json").write_text(json.dumps(log, indent=2))
    print(json.dumps(log, indent=2))
    sys.exit(0)

# ---------------------------------------------------------------- one link per drug-disease pair
ind = pd.read_csv(SRC)
pairs = (ind.groupby(["Drug_ID", "efo_id_graph"])
         .agg(Drug=("Drug", "first"), max_phase_for_ind=("max_phase_for_ind", "max"),
              efo_terms=("efo_term", lambda s: "; ".join(sorted(set(s.dropna())))),
              mesh_headings=("mesh_heading", lambda s: "; ".join(sorted(set(s.dropna())))),
              drugind_ids=("drugind_id", lambda s: ";".join(str(int(x)) for x in sorted(set(s.dropna())))),
              chembl_release=("chembl_release", "first"), chembl_release_date=("chembl_release_date", "first"),
              retrieved=("retrieved", "first"))
         .reset_index())

# which drugs / diseases exist, and which pairs already have TREATS
found = cy("""UNWIND $rows AS r
OPTIONAL MATCH (d:Drug {node_id: r.drug}) OPTIONAL MATCH (x:Disease {node_id: r.dis})
OPTIONAL MATCH (d)-[t:TREATS]->(x)
RETURN r.drug AS Drug_ID, r.dis AS efo_id_graph, d IS NOT NULL AS drug_in_graph,
       x IS NOT NULL AS disease_in_graph, count(t) AS existing_treats""",
           {"rows": [{"drug": a, "dis": b} for a, b in zip(pairs.Drug_ID, pairs.efo_id_graph)]})
pairs = pairs.merge(found, on=["Drug_ID", "efo_id_graph"], how="left")
pairs["status"] = "write"
pairs.loc[~pairs.disease_in_graph.astype(bool), "status"] = "skip: disease not in graph"
pairs.loc[~pairs.drug_in_graph.astype(bool), "status"] = "skip: drug not in graph"
pairs.loc[pairs.existing_treats.astype(int) > 0, "status"] = "skip: TREATS already present"
pairs.to_csv(OUT / "indication_pairs_plan.csv", index=False)
todo = pairs[pairs.status == "write"]
print(pairs.status.value_counts().to_string())
print(f"links to write: {len(todo):,} (drugs {todo.Drug_ID.nunique():,}, diseases {todo.efo_id_graph.nunique():,})")

if MODE == "dry run":
    print("DRY RUN: nothing written. Run with --apply to write.")
    sys.exit(0)

# ---------------------------------------------------------------- write
before = snapshot()
if before["treats_added"]:
    sys.exit(f"{before['treats_added']} links from this layer are already in the graph; run --undo first.")
today = date.today().isoformat()
rows = [{"drug": r.Drug_ID, "dis": r.efo_id_graph,
         "p": {"source": "ChEMBL indication", "layer": LAYER,
               "max_phase_for_ind": None if pd.isna(r.max_phase_for_ind) else float(r.max_phase_for_ind),
               "efo_terms": r.efo_terms, "mesh_headings": r.mesh_headings, "drugind_ids": r.drugind_ids,
               "chembl_release": r.chembl_release, "chembl_release_date": r.chembl_release_date,
               "retrieved": r.retrieved, "added": today}} for r in todo.itertuples()]
written = 0
for i in range(0, len(rows), 1000):
    n = cy("""UNWIND $rows AS r
MATCH (d:Drug {node_id: r.drug}) MATCH (x:Disease {node_id: r.dis})
CREATE (d)-[t:TREATS]->(x) SET t = r.p RETURN count(t) AS n""", {"rows": rows[i:i + 1000]}, write=True)
    written += int(n.n.iloc[0])
after = snapshot()
chk = cy("MATCH (d:Drug)-[t:TREATS {layer: $layer}]->(x:Disease) RETURN d.node_id AS d, x.node_id AS x",
         {"layer": LAYER})
log = {"date": today, "mode": "apply", "links_written": written, "planned": len(todo),
       "drugs": int(todo.Drug_ID.nunique()), "diseases": int(todo.efo_id_graph.nunique()),
       "skipped": pairs.status.value_counts().drop("write", errors="ignore").to_dict(),
       "before": before, "after": after,
       "links_match_plan": set(zip(chk.d, chk.x)) == set(zip(todo.Drug_ID, todo.efo_id_graph)),
       "relationship_types_unchanged": "TREATS already existed; no new relationship type",
       "chembl_release": str(todo.chembl_release.iloc[0]) if len(todo) else None}
(OUT / "write_log.json").write_text(json.dumps(log, indent=2))
(OUT / "undo_indications.cypher").write_text(f"MATCH ()-[r:TREATS {{layer: '{LAYER}'}}]->() DELETE r;\n")
print(json.dumps(log, indent=2))
