"""Write the 297 external drug-action directions into Neo4j as NEW properties on existing
Drug-[:ACTS_ON]->Gene relationships (requested 21 Oct 2026).

Only adds properties whose names start with `external_`; never changes `action` or any other
existing property, never creates or deletes nodes or relationships.

  external_action_direction   'down' | 'up'   (direction the drug pushes the gene product)
  external_action_label       the labels found, e.g. 'Antagonist / Antagonist | inhibitor'
  external_action_sources     'GtoPdb;DGIdb' etc.
  external_action_release     source versions used
  external_action_added       date written

Steps (default is a DRY RUN that changes nothing):
  1. match each (drug ChEMBL ID, gene Ensembl ID) to exactly one ACTS_ON relationship
  2. --apply: save a backup of the 297 relationships' current properties and an undo Cypher
     file, write the properties in one transaction, then verify

    python Paper2_Results/code/write_external_directions_to_graph.py            # dry run
    python Paper2_Results/code/write_external_directions_to_graph.py --apply    # write
    python Paper2_Results/code/write_external_directions_to_graph.py --undo     # remove them again
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
RES = ROOT / "Paper2_Results"
OUT = RES / "External_action_check" / "graph_write"
OUT.mkdir(parents=True, exist_ok=True)
PW = os.environ.get("NEO4J_PASSWORD") or getpass.getpass("Neo4j password: ")
RELEASE = "Guide to Pharmacology 2026.3.1; DGIdb 5 (API, Oct 2026); ChEMBL REST API (Oct 2026)"
PROPS = ["external_action_direction", "external_action_label", "external_action_sources",
         "external_action_release", "external_action_added"]


def cy(stmt, params=None, write=False):
    h = {"Authorization": "Basic " + base64.b64encode(f"neo4j:{PW}".encode()).decode(),
         "Content-Type": "application/json", "Accept": "application/json"}
    body = json.dumps({"statement": stmt, "parameters": params or {}, "accessMode": "Write" if write else "Read"}).encode()
    r = json.load(urllib.request.urlopen(urllib.request.Request(
        "http://127.0.0.1:7474/db/neo4j/query/v2", data=body, headers=h, method="POST"), timeout=600))
    if r.get("errors"):
        raise RuntimeError(r["errors"])
    return pd.DataFrame(r["data"]["values"], columns=r["data"]["fields"])


x = pd.read_csv(RES / "Task7_external_layer/external_layer_links.csv")
assert len(x) == 297 and not x.duplicated(["Drug_ID", "Gene_ID"]).any()
rows = [{"drug": d, "gene": g, "direction": p, "label": str(l)[:300], "sources": s.replace("external: ", "")}
        for d, g, p, l, s in zip(x.Drug_ID, x.Gene_ID, x.push, x.action_label, x.direction_source)]

if "--undo" in sys.argv:
    res = cy("MATCH ()-[r:ACTS_ON]->() WHERE r.external_action_direction IS NOT NULL "
             "REMOVE " + ", ".join(f"r.{p}" for p in PROPS) + " RETURN count(r) AS removed", write=True)
    print("removed external properties from", int(res.removed.iloc[0]), "relationships")
    sys.exit(0)

# 1. exact matching
m = cy("""UNWIND $rows AS row
MATCH (d:Drug {node_id: row.drug})-[r:ACTS_ON]->(g:Gene {ensembl_id: row.gene})
RETURN row.drug AS drug, row.gene AS gene, count(r) AS n, collect(r.action)[0] AS action,
       collect(r.external_action_direction)[0] AS already""", {"rows": rows})
print(f"pairs matched: {len(m)} / {len(rows)}; exactly one relationship each: {(m.n == 1).all()}; "
      f"already carry external properties: {int(m.already.notna().sum())}")
print("current graph action labels of these links:", m.action.value_counts().to_dict())
assert len(m) == 297 and (m.n == 1).all(), "every pair must match exactly one ACTS_ON relationship"
before = cy("MATCH ()-[r:ACTS_ON]->() RETURN count(r) AS rels, count(DISTINCT r.action) AS labels")
print("ACTS_ON relationships in graph:", int(before.rels.iloc[0]))
if "--apply" not in sys.argv:
    print("DRY RUN: nothing written. Re-run with --apply to write.")
    sys.exit(0)

# 2. backup + undo file
bk = cy("""UNWIND $rows AS row
MATCH (d:Drug {node_id: row.drug})-[r:ACTS_ON]->(g:Gene {ensembl_id: row.gene})
RETURN row.drug AS drug, row.gene AS gene, properties(r) AS properties""", {"rows": rows})
bk["properties"] = bk.properties.map(json.dumps)
bk.to_csv(OUT / "backup_ACTS_ON_properties_before_write.csv", index=False)
(OUT / "undo_external_properties.cypher").write_text(
    "// Removes ONLY the properties added by write_external_directions_to_graph.py\n"
    "MATCH ()-[r:ACTS_ON]->() WHERE r.external_action_direction IS NOT NULL\n"
    "REMOVE " + ", ".join(f"r.{p}" for p in PROPS) + "\nRETURN count(r) AS removed;\n", encoding="utf-8")

# 3. write (one transaction)
today = date.today().isoformat()
w = cy("""UNWIND $rows AS row
MATCH (d:Drug {node_id: row.drug})-[r:ACTS_ON]->(g:Gene {ensembl_id: row.gene})
SET r.external_action_direction = row.direction,
    r.external_action_label = row.label,
    r.external_action_sources = row.sources,
    r.external_action_release = $release,
    r.external_action_added = $today
RETURN count(r) AS written""", {"rows": rows, "release": RELEASE, "today": today}, write=True)
print("written:", int(w.written.iloc[0]))

# 4. verify
after = cy("MATCH ()-[r:ACTS_ON]->() RETURN count(r) AS rels, count(DISTINCT r.action) AS labels")
v = cy("""MATCH ()-[r:ACTS_ON]->() WHERE r.external_action_direction IS NOT NULL
RETURN r.external_action_direction AS dir, count(*) AS n""")
chk = cy("""UNWIND $rows AS row
MATCH (d:Drug {node_id: row.drug})-[r:ACTS_ON]->(g:Gene {ensembl_id: row.gene})
RETURN row.drug AS drug, row.gene AS gene, r.action AS action, r.external_action_direction AS dir""", {"rows": rows})
orig = dict(zip(zip(m.drug, m.gene), m.action))
log = {
    "date": today, "relationships_written": int(w.written.iloc[0]),
    "ACTS_ON_before": int(before.rels.iloc[0]), "ACTS_ON_after": int(after.rels.iloc[0]),
    "action_labels_unchanged": all(orig[(a, b)] == c for a, b, c in zip(chk.drug, chk.gene, chk.action)),
    "directions_match_file": all(d == r["direction"] for d, r in zip(chk.set_index(["drug", "gene"]).loc[[(r["drug"], r["gene"]) for r in rows]].dir, rows)),
    "external_by_direction": dict(zip(v.dir, v.n.astype(int))), "properties_added": PROPS, "release": RELEASE,
    "backup": "backup_ACTS_ON_properties_before_write.csv", "undo": "undo_external_properties.cypher",
}
(OUT / "write_log.json").write_text(json.dumps(log, indent=2), encoding="utf-8")
print(json.dumps(log, indent=2))
assert log["relationships_written"] == 297 and log["ACTS_ON_before"] == log["ACTS_ON_after"]
assert log["action_labels_unchanged"] and log["directions_match_file"]
