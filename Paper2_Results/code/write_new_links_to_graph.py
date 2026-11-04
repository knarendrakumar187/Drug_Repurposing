"""Graph preparation, step 0.4: write the new drug–gene links (from ChEMBL, DGIdb, GtoPdb, Open Targets) into Neo4j as a
separate, clearly labelled layer (requested 21 Oct 2026).

What is added (nothing existing is changed or deleted):
  * relationship type  (:Drug|ExternalDrug)-[:EXTERNAL_ACTS_ON]->(:Gene)   one per link
      properties: action_labels, push, sources, source_detail, evidence, n_sources, drug_stage,
                  mechanism, release, added, layer = 'external (Idea 4)'
  * new nodes for drugs that were not in the graph, labelled :Drug:ExternalDrug (normal Drug nodes;
    the extra ExternalDrug label and source = 'external (Idea 4)' mark them as outside data).
    They have no ACTS_ON links, so every graph-only result (which starts from ACTS_ON) is unchanged;
    only the total Drug node count rises (12,025 -> 12,544).

Graph-only queries use ACTS_ON and :Drug, so all graph-only results stay identical.

Default is a DRY RUN that changes nothing.
    python Paper2_Results/code/write_new_links_to_graph.py            # dry run
    python Paper2_Results/code/write_new_links_to_graph.py --apply    # write
    python Paper2_Results/code/write_new_links_to_graph.py --undo     # remove everything this script added
    python Paper2_Results/code/write_new_links_to_graph.py --as-drug  # give ExternalDrug-only nodes the Drug label too
Needs NEO4J_PASSWORD in the environment.
"""
import base64
import json
import os
import sys
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "Paper2_Results"
NL = RES / "New_drug_gene_links"
OUT = NL / "graph_write"
OUT.mkdir(parents=True, exist_ok=True)
PW = os.environ["NEO4J_PASSWORD"]
REL, XLAB = "EXTERNAL_ACTS_ON", "ExternalDrug"


def cy(stmt, params=None, write=False):
    h = {"Authorization": "Basic " + base64.b64encode(f"neo4j:{PW}".encode()).decode(),
         "Content-Type": "application/json", "Accept": "application/json"}
    body = json.dumps({"statement": stmt, "parameters": params or {}, "accessMode": "Write" if write else "Read"}).encode()
    r = json.load(urllib.request.urlopen(urllib.request.Request(
        "http://127.0.0.1:7474/db/neo4j/query/v2", data=body, headers=h, method="POST"), timeout=900))
    if r.get("errors"):
        raise RuntimeError(r["errors"])
    return pd.DataFrame(r["data"]["values"], columns=r["data"]["fields"])


def snapshot():
    s = cy("""CALL () { MATCH (d:Drug) RETURN count(d) AS drugs }
              CALL () { MATCH ()-[r:ACTS_ON]->() RETURN count(r) AS acts_on }
              CALL () { MATCH (n) RETURN count(n) AS nodes }
              CALL () { MATCH ()-[r]->() RETURN count(r) AS rels }
              CALL () { MATCH (x:ExternalDrug) RETURN count(x) AS external_drugs }
              CALL () { MATCH ()-[r:EXTERNAL_ACTS_ON]->() RETURN count(r) AS external_links }
              RETURN drugs, acts_on, nodes, rels, external_drugs, external_links""")
    return {k: int(v) for k, v in s.iloc[0].items()}


if "--undo" in sys.argv:
    r1 = cy(f"MATCH ()-[r:{REL}]->() DELETE r RETURN count(r) AS removed", write=True)
    r2 = cy(f"MATCH (x:{XLAB}) WHERE x.source = 'external (Idea 4)' AND NOT (x)--() DELETE x RETURN count(x) AS removed", write=True)
    print("removed links:", int(r1.removed.iloc[0]), "| removed ExternalDrug nodes:", int(r2.removed.iloc[0]))
    print(snapshot())
    sys.exit(0)

if "--as-drug" in sys.argv:
    b = snapshot()
    r = cy(f"MATCH (x:{XLAB}) WHERE NOT x:Drug AND x.source = 'external (Idea 4)' SET x:Drug RETURN count(x) AS n", write=True)
    a = snapshot()
    log = {"date": date.today().isoformat(), "relabelled_as_Drug": int(r.n.iloc[0]), "before": b, "after": a}
    (OUT / "relabel_log.json").write_text(json.dumps(log, indent=2), encoding="utf-8")
    print(json.dumps(log, indent=2))
    assert a["acts_on"] == b["acts_on"] and a["drugs"] == b["drugs"] + log["relabelled_as_Drug"]
    sys.exit(0)

# ------------------------------------------------------------------ the links
m = pd.read_csv(NL / "new_links_medicines_curated.csv")
assert (m.evidence == "curated drug-target").all() and not m.duplicated(["Drug_ID", "Gene_ID"]).any()
raw = pd.read_csv(NL / "raw_rows_all_sources.csv.gz", low_memory=False, usecols=["Drug_key", "chembl"])
aliases = raw.dropna().groupby("Drug_key").chembl.agg(lambda s: sorted(set(s))).to_dict()   # parent → IDs seen (salts)
nj = json.loads((NL / "numbers.json").read_text(encoding="utf-8"))
release = "; ".join(f"{k} {v}" for k, v in nj["release"].items())

# which Drug node (if any) does each drug map to? parent ID first, then its salt / child IDs
drug_ids = sorted(m.Drug_ID.unique())
cands = [{"key": d, "ids": [d] + [a for a in aliases.get(d, []) if a != d]} for d in drug_ids]
hit = cy("""UNWIND $c AS c
UNWIND range(0, size(c.ids) - 1) AS i
MATCH (d:Drug {node_id: c.ids[i]}) WHERE NOT d:ExternalDrug
WITH c.key AS key, i, d.node_id AS node ORDER BY i
RETURN key, collect(node)[0] AS node""", {"c": cands})
node_of = dict(zip(hit.key, hit.node))
m["graph_node"] = m.Drug_ID.map(node_of)
m["node_kind"] = m.graph_node.notna().map({True: "Drug (existing)", False: "ExternalDrug (new node)"})
genes_ok = cy("UNWIND $g AS g MATCH (x:Gene {ensembl_id: g}) RETURN count(DISTINCT x) AS n", {"g": sorted(m.Gene_ID.unique())})
dup = cy(f"""UNWIND $rows AS r MATCH (d:Drug {{node_id: r.node}})-[x:ACTS_ON]->(:Gene {{ensembl_id: r.gene}}) RETURN count(x) AS n""",
         {"rows": [{"node": a, "gene": b} for a, b in zip(m.graph_node, m.Gene_ID) if isinstance(a, str)]})
before = snapshot()
print(f"links to write: {len(m):,}; drugs: {m.Drug_ID.nunique():,} "
      f"(existing Drug nodes: {m.loc[m.graph_node.notna(), 'Drug_ID'].nunique():,}, new ExternalDrug nodes: {m.loc[m.graph_node.isna(), 'Drug_ID'].nunique():,})")
print(f"genes found in graph: {int(genes_ok.n.iloc[0]):,} / {m.Gene_ID.nunique():,}; links that already exist as ACTS_ON: {int(dup.n.iloc[0])}")
print("graph now:", before)
assert int(genes_ok.n.iloc[0]) == m.Gene_ID.nunique() and int(dup.n.iloc[0]) == 0
if before["external_links"] or before["external_drugs"]:
    sys.exit("Graph already has EXTERNAL_ACTS_ON / ExternalDrug — run --undo first.")
if "--apply" not in sys.argv:
    print("DRY RUN: nothing written. Re-run with --apply to write.")
    sys.exit(0)

# ------------------------------------------------------------------ write
today = date.today().isoformat()
(OUT / "undo_new_links.cypher").write_text(
    f"// Removes ONLY what write_new_links_to_graph.py added\nMATCH ()-[r:{REL}]->() DELETE r;\n"
    f"MATCH (x:{XLAB}) WHERE NOT x:Drug AND NOT (x)--() DELETE x;\n", encoding="utf-8")
newd = m[m.graph_node.isna()].drop_duplicates("Drug_ID")
c1 = cy(f"""UNWIND $d AS d CREATE (x:Drug:{XLAB} {{node_id: d.id, name: d.name, drug_stage: d.stage, max_phase: d.phase,
          molecule_type: d.mtype, source: 'external (Idea 4)', release: $release, added: $today}}) RETURN count(x) AS n""",
        {"d": [{"id": r.Drug_ID, "name": r.Drug, "stage": r.drug_stage, "phase": None if pd.isna(r.max_phase) else float(r.max_phase),
                "mtype": None if pd.isna(r.molecule_type) else str(r.molecule_type)} for r in newd.itertuples()],
         "release": release, "today": today}, write=True)


def props(r):
    return {"action_labels": "" if pd.isna(r.action_labels) else str(r.action_labels), "push": r.push, "sources": r.sources,
            "source_detail": r.source_detail, "evidence": r.evidence, "n_sources": int(r.n_sources), "drug_stage": r.drug_stage,
            "mechanism": "" if pd.isna(r.mechanism_text) else str(r.mechanism_text)}


rows = [{"drug": r.graph_node if isinstance(r.graph_node, str) else r.Drug_ID, "existing": isinstance(r.graph_node, str),
         "gene": r.Gene_ID, "p": props(r)} for r in m.itertuples()]
written = 0
for i in range(0, len(rows), 1000):
    w = cy(f"""UNWIND $rows AS r
MATCH (g:Gene {{ensembl_id: r.gene}})
CALL (r) {{ WITH r WHERE r.existing MATCH (d:Drug {{node_id: r.drug}}) RETURN d
           UNION WITH r WHERE NOT r.existing MATCH (d:{XLAB} {{node_id: r.drug}}) RETURN d }}
CREATE (d)-[x:{REL}]->(g)
SET x += r.p, x.layer = 'external (Idea 4)', x.release = $release, x.added = $today
RETURN count(x) AS n""", {"rows": rows[i:i + 1000], "release": release, "today": today}, write=True)
    written += int(w.n.iloc[0])

# ------------------------------------------------------------------ verify
after = snapshot()
chk = cy(f"MATCH (d)-[r:{REL}]->(g:Gene) RETURN d.node_id AS d, g.ensembl_id AS g, r.push AS push")
want = {(r["drug"], r["gene"], r["p"]["push"]) for r in rows}
log = {"date": today, "links_written": written, "external_drug_nodes_created": int(c1.n.iloc[0]),
       "before": before, "after": after,
       "ACTS_ON_unchanged": before["acts_on"] == after["acts_on"], "Drug_count_change": after["drugs"] - before["drugs"],
       "links_match_file": set(zip(chk.d, chk.g, chk.push)) == want,
       "relationship_type": REL, "new_node_label": XLAB, "release": release,
       "undo": "undo_new_links.cypher  (or: python Paper2_Results/code/write_new_links_to_graph.py --undo)"}
(OUT / "write_log.json").write_text(json.dumps(log, indent=2), encoding="utf-8")
m.drop(columns=["graph_node"]).to_csv(OUT / "links_written.csv", index=False)
print(json.dumps(log, indent=2))
assert written == len(m) and log["ACTS_ON_unchanged"] and log["Drug_count_change"] == len(newd) and log["links_match_file"]
assert after["external_drugs"] == len(newd) and after["nodes"] == before["nodes"] + len(newd)
