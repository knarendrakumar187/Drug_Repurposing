"""Paper 2 Task 2: extract every drug-gene (ACTS_ON) edge for the progression genes.

Source: the live Neo4j liver_kg graph only, read-only (HTTP query API). There is
no CSV fallback: if Neo4j cannot be reached the script stops, so every number
comes from the present graph.

The graph holds two kinds of ACTS_ON edges:
  * OptimusKG edges (action in `action`, per-source assertions, sources_direct, ...)
  * edges added on 21 Oct 2026 from ChEMBL, DGIdb, GtoPdb and Open Targets
    (layer = 'external (Idea 4)'; action labels in `action_labels`, direction in
    `push`, databases in `sources`; see Paper2_Results/New_drug_gene_links/).
Both are read and their own fields are kept. Where an OptimusKG edge gives no
direction (TARGET, ENZYME, TRANSPORTER, CARRIER or blank), the direction stored
on the same edge from GtoPdb/DGIdb/ChEMBL (`external_action_direction`, 297
edges; Paper2_Results/External_action_check/) is used. It never overrides a
direction OptimusKG already gives.

Password is read from the NEO4J_PASSWORD environment variable (prompted if unset).
Run from the repository root:

    python Paper2_Results/code/task2_drug_gene_edges.py
"""
import base64
import getpass
import json
import os
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Paper2_Results" / "Task2"
OUT.mkdir(parents=True, exist_ok=True)

GENES = ROOT / "Paper2_Results/Task1/Paper2_progression_genes.csv"
NEW_LINKS_LOG = ROOT / "Paper2_Results/New_drug_gene_links/graph_write/write_log.json"

NEO4J_URL = os.environ.get("NEO4J_HTTP", "http://127.0.0.1:7474")
NEO4J_DB = os.environ.get("NEO4J_DATABASE", "neo4j")
NEO4J_USER = os.environ.get("NEO4J_USERNAME", "neo4j")
EXTERNAL_LAYER = "external (Idea 4)"

# Drug action -> effect on the gene product ("switch down" / "switch up").
DOWN = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE_MODULATOR", "NEGATIVE_ALLOSTERIC_MODULATOR",
        "ALLOSTERIC_ANTAGONIST", "INVERSE_AGONIST", "DEGRADER"}
UP = {"AGONIST", "PARTIAL_AGONIST", "ACTIVATOR", "POSITIVE_MODULATOR", "POSITIVE_ALLOSTERIC_MODULATOR",
      "OPENER", "STABILISER"}
MIXED = {"MODULATOR", "BINDING_AGENT", "SUBSTRATE", "RELEASING_AGENT"}
UNKNOWN = {"TARGET", "ENZYME", "TRANSPORTER", "CARRIER"}
# Direction recorded on the added edges (`push`), as written by write_new_links_to_graph.py.
PUSH = {"down": "down", "up": "up", "conflicting": "mixed", "unknown": "unknown"}


def action_direction(action):
    a = str(action).upper()
    if a in DOWN:
        return "down"
    if a in UP:
        return "up"
    if a in MIXED:
        return "mixed"
    return "unknown"


def edge_direction(action, assertions):
    """Direction of an OptimusKG edge's consolidated action, demoted to 'mixed' when the
    individual source assertions point both down and up."""
    d = action_direction(action)
    if isinstance(assertions, str) and assertions:
        dirs = {action_direction(a.split(":", 1)[-1]) for a in assertions.split(";")}
        if {"down", "up"} <= dirs:
            return "mixed"
    return d


def blank(x):
    return x is None or (isinstance(x, float) and pd.isna(x)) or str(x).strip() == ""


# ------------------------------------------------------------------ inputs
genes = pd.read_csv(GENES)
ids = genes.ensembl_id.tolist()


# ------------------------------------------------------------------ Neo4j (no fallback)
def neo4j_query(statement, params, password):
    headers = {"Authorization": "Basic " + base64.b64encode(f"{NEO4J_USER}:{password}".encode()).decode(),
               "Content-Type": "application/json", "Accept": "application/json"}
    body = json.dumps({"statement": statement, "parameters": params, "accessMode": "Read"}).encode()
    req = urllib.request.Request(f"{NEO4J_URL}/db/{NEO4J_DB}/query/v2", data=body, headers=headers, method="POST")
    res = json.load(urllib.request.urlopen(req, timeout=300))
    if res.get("errors"):
        raise RuntimeError(res["errors"])
    return pd.DataFrame(res["data"]["values"], columns=res["data"]["fields"])


QUERY = """
MATCH (d:Drug)-[r:ACTS_ON]->(g:Gene)
WHERE g.ensembl_id IN $ids
RETURN d.node_id AS Drug_ID, d.name AS Drug, d.molecule_type AS molecule_type,
       d.max_phase AS max_phase, d.is_approved AS is_approved, d.has_been_withdrawn AS withdrawn,
       d.source AS drug_node_source, d.drug_stage AS drug_stage,
       g.ensembl_id AS Gene_ID, g.symbol AS graph_symbol,
       type(r) AS Relationship, r.layer AS layer,
       r.action AS okg_action, r.relation_assertions AS relation_assertions, r.relation_conflict AS relation_conflict,
       r.sources_direct AS sources_direct, r.sources_indirect AS sources_indirect,
       r.evidence_level AS evidence_level, r.mechanisms_of_action AS okg_mechanism, r.source_urls AS source_urls,
       r.action_labels AS ext_action_labels, r.push AS ext_push, r.sources AS ext_sources,
       r.source_detail AS ext_source_detail, r.evidence AS ext_evidence, r.mechanism AS ext_mechanism,
       r.release AS ext_release, r.added AS ext_added, r.n_sources AS n_sources,
       r.external_action_direction AS external_action_direction, r.external_action_label AS external_action_label,
       r.external_action_sources AS external_action_sources
"""
GRAPH_TOTALS = """
MATCH (d:Drug) WITH count(d) AS drugs
MATCH ()-[r:ACTS_ON]->(:Gene)
RETURN drugs, count(r) AS acts_on, sum(CASE WHEN r.layer = $layer THEN 1 ELSE 0 END) AS acts_on_added
"""
password = os.environ.get("NEO4J_PASSWORD") or getpass.getpass("Neo4j password: ")
try:
    edges = neo4j_query(QUERY, {"ids": ids}, password)
    drug_degree = neo4j_query("MATCH (d:Drug)-[r:ACTS_ON]->(:Gene) RETURN d.node_id AS Drug_ID, count(r) AS drug_total_targets",
                              {}, password)
    totals = neo4j_query(GRAPH_TOTALS, {"layer": EXTERNAL_LAYER}, password).iloc[0]
except Exception as exc:
    raise SystemExit(f"Neo4j query failed ({exc}). Task 2 uses the live graph only; start Neo4j and set "
                     f"NEO4J_PASSWORD, then rerun.")

# ------------------------------------------------------------------ build table
added = edges.layer.eq(EXTERNAL_LAYER)
edges["Link_origin"] = added.map({True: "ChEMBL/DGIdb/GtoPdb/OpenTargets (added 21 Oct 2026)", False: "OptimusKG"})

# Action exactly as recorded on the edge; blank -> "unknown".
edges["Action_raw"] = [
    (lab if not blank(lab) else "unknown") if is_add else (act if not blank(act) else "unknown")
    for is_add, lab, act in zip(added, edges.ext_action_labels, edges.okg_action)]
edges["Action_direction"] = [
    PUSH.get(str(push), "unknown") if is_add else edge_direction(act, ass)
    for is_add, push, act, ass in zip(added, edges.ext_push, edges.okg_action, edges.relation_assertions)]
# OptimusKG edges with no direction: use the GtoPdb/DGIdb/ChEMBL direction stored on the same edge.
fill = (~added) & edges.Action_direction.eq("unknown") & edges.external_action_direction.isin(["down", "up"])
edges.loc[fill, "Action_direction"] = edges.loc[fill, "external_action_direction"]
edges["Direction_source"] = "none"
edges.loc[added & edges.Action_direction.ne("unknown"), "Direction_source"] = "added edge (push)"
edges.loc[(~added) & edges.Action_direction.ne("unknown") & ~fill, "Direction_source"] = "OptimusKG action"
edges.loc[fill, "Direction_source"] = "GtoPdb/DGIdb/ChEMBL on OptimusKG edge"
edges["Evidence_source"] = [
    "External/" + str(src).replace(";", "+") if is_add else "OptimusKG/" + ("" if blank(sd) else str(sd).replace(";", "+"))
    for is_add, src, sd in zip(added, edges.ext_sources, edges.sources_direct)]
edges["mechanism_of_action"] = edges.okg_mechanism.where(~added, edges.ext_mechanism)
edges["evidence_level"] = edges.evidence_level.where(~added, edges.ext_evidence)
edges["Database_ID"] = edges.Drug_ID + "|" + edges.Gene_ID

ann = genes[["ensembl_id", "Gene", "progression_class", "paper1_group", "fibrosis_direction",
             "inflammation_direction", "trajectory_cluster", "Tier", "in_shortlist"]]
ann = ann.rename(columns={"ensembl_id": "Gene_ID", "trajectory_cluster": "cluster"})
df = edges.merge(ann, on="Gene_ID", how="left").merge(drug_degree, on="Drug_ID", how="left")
assert len(df) == len(edges), "annotation merge changed the edge count"

# Gene direction used for comparison: fibrosis direction where available, otherwise inflammation (disease-group).
df["gene_direction"] = df.fibrosis_direction.fillna(df.inflammation_direction)
df["gene_direction_axis"] = df.fibrosis_direction.notna().map({True: "fibrosis", False: "inflammation"})
df["direction_comparable"] = df.Action_direction.isin(["up", "down"]) & df.gene_direction.notna()
df["cluster_comparable"] = df.Action_direction.isin(["up", "down"]) & df.cluster.notna()

cols = ["Drug", "Drug_ID", "Gene", "Gene_ID", "Relationship", "Action_raw", "Action_direction", "Direction_source",
        "Evidence_source", "Database_ID", "Link_origin", "evidence_level", "n_sources", "relation_assertions",
        "relation_conflict", "mechanism_of_action", "sources_indirect", "source_urls",
        "ext_source_detail", "ext_release", "ext_added",
        "external_action_direction", "external_action_label", "external_action_sources",
        "molecule_type", "max_phase", "is_approved", "withdrawn", "drug_stage", "drug_total_targets",
        "progression_class", "paper1_group", "fibrosis_direction", "inflammation_direction",
        "gene_direction", "gene_direction_axis", "cluster", "Tier", "in_shortlist",
        "direction_comparable", "cluster_comparable"]
df = df[[c for c in cols if c in df.columns]].sort_values(["Gene", "Drug"])
df.to_csv(OUT / "drug_gene_edges.csv", index=False)

# ------------------------------------------------------------------ graph check
log = json.loads(NEW_LINKS_LOG.read_text(encoding="utf-8"))
gcheck = {
    "source_used": "Neo4j (live liver_kg, read-only)",
    "graph_drug_nodes": int(totals.drugs),
    "graph_acts_on_total": int(totals.acts_on),
    "graph_acts_on_added_21oct": int(totals.acts_on_added),
    "write_log_acts_on_after": int(log["after"]["acts_on"]) + int(log["after"]["external_links"]),
    "write_log_drugs_after": int(log["after"]["drugs"]) + int(log["after"]["external_drugs"]),
    "edges_on_progression_genes": len(df),
    "edges_optimuskg": int((df.Link_origin == "OptimusKG").sum()),
    "edges_added_21oct": int((df.Link_origin != "OptimusKG").sum()),
    "duplicate_drug_gene_pairs": int(df.duplicated(["Drug_ID", "Gene_ID"]).sum()),
    "edges_with_external_action_props": int(df.external_action_direction.notna().sum()) if "external_action_direction" in df else 0,
}
gcheck["graph_matches_write_log"] = (gcheck["graph_acts_on_total"] == gcheck["write_log_acts_on_after"]
                                     and gcheck["graph_drug_nodes"] == gcheck["write_log_drugs_after"])
(OUT / "graph_check.json").write_text(json.dumps(gcheck, indent=2))

# ------------------------------------------------------------------ summary numbers
n = len(df)
g_conn = df.Gene_ID.nunique()
known = int(df.Action_direction.isin(["up", "down"]).sum())
dir_counts = df.Action_direction.value_counts().reindex(["down", "up", "mixed", "unknown"]).fillna(0).astype(int)
cls_tot = genes.progression_class.value_counts()
cls_conn = df.groupby("progression_class").Gene_ID.nunique()
cls_drugs = df.groupby("progression_class").Drug_ID.nunique()
cls_edges = df.groupby("progression_class").size()
conflicts = int(df.relation_conflict.astype("boolean").fillna(False).sum())
approved = df.loc[df.is_approved == True, "Drug_ID"].nunique()
shortlist_conn = df.loc[df.in_shortlist == True, "Gene_ID"].nunique()
comp_dir = int(df.direction_comparable.sum())
comp_clu = int(df.cluster_comparable.sum())
by_origin = df.groupby("Link_origin").agg(edges=("Database_ID", "size"), genes=("Gene_ID", "nunique"),
                                          drugs=("Drug_ID", "nunique"),
                                          known=("Action_direction", lambda s: int(s.isin(["up", "down"]).sum())))

summary_rows = {
    "progression_genes": len(genes), "genes_with_drug": g_conn, "unique_drugs": df.Drug_ID.nunique(),
    "edges": n, "edges_optimuskg": gcheck["edges_optimuskg"], "edges_added_21oct": gcheck["edges_added_21oct"],
    "edges_known_direction": known, "edges_direction_comparable": comp_dir,
    "edges_cluster_comparable": comp_clu, "relation_conflict_edges": conflicts,
    "shortlist_genes_with_drug": shortlist_conn,
}
pd.Series(summary_rows).to_csv(OUT / "summary_numbers.csv", header=["value"])


def pct(a, b):
    return f"{100 * a / b:.1f}%" if b else "–"


def md(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out)


summary = f"""# Paper 2 — Task 2: drug–gene relationships

Gene set: the corrected Paper 1 union of **{len(genes):,} progression genes** (Task 1). "Inflammation" columns hold Paper 1's ordered disease-group axis.

**Source:** the present Neo4j graph only (live, read-only; no CSV fallback). Query: `MATCH (d:Drug)-[r:ACTS_ON]->(g:Gene) WHERE g.ensembl_id IN $ids`.
Script: `Paper2_Results/code/task2_drug_gene_edges.py`. Output: `drug_gene_edges.csv` ({n:,} rows).

## Graph check

{md(["Check", "Result"], [
    ("Drug nodes in graph", f"{gcheck['graph_drug_nodes']:,} (write log after 21 Oct: {gcheck['write_log_drugs_after']:,})"),
    ("ACTS_ON relationships in graph", f"{gcheck['graph_acts_on_total']:,} (write log: {gcheck['write_log_acts_on_after']:,}; of which added 21 Oct: {gcheck['graph_acts_on_added_21oct']:,})"),
    ("Graph matches the 21 Oct write log", "YES" if gcheck['graph_matches_write_log'] else "NO — investigate"),
    ("Duplicate drug–gene pairs among progression-gene edges", gcheck['duplicate_drug_gene_pairs']),
])}

## Headline numbers

{md(["Measure", "Value"], [
    ("Progression genes with ≥1 drug", f"{g_conn:,} / {len(genes):,} ({pct(g_conn, len(genes))})"),
    ("Unique drugs", f"{df.Drug_ID.nunique():,} ({approved:,} flagged approved)"),
    ("Total drug–gene relationships", f"{n:,}"),
    ("Relationships with known direction (down/up)", f"{known:,} ({pct(known, n)})"),
    ("…comparable with Paper 1 gene direction (positive/negative)", f"{comp_dir:,} ({pct(comp_dir, n)})"),
    ("…comparable with an F0–F4 trajectory cluster", f"{comp_clu:,} ({pct(comp_clu, n)})"),
    ("23-gene shortlist genes with ≥1 drug", f"{shortlist_conn} / 23"),
])}

"Comparable" means the drug pushes the gene product down or up **and** Paper 1 gives the gene a direction (fibrosis direction; disease-group direction for disease-group-only genes) or an F0–F4 trajectory cluster. Whether each comparison is direction-consistent is Task 7. This is an early mechanism-consistency check, not therapeutic proof.

## By link origin

{md(["Origin", "Edges", "Genes", "Drugs", "Known direction", "%"], [
    (o, f"{r.edges:,}", f"{r.genes:,}", f"{r.drugs:,}", f"{r.known:,}", pct(r.known, r.edges)) for o, r in by_origin.iterrows()])}

## Drug action direction

{md(["Direction", "Edges", "%"], [(k, f"{v:,}", pct(v, n)) for k, v in dir_counts.items()])}

- OptimusKG edges: direction from the consolidated `action` (down: {", ".join(sorted(DOWN))}; up: {", ".join(sorted(UP))}; mixed: {", ".join(sorted(MIXED))} or sources disagree; unknown: {", ".join(sorted(UNKNOWN))} or blank). OptimusKG `relation_conflict` edges: {conflicts:,}.
- Added edges: direction from the edge's `push` field (down / up; "conflicting" counted as mixed; unknown when no action label).
- `Action_raw` keeps the action text exactly as stored on the edge; blank actions are written as "unknown".
- OptimusKG edges with no direction take the GtoPdb/DGIdb/ChEMBL direction stored on the same edge (`external_action_direction`): {int(fill.sum()):,} edges. It never overrides a direction OptimusKG gives.

{md(["Direction source", "Edges"], [(k, f"{v:,}") for k, v in df.Direction_source.value_counts().items()])}

## By progression class

{md(["Class", "Genes", "Genes with drug", "%", "Unique drugs", "Edges"], [
    (c, f"{cls_tot.get(c, 0):,}", f"{cls_conn.get(c, 0):,}", pct(cls_conn.get(c, 0), cls_tot.get(c, 0)),
     f"{cls_drugs.get(c, 0):,}", f"{cls_edges.get(c, 0):,}") for c in ["fibrosis_only", "shared", "inflammation_only"]
])}

(Statistics comparing these groups are Task 3.)

## Notes

- `Database_ID` is the stable pair `drug ID | Ensembl gene ID`; the graph holds no per-assertion accession IDs.
- No binding strength, dose or tissue context is stored on the edges.
- Action labels such as TARGET, ENZYME, TRANSPORTER and CARRIER say what the gene product is for the drug, not whether the drug switches it down or up; they count as unknown direction.
"""
(OUT / "summary.md").write_text(summary, encoding="utf-8")
print(summary)
print(json.dumps(gcheck, indent=2))
