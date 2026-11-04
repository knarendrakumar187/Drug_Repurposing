"""
paper2_write_recovered_to_graph.py
----------------------------------
Write the two recoveries back onto the graph: drug properties on Drug nodes,
per-edge evidence on ACTS_ON relationships.

WHY HERE AND NOT ON GENE NODES
Nothing recovered is gene-level. It is drug-level (structure, molecule type,
clinical phase) or edge-level (asserting databases, reference URLs, written
mechanisms). Gene nodes are also the one place in this graph where duplicating
derived data has already gone wrong: `fibrosis_rho` sits on 3,645 Gene nodes
while 3,845 carry a TRACKS_FIBROSIS edge, and 191 have the property with no
edge. Two copies of one number, drifted apart.

Drug nodes carry `name` and `node_id` and nothing else; ACTS_ON carries
`action` and nothing else. There is no second copy to drift from, which is what
makes these two safe targets and Gene nodes an unsafe one.

SAFETY
This script is additive and reversible by construction:

  * It only ever SETs keys that do not already exist. Before writing it counts
    every target key across the graph and REFUSES to run if any is present,
    so it can never overwrite Team 1's data or its own earlier run.
  * It writes nothing without --apply. The default is a dry run that reports
    exactly what would change.
  * --rollback removes exactly the keys this script sets, and nothing else.
  * Values come from paper2_integrate_recovered_evidence.load_drug_metadata(),
    the same function that produced the delivered CSVs, so the graph cannot
    disagree with the tables.

Properties carry their provenance. Every recovered field that can come from
either recovery gets a companion *_source property naming which one supplied
it, exactly as the CSV columns do -- the mentor's instruction of 13 October was
to keep the two separate, and that has to survive into the graph.

`name` is never touched. Where OptimusKG stored an accession instead of a name,
the recovered name is written to `recovered_name` beside it, not over it.

Usage:
    uv run --with neo4j --with pandas python scripts/paper2_write_recovered_to_graph.py
    uv run --with neo4j --with pandas python scripts/paper2_write_recovered_to_graph.py --apply
    uv run --with neo4j --with pandas python scripts/paper2_write_recovered_to_graph.py --rollback
"""

import json
import math
import os
import sys

import pandas as pd
from neo4j import GraphDatabase

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paper2_integrate_recovered_evidence import load_drug_metadata

EVIDENCE = "data/optimuskg_recovered/drug_gene_evidence.csv"
BATCH = 2000

# Drug node: recovered field -> property name on the node.
DRUG_PROPS = [
    "smiles", "smiles_source",
    "molecule_type", "molecule_type_source",
    "max_phase", "max_phase_source",
    "recovered_name", "recovered_name_source",
    "is_approved", "has_been_withdrawn",
    "inchi_key", "cd_formula", "cd_mol_weight",
]

# ACTS_ON relationship properties.
EDGE_PROPS = [
    "sources_direct", "sources_indirect", "relation_assertions",
    "relation_conflict", "evidence_level", "source_urls",
    "mechanisms_of_action", "n_sources",
]


def driver():
    """Credentials come from the environment, never from source.
    NEO4J_PASSWORD (prompted if unset); NEO4J_URI, NEO4J_USERNAME and
    NEO4J_DATABASE default to bolt://localhost:7687, neo4j and neo4j."""
    import getpass
    pw = os.environ.get("NEO4J_PASSWORD") or getpass.getpass("Neo4j password: ")
    return (GraphDatabase.driver(os.environ.get("NEO4J_URI", "bolt://localhost:7687"),
                                 auth=(os.environ.get("NEO4J_USERNAME", "neo4j"), pw)),
            os.environ.get("NEO4J_DATABASE", "neo4j"))


def clean(v):
    """Neo4j takes no NaN. A missing value becomes an absent property."""
    if v is None:
        return None
    if isinstance(v, float) and math.isnan(v):
        return None
    if isinstance(v, (pd.Timestamp,)):
        return str(v)
    if hasattr(v, "item"):          # numpy scalar
        v = v.item()
    if isinstance(v, float) and math.isnan(v):
        return None
    return v


def existing_keys(sess):
    """Which target keys are already present, and on how many nodes/edges."""
    found = {}
    for p in DRUG_PROPS:
        n = sess.run(f"MATCH (d:Drug) WHERE d.`{p}` IS NOT NULL "
                     "RETURN count(d) AS n").single()["n"]
        if n:
            found[f"Drug.{p}"] = n
    for p in EDGE_PROPS:
        n = sess.run(f"MATCH ()-[r:ACTS_ON]->() WHERE r.`{p}` IS NOT NULL "
                     "RETURN count(r) AS n").single()["n"]
        if n:
            found[f"ACTS_ON.{p}"] = n
    return found


def build_drug_rows(md, graph_ids):
    md = md.copy()
    # `name` is never overwritten. The recovered name lives beside it.
    md["recovered_name"] = md["drug_name"]
    md["recovered_name_source"] = md["drug_name_source"]
    md = md[md.Drug_ID.isin(graph_ids)]

    rows = []
    for rec in md.to_dict("records"):
        props = {p: clean(rec.get(p)) for p in DRUG_PROPS}
        props = {k: v for k, v in props.items() if v is not None}
        if props:
            rows.append({"id": rec["Drug_ID"], "props": props})
    return rows


def build_edge_rows(ev):
    rows = []
    for rec in ev.to_dict("records"):
        props = {p: clean(rec.get(p)) for p in EDGE_PROPS}
        props = {k: v for k, v in props.items() if v is not None}
        if props:
            rows.append({"drug": rec["Drug_ID"], "gene": rec["Gene_ID"],
                         "props": props})
    return rows


def write_drugs(sess, rows):
    q = ("UNWIND $rows AS row "
         "MATCH (d:Drug {node_id: row.id}) "
         "SET d += row.props "
         "RETURN count(d) AS n")
    done = 0
    for i in range(0, len(rows), BATCH):
        done += sess.run(q, rows=rows[i:i + BATCH]).single()["n"]
    return done


def write_edges(sess, rows):
    q = ("UNWIND $rows AS row "
         "MATCH (d:Drug {node_id: row.drug})-[r:ACTS_ON]->"
         "(g:Gene {ensembl_id: row.gene}) "
         "SET r += row.props "
         "RETURN count(r) AS n")
    done = 0
    for i in range(0, len(rows), BATCH):
        done += sess.run(q, rows=rows[i:i + BATCH]).single()["n"]
    return done


def rollback(sess):
    print("ROLLBACK -- removing exactly the keys this script sets\n")
    for p in DRUG_PROPS:
        n = sess.run(f"MATCH (d:Drug) WHERE d.`{p}` IS NOT NULL "
                     f"REMOVE d.`{p}` RETURN count(d) AS n").single()["n"]
        if n:
            print(f"   Drug.{p:<24} removed from {n:,}")
    for p in EDGE_PROPS:
        n = sess.run(f"MATCH ()-[r:ACTS_ON]->() WHERE r.`{p}` IS NOT NULL "
                     f"REMOVE r.`{p}` RETURN count(r) AS n").single()["n"]
        if n:
            print(f"   ACTS_ON.{p:<21} removed from {n:,}")
    left = existing_keys(sess)
    print(f"\nremaining target keys: {left if left else 'none -- clean'}")


def main():
    apply_ = "--apply" in sys.argv
    roll = "--rollback" in sys.argv

    drv, db = driver()
    try:
        with drv.session(database=db) as sess:
            if roll:
                rollback(sess)
                return

            graph_ids = {r["id"] for r in sess.run(
                "MATCH (d:Drug) RETURN d.node_id AS id").data()}
            n_edges = sess.run("MATCH ()-[r:ACTS_ON]->() "
                               "RETURN count(r) AS n").single()["n"]
            print(f"graph: {len(graph_ids):,} Drug nodes, {n_edges:,} ACTS_ON edges")

            # --- refuse to overwrite anything, ever ------------------------
            clash = existing_keys(sess)
            if clash:
                print("\nREFUSING TO WRITE -- these target keys already exist:")
                for k, v in clash.items():
                    print(f"   {k}: {v:,}")
                print("\nThis script only ever adds. Run --rollback first if "
                      "you intend to rewrite them.")
                raise SystemExit(1)
            print("no target key exists yet -- safe to add")

            md, chembl_added, _ = load_drug_metadata()
            ev = pd.read_csv(EVIDENCE, low_memory=False)

            drug_rows = build_drug_rows(md, graph_ids)
            edge_rows = build_edge_rows(ev)
            print(f"\nprepared {len(drug_rows):,} Drug updates, "
                  f"{len(edge_rows):,} candidate edge updates")

            if not apply_:
                print("\nDRY RUN -- nothing written. Re-run with --apply.")
                show_plan(drug_rows, edge_rows)
                return

            n_d = write_drugs(sess, drug_rows)
            n_e = write_edges(sess, edge_rows)
            print(f"\nwrote: {n_d:,} Drug nodes, {n_e:,} ACTS_ON edges")

            verify(sess, n_edges)
    finally:
        drv.close()


def show_plan(drug_rows, edge_rows):
    from collections import Counter
    dc, ec = Counter(), Counter()
    for r in drug_rows:
        dc.update(r["props"].keys())
    for r in edge_rows:
        ec.update(r["props"].keys())
    print("\n   Drug node properties to be set")
    for p in DRUG_PROPS:
        if dc[p]:
            print(f"      {p:<24} {dc[p]:>7,}")
    print("\n   ACTS_ON properties to be set (where the edge exists)")
    for p in EDGE_PROPS:
        if ec[p]:
            print(f"      {p:<24} {ec[p]:>7,}")


def verify(sess, n_edges):
    """Confirm what landed, and that nothing pre-existing moved."""
    print("\nverification")
    for p in ["smiles", "molecule_type", "max_phase", "recovered_name"]:
        n = sess.run(f"MATCH (d:Drug) WHERE d.`{p}` IS NOT NULL "
                     "RETURN count(d) AS n").single()["n"]
        print(f"   Drug.{p:<18} {n:>7,}")
    for p in ["sources_direct", "evidence_level", "source_urls",
              "mechanisms_of_action"]:
        n = sess.run(f"MATCH ()-[r:ACTS_ON]->() WHERE r.`{p}` IS NOT NULL "
                     "RETURN count(r) AS n").single()["n"]
        print(f"   ACTS_ON.{p:<15} {n:>7,}")

    # nothing pre-existing may have changed
    now = sess.run("MATCH ()-[r:ACTS_ON]->() RETURN count(r) AS n").single()["n"]
    acts = sess.run("MATCH ()-[r:ACTS_ON]->() WHERE r.action IS NOT NULL "
                    "RETURN count(r) AS n").single()["n"]
    names = sess.run("MATCH (d:Drug) WHERE d.name IS NOT NULL "
                     "RETURN count(d) AS n").single()["n"]
    print(f"\n   ACTS_ON edges      {now:,}  (was {n_edges:,})")
    print(f"   r.action intact    {acts:,}")
    print(f"   d.name intact      {names:,}")
    assert now == n_edges, "edge count changed -- this script must not create edges"
    assert acts == 20674, f"action property count moved: {acts}"
    assert names == 12025, f"name property count moved: {names}"
    print("\n   nothing pre-existing was modified")


if __name__ == "__main__":
    main()
