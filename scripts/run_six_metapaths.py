"""Run the 6 metapaths for one drug set x one gene set, and write one Excel file.

Gene* = a gene from this run's gene set (it gives the stage); "Gene" = any
other gene, a direct interaction partner. All paths end at NAFLD (EFO_0003095).

 With pathway (needs Pathway -IMPLICATED_IN-> NAFLD, see build_pathway_nafld_links.py):
  A : Drug -ACTS_ON-> Gene* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD
  B : Drug -ACTS_ON-> Gene -INTERACTS_WITH- Gene* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD
  C : Drug -ACTS_ON-> Gene* -INTERACTS_WITH- Gene -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD
 Without pathway:
  A0: Drug -ACTS_ON-> Gene* -TRACKS_FIBROSIS-> NAFLD
  B0: Drug -ACTS_ON-> Gene -INTERACTS_WITH- Gene* -TRACKS_FIBROSIS-> NAFLD
  C0: Drug -ACTS_ON-> Gene* -INTERACTS_WITH- Gene -[TRACKS_FIBROSIS|ASSOCIATED_WITH]-> NAFLD

A route is kept only if it passes all rules:
  1. drug -> target direction known and correct (direction_check == "direction-consistent")
  2. the in-between gene (B/B0: drug target, C/C0: partner) has <= HUB_MAX partners
  3. pathway size 5-200 genes (A, B, C)
  4. the drug has no TREATS / OFF_LABEL_FOR to NAFLD, NASH or cirrhosis of liver

Output (in --out):
  <run>.xlsx           6 sheets A, B, C, A0, B0, C0 with the kept routes
  <run>_all_routes.csv every route found, result = KEEP or the drop reason(s)

route_id = metapath|drug_id|target|partner|gene_star|pathway_id identifies a
route across runs, so a later step can remove duplicates.

Example (from the project root; Paper2_Results/code/metapath_pipeline.py runs this for each gene set):
  python scripts/run_six_metapaths.py --run Run4_onlyconsistent_tier3 \
      --drugs Paper2_Results/Metapath_runs/only_consistent_drugs.xlsx \
      --genes Paper1_HighConfidence/11_final_all_3845_genes.csv --tier "Tier 3" \
      --metapaths B --out Paper2_Results/Metapath_runs
Gene files are recognised by their columns: Gene/Direction/Cluster/Largest_transition
(shortlist) or symbol/direction/cluster/largest_transition (Tier 1).
Password: NEO4J_PASSWORD, or else the neo4j-liver-kg entry in .mcp.json.
"""
import argparse
import json
import os
from pathlib import Path

import pandas as pd
from neo4j import GraphDatabase

ROOT = Path(__file__).resolve().parent.parent
NAFLD = "EFO_0003095"
KNOWN_LIVER = ["EFO_0003095", "EFO_1001249", "EFO_0001422"]  # NAFLD, NASH, cirrhosis
HUB_MAX = 500
SIZE_MIN, SIZE_MAX = 5, 200
SHEETS = ["A", "B", "C", "A0", "B0", "C0"]

# Each query binds: t = drug target, g = Gene*, mid = in-between gene (or null),
# s = the gene the last step leaves from. Tails then add the pathway or disease step.
HEADS = {
    "A": "MATCH (t) WHERE t.symbol IN $genes WITH pr, t, t AS g, null AS mid, t AS s",
    "B": """MATCH (t)-[:INTERACTS_WITH]-(g:Gene) WHERE g.symbol IN $genes AND g <> t
            WITH DISTINCT pr, t, g, t AS mid, g AS s""",
    "C": """MATCH (t) WHERE t.symbol IN $genes
            MATCH (t)-[:INTERACTS_WITH]-(h:Gene) WHERE h <> t
            WITH DISTINCT pr, t, t AS g, h AS mid, h AS s""",
}
START = """UNWIND $pairs AS pr
MATCH (:Drug {node_id: pr.drug_id})-[:ACTS_ON]->(t:Gene {symbol: pr.gene})
"""
RET = """pr.drug AS drug, pr.drug_id AS drug_id, t.symbol AS target,
         mid.symbol AS partner_gene,
         CASE WHEN mid IS NULL THEN null ELSE COUNT { (mid)-[:INTERACTS_WITH]-() } END AS partner_n_partners,
         g.symbol AS gene_star"""
PATHWAY_TAIL = f"""
MATCH (s)-[:IN_PATHWAY]->(p:Pathway)-[ip:IMPLICATED_IN]->(:Disease {{node_id: '{NAFLD}'}})
RETURN {RET}, p.node_id AS pathway_id, p.name AS pathway,
       coalesce(p.source, 'Reactome') AS pathway_library,
       COUNT {{ (p)<-[:IN_PATHWAY]-(:Gene) }} AS pathway_size, ip.clusters AS pathway_clusters,
       'IMPLICATED_IN' AS last_link"""
DISEASE_TAIL = {
    "A0": f"MATCH (s)-[:TRACKS_FIBROSIS]->(:Disease {{node_id: '{NAFLD}'}}) RETURN DISTINCT {RET}, 'TRACKS_FIBROSIS' AS last_link",
    "B0": f"MATCH (s)-[:TRACKS_FIBROSIS]->(:Disease {{node_id: '{NAFLD}'}}) RETURN DISTINCT {RET}, 'TRACKS_FIBROSIS' AS last_link",
    "C0": f"""MATCH (s)-[r:TRACKS_FIBROSIS|ASSOCIATED_WITH]->(:Disease {{node_id: '{NAFLD}'}})
              RETURN {RET}, collect(DISTINCT type(r)) AS links""",
}
QUERIES = {
    "A": START + HEADS["A"] + PATHWAY_TAIL,
    "B": START + HEADS["B"] + PATHWAY_TAIL,
    "C": START + HEADS["C"] + PATHWAY_TAIL,
    "A0": START + HEADS["A"] + "\n" + DISEASE_TAIL["A0"],
    "B0": START + HEADS["B"] + "\n" + DISEASE_TAIL["B0"],
    "C0": START + HEADS["C"] + "\n" + DISEASE_TAIL["C0"],
}

COLUMNS = ["run", "route_id", "metapath", "drug", "drug_id", "target", "drug_action",
           "direction_check", "drug_mixed_direction", "stable_drug", "partner_gene", "partner_n_partners", "gene_star",
           "gene_star_direction", "gene_star_cluster", "stage", "pathway", "pathway_id",
           "pathway_library", "pathway_size", "pathway_clusters", "stage_match", "last_link"]


def password():
    if os.environ.get("NEO4J_PASSWORD"):
        return os.environ["NEO4J_PASSWORD"]
    return json.load(open(ROOT / ".mcp.json"))["mcpServers"]["neo4j-liver-kg"]["env"]["NEO4J_PASSWORD"]


def load_drugs(path):
    """Drug -> target pairs from either drug file format.

    .csv  : one row per drug-gene pair, with link_type, drug_action and
        direction_check columns.
    .xlsx : one row per drug (e.g. only_consistent_drugs.xlsx), with the
        correct-direction targets in consistent_targets ("PTGS1;PTGS2") and
        direction = "only consistent" / "mixed". Only the listed targets are
        used; drug_action is filled from the graph later.
    drug_mixed_direction = True when the drug also pushes some other target
    the wrong way (flagged, not removed).
    """
    if str(path).lower().endswith((".xlsx", ".xls")):
        d = pd.read_excel(path)
        rows = [dict(drug=r.Drug, drug_id=r.Drug_ID, gene=g.strip(), drug_action=None,
                     direction_check="direction-consistent",
                     mixed_direction=r.direction == "mixed", stable_drug=bool(r.stable_drug))
                for r in d.itertuples() for g in str(r.consistent_targets).split(";") if g.strip()]
        return pd.DataFrame(rows)
    d = pd.read_csv(path)
    d = d[d.link_type == "pharmacological target"].copy()
    wrong = d[d.direction_check.str.contains("inconsistent", na=False)].drug_id.unique()
    return d.assign(mixed_direction=d.drug_id.isin(wrong), stable_drug=True)[
        ["drug", "drug_id", "gene", "drug_action", "direction_check", "mixed_direction", "stable_drug"]]


def load_genes(path, tier=None):
    g = pd.read_csv(path)
    if tier:
        col = "Tier" if "Tier" in g.columns else "tier"
        g = g[g[col].str.startswith(tier)]
    if "Gene" in g.columns:
        g = g.rename(columns={"Gene": "symbol", "Direction": "direction", "Cluster": "cluster",
                              "Largest_transition": "largest_transition"})
    g["direction"] = g.direction.map({"increasing": "rises", "decreasing": "falls"}).fillna(g.direction)
    return g.set_index("symbol")[["direction", "cluster", "largest_transition"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="run name, e.g. Run4_onlyconsistent_tier3")
    ap.add_argument("--drugs", required=True, help="drug file (.xlsx one row per drug, or .csv one row per drug-gene pair)")
    ap.add_argument("--genes", required=True, help="gene set file (shortlist, Tier 1, or all 3,845)")
    ap.add_argument("--tier", default=None,
                    help='keep only genes whose Tier starts with this, e.g. "Tier 2"')
    ap.add_argument("--out", default=str(ROOT / "Paper2_Metapath_runs"))
    ap.add_argument("--metapaths", default=",".join(SHEETS),
                    help='comma-separated subset of A,B,C,A0,B0,C0 to run, e.g. "B" (default: all six)')
    args = ap.parse_args()
    run_mps = [m.strip() for m in args.metapaths.split(",") if m.strip()]
    assert set(run_mps) <= set(SHEETS), f"unknown metapath in {run_mps}"

    tg = load_drugs(args.drugs)
    pairs = tg[["drug", "drug_id", "gene"]].drop_duplicates().to_dict("records")
    action = {(r.drug_id, r.gene): r.drug_action for r in tg.itertuples()}
    direction = {(r.drug_id, r.gene): r.direction_check for r in tg.itertuples()}
    mixed = dict(zip(tg.drug_id, tg.mixed_direction))
    stable = dict(zip(tg.drug_id, tg.stable_drug))
    genes = load_genes(args.genes, args.tier)
    print(f"{args.run}: drugs {tg.drug_id.nunique()} | drug-target pairs {len(pairs)} | Gene* {len(genes)}")

    drv = GraphDatabase.driver(os.environ.get("NEO4J_URI", "bolt://localhost:7687"),
                               auth=(os.environ.get("NEO4J_USERNAME", "neo4j"), password()))
    rows = []
    with drv.session(database="neo4j") as ses:
        if not ses.run("MATCH ()-[r:IMPLICATED_IN]->(:Disease {node_id:$d}) RETURN count(r) > 0 AS ok",
                       d=NAFLD).single()["ok"]:
            raise SystemExit("No Pathway -IMPLICATED_IN-> NAFLD links. Run build_pathway_nafld_links.py --load first.")
        known = {r["id"]: r["n"] for r in ses.run(
            """UNWIND $ids AS id MATCH (d:Drug {node_id: id})
               OPTIONAL MATCH (d)-[:TREATS|OFF_LABEL_FOR]->(x:Disease) WHERE x.node_id IN $liver
               RETURN id, count(x) AS n""",
            ids=list({p["drug_id"] for p in pairs}), liver=KNOWN_LIVER)}
        missing = [p for p in pairs if pd.isna(action.get((p["drug_id"], p["gene"])))]
        for r in ses.run("""UNWIND $pairs AS pr
                            MATCH (:Drug {node_id: pr.drug_id})-[a:ACTS_ON]->(:Gene {symbol: pr.gene})
                            RETURN pr.drug_id AS id, pr.gene AS gene, a.unified_actions AS acts""",
                         pairs=missing):
            acts = [x for x in (r["acts"] or []) if x != "TARGET"] or (r["acts"] or [])
            action[(r["id"], r["gene"])] = ";".join(acts)
        for mp in run_mps:
            for r in ses.run(QUERIES[mp], pairs=pairs, genes=list(genes.index)).data():
                if mp == "C0":
                    r["last_link"] = "+".join(sorted(r.pop("links")))
                r["metapath"] = mp
                rows.append(r)
    drv.close()

    d = pd.DataFrame(rows)
    for c in ["pathway", "pathway_id", "pathway_library", "pathway_size", "pathway_clusters"]:
        if c not in d:
            d[c] = None
    d["run"] = args.run
    d["drug_action"] = [action[(a, b)] for a, b in zip(d.drug_id, d.target)]
    d["direction_check"] = [direction[(a, b)] for a, b in zip(d.drug_id, d.target)]
    d["drug_mixed_direction"] = d.drug_id.map(mixed)
    d["stable_drug"] = d.drug_id.map(stable)
    d["gene_star_direction"] = d.gene_star.map(genes.direction)
    d["gene_star_cluster"] = d.gene_star.map(genes.cluster)
    d["stage"] = d.gene_star.map(genes.largest_transition)
    d["stage_match"] = [(c in pc.split(";")) if isinstance(pc, str) else None
                        for c, pc in zip(d.gene_star_cluster, d.pathway_clusters)]
    d["route_id"] = ["|".join("" if pd.isna(x) else str(x) for x in t) for t in
                     zip(d.metapath, d.drug_id, d.target, d.partner_gene, d.gene_star, d.pathway_id)]

    reasons = []
    for r in d.itertuples():
        why = []
        if r.direction_check != "direction-consistent":
            why.append("drug direction not known/correct")
        if pd.notna(r.partner_n_partners) and r.partner_n_partners > HUB_MAX:
            why.append(f"hub in-between gene (>{HUB_MAX} partners)")
        if r.metapath in ("A", "B", "C") and not SIZE_MIN <= r.pathway_size <= SIZE_MAX:
            why.append(f"pathway size outside {SIZE_MIN}-{SIZE_MAX}")
        if known.get(r.drug_id, 0):
            why.append("drug already treats a liver disease")
        reasons.append("; ".join(why) or "KEEP")
    d["result"] = reasons
    d = d.drop_duplicates("route_id")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    d.sort_values(["metapath", "result", "drug"]).to_csv(out / f"{args.run}_all_routes.csv", index=False)
    with pd.ExcelWriter(out / f"{args.run}.xlsx") as xw:
        for mp in run_mps:
            k = d[(d.metapath == mp) & (d.result == "KEEP")][COLUMNS].sort_values(["drug", "gene_star", "pathway"])
            if k.empty:
                k = pd.DataFrame([{"run": args.run, "metapath": mp, "route_id": "no kept routes"}], columns=COLUMNS)
            k.to_excel(xw, sheet_name=mp, index=False)

    print(f"routes found {len(d)} | kept {int((d.result == 'KEEP').sum())}")
    s = d.groupby("metapath").agg(found=("route_id", "size"), kept=("result", lambda x: int((x == "KEEP").sum())),
                                  kept_drugs=("drug", lambda x: x[d.loc[x.index, "result"] == "KEEP"].nunique()))
    print(s.reindex(run_mps).fillna(0).astype(int).to_string())
    k = d[d.result == "KEEP"]
    print("kept drugs:", sorted(k.drug.unique()))
    print(f"written: {out / (args.run + '.xlsx')}")


if __name__ == "__main__":
    main()
