"""Longer version of metapath B: two genes between the drug and Gene*.

  B  : Drug -ACTS_ON-> Gene1 -INTERACTS_WITH- Gene* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD
  BB : Drug -ACTS_ON-> Gene1 -INTERACTS_WITH- Gene2 -INTERACTS_WITH- Gene* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD

Run only for the drugs that already have metapath B routes in a run file
(made by run_six_metapaths.py --metapaths B). Same rules as run_six_metapaths.py:
  1. only the drug's direction-consistent targets (from the drug file)
  2. no hub in-between gene: Gene1 and Gene2 each have <= HUB_MAX partners
  3. pathway size 5-200 genes
  4. Gene1, Gene2 and Gene* are three different genes
(The drugs already passed "no TREATS/OFF_LABEL_FOR to a liver disease" in the run.)

Output (in --out):
  <name>.xlsx  sheet Summary (drug -> NAFLD, number of routes) and sheet BB (one row per route)
  <name>.csv   the same routes with extra columns (partner counts, stage, pathway id, also_in_B)

Example (from the project root):
  python scripts/run_two_gene_metapath.py --run-file Paper2_Results/Metapath_runs/Run4_onlyconsistent_tier3.xlsx \
      --drugs Paper2_Results/Metapath_runs/only_consistent_drugs.xlsx \
      --genes Paper1_HighConfidence/11_final_all_3845_genes.csv --tier "Tier 3" \
      --name Run4_onlyconsistent_tier3_BB --out Paper2_Results/Metapath_runs
"""
import argparse
import json
import os
from pathlib import Path

import pandas as pd
from neo4j import GraphDatabase

ROOT = Path(__file__).resolve().parent.parent
NAFLD = "EFO_0003095"
NAFLD_NAME = "non-alcoholic fatty liver disease (NAFLD)"
HUB_MAX = 500
SIZE_MIN, SIZE_MAX = 5, 200

PATHWAYS = """
MATCH (p:Pathway)-[:IMPLICATED_IN]->(:Disease {node_id: $nafld})
WITH p, COUNT { (p)<-[:IN_PATHWAY]-(:Gene) } AS size
WHERE $smin <= size <= $smax
RETURN p.node_id AS id, p.name AS name, size
"""
# Gene1 = drug target, Gene2 = its partner, Gene* = a gene of this run's gene set.
QUERY = """
UNWIND $pairs AS pr
MATCH (:Drug {node_id: pr.drug_id})-[:ACTS_ON]->(t:Gene {symbol: pr.gene})
WITH pr, t, COUNT { (t)-[:INTERACTS_WITH]-() } AS t_n WHERE t_n <= $hub
MATCH (t)-[:INTERACTS_WITH]-(h:Gene) WHERE h <> t
WITH DISTINCT pr, t, t_n, h
WITH pr, t, t_n, h, COUNT { (h)-[:INTERACTS_WITH]-() } AS h_n WHERE h_n <= $hub
MATCH (h)-[:INTERACTS_WITH]-(g:Gene) WHERE g.symbol IN $genes AND g <> t AND g <> h
WITH DISTINCT pr, t, t_n, h, h_n, g
MATCH (g)-[:IN_PATHWAY]->(p:Pathway) WHERE p.node_id IN $pathways
RETURN pr.drug AS drug, pr.drug_id AS drug_id, t.symbol AS gene1, t_n AS gene1_n_partners,
       h.symbol AS gene2, h_n AS gene2_n_partners, g.symbol AS gene_star,
       EXISTS { (t)-[:INTERACTS_WITH]-(g) } AS gene1_touches_gene_star,
       p.node_id AS pathway_id
"""


def password():
    if os.environ.get("NEO4J_PASSWORD"):
        return os.environ["NEO4J_PASSWORD"]
    return json.load(open(ROOT / ".mcp.json"))["mcpServers"]["neo4j-liver-kg"]["env"]["NEO4J_PASSWORD"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-file", required=True, help="run workbook with a B sheet (gives the drugs)")
    ap.add_argument("--drugs", required=True, help="drug file, e.g. only_consistent_drugs.xlsx")
    ap.add_argument("--genes", required=True, help="gene set file, e.g. 11_final_all_3845_genes.csv")
    ap.add_argument("--tier", default=None, help='keep genes whose Tier starts with this, e.g. "Tier 3"')
    ap.add_argument("--name", required=True, help="output name, e.g. Run4_onlyconsistent_tier3_BB")
    ap.add_argument("--out", default=str(ROOT / "Paper2_Metapath_runs"))
    ap.add_argument("--count-only", action="store_true", help="only print route counts")
    args = ap.parse_args()

    b = pd.read_excel(args.run_file, sheet_name="B")
    b = b.rename(columns=str.lower)
    b_drugs = set(b.drug.dropna())
    b_routes = b.groupby("drug").size()

    d = pd.read_excel(args.drugs)
    pairs = [dict(drug=r.Drug, drug_id=r.Drug_ID, gene=g.strip())
             for r in d.itertuples() if r.Drug in b_drugs
             for g in str(r.consistent_targets).split(";") if g.strip()]
    missing = b_drugs - {p["drug"] for p in pairs}
    if missing:
        raise SystemExit(f"B drugs not found in the drug file: {sorted(missing)}")

    g = pd.read_csv(args.genes)
    if args.tier:
        g = g[g.Tier.str.startswith(args.tier)]
    genes = g.set_index("Gene")
    print(f"drugs {len(b_drugs)} | drug-target pairs {len(pairs)} | Gene* {len(genes)}")

    drv = GraphDatabase.driver(os.environ.get("NEO4J_URI", "bolt://localhost:7687"),
                               auth=(os.environ.get("NEO4J_USERNAME", "neo4j"), password()))
    with drv.session(database=os.environ.get("NEO4J_DATABASE", "neo4j")) as ses:
        pw = pd.DataFrame(ses.run(PATHWAYS, nafld=NAFLD, smin=SIZE_MIN, smax=SIZE_MAX).data()).set_index("id")
        print(f"NAFLD pathways of size {SIZE_MIN}-{SIZE_MAX}: {len(pw)}")
        rows = []
        for i in range(0, len(pairs), 20):
            rows += ses.run(QUERY, pairs=pairs[i:i + 20], genes=list(genes.index), hub=HUB_MAX,
                            pathways=list(pw.index)).data()
    drv.close()

    r = pd.DataFrame(rows).drop_duplicates(["drug_id", "gene1", "gene2", "gene_star", "pathway_id"])
    print(f"routes {len(r)} | drugs with routes {r.drug.nunique() if len(r) else 0} of {len(b_drugs)}")
    if args.count_only or r.empty:
        return
    r["pathway"] = r.pathway_id.map(pw.name)
    r["pathway_size"] = r.pathway_id.map(pw["size"])
    r["gene_star_direction"] = r.gene_star.map(genes.Direction)
    r["stage"] = r.gene_star.map(genes.Largest_transition)
    r["disease"] = NAFLD_NAME
    r = r.sort_values(["drug", "gene1", "gene2", "gene_star", "pathway"])

    s = r.groupby("drug").size().rename("BB routes").reset_index()
    s["B routes"] = s.drug.map(b_routes).fillna(0).astype(int)
    s = s.sort_values(["BB routes", "drug"], ascending=[False, True]).reset_index(drop=True)
    s.insert(0, "No.", s.index + 1)
    s.insert(1, "Drug -> Disease", s.drug + " -> NAFLD")
    s.insert(3, "Disease", NAFLD_NAME)
    s = s.rename(columns={"drug": "Drug"})
    no = sorted(b_drugs - set(r.drug))
    if no:
        s = pd.concat([s, pd.DataFrame({"Drug -> Disease": [x + " -> NAFLD" for x in no], "Drug": no,
                                        "Disease": NAFLD_NAME, "BB routes": 0,
                                        "B routes": [int(b_routes.get(x, 0)) for x in no]})])
        s["No."] = range(1, len(s) + 1)

    simple = r[["drug", "gene1", "gene2", "gene_star", "pathway", "disease"]]
    simple.columns = ["Drug", "Gene 1", "Gene 2", "Gene*", "Pathway", "Disease"]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    r.to_csv(out / f"{args.name}.csv", index=False)
    with pd.ExcelWriter(out / f"{args.name}.xlsx") as xw:
        s.to_excel(xw, sheet_name="Summary", index=False)
        # an Excel sheet holds at most 1,048,576 rows: split long route lists over BB, BB (2), ...
        # (the CSV written above always holds every route in one file)
        step = 1_000_000
        for k, i in enumerate(range(0, max(len(simple), 1), step), 1):
            simple.iloc[i:i + step].to_excel(xw, sheet_name="BB" if k == 1 else f"BB ({k})", index=False)
    print(f"written: {out / (args.name + '.xlsx')}")


if __name__ == "__main__":
    main()
