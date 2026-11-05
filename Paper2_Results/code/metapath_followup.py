"""Paper 2 metapath follow-up: direct diseases, liver-treating drugs ranked by BB routes, stages.

For one metapath run (the BB workbook and route table made by scripts/run_two_gene_metapath.py),
makes the follow-up tables:

  <name>_direct_diseases.xlsx
      Summary      one row per drug with BB routes: indirect (BB) routes to NAFLD and every direct
                   TREATS / OFF_LABEL_FOR / CONTRAINDICATED_IN disease in the graph
      Direct links one row per drug-disease link
      Liver links  only links to liver diseases (Task 4 liver set: EFO "liver disease" subtree
                   minus non-hepatic-tissue entries)
  <name>_liver_treating_drugs_final.xlsx
      Final drugs  drugs with >=1 TREATS link to a liver disease and no CONTRAINDICATED_IN link to a
                   liver disease, ranked by BB routes (ties share the lowest rank)
      Removed drugs drugs that treat a liver disease but are contraindicated in one
  <name>_top_drugs_stages.xlsx  (final drugs ranked 1-25, ties included)
      Drug stages, Drug - Gene star detail, Stage totals: BB routes by the Paper 1 stage step
      (largest adjacent-stage change) of each route's Gene*.

Off-label use alone does not make a drug "liver-treating".
Read-only Neo4j. Password from NEO4J_PASSWORD (prompted if unset). Run from the repository root:

    python Paper2_Results/code/metapath_followup.py --bb Paper2_Results/Metapath_runs/<name>_BB --genes <gene file>
"""
import argparse
import base64
import getpass
import json
import os
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
LIVER_SET = ROOT / "Paper2_Results/Task4/liver_disease_set.csv"
STAGES = ["F0→F1", "F1→F2", "F2→F3", "F3→F4"]
TOP_RANK = 25
LINK = {"TREATS": "treats", "OFF_LABEL_FOR": "off-label for", "CONTRAINDICATED_IN": "contraindicated in"}


def cy(stmt, params, pw):
    h = {"Authorization": "Basic " + base64.b64encode(f"{os.environ.get('NEO4J_USERNAME', 'neo4j')}:{pw}".encode()).decode(),
         "Content-Type": "application/json", "Accept": "application/json"}
    body = json.dumps({"statement": stmt, "parameters": params, "accessMode": "Read"}).encode()
    r = json.load(urllib.request.urlopen(urllib.request.Request(
        os.environ.get("NEO4J_HTTP", "http://127.0.0.1:7474") + "/db/" + os.environ.get("NEO4J_DATABASE", "neo4j") + "/query/v2", data=body, headers=h,
        method="POST"), timeout=600))
    if r.get("errors"):
        raise RuntimeError(r["errors"])
    return pd.DataFrame(r["data"]["values"], columns=r["data"]["fields"])


def join(s):
    return "; ".join(sorted(set(s.dropna())))


ap = argparse.ArgumentParser()
ap.add_argument("--bb", required=True, help="path stem of the BB output, without .xlsx/.csv")
ap.add_argument("--genes", required=True, help="gene file used for the run (Gene, Direction, F0-F4, Largest_transition)")
args = ap.parse_args()
stem = Path(args.bb)
name = stem.name

summ = pd.read_excel(f"{stem}.xlsx", sheet_name="Summary")
bb = summ[summ["BB routes"] > 0][["Drug", "BB routes", "B routes"]]
routes = pd.read_csv(f"{stem}.csv")
ids = routes.drop_duplicates("drug")[["drug", "drug_id"]].rename(columns={"drug": "Drug", "drug_id": "Drug_ID"})
bb = bb.merge(ids, on="Drug", how="left")
assert bb.Drug_ID.notna().all(), "drug IDs missing for some BB drugs"

pw = os.environ.get("NEO4J_PASSWORD") or getpass.getpass("Neo4j password: ")
links = cy("""UNWIND $ids AS id MATCH (d:Drug {node_id: id})-[r:TREATS|OFF_LABEL_FOR|CONTRAINDICATED_IN]->(x:Disease)
RETURN id AS Drug_ID, type(r) AS rel, x.node_id AS disease_id, coalesce(x.name, x.node_id) AS Disease""",
           {"ids": list(bb.Drug_ID)}, pw)
liver = pd.read_csv(LIVER_SET)
liver_ids = set(liver.loc[~liver.excluded_non_hepatic, "disease_id"])
links = links.drop_duplicates(["Drug_ID", "rel", "disease_id"]).merge(bb[["Drug_ID", "Drug"]], on="Drug_ID")
links["Link"] = links.rel.map(LINK)
links["is_liver"] = links.disease_id.isin(liver_ids)

# ---------------------------------------------------------------- direct diseases workbook
g = links.groupby("Drug")
S = bb.rename(columns={"BB routes": "Indirect routes to NAFLD"}).copy()
for rel, col in [("TREATS", "Treats"), ("OFF_LABEL_FOR", "Off-label for"), ("CONTRAINDICATED_IN", "Contraindicated in")]:
    sub = links[links.rel == rel].groupby("Drug").Disease
    S[f"{col.split()[0]} (n)"] = S.Drug.map(sub.nunique()).fillna(0).astype(int)
    S[col] = S.Drug.map(sub.agg(join)).fillna("")
lv = links[links.is_liver]
S["Liver-related direct links"] = S.Drug.map(lv.assign(t=lv.Link + " " + lv.Disease).groupby("Drug").t.agg(join)).fillna("")
S = S.sort_values(["Indirect routes to NAFLD", "Drug"], ascending=[False, True]).reset_index(drop=True)
S.insert(0, "No.", S.index + 1)
S = S[["No.", "Drug", "Indirect routes to NAFLD", "Treats (n)", "Off-label (n)", "Contraindicated (n)",
       "Treats", "Off-label for", "Contraindicated in", "Liver-related direct links", "B routes"]]
direct = links.sort_values(["Drug", "Link", "Disease"])[["Drug", "Link", "Disease"]]
liver_links = lv.sort_values(["Drug", "Link", "Disease"])[["Drug", "Link", "Disease"]]
with pd.ExcelWriter(stem.parent / f"{name}_direct_diseases.xlsx") as xw:
    S.to_excel(xw, sheet_name="Summary", index=False)
    direct.to_excel(xw, sheet_name="Direct links", index=False)
    liver_links.to_excel(xw, sheet_name="Liver links", index=False)

# ---------------------------------------------------------------- liver-treating drugs, ranked
treats_liver = lv[lv.rel == "TREATS"].groupby("Drug").Disease.agg(join)
ci_liver = lv[lv.rel == "CONTRAINDICATED_IN"].groupby("Drug").Disease.agg(join)
T = bb[bb.Drug.isin(treats_liver.index)].copy()
T["Liver diseases it treats"] = T.Drug.map(treats_liver)
T["Liver contraindications"] = T.Drug.map(ci_liver)
final = T[T["Liver contraindications"].isna()].sort_values(["BB routes", "Drug"], ascending=[False, True])
final.insert(0, "Rank", final["BB routes"].rank(method="min", ascending=False).astype(int))
final = final.rename(columns={"BB routes": "BB routes to NAFLD"})[["Rank", "Drug", "BB routes to NAFLD",
                                                                     "Liver diseases it treats"]]
removed = T[T["Liver contraindications"].notna()].sort_values("BB routes", ascending=False).rename(
    columns={"BB routes": "BB routes to NAFLD"})[["Drug", "BB routes to NAFLD", "Liver diseases it treats",
                                                   "Liver contraindications"]]
with pd.ExcelWriter(stem.parent / f"{name}_liver_treating_drugs_final.xlsx") as xw:
    final.to_excel(xw, sheet_name=f"Final {len(final)} drugs", index=False)
    removed.to_excel(xw, sheet_name=f"Removed {len(removed)} drugs", index=False)

# ---------------------------------------------------------------- stages of the top drugs
genes = pd.read_csv(args.genes).set_index("Gene")
top = final[final.Rank <= TOP_RANK]
r = routes[routes.drug.isin(top.Drug)].copy()
r["stage"] = r.gene_star.map(genes.Largest_transition)
st = r.groupby(["drug", "stage"]).size().unstack(fill_value=0).reindex(columns=STAGES, fill_value=0)
D = top[["Rank", "Drug", "BB routes to NAFLD"]].rename(columns={"BB routes to NAFLD": "BB routes"}).set_index("Drug")
D["Gene* (n)"] = r.groupby("drug").gene_star.nunique()
for s_ in STAGES:
    D[f"{s_} routes"] = st[s_]
D["Main stage (most routes)"] = st.idxmax(axis=1)
D["% routes in main stage"] = (100 * st.max(axis=1) / st.sum(axis=1)).round(1)
D["All stages reached"] = st.apply(lambda row: ", ".join(c for c in STAGES if row[c] > 0), axis=1)
D["Liver disease it treats"] = top.set_index("Drug")["Liver diseases it treats"]
D = D.reset_index()[["Rank", "Drug", "BB routes", "Gene* (n)"] + [f"{s_} routes" for s_ in STAGES]
                    + ["Main stage (most routes)", "% routes in main stage", "All stages reached",
                       "Liver disease it treats"]]
gs = r.groupby(["drug", "gene_star"]).size().rename("BB routes").reset_index()
gs["Stage (biggest change)"] = gs.gene_star.map(genes.Largest_transition)
gs["Direction"] = gs.gene_star.map(genes.Direction).map({"increasing": "rises with fibrosis",
                                                         "decreasing": "falls with fibrosis"})
for c in ["F0", "F1", "F2", "F3", "F4"]:
    gs[c] = gs.gene_star.map(genes[c])
gs = gs.rename(columns={"drug": "Drug", "gene_star": "Gene*"})
gs["_rank"] = gs.Drug.map(top.set_index("Drug").Rank)
gs = gs.sort_values(["_rank", "Drug", "BB routes"], ascending=[True, True, False]).drop(columns="_rank")
tot = pd.DataFrame({"Stage": STAGES,
                    "Drugs reaching it": [int((st[s_] > 0).sum()) for s_ in STAGES],
                    "Routes": [int(st[s_].sum()) for s_ in STAGES],
                    "Main stage for (n drugs)": [int((D["Main stage (most routes)"] == s_).sum()) for s_ in STAGES]})
with pd.ExcelWriter(stem.parent / f"{name}_top_drugs_stages.xlsx") as xw:
    D.to_excel(xw, sheet_name="Drug stages", index=False)
    gs.to_excel(xw, sheet_name="Drug - Gene star detail", index=False)
    tot.to_excel(xw, sheet_name="Stage totals", index=False)

info = {"run": name, "drugs_with_BB_routes": len(bb), "drugs_treating_a_liver_disease": len(T),
        "removed_liver_contraindicated": len(removed), "final_drugs": len(final), "stage_table_drugs": len(top)}
(stem.parent / f"{name}_followup.json").write_text(json.dumps(info, indent=2))
print(json.dumps(info, indent=2))
print(final.head(15).to_string(index=False))
