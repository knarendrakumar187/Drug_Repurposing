"""Paper 2 Task 4: drug coverage of progression genes and repurposing candidate generation.

Main analysis uses pharmacological-target links only (same definition as Task 3: OptimusKG
actions other than ENZYME / TRANSPORTER / CARRIER, curated links added from ChEMBL, DGIdb,
GtoPdb and Open Targets, and OptimusKG links carrying a GtoPdb/DGIdb/ChEMBL action). Any-link counts are kept as secondary
columns. Repurposing information (approval, phase, withdrawal, existing disease uses)
comes from the knowledge graph only (project scope decision).

Outputs (Paper2_Results/Task4/):
  drug_coverage.csv, top50_by_coverage.csv, repurposing_candidates.csv,
  network_nodes.csv, network_edges.csv, network_top30_candidates.png, summary.md

This is candidate generation, not ranking: no drug is called the best candidate.
Present Neo4j graph only (read-only, no CSV fallback); password from NEO4J_PASSWORD (prompted if unset).
Run from the repository root:

    python Paper2_Results/code/task4_drug_coverage.py
"""
import base64
import getpass
import json
import os
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Paper2_Results" / "Task4"
OUT.mkdir(parents=True, exist_ok=True)

GENES = ROOT / "Paper2_Results/Task1/Paper2_progression_genes.csv"
EDGES = ROOT / "Paper2_Results/Task2/drug_gene_edges.csv"
TERMS = ROOT / "Paper1_Results/Task4/task4_all_significant.csv"
# Indications of the drugs added on 21 Oct (they have no TREATS links): made by task4_added_drug_indications.py
ADDED_IND = ROOT / "Paper2_Results/Task4/added_drug_indications_chembl.csv"
ADDED_IND_META = ROOT / "Paper2_Results/Task4/added_drug_indications_chembl.json"
# Name fallback for ChEMBL indications whose EFO/MONDO ID is not a graph liver disease
FIB_LIVER_NAME = (r"fatty liver|steatohepatitis|cirrhosis|(liver|hepatic) fibrosis|portal hypertension|"
                  r"(liver|hepatic) failure|NAFLD|NASH|MASLD|MASH")

PK_ACTIONS = ["ENZYME", "TRANSPORTER", "CARRIER"]
EXTERNAL_LAYER = "external (Idea 4)"   # links added from ChEMBL/DGIdb/GtoPdb/Open Targets
LIVER_ROOT = "EFO_0001421"          # "liver disease"
# The EFO "liver disease" subtree also holds non-hepatic-tissue diseases (malaria with a liver
# stage, gallbladder / ampulla / extrahepatic bile-duct tumours, erythropoietic porphyrias).
# They are removed so that e.g. antimalarials are not mistaken for existing liver drugs.
LIVER_EXCLUDE = r"malaria|blackwater|gallbladder|gall bladder|ampulla|extrahepatic|bladder hepatoid|(?<!hepato)erythropoietic"
# Only drugs already used for a FIBROSIS-TYPE liver disease are excluded from the candidates.
# Drugs used only for liver cancer, viral hepatitis, bile-duct cancer, etc. are kept, because
# using them for fibrosis would still be repurposing. Same fibrosis-spectrum definition as Task 5.
FIB_SPECTRUM = r"fatty liver|steatohepatitis|cirrhosis|fibrosis|portal hypertension|liver failure"
MIN_PHASE = 2
TOP_N_TABLE = 50
TOP_N_NETWORK = 30

NEO4J_URL = os.environ.get("NEO4J_HTTP", "http://127.0.0.1:7474")
NEO4J_DB = os.environ.get("NEO4J_DATABASE", "neo4j")
NEO4J_USER = os.environ.get("NEO4J_USERNAME", "neo4j")


def neo4j_query(statement, params, password):
    headers = {"Authorization": "Basic " + base64.b64encode(f"{NEO4J_USER}:{password}".encode()).decode(),
               "Content-Type": "application/json", "Accept": "application/json"}
    body = json.dumps({"statement": statement, "parameters": params, "accessMode": "Read"}).encode()
    req = urllib.request.Request(f"{NEO4J_URL}/db/{NEO4J_DB}/query/v2", data=body, headers=headers, method="POST")
    res = json.load(urllib.request.urlopen(req, timeout=600))
    if res.get("errors"):
        raise RuntimeError(res["errors"])
    return pd.DataFrame(res["data"]["values"], columns=res["data"]["fields"])


# ------------------------------------------------------------------ inputs
genes = pd.read_csv(GENES)
edges = pd.read_csv(EDGES, low_memory=False)
edges["pharm"] = ((edges.Link_origin != "OptimusKG") | edges.external_action_direction.notna()
                  | (~edges.Action_raw.isin(PK_ACTIONS) & (edges.Action_raw != "unknown")))
pe = edges[edges.pharm]
drug_ids = sorted(edges.Drug_ID.unique())

pw = os.environ.get("NEO4J_PASSWORD") or getpass.getpass("Neo4j password: ")

totals = neo4j_query("""
UNWIND $ids AS id
MATCH (d:Drug {node_id: id})-[r:ACTS_ON]->(g:Gene)
RETURN id AS Drug_ID, count(DISTINCT g) AS all_targets_graph,
       count(DISTINCT CASE WHEN r.layer = $layer OR r.external_action_direction IS NOT NULL
                                OR (r.action IS NOT NULL AND NOT r.action IN $pk) THEN g END) AS pharm_targets_graph
""", {"ids": drug_ids, "pk": PK_ACTIONS, "layer": EXTERNAL_LAYER}, pw)

liver = neo4j_query("""
MATCH (:Disease {node_id: $root})-[:PARENT_OF*0..20]->(c:Disease)
RETURN DISTINCT c.node_id AS disease_id, c.name AS disease_name
""", {"root": LIVER_ROOT}, pw)
liver["disease_name"] = liver.disease_name.fillna(liver.disease_id)
liver["excluded_non_hepatic"] = liver.disease_name.str.contains(LIVER_EXCLUDE, case=False, regex=True)
liver["fibrosis_type"] = ~liver.excluded_non_hepatic & liver.disease_name.str.contains(FIB_SPECTRUM, case=False, regex=True)
liver.sort_values(["excluded_non_hepatic", "disease_name"]).to_csv(OUT / "liver_disease_set.csv", index=False)
liver_ids = set(liver.loc[~liver.excluded_non_hepatic, "disease_id"])
fib_liver_ids = set(liver.loc[liver.fibrosis_type, "disease_id"])

uses = neo4j_query("""
UNWIND $ids AS id
MATCH (d:Drug {node_id: id})-[r:TREATS|OFF_LABEL_FOR|CONTRAINDICATED_IN]->(x:Disease)
RETURN id AS Drug_ID, type(r) AS use_type, x.node_id AS disease_id, x.name AS disease_name
""", {"ids": drug_ids}, pw)
uses["use_source"] = "graph"
# ChEMBL indications for the drugs added on 21 Oct (any clinical phase counts as an existing use)
ind = pd.read_csv(ADDED_IND)
ind_meta = json.loads(ADDED_IND_META.read_text(encoding="utf-8"))
ind = ind[ind.Drug_ID.isin(drug_ids)]
uses = pd.concat([uses, pd.DataFrame({
    "Drug_ID": ind.Drug_ID, "use_type": "TREATS", "disease_id": ind.efo_id_graph,
    "disease_name": ind.efo_term.fillna(ind.mesh_heading),
    "use_source": "ChEMBL indication (phase " + ind.max_phase_for_ind.map(lambda x: f"{x:g}") + ")"})],
    ignore_index=True)
uses["disease_name"] = uses.disease_name.fillna(uses.disease_id)
chembl_use = uses.use_source.ne("graph")
name_fib = uses.disease_name.str.contains(FIB_LIVER_NAME, case=False, regex=True)
uses["is_liver"] = uses.disease_id.isin(liver_ids) | (chembl_use & name_fib)
uses["is_fib_liver"] = uses.disease_id.isin(fib_liver_ids) | (chembl_use & name_fib)

# ------------------------------------------------------------------ main Paper 1 pathway per gene
DB_PREF = {"Hallmark": 0, "Reactome": 1, "KEGG": 2, "GO_BP": 3}
t = pd.read_csv(TERMS, encoding="utf-8-sig")
tl = t.assign(gene=t.genes.str.split(r",\s*")).explode("gene")
tl["pref"] = tl.database.map(DB_PREF)
tl = tl.sort_values(["pref", "q_value"]).drop_duplicates("gene")
main_pw = dict(zip(tl.gene, tl.database + ": " + tl.term_label))
genes["main_paper1_pathway"] = genes.Gene.map(main_pw)

# ------------------------------------------------------------------ per-drug coverage
GROUP_COLS = {
    "fibrosis_only_positive": "fib_only_pos", "fibrosis_only_negative": "fib_only_neg",
    "shared_positive": "shared_pos", "shared_negative": "shared_neg",
    "shared_discordant_fib_negative": "shared_discordant",
}


def per_drug(e):
    g = e.groupby("Drug_ID")
    out = pd.DataFrame({
        "n_progression_genes": g.Gene_ID.nunique(),
        "n_fibrosis_genes": e[e.progression_class != "inflammation_only"].groupby("Drug_ID").Gene_ID.nunique(),
        "fibrosis_only_targets": e[e.progression_class == "fibrosis_only"].groupby("Drug_ID").Gene_ID.nunique(),
        "shared_targets": e[e.progression_class == "shared"].groupby("Drug_ID").Gene_ID.nunique(),
        "inflammation_only_targets": e[e.progression_class == "inflammation_only"].groupby("Drug_ID").Gene_ID.nunique(),
    })
    for grp, col in GROUP_COLS.items():
        out[col] = e[e.paper1_group == grp].groupby("Drug_ID").Gene_ID.nunique()
    for c in ["C1", "C2", "C3", "C4", "C5"]:
        out[c] = e[e.cluster == c].groupby("Drug_ID").Gene_ID.nunique()
    for d in ["down", "up", "mixed", "unknown"]:
        out[f"n_{d}_actions"] = e[e.Action_direction == d].groupby("Drug_ID").size()
    out["mechanism_annotated"] = g.apply(lambda s: bool((s.Action_direction != "unknown").any()
                                                        or s.mechanism_of_action.notna().any()), include_groups=False)
    out["shortlist_targets"] = e[e.in_shortlist].groupby("Drug_ID").Gene.agg(lambda s: ";".join(sorted(set(s))))
    out["target_genes"] = g.Gene.agg(lambda s: ";".join(sorted(set(s))))
    return out


main = per_drug(pe)
anyl = edges.groupby("Drug_ID").Gene_ID.nunique().rename("n_progression_genes_any_link")
info = (edges.sort_values("Drug_ID").drop_duplicates("Drug_ID")
        .set_index("Drug_ID")[["Drug", "molecule_type", "max_phase", "is_approved", "withdrawn", "drug_stage"]])
# Drugs added on 21 Oct carry drug_stage ("approved" / "clinical trial") instead of is_approved / has_been_withdrawn.
info["added_drug_node"] = info.is_approved.isna() & info.drug_stage.notna()
info["is_approved"] = info.is_approved.astype("boolean").fillna(info.drug_stage.eq("approved")).astype(bool)
# "approved" means phase 4; added drugs recorded only as approved (RxNorm/GtoPdb/NCIt IDs) have no max_phase.
info["phase_from_stage"] = info.added_drug_node & info.max_phase.isna() & info.drug_stage.eq("approved")
info.loc[info.phase_from_stage, "max_phase"] = 4.0
info = info.drop(columns="drug_stage")

cov = info.join(main, how="left").join(anyl).join(totals.set_index("Drug_ID"))
num = [c for c in cov.columns if c not in ("Drug", "molecule_type", "max_phase", "is_approved", "withdrawn", "added_drug_node", "phase_from_stage",
                                           "shortlist_targets", "target_genes", "mechanism_annotated")]
cov["mechanism_annotated"] = cov.mechanism_annotated.fillna(False).astype(bool)
cov[num] = cov[num].fillna(0).astype(int)
cov["progression_target_fraction"] = np.where(cov.pharm_targets_graph > 0,
                                              cov.n_progression_genes / cov.pharm_targets_graph.replace(0, np.nan), np.nan)
cov["progression_target_fraction"] = cov.progression_target_fraction.round(3)
cov["n_shortlist_targets"] = cov.shortlist_targets.fillna("").map(lambda s: len([x for x in s.split(";") if x]))

# disease uses
u = uses.groupby("Drug_ID")
cov["n_treats"] = uses[uses.use_type == "TREATS"].groupby("Drug_ID").disease_id.nunique()
cov["has_liver_indication"] = uses[uses.is_liver & uses.use_type.isin(["TREATS", "OFF_LABEL_FOR"])].groupby("Drug_ID").size() > 0
cov["liver_indications"] = (uses[uses.is_liver & uses.use_type.isin(["TREATS", "OFF_LABEL_FOR"])]
                            .groupby("Drug_ID").disease_name.agg(lambda s: "; ".join(sorted(set(s)))))
use_ok = uses.use_type.isin(["TREATS", "OFF_LABEL_FOR"])
cov["has_fibrosis_type_liver_use"] = uses[uses.is_fib_liver & use_ok].groupby("Drug_ID").size() > 0
cov["fibrosis_type_liver_uses"] = (uses[uses.is_fib_liver & use_ok].groupby("Drug_ID").disease_name
                                   .agg(lambda s: "; ".join(sorted(set(s)))))
cov["liver_contraindicated"] = uses[uses.is_liver & (uses.use_type == "CONTRAINDICATED_IN")].groupby("Drug_ID").size() > 0
cov["liver_contraindications"] = (uses[uses.is_liver & (uses.use_type == "CONTRAINDICATED_IN")]
                                  .groupby("Drug_ID").disease_name.agg(lambda s: "; ".join(sorted(set(s)))))
cov["current_indications_examples"] = (uses[uses.use_type == "TREATS"].groupby("Drug_ID").disease_name
                                       .agg(lambda s: "; ".join(sorted(set(s))[:5])))
cov["n_treats"] = cov.n_treats.fillna(0).astype(int)
for c in ["has_liver_indication", "has_fibrosis_type_liver_use", "liver_contraindicated"]:
    cov[c] = cov[c].astype("boolean").fillna(False).astype(bool)

cov["kept_with_other_liver_use"] = cov.has_liver_indication & ~cov.has_fibrosis_type_liver_use
cov = cov.reset_index().rename(columns={"index": "Drug_ID"})
order = ["Drug", "Drug_ID", "molecule_type", "max_phase", "is_approved", "withdrawn", "added_drug_node", "phase_from_stage",
         "n_progression_genes", "n_fibrosis_genes", "fibrosis_only_targets", "shared_targets", "inflammation_only_targets",
         "fib_only_pos", "fib_only_neg", "shared_pos", "shared_neg", "shared_discordant",
         "C1", "C2", "C3", "C4", "C5", "mechanism_annotated", "n_down_actions", "n_up_actions", "n_mixed_actions", "n_unknown_actions",
         "n_shortlist_targets", "shortlist_targets", "pharm_targets_graph", "all_targets_graph",
         "progression_target_fraction", "n_progression_genes_any_link",
         "n_treats", "has_liver_indication", "liver_indications", "has_fibrosis_type_liver_use", "fibrosis_type_liver_uses",
         "kept_with_other_liver_use", "liver_contraindicated", "liver_contraindications",
         "current_indications_examples", "target_genes"]
cov = cov[order].sort_values(["n_progression_genes", "n_fibrosis_genes", "progression_target_fraction"],
                             ascending=False)
cov.to_csv(OUT / "drug_coverage.csv", index=False)
top50 = cov[cov.n_progression_genes > 0].head(TOP_N_TABLE)
top50.to_csv(OUT / "top50_by_coverage.csv", index=False)

# ------------------------------------------------------------------ repurposing candidates
cand = cov[(cov.max_phase >= MIN_PHASE) & (~cov.withdrawn.fillna(False).astype(bool)) &
           (cov.n_progression_genes >= 1) & (~cov.has_fibrosis_type_liver_use)].copy()
cand = cand.sort_values(["n_fibrosis_genes", "progression_target_fraction", "n_progression_genes"], ascending=False)
cand.insert(0, "candidate_order", range(1, len(cand) + 1))
cand.to_csv(OUT / "repurposing_candidates.csv", index=False)
mech = cand[cand.mechanism_annotated].copy()
mech.insert(1, "mechanism_order", range(1, len(mech) + 1))
mech.to_csv(OUT / "repurposing_candidates_mechanism_annotated.csv", index=False)

# ------------------------------------------------------------------ network (top 30 candidates)
top = mech.head(TOP_N_NETWORK)   # mechanism-annotated candidates: interpretable drug-target biology
ne = pe[pe.Drug_ID.isin(top.Drug_ID)]
gann = genes.set_index("ensembl_id")
edges_out = pd.DataFrame({
    "source": ne.Drug_ID, "target": ne.Gene_ID, "source_name": ne.Drug, "target_name": ne.Gene,
    "action": ne.Action_raw, "action_direction": ne.Action_direction, "evidence_level": ne.evidence_level,
    "link_origin": ne.Link_origin, "gene_class": ne.progression_class, "gene_direction": ne.gene_direction,
    "gene_cluster": ne.cluster,
})
edges_out.to_csv(OUT / "network_edges.csv", index=False)
gn = gann.loc[ne.Gene_ID.unique()]
drug_dir = gn.reset_index()
node_genes = pd.DataFrame({
    "id": gn.index, "label": gn.Gene, "node_type": "gene", "class": gn.progression_class,
    "direction": gn.fibrosis_direction.fillna(gn.inflammation_direction), "cluster": gn.trajectory_cluster,
    "main_paper1_pathway": gn.main_paper1_pathway, "in_shortlist": gn.in_shortlist, "is_approved": np.nan,
    "n_progression_genes": np.nan,
})
node_drugs = pd.DataFrame({
    "id": top.Drug_ID, "label": top.Drug, "node_type": "drug", "class": "drug", "direction": np.nan, "cluster": np.nan,
    "main_paper1_pathway": np.nan, "in_shortlist": False, "is_approved": top.is_approved,
    "n_progression_genes": top.n_progression_genes,
})
nodes_out = pd.concat([node_drugs, node_genes], ignore_index=True)
nodes_out.to_csv(OUT / "network_nodes.csv", index=False)

COL = {"fibrosis_only": "#c0392b", "shared": "#7d3c98", "inflammation_only": "#2471a3", "drug": "#d4ac0d"}
ECOL = {"down": "#2471a3", "up": "#c0392b", "mixed": "#d4ac0d", "unknown": "#c5ccd3"}
# Picture only (CSV files stay complete): fibrosis-associated genes, and drugs with > PIC_MAX_GENES
# such genes are listed in a note instead of drawn, so the remaining structure is readable.
PIC_MAX_GENES = 40
fib_gene_ids = set(genes.loc[genes.progression_class != "inflammation_only", "ensembl_id"])
pic_edges = edges_out[edges_out.target.isin(fib_gene_ids)]
big = pic_edges.groupby("source").target.nunique()
big = big[big > PIC_MAX_GENES]
pic_edges = pic_edges[~pic_edges.source.isin(big.index)]
pic_nodes = nodes_out[nodes_out.id.isin(set(pic_edges.source) | set(pic_edges.target))]
G = nx.Graph()
for _, r in pic_nodes.iterrows():
    G.add_node(r.id, **{k: r[k] for k in ["label", "node_type", "class", "in_shortlist"]})
for _, r in pic_edges.iterrows():
    G.add_edge(r.source, r.target, direction=r.action_direction)
pos = nx.kamada_kawai_layout(G, pos=nx.spring_layout(G, k=1.0, seed=42, iterations=300))
fig, ax = plt.subplots(figsize=(17, 14))
for d, c in ECOL.items():
    es = [(u, v) for u, v, a in G.edges(data=True) if a["direction"] == d]
    nx.draw_networkx_edges(G, pos, edgelist=es, edge_color=c, width=1.6 if d in ("down", "up") else 0.6,
                           alpha=0.85 if d in ("down", "up") else 0.5, ax=ax)
gnodes = [n for n, a in G.nodes(data=True) if a["node_type"] == "gene"]
dnodes = [n for n, a in G.nodes(data=True) if a["node_type"] == "drug"]
nx.draw_networkx_nodes(G, pos, nodelist=gnodes, node_size=[180 if G.nodes[n]["in_shortlist"] else 60 for n in gnodes],
                       node_color=[COL[G.nodes[n]["class"]] for n in gnodes],
                       edgecolors=["#000000" if G.nodes[n]["in_shortlist"] else "none" for n in gnodes],
                       linewidths=1.5, ax=ax)
nx.draw_networkx_nodes(G, pos, nodelist=dnodes, node_shape="s", node_size=420, node_color=COL["drug"],
                       edgecolors="#7d6608", ax=ax)
nx.draw_networkx_labels(G, pos, labels={n: G.nodes[n]["label"] for n in dnodes}, font_size=8, font_weight="bold", ax=ax)
hub = {n for n in gnodes if G.degree(n) >= 3 or G.nodes[n]["in_shortlist"]}
nx.draw_networkx_labels(G, pos, labels={n: G.nodes[n]["label"] for n in hub}, font_size=7, font_color="#1f2933", ax=ax)
handles = [plt.Line2D([], [], marker="s", ls="", color=COL["drug"], ms=10, label="Drug (candidate)"),
           plt.Line2D([], [], marker="o", ls="", color=COL["fibrosis_only"], ms=8, label="Fibrosis-only gene"),
           plt.Line2D([], [], marker="o", ls="", color=COL["shared"], ms=8, label="Shared gene"),
           plt.Line2D([], [], marker="o", ls="", color=COL["inflammation_only"], ms=8, label="Disease-group-only gene"),
           plt.Line2D([], [], marker="o", ls="", color="white", mec="black", ms=9, label="Paper 1 shortlist gene"),
           plt.Line2D([], [], color=ECOL["down"], lw=2, label="Drug pushes gene ⬇"),
           plt.Line2D([], [], color=ECOL["up"], lw=2, label="Drug pushes gene ⬆"),
           plt.Line2D([], [], color=ECOL["unknown"], lw=2, label="Direction unknown")]
ax.legend(handles=handles, loc="lower left", fontsize=9, frameon=False)
ax.set_title(f"Top {TOP_N_NETWORK} mechanism-annotated repurposing candidates → fibrosis-associated gene targets "
             "(pharmacological links; present knowledge graph)", fontsize=13, fontweight="bold", loc="left")
if len(big):
    names = top.set_index("Drug_ID").Drug
    ax.text(0.99, 0.01, "Not drawn (too many targets to show): " + ", ".join(
        f"{names[d]} ({n} fibrosis genes)" for d, n in big.items()) + ".\nFull network incl. disease-group-only genes: "
        "network_nodes.csv / network_edges.csv.", transform=ax.transAxes, ha="right", fontsize=9, color="#616e7c")
ax.axis("off")
fig.savefig(OUT / "network_top30_candidates.png", dpi=150, bbox_inches="tight")
plt.close(fig)


# ------------------------------------------------------------------ summary
def md(df, cols, fmt=None):
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join("" if pd.isna(r[c]) else str(r[c]) for c in cols) + " |")
    return "\n".join(out)


drugs_hit = cov[cov.n_progression_genes > 0]
phase_tab = drugs_hit.max_phase.fillna(-9).value_counts().sort_index(ascending=False)
filters = [
    ("Drugs with ≥1 pharmacological link to a progression gene", len(drugs_hit)),
    (f"…with max_phase ≥ {MIN_PHASE}", int((drugs_hit.max_phase >= MIN_PHASE).sum())),
    ("…and not withdrawn", int(((drugs_hit.max_phase >= MIN_PHASE) & ~drugs_hit.withdrawn.fillna(False).astype(bool)).sum())),
    ("…and no existing FIBROSIS-TYPE liver-disease TREATS/OFF_LABEL_FOR use → **repurposing candidates**", len(cand)),
    ("…of which mechanism-annotated (≥1 link with known action or mechanism text) → **main list**", len(mech)),
]
n_liver_ind = int(drugs_hit.has_liver_indication.sum())
n_fib_ind = int(drugs_hit.has_fibrosis_type_liver_use.sum())
elig = (drugs_hit.max_phase >= MIN_PHASE) & ~drugs_hit.withdrawn.fillna(False).astype(bool)
n_old_excl = int((elig & drugs_hit.has_liver_indication).sum())
n_back = int(cand.kept_with_other_liver_use.sum())
back_show = cand[cand.kept_with_other_liver_use].head(20).copy()
back_show["phase"] = back_show.max_phase.map(lambda x: f"{x:g}")
mech_show = mech.head(20).copy()
mech_show["phase"] = mech_show.max_phase.map(lambda x: f"{x:g}")
mech_show["frac"] = mech_show.progression_target_fraction.map(lambda x: f"{x:.2f}")
mech_show["liver CI"] = mech_show.liver_contraindicated.map({True: "⚠️ yes", False: ""})
mech_show["now used for"] = mech_show.current_indications_examples.fillna("").map(lambda x: "; ".join(x.split("; ")[:3]))
top_show = cand.head(20).copy()
top_show["phase"] = top_show.max_phase.map(lambda x: f"{x:g}")
top_show["frac"] = top_show.progression_target_fraction.map(lambda x: f"{x:.2f}")
top_show["liver CI"] = top_show.liver_contraindicated.map({True: "⚠️ yes", False: ""})
t50 = top50.head(15).copy()
t50["frac"] = t50.progression_target_fraction.map(lambda x: f"{x:.2f}")
t50["liver use"] = t50.has_liver_indication.map({True: "yes", False: ""})
with_short = cand[cand.n_shortlist_targets > 0]
n_added_cand = int(cand.added_drug_node.sum())
added_all = cov[cov.added_drug_node & (cov.n_progression_genes > 0)]
added_rows = [
    ("Added drugs with ≥1 pharmacological progression-gene link", len(added_all)),
    ("…with a ChEMBL ID (graph ID or exact name match)", int(added_all.Drug_ID.isin(pd.read_csv(ADDED_IND.with_name("added_drug_chembl_ids.csv")).dropna(subset=["chembl_id"]).Drug_ID).sum())),
    ("…with ≥1 ChEMBL indication", int(added_all.Drug_ID.isin(ind.Drug_ID).sum())),
    ("…with a liver-disease indication", int(added_all.has_liver_indication.sum())),
    ("…with a fibrosis-type liver indication → excluded", int(added_all.has_fibrosis_type_liver_use.sum())),
    ("…approved with no phase number, treated as phase 4", int(added_all.phase_from_stage.sum())),
    ("…clinical-trial with no phase number (cannot confirm phase ≥ 2, excluded)", int((added_all.max_phase.isna()).sum())),
    ("Added drugs among the final candidates", n_added_cand),
]
known_liver = drugs_hit[drugs_hit.has_fibrosis_type_liver_use].sort_values("n_fibrosis_genes", ascending=False).head(10)

summary = f"""# Paper 2 — Task 4: drug coverage and repurposing candidates

Gene set: the corrected Paper 1 union of **{len(genes):,} genes** (Task 1). "Inflammation" columns hold Paper 1's ordered disease-group axis.
**Scope:** knowledge graph only (project scope decision: repurposing, not reversal). Main analysis uses **pharmacological-target links** (same definition as Task 3); any-link gene counts are in `n_progression_genes_any_link`.
Script: `Paper2_Results/code/task4_drug_coverage.py`.

**This is candidate generation only.** Coverage is not a therapeutic score. Broad, heavily annotated drugs (e.g. multi-kinase inhibitors) touch many genes, so `progression_target_fraction` (progression targets / all pharmacological targets in the graph) is reported alongside coverage. No drug here is "the best candidate".

**Exclusion rule:** only drugs already used for a fibrosis-type liver disease are excluded. Fibrosis-type = liver diseases whose name matches *fatty liver, steatohepatitis, cirrhosis, fibrosis, portal hypertension, liver failure* (the Task 5 fibrosis-spectrum definition; {len(fib_liver_ids)} of {len(liver_ids)} liver diseases, flagged in `liver_disease_set.csv`). Drugs used only for liver cancer, viral hepatitis, bile-duct cancer or other liver diseases are now **kept**: using them for fibrosis is still repurposing.

## From drugs to repurposing candidates

{md(pd.DataFrame(filters, columns=["Step", "Drugs"]), ["Step", "Drugs"])}

- Liver disease = `{LIVER_ROOT}` "liver disease" + all `PARENT_OF` descendants, **minus {int(liver.excluded_non_hepatic.sum())} non-hepatic-tissue entries** that EFO files under it (malaria, gallbladder / ampulla / extrahepatic bile-duct tumours, erythropoietic porphyrias) = **{len(liver_ids)} liver diseases**. Both kept and removed entries are listed in `liver_disease_set.csv`. Without this step, antimalarials such as those for *Plasmodium* malaria would be wrongly treated as existing liver drugs.
- {n_liver_ind:,} drugs hitting progression genes have some liver-disease TREATS/OFF_LABEL_FOR use; only the **{n_fib_ind:,}** with a fibrosis-type use are excluded now (reference set in `drug_coverage.csv`).
- {int(cand.liver_contraindicated.sum()):,} candidates are **contraindicated** in a liver disease. They are kept and flagged (`liver_contraindicated`), as requested.
- {len(with_short):,} candidates hit at least one gene in the Paper 1 23-gene shortlist.

## Drugs added to the graph on 21 Oct 2026

These {len(added_all):,} drugs came with drug–gene links only (no TREATS / OFF_LABEL_FOR links). Their existing uses were taken from ChEMBL indications ({ind_meta["chembl_release"]}, released {ind_meta["chembl_release_date"]}, retrieved {ind_meta["retrieved"]}; `added_drug_indications_chembl.csv`), matched by ChEMBL ID or, for RxNorm/GtoPdb/NCIt IDs, by exact name (`added_drug_chembl_ids.csv`). An indication at any clinical phase counts as an existing use. Liver diseases are matched by graph disease ID, or by name for IDs outside the graph's liver set.

{md(pd.DataFrame(added_rows, columns=["Step", "Drugs"]), ["Step", "Drugs"])}

## Kept although they have a non-fibrosis liver use (top 20 by fibrosis-gene coverage)

{md(back_show, ["candidate_order", "Drug", "phase", "n_fibrosis_genes", "n_progression_genes", "mechanism_annotated", "liver_indications"])}

## Main list: top 20 mechanism-annotated repurposing candidates

Ordered by fibrosis-gene coverage, then progression-target fraction. "Mechanism-annotated" means the graph says *how* the drug acts on at least one of its progression-gene targets (inhibitor, agonist, … or a mechanism-of-action text). These are the interpretable candidates.

{md(mech_show, ["mechanism_order", "Drug", "phase", "n_fibrosis_genes", "n_progression_genes", "pharm_targets_graph", "frac", "n_down_actions", "n_up_actions", "shortlist_targets", "liver CI", "now used for"])}

Full list: `repurposing_candidates_mechanism_annotated.csv` ({len(mech):,} drugs).

## All candidates, including drugs with no mechanism annotation

Drugs without any mechanism annotation (often nutrients, metal ions or cofactors that OptimusKG lists as binding many enzymes with the generic label "TARGET") are kept for completeness, but are not meaningful repurposing hypotheses without mechanism evidence.

{md(top_show, ["candidate_order", "Drug", "phase", "n_fibrosis_genes", "fib_only_pos", "fib_only_neg", "shared_pos", "shared_neg", "n_progression_genes", "pharm_targets_graph", "frac", "n_down_actions", "n_up_actions", "shortlist_targets", "liver CI"])}

Full list: `repurposing_candidates.csv` ({len(cand):,} drugs, with current indications, target genes, cluster counts and action counts).

## Top 15 drugs by raw progression-gene coverage (all drugs, incl. existing liver drugs)

{md(t50, ["Drug", "max_phase", "n_progression_genes", "n_fibrosis_genes", "inflammation_only_targets", "pharm_targets_graph", "frac", "liver use"])}

Full top 50: `top50_by_coverage.csv`.

## Drugs that already have a fibrosis-type liver use (reference, not candidates)

{md(known_liver, ["Drug", "max_phase", "n_fibrosis_genes", "n_progression_genes", "fibrosis_type_liver_uses"])}

## Network
`network_nodes.csv` / `network_edges.csv` (Cytoscape-ready): top {TOP_N_NETWORK} candidates → their progression-gene targets ({len(nodes_out) - len(top)} genes, {len(edges_out)} links). Gene attributes: class, direction, trajectory cluster, main Paper 1 pathway (Hallmark > Reactome > KEGG > GO, most significant), shortlist flag. Drug attribute: is_approved. Edge attribute: action direction.

![Network](network_top30_candidates.png)

## Caveats
- Coverage counts known drug–target annotations, which favour well-studied drugs and genes.
- Most links have unknown push direction (Task 2). Whether a candidate's action opposes the fibrosis direction is checked in Task 7.
- "No fibrosis-type liver use" means none is recorded in the graph (or, for added drugs, in ChEMBL), not necessarily none in practice.
- Candidates are **repurposing hypotheses**, not treatments.
"""
(OUT / "summary.md").write_text(summary, encoding="utf-8")
print(summary)
