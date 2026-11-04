"""Graph preparation, step 0.5: drug action directions for OptimusKG links that have none.

Run on 21 Oct 2026. The OptimusKG layer gives no push direction for 7,919 of its 8,948 drug–gene
links on the progression genes (TARGET / ENZYME / TRANSPORTER / CARRIER). This looks each of those links
up in three free public resources and measures how many get a direction:

  1. ChEMBL drug mechanisms (REST API)       matched by drug ChEMBL ID (and its parent
                                            compound) + target gene symbol
  2. IUPHAR/BPS Guide to Pharmacology        matched by ligand ChEMBL ID + target Ensembl ID
     (public interactions.csv download)
  3. DGIdb 5 (GraphQL API; aggregates GtoPdb, TTD, Drug Repurposing Hub, DrugBank-free
     sources, ChEMBL …)                     matched by gene symbol + drug ChEMBL ID or name

The 297 links that get a clear, agreeing direction are written to the graph as external_action_*
properties by write_external_directions_to_graph.py. Outputs: Paper2_Results/External_action_check/.
Raw downloads are cached in that folder's raw/ subfolder and only read as data.

    python Paper2_Results/code/external_action_check.py
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "Paper2_Results"
OUT = RES / "External_action_check"
RAW = OUT / "raw"
RAW.mkdir(parents=True, exist_ok=True)

DOWN_WORDS = ("INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE", "INVERSE AGONIST", "INVERSE_AGONIST", "DEGRADER",
              "INHIBITION", "INHIBITORY", "ANTISENSE", "RNAI", "DISRUPTING", "SUPPRESSOR", "CHANNEL BLOCKER",
              "GATING INHIBITOR", "ALLOSTERIC ANTAGONIST")
UP_WORDS = ("AGONIST", "ACTIVATOR", "POSITIVE", "OPENER", "STABILISER", "STABILIZER", "INDUCER", "ACTIVATION",
            "ACTIVATING", "POTENTIATOR", "ENHANCER")


def direction(label):
    s = str(label).upper().replace("_", " ")
    if not s or s in ("NAN", "NONE"):
        return None
    if any(w in s for w in DOWN_WORDS):          # checked first: 'inverse agonist', 'negative modulator'
        return "down"
    if any(w in s for w in UP_WORDS):
        return "up"
    return None                                  # binding, substrate, modulator, unknown …


def get_json(url, cache_name, data=None):
    f = RAW / cache_name
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))
    req = urllib.request.Request(url, data=data, headers={"Accept": "application/json", "Content-Type": "application/json"})
    for attempt in range(4):
        try:
            js = json.load(urllib.request.urlopen(req, timeout=120))
            f.write_text(json.dumps(js), encoding="utf-8")
            time.sleep(0.2)
            return js
        except Exception:
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"failed: {url[:120]}")


# ------------------------------------------------------------------ the unknown links
e = pd.read_csv(RES / "Task2/drug_gene_edges.csv")
unk = e[e.Action_direction == "unknown"].copy()
drugs = sorted(e.Drug_ID.unique())          # all drugs: unknown links + the positive control
gene_sym = dict(zip(e.Gene_ID, e.Gene))
print(f"unknown links {len(unk):,}; drugs {len(drugs):,}; genes {unk.Gene_ID.nunique():,}")
batches = [drugs[i:i + 40] for i in range(0, len(drugs), 40)]

# ------------------------------------------------------------------ 1. ChEMBL mechanisms
parent = {}
for i, b in enumerate(batches):
    js = get_json("https://www.ebi.ac.uk/chembl/api/data/molecule.json?limit=1000&only=molecule_chembl_id,molecule_hierarchy"
                  "&molecule_chembl_id__in=" + ",".join(b), f"chembl_molecule_{i}.json")
    for m in js["molecules"]:
        h = m.get("molecule_hierarchy") or {}
        parent[m["molecule_chembl_id"]] = h.get("parent_chembl_id") or m["molecule_chembl_id"]
query_ids = sorted(set(drugs) | set(parent.values()))
mech = []
for i in range(0, len(query_ids), 40):
    b = query_ids[i:i + 40]
    offset = 0
    while True:
        js = get_json("https://www.ebi.ac.uk/chembl/api/data/mechanism.json?limit=1000&offset=%d&molecule_chembl_id__in=%s"
                      % (offset, ",".join(b)), f"chembl_mech_{i}_{offset}.json")
        mech += js["mechanisms"]
        if not js["page_meta"]["next"]:
            break
        offset += 1000
mech = pd.DataFrame(mech)
targets = sorted(mech.target_chembl_id.dropna().unique())
tsym = {}
for i in range(0, len(targets), 40):
    js = get_json("https://www.ebi.ac.uk/chembl/api/data/target.json?limit=1000&only=target_chembl_id,target_components"
                  "&target_chembl_id__in=" + ",".join(targets[i:i + 40]), f"chembl_target_{i}.json")
    for t in js["targets"]:
        syms = set()
        for c in t.get("target_components") or []:
            for s in c.get("target_component_synonyms") or []:
                if s.get("syn_type") == "GENE_SYMBOL":
                    syms.add(s["component_synonym"].upper())
        tsym[t["target_chembl_id"]] = syms
rows = []
child_of = {}
for d, p in parent.items():
    child_of.setdefault(p, set()).add(d)
for m in mech.itertuples():
    for sym in tsym.get(m.target_chembl_id, ()):
        for d in {m.molecule_chembl_id} | child_of.get(m.molecule_chembl_id, set()):
            rows.append((d, sym, m.action_type, m.mechanism_of_action))
chembl = pd.DataFrame(rows, columns=["Drug_ID", "symbol", "action", "mechanism"]).drop_duplicates()
chembl["dir"] = chembl.action.map(direction)
print(f"ChEMBL: {len(mech):,} mechanisms, {len(chembl):,} drug–gene rows")

# ------------------------------------------------------------------ 2. Guide to Pharmacology
gi = pd.read_csv(RAW / "gtopdb_interactions.csv", skiprows=1, low_memory=False)
gm = pd.read_csv(RAW / "gtopdb_ligand_id_mapping.csv", skiprows=1, low_memory=False)
lig2chembl = dict(zip(gm["Ligand id"].astype(str), gm["ChEMBl ID"]))
gi = gi[gi["Target Species"].fillna("").str.contains("Human", case=False)]
gi["Drug_ID"] = gi["Ligand ID"].astype(str).map(lig2chembl)
gi = gi.dropna(subset=["Drug_ID", "Target Ensembl Gene ID"])
gi = gi.assign(Gene_ID=gi["Target Ensembl Gene ID"].astype(str).str.split("|")).explode("Gene_ID")
gtop = pd.DataFrame({"Drug_ID": gi.Drug_ID, "Gene_ID": gi.Gene_ID.str.strip(),
                     "action": gi["Type"].fillna("") + " / " + gi["Action"].fillna(""),
                     "affinity": gi["Affinity Units"].fillna("") + " " + gi["Affinity Median"].astype(str)})
gtop["dir"] = [direction(t) or direction(a) for t, a in zip(gi["Type"], gi["Action"])]
print(f"GtoPdb: {len(gtop):,} human ligand–target rows with ChEMBL + Ensembl IDs")

# ------------------------------------------------------------------ 3. DGIdb
genes_u = sorted({gene_sym[g] for g in e.Gene_ID})
drows = []
for i in range(0, len(genes_u), 25):
    names = json.dumps(genes_u[i:i + 25])
    q = '{ genes(names: %s) { nodes { name interactions { drug { name conceptId } interactionTypes { type directionality } sources { sourceDbName } } } } }' % names
    js = get_json("https://dgidb.org/api/graphql", f"dgidb_{i}.json", data=json.dumps({"query": q}).encode())
    for gnode in js["data"]["genes"]["nodes"]:
        for it in gnode["interactions"]:
            types = it["interactionTypes"] or []
            dirs = {("down" if t["directionality"] == "INHIBITORY" else "up" if t["directionality"] == "ACTIVATING" else None)
                    for t in types} - {None}
            drows.append((gnode["name"].upper(), it["drug"]["name"].upper(), it["drug"]["conceptId"],
                          ";".join(sorted(t["type"] for t in types)), ";".join(sorted(dirs)),
                          ";".join(sorted(s["sourceDbName"] for s in it["sources"]))))
dg = pd.DataFrame(drows, columns=["symbol", "drug_name", "concept", "action", "dirs", "sources"])
dg["dir"] = dg.dirs.map(lambda s: None if not s else ("mixed" if ";" in s else s))
dg["Drug_ID"] = dg.concept.str.extract(r"chembl:(CHEMBL\d+)", expand=False)
print(f"DGIdb: {len(dg):,} gene–drug interactions for {len(genes_u):,} genes")

# ------------------------------------------------------------------ match to our unknown links
u = unk[["Drug_ID", "Drug", "Gene_ID", "Gene", "Action_raw", "progression_class", "fibrosis_direction",
         "inflammation_direction", "paper1_group", "cluster"]].copy()
u["symbol"] = u.Gene.str.upper()
u["drug_upper"] = u.Drug.str.upper()


def pick(df, keys, label):
    g = df.dropna(subset=["dir"]).groupby(keys).agg(dirs=("dir", lambda s: ";".join(sorted(set(s)))),
                                                  actions=("action", lambda s: "; ".join(sorted(set(map(str, s))))[:200]))
    g[f"{label}_dir"] = g.dirs.map(lambda s: "mixed" if ";" in s else s)
    return g[[f"{label}_dir"]].join(g.actions.rename(f"{label}_action"))


u = u.join(pick(chembl, ["Drug_ID", "symbol"], "chembl"), on=["Drug_ID", "symbol"])
u = u.join(pick(gtop, ["Drug_ID", "Gene_ID"], "gtopdb"), on=["Drug_ID", "Gene_ID"])
dg_id = pick(dg.dropna(subset=["Drug_ID"]), ["Drug_ID", "symbol"], "dgidb_id")
dg_nm = pick(dg, ["drug_name", "symbol"], "dgidb_nm").rename_axis(["drug_upper", "symbol"])
u = u.join(dg_id, on=["Drug_ID", "symbol"]).join(dg_nm, on=["drug_upper", "symbol"])
u["dgidb_dir"] = u.dgidb_id_dir.fillna(u.dgidb_nm_dir)
u["dgidb_action"] = u.dgidb_id_action.fillna(u.dgidb_nm_action)


def combine(r):
    ds = {d for d in (r.chembl_dir, r.gtopdb_dir, r.dgidb_dir) if isinstance(d, str)}
    if not ds:
        return "still unknown"
    if ds == {"down"}:
        return "down"
    if ds == {"up"}:
        return "up"
    return "conflicting"


u["external_direction"] = u.apply(combine, axis=1)
u["n_sources_agreeing"] = u.apply(lambda r: sum(isinstance(x, str) and x == r.external_direction
                                                for x in (r.chembl_dir, r.gtopdb_dir, r.dgidb_dir)), axis=1)
gdir = np.where(u.progression_class == "inflammation_only", u.inflammation_direction, u.fibrosis_direction)
u["gene_direction"] = gdir
u["graph_action_type"] = np.where(u.Action_raw == "TARGET", "TARGET (binding, no direction)", "ENZYME/TRANSPORTER/CARRIER")
keep = ["Drug_ID", "Drug", "Gene_ID", "Gene", "Action_raw", "graph_action_type", "progression_class", "gene_direction",
        "chembl_dir", "chembl_action", "gtopdb_dir", "gtopdb_action", "dgidb_dir", "dgidb_action",
        "external_direction", "n_sources_agreeing"]
u[keep].to_csv(OUT / "unknown_links_external_lookup.csv", index=False)

# ------------------------------------------------------------------ positive control: links whose direction the graph knows
k = e[e.Action_direction.isin(["up", "down"])][["Drug_ID", "Drug", "Gene_ID", "Gene", "Action_direction"]].copy()
k["symbol"] = k.Gene.str.upper()
k["drug_upper"] = k.Drug.str.upper()
k = k.join(pick(chembl, ["Drug_ID", "symbol"], "chembl"), on=["Drug_ID", "symbol"])
k = k.join(pick(gtop, ["Drug_ID", "Gene_ID"], "gtopdb"), on=["Drug_ID", "Gene_ID"])
k = k.join(dg_id, on=["Drug_ID", "symbol"]).join(dg_nm, on=["drug_upper", "symbol"])
k["dgidb_dir"] = k.dgidb_id_dir.fillna(k.dgidb_nm_dir)
ctrl = []
for lab, col in [("ChEMBL mechanisms", "chembl_dir"), ("Guide to Pharmacology", "gtopdb_dir"), ("DGIdb", "dgidb_dir")]:
    f_ = k[col].isin(["up", "down"])
    a_ = f_ & (k[col] == k.Action_direction)
    ctrl.append((lab, int(f_.sum()), f"{100 * f_.mean():.0f}%", int(a_.sum()), f"{100 * a_.sum() / max(f_.sum(), 1):.0f}%"))
ctrl = pd.DataFrame(ctrl, columns=["source", "known links found", "% found", "agree with graph", "% agree"]).set_index("source")
ctrl.to_csv(OUT / "positive_control.csv")


# ------------------------------------------------------------------ summary
def md(df, idx=True):
    x = df.reset_index() if idx else df
    lines = ["| " + " | ".join(map(str, x.columns)) + " |", "|" + "---|" * x.shape[1]]
    for row in x.values:
        lines.append("| " + " | ".join(f"{v:,}" if isinstance(v, (int, np.integer)) else str(v) for v in row) + " |")
    return "\n".join(lines)


n = len(u)
by_src = pd.DataFrame({
    "links given a direction": [int(u.chembl_dir.isin(["down", "up"]).sum()), int(u.gtopdb_dir.isin(["down", "up"]).sum()),
                                int(u.dgidb_dir.isin(["down", "up"]).sum()), int(u.external_direction.isin(["down", "up"]).sum())],
}, index=["ChEMBL mechanisms", "Guide to Pharmacology", "DGIdb", "Any source (agreeing)"])
by_src["% of 7,919 unknown links"] = (100 * by_src["links given a direction"] / n).round(1)
by_type = pd.crosstab(u.graph_action_type, u.external_direction).reindex(columns=["down", "up", "conflicting", "still unknown"], fill_value=0)
resolved = u[u.external_direction.isin(["down", "up"])]
known_before = int((e.Action_direction.isin(["up", "down"])).sum())
after = known_before + len(resolved)
summary = f"""# Drug action directions for OptimusKG links without one (graph preparation, step 0.5)

**Question:** the knowledge graph gives no push direction for {n:,} of {len(e):,} drug–gene links (TARGET {int((u.Action_raw == 'TARGET').sum()):,}; ENZYME / TRANSPORTER / CARRIER {int((u.Action_raw != 'TARGET').sum()):,}). Can other resources supply it?
**Answer: partly.** {len(resolved):,} of the {n:,} unknown links ({100 * len(resolved) / n:.1f}%) get a clear direction from at least one source with no disagreement. A further {int((u.external_direction == 'conflicting').sum()):,} get conflicting labels. That would raise known-direction links from **{known_before:,} ({100 * known_before / len(e):.1f}%) to {after:,} ({100 * after / len(e):.1f}%)**.

## Positive control (does the matching work?)
The same lookup was run on the {len(k):,} links whose direction the graph already knows:

{md(ctrl)}

The sources find many known links and agree with the graph almost always, so the matching works. **The low fill rate below is real:** most unknown links have no curated direction anywhere.

## Coverage of the unknown links, by source
{md(by_src)}

Sources (all free): **ChEMBL** drug mechanisms (EBI REST API; matched by drug ChEMBL ID, including the parent compound of salts, and target gene symbol); **IUPHAR/BPS Guide to Pharmacology** 2026.3 (public interactions.csv; matched by ligand ChEMBL ID and target Ensembl ID, also gives affinity values); **DGIdb 5** (GraphQL API; combines Guide to Pharmacology, TTD, Drug Repurposing Hub, ChEMBL and others; matched by gene symbol and drug ChEMBL ID or name).

## By graph label
{md(by_type)}

ChEMBL adds almost nothing because the knowledge graph already contains the ChEMBL / Open Targets mechanisms: those links are the ones that already have INHIBITOR, AGONIST and similar labels. The graph's TARGET links are mostly DrugBank "binds to" entries for which no curated database states a direction. TARGET links (binding with no stated direction) are the ones where external sources help most. ENZYME / TRANSPORTER / CARRIER links mostly describe drug metabolism; an external "inhibitor" label there (e.g. a CYP inhibitor) is real but rarely therapeutic.

## Other resources considered
- **DrugBank** (full action fields) needs an academic licence; the graph already uses its public part.
- **CTD** (Comparative Toxicogenomics Database) gives curated chemical → gene "increases / decreases" statements, but mostly for expression, not protein activity, and is a large download.
- **TTD** is included through DGIdb.
- **ChEMBL activities / BindingDB** give affinity (Ki, IC50), not direction; useful later for potency.

## Caveats
- Directions come from GtoPdb, DGIdb and ChEMBL; they are stored on the graph links as separate external_action_* properties.
- Gene matching in ChEMBL and DGIdb uses gene symbols; Guide to Pharmacology uses Ensembl IDs.
- Labels like "modulator" or "binder" still give no direction and stay unknown.

Per-link detail: `unknown_links_external_lookup.csv`. Raw downloads: `raw/`. Script: `code/external_action_check.py`.
"""
(OUT / "summary.md").write_text(summary, encoding="utf-8")
print(summary)
