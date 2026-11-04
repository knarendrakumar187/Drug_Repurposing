"""Graph preparation, step 0.4: find drug–gene links that the OptimusKG layer does not have.

Run on 21 Oct 2026, before the links were added to the graph. For all 6,222 progression genes this
collects every drug / compound that acts on the gene in four free public resources:

  1. IUPHAR/BPS Guide to Pharmacology 2026.3.1  (public interactions.csv, Ensembl IDs)
  2. DGIdb 5                                    (public interactions.tsv download, gene symbols)
  3. ChEMBL drug mechanisms                     (REST API, every mechanism, matched by target gene symbol)
  4. Open Targets Platform, latest release      (GraphQL API, drugAndClinicalCandidates per Ensembl ID)

Each drug is reduced to its ChEMBL *parent* compound (so salts and the free base count once) and
each pair is checked against every Drug-[:ACTS_ON]->Gene link the graph had for these genes at that
time (8,948 OptimusKG links, read from the drug-gene table in Task2/drug_gene_edges.csv). Pairs not in
the graph are the NEW links; write_new_links_to_graph.py then adds them to Neo4j. Run against the
final graph (which already holds them) this search would find no new links.

Outputs: Paper2_Results/New_drug_gene_links/.
Downloads are cached in its raw/ folder and only read as data.

    python Paper2_Results/code/new_links_search.py
"""
import json
from concurrent.futures import ThreadPoolExecutor
import re
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "Paper2_Results"
OUT = RES / "New_drug_gene_links"
RAW = OUT / "raw"
RAW.mkdir(parents=True, exist_ok=True)
OLD_RAW = RES / "External_action_check" / "raw"      # Guide to Pharmacology download from the 21 Oct check

DOWN_WORDS = ("INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE", "INVERSE AGONIST", "DEGRADER", "INHIBITION",
              "INHIBITORY", "ANTISENSE", "RNAI", "DISRUPTING", "SUPPRESSOR", "ALLOSTERIC ANTAGONIST")
UP_WORDS = ("AGONIST", "ACTIVATOR", "POSITIVE", "OPENER", "STABILISER", "STABILIZER", "INDUCER", "ACTIVATION",
            "ACTIVATING", "POTENTIATOR", "ENHANCER")


def direction(label):
    s = str(label).upper().replace("_", " ")
    if not s.strip() or s in ("NAN", "NONE"):
        return None
    if any(w in s for w in DOWN_WORDS):          # checked first: 'inverse agonist', 'negative modulator'
        return "down"
    if any(w in s for w in UP_WORDS):
        return "up"
    return None


def get_json(url, cache_name, data=None):
    f = RAW / cache_name
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))
    req = urllib.request.Request(url, data=data, headers={"Accept": "application/json", "Content-Type": "application/json"})
    for attempt in range(5):
        try:
            js = json.load(urllib.request.urlopen(req, timeout=180))
            if isinstance(js, dict) and js.get("errors") and not js.get("data"):
                raise RuntimeError(str(js["errors"])[:300])
            f.write_text(json.dumps(js), encoding="utf-8")
            time.sleep(0.15)
            return js
        except Exception as ex:                  # noqa: BLE001  retry network hiccups
            last = ex
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"failed: {url[:120]} {last}")


def norm_name(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


# ------------------------------------------------------------------ our genes and the graph's links
genes = pd.read_csv(RES / "Task1/Paper2_progression_genes.csv")
assert len(genes) == 6222
sym2ens = dict(zip(genes.Gene.str.upper(), genes.ensembl_id))
ens2sym = dict(zip(genes.ensembl_id, genes.Gene))
e = pd.read_csv(RES / "Task2/drug_gene_edges.csv")
print(f"genes {len(genes):,}; graph drug–gene links {len(e):,} ({e.Drug_ID.nunique():,} drugs, {e.Gene_ID.nunique():,} genes)")

# ------------------------------------------------------------------ 1. Guide to Pharmacology
gi = pd.read_csv(OLD_RAW / "gtopdb_interactions.csv", skiprows=1, low_memory=False)
gm = pd.read_csv(OLD_RAW / "gtopdb_ligand_id_mapping.csv", skiprows=1, low_memory=False)
lig2chembl = dict(zip(gm["Ligand id"].astype(str), gm["ChEMBl ID"]))
gi = gi[gi["Target Species"].fillna("").str.contains("Human", case=False)]
gi = gi[gi["Endogenous"].astype(str).str.lower() != "true"]          # body's own ligands are not drugs
gi = gi.assign(Gene_ID=gi["Target Ensembl Gene ID"].astype(str).str.split("|")).explode("Gene_ID")
gi["Gene_ID"] = gi.Gene_ID.str.strip()
gi = gi[gi.Gene_ID.isin(ens2sym)]
gtop = pd.DataFrame({
    "source": "GtoPdb", "chembl": gi["Ligand ID"].astype(str).map(lig2chembl).values,
    "drug_name": gi["Ligand"].values, "Gene_ID": gi.Gene_ID.values,
    "action": (gi["Type"].fillna("") + " / " + gi["Action"].fillna("")).values,
    "approved_flag": gi["Approved"].astype(str).str.lower().eq("true").values,
    "phase": np.nan, "ligand_type": gi["Ligand Type"].values,
    "alt_id": ("GtoPdb:" + gi["Ligand ID"].astype(str)).values})
gtop["dir"] = [direction(t) or direction(a) for t, a in zip(gi["Type"], gi["Action"])]
print(f"GtoPdb: {len(gtop):,} human, non-endogenous ligand–gene rows on our genes")

# ------------------------------------------------------------------ 2. DGIdb (whole database, kept for our 6,222 genes)
dt = pd.read_csv(RAW / "dgidb_interactions.tsv", sep="	", low_memory=False, na_values=["NULL"], keep_default_na=False)
dt["Gene_ID"] = dt.gene_name.astype(str).str.upper().map(sym2ens)
dt = dt.dropna(subset=["Gene_ID"])
# DGIdb's copy of Guide to Pharmacology is broken: GtoPdb *target* numbers were read as NCBI gene numbers
# (e.g. GtoPdb target 215 = dopamine D2 receptor, NCBI gene 215 = ABCD1), so every such row sits on the
# wrong gene. Guide to Pharmacology is used directly above, so these rows are dropped.
# DGIdb's copy of ChEMBL also keeps virus / bacteria targets matched by short symbols (e.g. hepatitis B
# polymerase "P" → human OCA2). ChEMBL is used directly below with human targets only, so these go too.
n_bad = int((dt.interaction_source_db_name == "GuideToPharmacology").sum())
n_bad2 = int((dt.interaction_source_db_name == "ChEMBL").sum())
dt = dt[~dt.interaction_source_db_name.isin(["GuideToPharmacology", "ChEMBL"])]
print(f"DGIdb: dropped {n_bad:,} GuideToPharmacology rows (gene-ID mix-up) and {n_bad2:,} ChEMBL rows "
      f"(non-human targets; ChEMBL used directly) on our genes")
cid = dt.drug_concept_id.fillna("").astype(str)
drows = list(zip(["DGIdb"] * len(dt),
                 [c.split(":", 1)[1].upper() if c.lower().startswith("chembl:") else None for c in cid],
                 dt.drug_name.where(~dt.drug_name.astype(str).str.upper().str.startswith("CHEMBL:"), dt.drug_claim_name),
                 dt.Gene_ID, dt.interaction_type.fillna(""), dt.approved.astype(str).str.upper().eq("TRUE"),
                 [np.nan] * len(dt), [None] * len(dt), cid, dt.interaction_type.map(direction),
                 dt.interaction_source_db_name.astype(str)))
dg = pd.DataFrame(drows, columns=["source", "chembl", "drug_name", "Gene_ID", "action", "approved_flag", "phase",
                                  "ligand_type", "alt_id", "dir", "dgidb_sources"])
print(f"DGIdb: {len(dg):,} drug–gene rows on {dg.Gene_ID.nunique():,} of our genes")

# ------------------------------------------------------------------ 3. ChEMBL mechanisms (all of them)
mech, offset = [], 0
while True:
    js = get_json(f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?limit=1000&offset={offset}", f"chembl_allmech_{offset}.json")
    mech += js["mechanisms"]
    if not js["page_meta"]["next"]:
        break
    offset += 1000
mech = pd.DataFrame(mech)
targets = sorted(mech.target_chembl_id.dropna().unique())
tsym = {}
for i in range(0, len(targets), 50):
    js = get_json("https://www.ebi.ac.uk/chembl/api/data/target.json?limit=1000&only=target_chembl_id,organism,target_components"
                  "&target_chembl_id__in=" + ",".join(targets[i:i + 50]), f"chembl_alltarget_{i}.json")
    for t in js["targets"]:
        if t.get("organism") != "Homo sapiens":
            continue
        s = {syn["component_synonym"].upper() for c in t.get("target_components") or []
             for syn in c.get("target_component_synonyms") or [] if syn.get("syn_type") == "GENE_SYMBOL"}
        tsym[t["target_chembl_id"]] = {sym2ens[x] for x in s if x in sym2ens}
crows = []
for m in mech.itertuples():
    for gid in tsym.get(m.target_chembl_id, ()):
        crows.append(("ChEMBL", m.molecule_chembl_id, None, gid, m.action_type, False, m.max_phase, None,
                      m.mechanism_of_action, direction(m.action_type)))
chem = pd.DataFrame(crows, columns=["source", "chembl", "drug_name", "Gene_ID", "action", "approved_flag", "phase",
                                    "ligand_type", "alt_id", "dir"])
print(f"ChEMBL: {len(mech):,} mechanisms; {len(chem):,} rows on our genes")

# ------------------------------------------------------------------ 4. Open Targets (latest)
meta = get_json("https://api.platform.opentargets.org/api/v4/graphql", "ot_meta.json",
                data=json.dumps({"query": "{ meta { apiVersion { x y z } dataVersion { year month } } }"}).encode())
ot_release = "Open Targets Platform {year}.{month}".format(**meta["data"]["meta"]["dataVersion"])
STAGE = {"APPROVAL": 4, "PHASE_4": 4, "PHASE_3": 3, "PHASE_2_3": 2.5, "PHASE_2": 2, "PHASE_1_2": 1.5,
         "PHASE_1": 1, "EARLY_PHASE_1": 0.5, "PRECLINICAL": 0}
ens = sorted(genes.ensembl_id)
orows = []


def ot_batch(i):
    parts = " ".join(
        f't{j}: target(ensemblId:"{g}") {{ id drugAndClinicalCandidates {{ rows {{ maxClinicalStage drug {{ id name drugType '
        f'maximumClinicalStage mechanismsOfAction {{ rows {{ actionType mechanismOfAction targets {{ id }} }} }} }} }} }} }}'
        for j, g in enumerate(ens[i:i + 15]))
    return get_json("https://api.platform.opentargets.org/api/v4/graphql", f"ot_{i}.json",
                    data=json.dumps({"query": "{ " + parts + " }"}).encode())


with ThreadPoolExecutor(8) as pool:                  # 8 requests at a time; results cached in raw/
    ot_pages = list(pool.map(ot_batch, range(0, len(ens), 15)))
for js in ot_pages:
    for t in (js["data"] or {}).values():
        if not t:
            continue
        for r in t["drugAndClinicalCandidates"]["rows"]:
            d = r["drug"] or {}
            moa = [x for x in ((d.get("mechanismsOfAction") or {}).get("rows") or [])
                   if t["id"] in {y["id"] for y in x["targets"] or []}]
            acts = sorted({x["actionType"] for x in moa if x["actionType"]})
            dirs = {direction(a) for a in acts} - {None}
            stage = d.get("maximumClinicalStage") or r.get("maxClinicalStage")
            orows.append(("OpenTargets", d.get("id"), d.get("name"), t["id"], ";".join(acts),
                          stage in ("APPROVAL", "PHASE_4"), STAGE.get(stage, np.nan), d.get("drugType"),
                          "; ".join(sorted({x["mechanismOfAction"] for x in moa}))[:200],
                          None if not dirs else ("mixed" if len(dirs) > 1 else dirs.pop())))
ot = pd.DataFrame(orows, columns=["source", "chembl", "drug_name", "Gene_ID", "action", "approved_flag", "phase",
                                  "ligand_type", "alt_id", "dir"])
print(f"{ot_release}: {len(ot):,} drug–gene rows on {ot.Gene_ID.nunique():,} of our genes")

allr = pd.concat([gtop, dg, chem, ot], ignore_index=True)
# What kind of evidence is each row? Only curated drug-target statements count as a real new link.
CURATED_DGIDB = {"TdgClinicalTrial", "TEND", "TALC", "FDA"}   # gene given by symbol or UniProt ID
NAME_MATCHED_DGIDB = {"TTD"}                # curated, but the gene is matched from a protein name (some errors, e.g. KIR→GEM)
SCREEN_DGIDB = {"DTC"}                      # Drug Target Commons: measured binding in lab assays
EVID_RANK = {"curated drug-target": 0, "curated, TTD only (gene matched by protein name)": 1,
             "lab binding screen": 2, "literature / biomarker only": 3}
allr["evidence"] = np.select(
    [allr.source.isin(["GtoPdb", "ChEMBL", "OpenTargets"]) | allr.dgidb_sources.isin(CURATED_DGIDB),
     allr.dgidb_sources.isin(NAME_MATCHED_DGIDB), allr.dgidb_sources.isin(SCREEN_DGIDB)],
    ["curated drug-target", "curated, TTD only (gene matched by protein name)", "lab binding screen"],
    "literature / biomarker only")
allr["source_detail"] = np.where(allr.source == "DGIdb", "DGIdb/" + allr.dgidb_sources.astype(str), allr.source)

# ------------------------------------------------------------------ drug identity: ChEMBL parent + names
ids = sorted(set(allr.chembl.dropna()) | set(e.Drug_ID))
mol = {}


def mol_batch(i):
    return get_json("https://www.ebi.ac.uk/chembl/api/data/molecule.json?limit=1000&only=molecule_chembl_id,pref_name,"
                    "max_phase,molecule_type,molecule_hierarchy&molecule_chembl_id__in=" + ",".join(ids[i:i + 50]),
                    f"chembl_allmol_{i}.json")


with ThreadPoolExecutor(6) as pool:
    mol_pages = list(pool.map(mol_batch, range(0, len(ids), 50)))
for js in mol_pages:
    for m in js["molecules"]:
        h = m.get("molecule_hierarchy") or {}
        mol[m["molecule_chembl_id"]] = (h.get("parent_chembl_id") or m["molecule_chembl_id"], m.get("pref_name"),
                                        m.get("max_phase"), m.get("molecule_type"))
parent = {k: v[0] for k, v in mol.items()}
pinfo = {}                                          # parent → (name, phase, type)
for k, (p, n, ph, ty) in mol.items():
    if k == p or p not in pinfo:
        pinfo[p] = (n, float(ph) if ph not in (None, "") else np.nan, ty)
# name → ChEMBL parent, for DGIdb / GtoPdb rows without a ChEMBL ID
name2id = {}
for d, n in zip(e.Drug_ID, e.Drug):
    name2id.setdefault(norm_name(n), parent.get(d, d))
for k, (p, n, ph, ty) in mol.items():
    if n:
        name2id.setdefault(norm_name(n), p)
for d, n in zip(allr.chembl, allr.drug_name):
    if isinstance(d, str) and isinstance(n, str):
        name2id.setdefault(norm_name(n), parent.get(d, d))
for n, inn, c in zip(gm.Name, gm.INN, gm["ChEMBl ID"]):
    if isinstance(c, str):
        for x in (n, inn):
            if isinstance(x, str):
                name2id.setdefault(norm_name(x), parent.get(c, c))
allr["Drug_key"] = [parent.get(c, c) if isinstance(c, str) else name2id.get(norm_name(n)) for c, n in zip(allr.chembl, allr.drug_name)]
allr["matched_by"] = np.where(allr.chembl.notna(), "ChEMBL ID", np.where(allr.Drug_key.notna(), "drug name", "no ChEMBL match"))
allr["Drug_key"] = allr.Drug_key.fillna(allr.alt_id.where(allr.source.isin(["GtoPdb", "DGIdb"])))
allr["Drug_key"] = allr.Drug_key.fillna("NAME:" + allr.drug_name.map(norm_name))

# ------------------------------------------------------------------ the graph's own links (parent level)
e["Drug_key"] = e.Drug_ID.map(lambda d: parent.get(d, d))
graph_pairs = set(zip(e.Drug_key, e.Gene_ID)) | set(zip(e.Drug_ID, e.Gene_ID))
graph_drugs = set(e.Drug_key) | set(e.Drug_ID)
dd = pd.read_excel(RES / "Graph_files/3_drug_disease.xlsx", sheet_name="drug_disease_links")
dcol = "drug_id"
graph_drugs_any = graph_drugs | {parent.get(d, d) for d in dd[dcol].astype(str)} | set(dd[dcol].astype(str))
allr["in_graph"] = [(k, g) in graph_pairs for k, g in zip(allr.Drug_key, allr.Gene_ID)]

# ------------------------------------------------------------------ positive control: does each source find the graph's links?
pharm = e[~e.Action_raw.isin(["ENZYME", "TRANSPORTER", "CARRIER"])]
ctrl = []
for src in ["GtoPdb", "DGIdb", "ChEMBL", "OpenTargets"]:
    s = allr[allr.source == src]
    sp = set(zip(s.Drug_key, s.Gene_ID))
    for lab, sub in [("all graph links", e), ("pharmacological graph links", pharm)]:
        hit = sum((k, g) in sp for k, g in zip(sub.Drug_key, sub.Gene_ID))
        ctrl.append((src, lab, len(sub), hit, round(100 * hit / len(sub), 1)))
ctrl = pd.DataFrame(ctrl, columns=["source", "graph links tested", "n", "found by source", "% found"])
ctrl.to_csv(OUT / "positive_control.csv", index=False)
print(ctrl.to_string(index=False))

# ------------------------------------------------------------------ one row per NEW pair
new = allr[~allr.in_graph].copy()


def agg_dir(s):
    d = set(x for x in s if isinstance(x, str))
    if not d:
        return "unknown"
    return d.pop() if len(d) == 1 and "mixed" not in d else "conflicting"


g = new.groupby(["Drug_key", "Gene_ID"])
pairs = g.agg(
    sources=("source", lambda s: ";".join(sorted(set(s)))),
    source_detail=("source_detail", lambda s: ";".join(sorted(set(s)))),
    evidence=("evidence", lambda s: min(s, key=EVID_RANK.get)),
    n_sources=("source", "nunique"),
    drug_name_src=("drug_name", lambda s: next((x for x in s if isinstance(x, str)), "")),
    approved_any=("approved_flag", "max"),
    phase_src=("phase", "max"),
    action_labels=("action", lambda s: "; ".join(sorted({str(x) for x in s if isinstance(x, str) and x.strip()}))[:300]),
    mechanism_text=("alt_id", lambda s: "; ".join(sorted({str(x) for x in s if isinstance(x, str) and not re.match(r"^(GtoPdb|chembl|rxcui|ncit|iuphar|drugbank|wikidata):", str(x), re.I)}))[:200]),
    push=("dir", agg_dir),
    matched_by=("matched_by", lambda s: ";".join(sorted(set(s)))),
    ligand_type=("ligand_type", lambda s: next((x for x in s if isinstance(x, str)), "")),
).reset_index()
pairs["Drug"] = [pinfo.get(k, (None,))[0] or n for k, n in zip(pairs.Drug_key, pairs.drug_name_src)]
pairs["chembl_max_phase"] = [pinfo.get(k, (None, np.nan))[1] for k in pairs.Drug_key]
pairs["molecule_type"] = [pinfo.get(k, (None, None, None))[2] or lt for k, lt in zip(pairs.Drug_key, pairs.ligand_type)]
pairs["max_phase"] = pairs[["chembl_max_phase", "phase_src"]].max(axis=1)
pairs["approved"] = pairs.approved_any.astype(bool) | (pairs.max_phase >= 4)
pairs["drug_stage"] = np.select([pairs.approved, pairs.max_phase >= 1], ["approved", "clinical trial"], "research compound")
pairs["drug_in_graph"] = np.select([pairs.Drug_key.isin(graph_drugs), pairs.Drug_key.isin(graph_drugs_any)],
                                   ["yes, acts on other progression genes", "yes, other links only"], "no (new drug)")
gi_ = genes.set_index("ensembl_id")
pairs["Gene"] = pairs.Gene_ID.map(ens2sym)
for c in ["progression_class", "paper1_group", "fibrosis_direction", "inflammation_direction", "trajectory_cluster", "Tier", "in_shortlist"]:
    pairs[c] = pairs.Gene_ID.map(gi_[c])
pairs["gene_had_drug_in_graph"] = pairs.Gene_ID.isin(set(e.Gene_ID))
pairs["gene_direction"] = np.where(pairs.progression_class == "inflammation_only", pairs.inflammation_direction, pairs.fibrosis_direction)
pairs["medicine"] = pairs.drug_stage.isin(["approved", "clinical trial"])
cols = ["Drug_key", "Drug", "drug_stage", "max_phase", "molecule_type", "drug_in_graph", "Gene_ID", "Gene", "progression_class", "paper1_group", "gene_direction", "trajectory_cluster",
        "Tier", "in_shortlist", "gene_had_drug_in_graph", "evidence", "sources", "source_detail", "n_sources", "action_labels", "mechanism_text", "push",
        "matched_by"]
pairs["evidence_rank"] = pairs.evidence.map(EVID_RANK)
pairs = pairs.sort_values(["medicine", "evidence_rank", "n_sources", "drug_stage", "Drug"], ascending=[False, True, False, True, True])
pairs[cols].rename(columns={"Drug_key": "Drug_ID"}).to_csv(OUT / "new_links_all.csv", index=False)
med = pairs[pairs.medicine]
med[cols].rename(columns={"Drug_key": "Drug_ID"}).to_csv(OUT / "new_links_medicines.csv", index=False)
allr.to_csv(OUT / "raw_rows_all_sources.csv.gz", index=False, compression="gzip")

# ------------------------------------------------------------------ numbers
def cnt(df):
    return dict(pairs=len(df), drugs=df.Drug_key.nunique(), genes=df.Gene_ID.nunique(),
                new_genes=df.loc[~df.gene_had_drug_in_graph, "Gene_ID"].nunique())


src_tab = []
for src in ["GtoPdb", "DGIdb", "ChEMBL", "OpenTargets"]:
    s = allr[allr.source == src]
    sn = s[~s.in_graph]
    sp = pairs[pairs.sources.str.contains(src)]
    src_tab.append((src, s.groupby(["Drug_key", "Gene_ID"]).ngroups, s[s.in_graph].groupby(["Drug_key", "Gene_ID"]).ngroups,
                    sn.groupby(["Drug_key", "Gene_ID"]).ngroups, int(sp.medicine.sum()), int((sp.n_sources == 1).sum())))
src_tab = pd.DataFrame(src_tab, columns=["source", "drug–gene pairs on our genes", "already in graph", "NEW pairs",
                                         "NEW pairs with a medicine", "found ONLY by this source"])
src_tab.to_csv(OUT / "by_source.csv", index=False)
stage_tab = pd.DataFrame({k: cnt(pairs[pairs.drug_stage == k]) for k in ["approved", "clinical trial", "research compound"]}).T
stage_tab.loc["all"] = list(cnt(pairs).values())
stage_tab.to_csv(OUT / "by_drug_stage.csv")
ev_tab = pd.crosstab([pairs.evidence], pairs.drug_stage, margins=True, margins_name="all")
ev_tab.to_csv(OUT / "by_evidence_and_stage.csv")
med_all = med
med = med[med.evidence == "curated drug-target"]          # headline set from here on
med[cols].rename(columns={"Drug_key": "Drug_ID"}).to_csv(OUT / "new_links_medicines_curated.csv", index=False)
gene_new = med[~med.gene_had_drug_in_graph].groupby(["Gene_ID", "Gene", "progression_class", "Tier", "in_shortlist"], dropna=False).agg(
    n_medicines=("Drug_key", "nunique"), approved=("drug_stage", lambda s: int((s == "approved").sum())),
    examples=("Drug", lambda s: "; ".join(sorted(set(map(str, s)))[:6]))).reset_index().sort_values(
    ["in_shortlist", "approved", "n_medicines"], ascending=False)
gene_new.to_csv(OUT / "genes_reached_for_first_time.csv", index=False)
drug_tab = med.groupby(["Drug_key", "Drug", "drug_stage", "drug_in_graph"], dropna=False).agg(
    new_genes=("Gene_ID", "nunique"),
    first_drug_for=("gene_had_drug_in_graph", lambda s: int((~s).sum())),
    genes=("Gene", lambda s: ";".join(sorted(set(map(str, s)))))).reset_index().sort_values(["new_genes"], ascending=False)
drug_tab.rename(columns={"Drug_key": "Drug_ID"}).to_csv(OUT / "medicines_by_drug.csv", index=False)

new_med_drugs = set(med.Drug_key)
newly_candidate = med[med.drug_in_graph.ne("yes, acts on other progression genes")].Drug_key.nunique()
numbers = {
    "new_pairs_medicines_any_evidence": cnt(med_all),
    "new_pairs_by_evidence": {k: cnt(pairs[pairs.evidence == k]) for k in EVID_RANK},
    "headline_set": "approved or clinical-trial drugs with a curated drug-target source",
    "release": {"GtoPdb": "2026.3.1", "DGIdb": "5 (interactions.tsv, downloaded Oct 2026)", "ChEMBL": "REST API (Oct 2026)", "OpenTargets": ot_release},
    "new_pairs_all": cnt(pairs), "new_pairs_medicines_curated": cnt(med),
    "new_pairs_approved_curated": cnt(med[med.drug_stage == "approved"]),
    "new_pairs_approved": cnt(pairs[pairs.drug_stage == "approved"]),
    "medicine_pairs_with_direction": int(med.push.isin(["up", "down"]).sum()),
    "medicine_pairs_2plus_sources": int((med.n_sources >= 2).sum()),
    "medicines_new_to_progression_genes": int(newly_candidate),
    "shortlist_genes_reached_first_time": int(gene_new.in_shortlist.sum()),
    "genes_with_drug_graph": int(e.Gene_ID.nunique()),
    "genes_with_drug_graph_plus_new": int(len(set(e.Gene_ID) | set(med.Gene_ID))),
    "neo4j_check": "21 Oct 2026: 120 random new pairs -> 0 have any Drug–Gene relationship in Neo4j; "
                   "150 random new-link drugs -> 121 exist as Drug nodes, 7 have any ACTS_ON",
}
(OUT / "numbers.json").write_text(json.dumps(numbers, indent=2, default=str, ensure_ascii=False), encoding="utf-8")
print(json.dumps(numbers, indent=2, default=str))
print(src_tab.to_string(index=False))
print(ev_tab)
print(stage_tab)
print(gene_new.head(15).to_string(index=False))
print(drug_tab.head(20).to_string(index=False))
