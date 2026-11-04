"""Paper 2 Task 4 (input step): known indications of the drugs added to the graph on 21 Oct 2026.

The 519 Drug nodes added from ChEMBL, DGIdb, GtoPdb and Open Targets carry drug-gene
links only, with no TREATS / OFF_LABEL_FOR links. Task 4 needs their existing disease
uses to decide whether a drug is already used for a fibrosis-type liver disease. This
script fetches their indications from ChEMBL (the source the added drugs' IDs come
from) and saves them with the ChEMBL release and retrieval date. Task 4 reads the file.

Indications are matched on the drug's ChEMBL ID as a molecule or as a parent molecule,
so salt forms recorded under the parent are included. Added drugs whose graph ID is not
a ChEMBL ID (RxNorm, GtoPdb, NCIt) are first matched to ChEMBL by exact name (preferred
name, then synonym); the match method is recorded and unmatched drugs are listed. Every
phase is kept (`max_phase_for_ind`); Task 4 decides how phases are used.

Needs internet access; no Neo4j. Run from the repository root:

    python Paper2_Results/code/task4_added_drug_indications.py
"""
import json
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Paper2_Results" / "Task4"
OUT.mkdir(parents=True, exist_ok=True)

EDGES = ROOT / "Paper2_Results/Task2/drug_gene_edges.csv"
API = "https://www.ebi.ac.uk/chembl/api/data"
BATCH = 40
FIELDS = ["molecule_chembl_id", "parent_molecule_chembl_id", "efo_id", "efo_term", "mesh_id", "mesh_heading",
          "max_phase_for_ind", "drugind_id"]


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}),
                                        timeout=120) as r:
                return json.load(r)
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))


def indications(ids, field):
    rows, url = [], f"{API}/drug_indication.json?" + urllib.parse.urlencode(
        {f"{field}__in": ",".join(ids), "limit": 1000})
    while url:
        d = get(url)
        rows += [{k: r.get(k) for k in FIELDS} for r in d["drug_indications"]]
        nxt = d["page_meta"]["next"]
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return rows


def chembl_by_name(name):
    """Exact (case-insensitive) ChEMBL match on preferred name, then on synonym."""
    for field, how in [("pref_name__iexact", "preferred name"),
                       ("molecule_synonyms__molecule_synonym__iexact", "synonym")]:
        d = get(f"{API}/molecule.json?" + urllib.parse.urlencode({field: name, "limit": 5}))
        hits = {m["molecule_chembl_id"] for m in d["molecules"]}
        if len(hits) == 1:
            return hits.pop(), how
        if len(hits) > 1:
            return None, f"ambiguous ({how}: {';'.join(sorted(hits))})"
    return None, "not found"


edges = pd.read_csv(EDGES, low_memory=False)
added = edges.loc[edges.Link_origin != "OptimusKG", ["Drug_ID", "Drug", "drug_stage", "is_approved"]]
added = added[added.is_approved.isna() & added.drug_stage.notna()].drop_duplicates("Drug_ID").copy()
print(f"added drug nodes on progression genes: {len(added)}")

# ChEMBL ID for every added drug: its own ID, or an exact name match
is_chembl = added.Drug_ID.str.upper().str.startswith("CHEMBL")
added["chembl_id"] = added.Drug_ID.where(is_chembl).str.upper()
added["chembl_match"] = is_chembl.map({True: "graph ID", False: None})
todo = list(added.index[~is_chembl])
with ThreadPoolExecutor(max_workers=8) as pool:   # ChEMBL name search takes ~15 s per query
    for n, (i, (cid, how)) in enumerate(zip(todo, pool.map(lambda j: chembl_by_name(added.at[j, "Drug"]), todo)), 1):
        added.at[i, "chembl_id"], added.at[i, "chembl_match"] = cid, how
        if n % 20 == 0:
            print(f"  name matches {n}/{len(todo)}", flush=True)
added.to_csv(OUT / "added_drug_chembl_ids.csv", index=False)
print(added.chembl_match.value_counts().to_string())

lookup = added.dropna(subset=["chembl_id"])
ids = sorted(set(lookup.chembl_id))
status = get(f"{API}/status.json")
rows = []
for i in range(0, len(ids), BATCH):
    chunk = ids[i:i + BATCH]
    rows += indications(chunk, "molecule_chembl_id") + indications(chunk, "parent_molecule_chembl_id")
    print(f"  {min(i + BATCH, len(ids))}/{len(ids)}", flush=True)

ind = pd.DataFrame(rows, columns=FIELDS).drop_duplicates("drugind_id")
# assign each indication to the ChEMBL ID it was found for (molecule or parent), then to the graph drug
ind["chembl_id"] = ind.molecule_chembl_id.where(ind.molecule_chembl_id.isin(ids), ind.parent_molecule_chembl_id)
ind = ind[ind.chembl_id.isin(ids)].merge(lookup[["chembl_id", "Drug_ID", "Drug", "chembl_match"]], on="chembl_id")
ind["max_phase_for_ind"] = pd.to_numeric(ind.max_phase_for_ind, errors="coerce")
ind["efo_id_graph"] = ind.efo_id.str.replace(":", "_", regex=False)   # graph Disease IDs use EFO_0001234 form
ind["chembl_release"] = status.get("chembl_db_version")
ind["chembl_release_date"] = status.get("chembl_release_date")
ind["retrieved"] = date.today().isoformat()
ind = ind.sort_values(["Drug", "max_phase_for_ind"], ascending=[True, False])
ind.to_csv(OUT / "added_drug_indications_chembl.csv", index=False)

meta = {"chembl_release": status.get("chembl_db_version"), "chembl_release_date": status.get("chembl_release_date"),
        "retrieved": date.today().isoformat(), "added_drugs": len(added),
        "chembl_match": added.chembl_match.fillna("none").value_counts().to_dict(),
        "added_drugs_with_chembl_id": int(added.chembl_id.notna().sum()),
        "added_drugs_with_indication": int(ind.Drug_ID.nunique()), "indication_rows": len(ind)}
(OUT / "added_drug_indications_chembl.json").write_text(json.dumps(meta, indent=2))
print(json.dumps(meta, indent=2))
