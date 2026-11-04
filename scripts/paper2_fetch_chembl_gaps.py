"""
paper2_fetch_chembl_gaps.py
---------------------------
Fill the drug-metadata gaps OptimusKG leaves, from ChEMBL.

WHY ONLY THE GAPS
OptimusKG is itself built from ChEMBL and DrugBank, so most of what ChEMBL holds
is already recovered by paper2_recover_optimuskg_evidence.py. Re-fetching all
12,025 drugs would be slow and would mostly re-download what we have. This asks
ChEMBL only where OptimusKG is short.

WHAT THIS ADDS TO PAPER 2, MEASURED ON THE RIGHT DENOMINATOR
Kept as a separate step and reported in its own folder
(Paper2_Results/Recovery2_ChEMBL) on the mentor's instruction of 13 October:
two jobs were asked for, so the two stay separate.

The contribution depends entirely on which drugs you count:

  all 18,640 OptimusKG raw records   +2,638 SMILES  +938 names  +688 phases
  the graph's 12,025 Drug nodes          +0 SMILES   +16 names   +28 phases
  the 2,404 drugs in Task 4              +0 SMILES    +8 names   +17 phases

Only the last line matters for the paper: 25 field values, no structures, and
no result changes. The drugs missing a structure are precisely the ones that
never reach a progression gene -- ChEMBL holds records for 2,654 of them and
has a structure for none, because they are biologics, antibodies and mixtures
that no SMILES can represent.

A 60-drug pilot first suggested ChEMBL would supply structures for 57 of 59.
It sampled all 18,640 records rather than the subset in use, so it drew almost
entirely from drugs Paper 2 never touches. The pilot was not wrong about
ChEMBL; it was wrong about which drugs mattered. An earlier four-drug test
suggested the opposite -- it drew cisplatin and bismuth subsalicylate,
inorganics with no standard SMILES. Both small samples misled, in opposite
directions, for the same reason.

NETWORK
curl does the fetching. The requests library cannot resolve ebi.ac.uk from
inside the uv subprocess while curl reaches it fine, and shelling out is less
trouble than diagnosing that.

Usage:
    uv run --with pandas python scripts/paper2_fetch_chembl_gaps.py
"""

import json
import os
import subprocess
import sys
import time

import pandas as pd

RECOVERED = "data/optimuskg_recovered/drug_metadata.csv"
OUT = "data/chembl_recovered/chembl_gap_fill.csv"
BATCH = 20
API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"


def fetch(ids):
    url = f"{API}?molecule_chembl_id__in={','.join(ids)}&limit=100"
    for attempt in range(3):
        r = subprocess.run(["curl", "-s", "--max-time", "90", url],
                           capture_output=True, text=True)
        try:
            return json.loads(r.stdout).get("molecules", [])
        except json.JSONDecodeError:
            time.sleep(2 * (attempt + 1))
    print(f"   batch failed after 3 attempts: {ids[0]}...")
    return []


def main():
    m = pd.read_csv(RECOVERED)

    # drugs the graph actually holds, where OptimusKG is missing something
    need = m[(m.canonical_smiles.isna())
             | (m.type.isna())
             | (m.maximum_clinical_trial_phase.isna())
             | (m.name == m.Drug_ID)]
    ids = need.Drug_ID.dropna().unique().tolist()
    print(f"OptimusKG drug records: {len(m):,}")
    print(f"records with a gap:     {len(ids):,}")
    print(f"batches of {BATCH}:          {(len(ids) + BATCH - 1) // BATCH:,}\n")

    rows, done = [], 0
    for i in range(0, len(ids), BATCH):
        for x in fetch(ids[i:i + BATCH]):
            struct = x.get("molecule_structures") or {}
            rows.append({
                "Drug_ID": x.get("molecule_chembl_id"),
                "chembl_pref_name": x.get("pref_name"),
                "chembl_molecule_type": x.get("molecule_type"),
                "chembl_max_phase": x.get("max_phase"),
                "chembl_first_approval": x.get("first_approval"),
                "chembl_withdrawn": x.get("withdrawn_flag"),
                "chembl_oral": x.get("oral"),
                "chembl_canonical_smiles": struct.get("canonical_smiles"),
                "chembl_inchi_key": struct.get("standard_inchi_key"),
            })
        done += BATCH
        if done % 400 == 0:
            print(f"   {min(done, len(ids)):,} / {len(ids):,}")

    out = pd.DataFrame(rows).drop_duplicates(subset="Drug_ID")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT, index=False)

    print(f"\nreturned: {len(out):,} drug records")
    for c in ["chembl_pref_name", "chembl_molecule_type", "chembl_max_phase",
              "chembl_canonical_smiles"]:
        print(f"   {c:<26} {out[c].notna().sum():>6,}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
