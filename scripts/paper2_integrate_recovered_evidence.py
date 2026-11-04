"""
paper2_integrate_recovered_evidence.py
--------------------------------------
Fold the recovered OptimusKG evidence and chemistry into Tasks 2, 4 and 7.

WHAT THIS FIXES
Three things the plan asked for and the tasks could not supply, because
build_optimuskg_layer.py read them and kept only `name` and `relation`:

  Task 2  Evidence/Source and Database ID      said "not available per edge"
  Task 7  Evidence source and evidence level   said "not recorded"
  Task 4  filter on molecule type and phase    impossible; name heuristic only

All three are now answerable from OptimusKG's own source files, with ChEMBL
filling the structure gap OptimusKG leaves.

USE IN THIS ANALYSIS
Only load_drug_metadata() is used: paper2_write_recovered_to_graph.py calls it to write the
recovered drug and drug-gene evidence onto the graph (graph-preparation steps 0.1-0.3). The
command-line part of this script is not used.

"""

import os

import pandas as pd

import sys, os as _os
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paper2_drug_filters import is_inorganic

# Two recoveries, two sources, two directories -- deliberately not pooled.
# Every field either source supplies carries a *_source column saying which one
# it came from, so a reader can always separate them again.
OKG = "data/optimuskg_recovered"
CHEMBL_DIR = "data/chembl_recovered"
EVIDENCE = f"{OKG}/drug_gene_evidence.csv"
METADATA = f"{OKG}/drug_metadata.csv"
CHEMBL = f"{CHEMBL_DIR}/chembl_gap_fill.csv"

T2 = "Paper2_Results/Task2/task2_drug_gene_edges.csv"
T4 = "Paper2_Results/Task4/task4_drug_coverage.csv"
T7 = "Paper2_Results/Task7/task7_action_direction_evidence.csv"

# Team 1's instruction was to filter on molecule type and clinical phase. Those
# fields are now available and they DO NOT WORK for this: all 70 name-flagged
# entries are "Small molecule", and ZINC ACETATE and ZINC CHLORIDE are
# max_phase 4.0 and is_approved True -- they genuinely are approved small
# molecules. The classification cannot separate a mineral salt from a drug
# because, as classified, there is nothing to separate.
#
# What does discriminate, for the inorganics at least, is chemistry: a molecule
# with no carbon in its formula is not an organic therapeutic. That rule is
# certain where it fires and catches ZnCl2, MgSO4, Na2SO4, AlCl3, CaCl2 and NO.
# It cannot catch the organic nutrients -- ATP, FAD, glutamic acid all contain
# carbon -- so the name heuristic stays as the broader, cruder net beside it.
#
# Both flags are reported. Neither replaces the other.
NON_THERAPEUTIC_TYPES = {"Unknown", "Inorganic", "Element"}


def coalesce(df, okg_col, chembl_col, out, label_okg="OptimusKG",
             label_chembl="ChEMBL"):
    """Prefer OptimusKG, fall back to ChEMBL, and RECORD WHICH.

    The mentor's instruction of 13 October was to keep the two recoveries
    separate. A plain fillna() would merge them into one column where nobody
    could tell them apart again, so every coalesced field gets a companion
    <out>_source column. Where ChEMBL is absent the column still exists and
    reads OptimusKG or nothing, so the schema does not change with the inputs.
    """
    okg = df[okg_col]
    df[out] = okg
    df[f"{out}_source"] = okg.notna().map({True: label_okg, False: None})
    if chembl_col is not None and chembl_col in df.columns:
        fill = okg.isna() & df[chembl_col].notna()
        df.loc[fill, out] = df.loc[fill, chembl_col]
        df.loc[fill, f"{out}_source"] = label_chembl
    return int((df[f"{out}_source"] == label_chembl).sum())


def load_drug_metadata():
    """OptimusKG drug properties, with ChEMBL filling only what OptimusKG lacks.

    Returns the table plus a per-field count of what ChEMBL contributed, so the
    caller can report the two recoveries separately rather than as one blended
    "recovered" figure.
    """
    m = pd.read_csv(METADATA, low_memory=False)
    m = m[["Drug_ID", "name", "type", "maximum_clinical_trial_phase",
           "is_approved", "has_been_withdrawn", "canonical_smiles",
           "inchi_key", "cd_formula", "cd_mol_weight"]]

    # OptimusKG stores an unnamed drug's name AS its ChEMBL accession, so an
    # accession in the name field means "no name", not a name.
    m["okg_name"] = m.name.where(m.name != m.Drug_ID)

    have_chembl = os.path.exists(CHEMBL)
    if have_chembl:
        c = pd.read_csv(CHEMBL, low_memory=False)
        m = m.merge(c, on="Drug_ID", how="left")
    else:
        print("   (no ChEMBL gap-fill file; using OptimusKG alone)")

    added = {
        "drug_name":     coalesce(m, "okg_name", "chembl_pref_name", "drug_name"),
        "molecule_type": coalesce(m, "type", "chembl_molecule_type", "molecule_type"),
        "max_phase":     coalesce(m, "maximum_clinical_trial_phase",
                                  "chembl_max_phase", "max_phase"),
        "smiles":        coalesce(m, "canonical_smiles",
                                  "chembl_canonical_smiles", "smiles"),
    }

    m["has_smiles"] = m.smiles.notna()
    m["named"] = m.drug_name.notna()
    return m, added, have_chembl


def main():
    ev = pd.read_csv(EVIDENCE, low_memory=False)
    md, chembl_added, have_chembl = load_drug_metadata()
    print(f"recovered evidence edges : {len(ev):,}")
    print(f"recovered drug records   : {len(md):,}")

    # Two recoveries, reported apart. Field by field: what OptimusKG gave, then
    # what ChEMBL added on top of it. Never one pooled "recovered" count.
    print("\n   field           OptimusKG    +ChEMBL      total")
    for f in ["drug_name", "molecule_type", "max_phase", "smiles"]:
        src = md[f"{f}_source"]
        okg = int((src == "OptimusKG").sum())
        ch = int((src == "ChEMBL").sum())
        print(f"   {f:<14} {okg:>9,} {ch:>10,} {okg + ch:>10,}")
    if have_chembl and sum(chembl_added.values()) == 0:
        print("\n   ChEMBL added nothing to these records.")

    ev_cols = ["Drug_ID", "Gene_ID", "sources_direct", "sources_indirect",
               "relation_assertions", "relation_conflict", "evidence_level",
               "source_urls", "mechanisms_of_action", "n_sources"]
    ev_small = ev[ev_cols]

    # ---- Task 2 ------------------------------------------------------------
    t2 = pd.read_csv(T2, low_memory=False)
    before = t2.Evidence_Source.iloc[0] if "Evidence_Source" in t2 else None
    t2 = t2.drop(columns=[c for c in ev_cols if c in t2.columns
                          and c not in ("Drug_ID", "Gene_ID")], errors="ignore")
    t2 = t2.merge(ev_small, on=["Drug_ID", "Gene_ID"], how="left")
    matched = t2.evidence_level.notna().sum()
    t2["Evidence_Source"] = (
        "OptimusKG: " + t2.sources_direct.fillna("unattributed")
        + t2.sources_indirect.fillna("").radd(" via ").where(
            t2.sources_indirect.notna(), ""))
    t2["Evidence_Level"] = t2.evidence_level.fillna("not matched to source")
    # evidence_level and Evidence_Level were identical on all 8,948 rows. Two
    # columns differing only in case are a hazard in Excel and in any
    # case-insensitive import, so only the Title_Case one -- matching Task 2's
    # own Drug / Gene / Evidence_Source convention -- is delivered.
    t2 = t2.drop(columns=["evidence_level"])
    t2.to_csv(T2, index=False)
    print(f"\nTask 2: {matched:,} of {len(t2):,} edges matched to source evidence")
    print(f"   Evidence_Source was: {str(before)[:60]!r}")
    print(f"   now e.g.           : {t2.Evidence_Source.iloc[0][:60]!r}")

    # ---- Task 7 ------------------------------------------------------------
    t7 = pd.read_csv(T7, low_memory=False)
    t7 = t7.drop(columns=[c for c in ev_cols if c in t7.columns
                          and c not in ("Drug_ID", "Gene_ID")], errors="ignore")
    t7 = t7.merge(ev_small, on=["Drug_ID", "Gene_ID"], how="left")
    t7["Evidence_Source"] = ("OptimusKG: "
                             + t7.sources_direct.fillna("unattributed"))
    t7["Evidence_Level"] = t7.evidence_level.fillna("not matched to source")
    n_t7_evidence = int(t7.evidence_level.notna().sum())
    t7 = t7.drop(columns=["evidence_level"])        # see the Task 2 note above
    t7.to_csv(T7, index=False)
    print(f"Task 7: {n_t7_evidence:,} of {len(t7):,} "
          "relationships now carry a real evidence level")
    print(t7.Evidence_Level.value_counts().to_string())

    # ---- Task 4 ------------------------------------------------------------
    t4 = pd.read_csv(T4, low_memory=False)
    # Every recovered field travels with its *_source companion, so OptimusKG
    # and ChEMBL stay separable in the delivered table and not just upstream.
    keep = ["Drug_ID",
            "drug_name", "drug_name_source",
            "molecule_type", "molecule_type_source",
            "max_phase", "max_phase_source",
            "smiles", "smiles_source",
            "is_approved", "has_been_withdrawn", "has_smiles",
            "inchi_key", "cd_formula", "cd_mol_weight"]
    # Drop the columns this script owns, INCLUDING ones it no longer produces.
    # The earlier version dropped only what was in `keep`, so a column written
    # by a superseded version of this script survived every rerun forever --
    # which is how flag_no_structure and likely_non_therapeutic_v2 came to sit
    # in the delivered table, written by nothing and flagging 1 drug where
    # likely_non_therapeutic flags 726. A "_v2" name on a dead column is worse
    # than no column: it reads as the newer, better one.
    RETIRED = ["flag_no_structure", "likely_non_therapeutic_v2"]
    t4 = t4.drop(columns=[c for c in keep + RETIRED if c in t4.columns
                          and c != "Drug_ID"], errors="ignore")
    t4 = t4.merge(md[keep], on="Drug_ID", how="left")

    # Team 1 asked to filter on molecule type and clinical phase. Both fields
    # are now available and NEITHER WORKS: all 70 name-flagged entries are
    # "Small molecule", and ZINC ACETATE and ZINC CHLORIDE carry max_phase 4.0
    # with is_approved True. They genuinely are approved small molecules, so
    # the classification has nothing to separate.
    #
    # Chemistry does discriminate, for the inorganics: a formula with no carbon
    # is not an organic therapeutic. Certain where it fires, and it catches
    # ZnCl2, MgSO4, Na2SO4, AlCl3, CaCl2 and NO. It cannot reach the organic
    # nutrients -- ATP, FAD and glutamic acid all contain carbon -- so the name
    # heuristic stays beside it as the broader, cruder net. Neither replaces
    # the other, and both are reported.
    t4["flag_molecule_type"] = t4.molecule_type.isin(NON_THERAPEUTIC_TYPES)
    # The same carbon test the shared filter uses, so Task 4 and Task 8 cannot
    # disagree about an entity again. Magnesium sulfate was flagged here by
    # formula and cleared by the name whitelist, so it reached Task 8 as
    # candidate 24; inorganics are now out of the whitelist and the test lives
    # in paper2_drug_filters.is_inorganic.
    t4["flag_inorganic"] = t4.cd_formula.map(is_inorganic)
    t4["likely_non_therapeutic"] = t4.likely_non_therapeutic | t4.flag_inorganic
    t4.to_csv(T4, index=False)

    print(f"\nTask 4: molecule type known for "
          f"{int(t4.molecule_type.notna().sum()):,} of {len(t4):,} drugs, "
          f"formula for {int(t4.cd_formula.notna().sum()):,}")
    print("   molecule_type does NOT separate salts from drugs -- the "
          "name-flagged entries are all:")
    print("      " + t4[t4.flag_name].molecule_type.value_counts().to_string())
    print(f"   flag_name (heuristic)     {int(t4.flag_name.sum()):>5,}")
    print(f"   flag_inorganic (formula)  {int(t4.flag_inorganic.sum()):>5,}"
          f"   newly caught: {int((t4.flag_inorganic & ~t4.flag_name).sum()):,}")
    print(f"   combined                  {int(t4.likely_non_therapeutic.sum()):>5,}")
    print(f"   SMILES now available      {int(t4.has_smiles.fillna(False).sum()):>5,}")


if __name__ == "__main__":
    main()
