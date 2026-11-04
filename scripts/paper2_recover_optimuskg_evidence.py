"""
paper2_recover_optimuskg_evidence.py
------------------------------------
Recover the drug evidence and chemistry that build_optimuskg_layer.py discarded.

WHY THIS EXISTS
The Paper 2 plan asks for an "Evidence/Source" and an "evidence level" on every
drug-gene relationship. Tasks 2 and 7 could not supply them: the graph's
ACTS_ON edges carry one property, `action`, and Drug nodes carry `node_id` and
`name`.

That is a property of the LOADER, not of the data. build_optimuskg_layer.py
reads `properties` from every parquet and then throws it away:

    cols = ["from", "to", "relation", "properties"]
    ...
    if fname == "disease_gene":   ...keep evidence_score...
    else:                          e["score"] = ""        # properties dropped

and for nodes it keeps only the display name:

    pd.DataFrame({"node_id": n["id"], "name": n["properties"].apply(name_of)})

Everything else OptimusKG carries was read into memory and discarded. This
recovers it from the same two source files, at the same DOI already recorded in
.gitignore -- the project's own knowledge layer, not a new resource.

WHAT COMES BACK

  drug_gene.parquet (36,611 edges), per edge:
    sources               which databases asserted it, direct and indirect
    relation_assertions   the per-source relation claim
    relation_conflict     whether the sources disagree
    source_urls           FDA label PDFs, on 38% of the 36,611 records
    mechanisms_of_action  free text, on 38% of the 36,611 records

    Those two percentages are over EVERY recovered record. On the 8,948
    edges Paper 2 actually uses the rate is 11.5%, because most records
    here are for edges outside the progression-gene set. Quoting the 38%
    against Paper 2's table is a denominator error, and was one until a
    reviewer caught it on 13 October.

  drug.parquet (18,640 nodes), per drug:
    type, maximum_clinical_trial_phase, is_approved, has_been_withdrawn
    canonical_smiles, inchi_key, cd_formula, cd_mol_weight
    synonyms, trade_names, description, year_of_first_approval

TWO FINDINGS

1. SMILES are available: canonical_smiles is present for 9,036 of our 12,025
   drugs in OptimusKG, but was not loaded when the graph was built.

2. The unnamed drugs are genuinely unnamed. All 1,163 carry a `name` in
   OptimusKG that IS their ChEMBL accession; 28 have synonyms, none have trade
   names. No lookup will name them. But `type` and
   `maximum_clinical_trial_phase` are present, which is what Team 1's
   "filter on molecule type and clinical phase" actually needed -- so the
   filtering problem is solved even though the naming problem is not.

SCOPE
Run on the user's instruction of 13 October 2026, ahead of the mentor's reply
to the question of whether "focus only on the knowledge graph" permits reading
OptimusKG's own discarded fields. Recorded here so the provenance is visible if
that answer comes back differently.

Usage:
    uv run --with pandas --with pyarrow python scripts/paper2_recover_optimuskg_evidence.py
"""

import os

import pandas as pd

OKG = "data/optimuskg"
OUT = "data/optimuskg_recovered"

DRUG_FIELDS = ["name", "type", "maximum_clinical_trial_phase", "is_approved",
               "has_been_withdrawn", "black_box_warning", "canonical_smiles",
               "inchi_key", "cd_formula", "cd_mol_weight", "calculated_log_p",
               "year_of_first_approval", "description"]
LIST_FIELDS = ["synonyms", "trade_names"]


def prop(p, key):
    return dict(p).get(key) if hasattr(p, "keys") else None


def as_text(v, sep="|"):
    """Arrays become a delimited string; everything else passes through."""
    if v is None:
        return None
    if hasattr(v, "tolist"):
        v = v.tolist()
    if isinstance(v, (list, tuple)):
        return sep.join(str(x) for x in v) if v else None
    return v


def main():
    os.makedirs(OUT, exist_ok=True)

    # ---- drug nodes --------------------------------------------------------
    d = pd.read_parquet(f"{OKG}/drug.parquet")
    out = pd.DataFrame({"Drug_ID": d["id"]})
    for f in DRUG_FIELDS:
        out[f] = d.properties.apply(lambda p, f=f: as_text(prop(p, f)))
    for f in LIST_FIELDS:
        out[f] = d.properties.apply(lambda p, f=f: as_text(prop(p, f)))
    out["has_smiles"] = out.canonical_smiles.notna()
    out.to_csv(f"{OUT}/drug_metadata.csv", index=False)
    print(f"drug nodes recovered: {len(out):,}")
    for c in ["name", "type", "maximum_clinical_trial_phase",
              "canonical_smiles", "is_approved"]:
        print(f"   {c:<30} {out[c].notna().sum():>7,}")

    # ---- drug-gene edge evidence ------------------------------------------
    e = pd.read_parquet(f"{OKG}/drug_gene.parquet")

    def sources(p, which):
        s = prop(p, "sources")
        if s is None:
            return None
        v = dict(s).get(which) if hasattr(s, "keys") else None
        return as_text(v, sep=";")

    def assertions(p):
        a = prop(p, "relation_assertions")
        if a is None:
            return None
        return ";".join(f"{dict(x).get('source')}:{dict(x).get('relation')}"
                        for x in a)

    ev = pd.DataFrame({
        "Drug_ID": e["from"], "Gene_ID": e["to"],
        "relation": e["relation"],
        "sources_direct": e.properties.apply(lambda p: sources(p, "direct")),
        "sources_indirect": e.properties.apply(lambda p: sources(p, "indirect")),
        "relation_assertions": e.properties.apply(assertions),
        "relation_conflict": e.properties.apply(
            lambda p: bool(prop(p, "relation_conflict"))),
        "source_urls": e.properties.apply(
            lambda p: as_text(prop(p, "source_urls"), sep=" ")),
        "mechanisms_of_action": e.properties.apply(
            lambda p: as_text(prop(p, "mechanisms_of_action"), sep="; ")),
    })

    # Evidence level, derived rather than asserted: how many independent
    # databases stand behind the edge, and do they agree.
    ev["n_sources"] = ev.relation_assertions.fillna("").apply(
        lambda s: len([x for x in s.split(";") if x]))
    ev["evidence_level"] = ev.apply(
        lambda r: ("conflicting" if r.relation_conflict else
                   "multi-source" if r.n_sources > 1 else
                   "single-source" if r.n_sources == 1 else "unattributed"),
        axis=1)
    ev["has_reference_url"] = ev.source_urls.notna()
    ev.to_csv(f"{OUT}/drug_gene_evidence.csv", index=False)

    print(f"\ndrug-gene edges recovered: {len(ev):,}")
    print(ev.evidence_level.value_counts().to_string())
    print(f"   with a reference URL      {int(ev.has_reference_url.sum()):,}")
    print(f"   with a stated mechanism   {int(ev.mechanisms_of_action.notna().sum()):,}")
    print(f"\nsources named on these edges:")
    sd = ev.sources_direct.fillna("").str.split(";").explode()
    print(sd[sd != ""].value_counts().to_string())
    print(f"\nwrote {OUT}/drug_metadata.csv")
    print(f"wrote {OUT}/drug_gene_evidence.csv")


if __name__ == "__main__":
    main()
