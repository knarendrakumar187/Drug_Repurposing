"""
paper2_drug_filters.py
----------------------
One definition of "is this entity actually a therapeutic?", shared by every
task that needs it.

WHY THIS MODULE EXISTS
The rule lived inside Task 4 alone, so Task 3 counted a gene as "druggable"
when its only drug was zinc chloride, ATP or an unnamed ChEMBL accession.

WHY IT WAS REWRITTEN, 13 OCTOBER
The first version was too greedy in two ways, both found in review:

  1. SALT AND ESTER WORDS ANYWHERE IN THE NAME.  "ACETATE" flagged ZINC
     ACETATE, which is right, and CORTISONE ACETATE, HYDROCORTISONE ACETATE,
     DEXAMETHASONE ACETATE, PREDNISOLONE ACETATE, MEDROXYPROGESTERONE ACETATE
     and MEGESTROL ACETATE, which is wrong. Those are steroid drugs whose
     ester form is part of the medicine, not a mineral counter-ion.

  2. THE DUPLICATE RULE APPLIED TO EVERYTHING.  Any drug sharing 95% of its
     targets with a flagged drug was flagged too. That caught OMEPRAZOLE,
     LANSOPRAZOLE, RABEPRAZOLE, ERYTHROMYCIN, WARFARIN, LACIDIPINE,
     FLUVOXAMINE, BROMAZEPAM and a dozen antivirals -- about thirty real
     medicines -- because sharing a CYP enzyme panel is not evidence of being
     the same substance. Worse, it cascaded: DEXAMETHASONE was flagged
     because DEXAMETHASONE ACETATE had been.

THE RULE NOW
  * A mineral or simple chemical is judged by the FIRST word of its name. The
    entity in "ZINC ACETATE" is zinc; the entity in "CORTISONE ACETATE" is
    cortisone. Salt and ester words no longer flag anything on their own.
  * Standalone nutrients, cofactors and solvents are matched whole.
  * The duplicate rule now applies ONLY to unnamed accessions, which is what
    it was built for: CHEMBL1201279 shares all 126 of its targets with ZINC
    ACETATE and no name rule can catch it. A named medicine is never flagged
    for the company it keeps.
  * A short whitelist overrides everything, for real medicines whose names
    begin with a mineral: lithium carbonate, magnesium sulfate and the rest.

It is still a heuristic and it still marks rather than deletes. Molecule type
and clinical phase cannot do this job -- ZINC ACETATE is a `Small molecule` at
`max_phase 4.0` because it genuinely is an approved small molecule -- so names
are what is left, and both filtered and unfiltered figures are reported.
"""

import re

# Edge-level: the action describes metabolism or transport, not treatment.
PHARMACOKINETIC = {"ENZYME", "TRANSPORTER", "CARRIER", "SUBSTRATE"}

# Real medicines whose names start with a mineral. Checked first, and they win.
# Real ORGANIC medicines whose names start with a mineral word. Checked before
# the name rule, and they win it.
#
# Carbon is the dividing line, and it is not arbitrary: an entity with no carbon
# is not an organic therapeutic, so it is not a repurposing candidate even when
# it is a genuine medicine. MAGNESIUM SULFATE, POTASSIUM IODIDE and SODIUM
# FLUORIDE are real drugs and are deliberately NOT here -- they were, and Task 4
# flagged them by formula while this whitelist cleared them by name, so Task 8
# listed magnesium sulfate as candidate 24. One rule now governs both.
WHITELIST = {
    "LITHIUM CARBONATE", "LITHIUM CITRATE", "LITHIUM CITRATE ANHYDROUS",
    "LITHIUM SUCCINATE", "MAGNESIUM SALICYLATE", "SODIUM OXYBATE",
    "SODIUM NITROPRUSSIDE", "SODIUM STIBOGLUCONATE", "SODIUM AUROTHIOMALATE",
    "SODIUM PHENYLBUTYRATE", "SODIUM BENZOATE", "SODIUM VALPROATE",
    "CALCIUM FOLINATE", "CALCIUM LEVOFOLINATE", "CALCIUM CARBIMIDE",
    "FERROUS FUMARATE", "BERAPROST SODIUM", "CHOLINE ALFOSCERATE",
    "CHOLINE SALICYLATE", "INOSITOL NIACINATE", "METOCURINE IODIDE",
    "EDETATE CALCIUM DISODIUM ANHYDROUS",
}

# The entity IS the mineral when the name begins with one.
MINERAL_FIRST = re.compile(
    r"^(ZINC|CALCIUM|MAGNESIUM|POTASSIUM|SODIUM|LITHIUM|FERROUS|FERRIC|IRON|"
    r"COPPER|MANGANESE|SELENIUM|CHROMIUM|ALUMINUM|ALUMINIUM|STRONTIUM|"
    r"AMMONIUM|BARIUM|BISMUTH|SILVER|STANNOUS|TITANIUM|BORIC|PHOSPHORIC)\b",
    re.I)

# Nutrients, cofactors, solvents and simple chemicals, matched as whole names
# or as the leading term, never as a suffix on somebody else's drug.
SIMPLE_ENTITY = re.compile(
    r"^(ALCOHOL|WATER|OXYGEN|NITROGEN|NITRIC OXIDE|GLYCEROL|UREA|GLUCOSE|"
    r"SUCROSE|DEXTROSE|LACTOSE|STARCH|MANNITOL|SORBITOL|XYLITOL)\b"
    r"|^(VITAMIN|ASCORBIC ACID|CHOLECALCIFEROL|ERGOCALCIFEROL|TOCOPHEROL|"
    r"ALPHA-TOCOPHEROL|RETINOL|THIAMINE|RIBOFLAVIN|NIACIN|PYRIDOXAL|"
    r"PYRIDOXAMINE|PYRIDOXINE|FOLIC ACID|CYANOCOBALAMIN|BIOTIN|PANTOTHEN)\b"
    r"|^(ADENOSINE TRIPHOSPHATE|ADENOSINE PHOSPHATE|GLUTAMIC ACID|"
    r"ASPARTIC ACID|GLYCINE|ALANINE|SERINE|CYSTEINE|TAURINE|CREATINE|"
    r"CHOLINE|INOSITOL|CARNITINE|ORNITHINE|CITRULLINE)\b"
    r"|NICOTINAMIDE ADENINE DINUCLEOTIDE|DINUCLEOTIDE|^FLAVIN|"
    r"^COENZYME|^HEME$|^HEMIN$", re.I)

UNNAMED = re.compile(r"^CHEMBL\d+$")

# Chemistry, certain where it fires: no carbon means not an organic
# therapeutic. Lives here rather than in Task 4 so every task applies the same
# test; Task 4 is the only one that has formulas, and it passes them in.
CARBON = re.compile(r"C(?![a-z])")


def is_inorganic(formula):
    """True when a formula is present and contains no carbon."""
    if formula is None or (isinstance(formula, float) and formula != formula):
        return False
    return not bool(CARBON.search(str(formula)))

DUPLICATE_OVERLAP = 0.95
DUPLICATE_MIN_TARGETS = 5


def is_non_therapeutic_name(name):
    """Name-based judgement for a single entity. Whitelist wins."""
    n = ("" if name is None else str(name)).strip()
    if n.upper() in WHITELIST:
        return False
    return bool(MINERAL_FIRST.match(n) or SIMPLE_ENTITY.match(n)
                or SIMPLE_ENTITY.search(n) and not MINERAL_FIRST.match(n)
                and SIMPLE_ENTITY.match(n))


def _name_flag(name):
    n = ("" if name is None else str(name)).strip()
    if n.upper() in WHITELIST:
        return False
    if MINERAL_FIRST.match(n):
        return True
    return bool(SIMPLE_ENTITY.match(n)) or bool(
        re.search(r"NICOTINAMIDE ADENINE DINUCLEOTIDE|DINUCLEOTIDE", n, re.I))


def flag_entities(names_by_id, targets_by_id, formulas_by_id=None):
    """Which drug ids are probably not therapeutics.

    names_by_id    {drug_id: name}
    targets_by_id  {drug_id: set of gene ids}
    formulas_by_id {drug_id: molecular formula}, optional. Where a formula is
                   present and carbon-free, the entity is inorganic and the
                   whitelist does not rescue it.

    Returns {drug_id: {"flag_name", "flag_unnamed", "flag_duplicate",
                       "non_therapeutic"}}.
    """
    by_name, by_unnamed = {}, {}
    for did, name in names_by_id.items():
        n = "" if name is None else str(name)
        by_name[did] = _name_flag(n)
        by_unnamed[did] = bool(UNNAMED.match(n.strip()))

    flagged_sets = [targets_by_id[d] for d, f in by_name.items()
                    if f and d in targets_by_id]

    out = {}
    for did in names_by_id:
        s = targets_by_id.get(did)
        dup = False
        # ONLY unnamed accessions. A named medicine is never flagged for
        # sharing a CYP panel with a mineral salt -- that mistake caught
        # omeprazole, erythromycin and warfarin.
        if (by_unnamed[did] and s and len(s) >= DUPLICATE_MIN_TARGETS):
            dup = any(len(s & f) / len(s) >= DUPLICATE_OVERLAP
                      for f in flagged_sets)
        inorg = is_inorganic((formulas_by_id or {}).get(did))
        out[did] = {"flag_name": by_name[did],
                    "flag_unnamed": by_unnamed[did],
                    "flag_duplicate": dup,
                    "flag_inorganic": inorg,
                    "non_therapeutic": (by_name[did] or by_unnamed[did]
                                        or dup or inorg)}
    return out
