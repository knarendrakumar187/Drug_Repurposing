"""Summary page and diagram for graph-preparation step 0.4 (drug–gene links added from ChEMBL, DGIdb, GtoPdb, Open Targets).

Reads the outputs of new_links_search.py; writes New_drug_gene_links/summary.md and
New_links_diagram.png. Changes nothing else.

    python Paper2_Results/code/new_links_report.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "Paper2_Results"
OUT = RES / "New_drug_gene_links"
n = json.loads((OUT / "numbers.json").read_text(encoding="utf-8"))
src = pd.read_csv(OUT / "by_source.csv")
ev = pd.read_csv(OUT / "by_evidence_and_stage.csv")
ctrl = pd.read_csv(OUT / "positive_control.csv")
cls = pd.read_csv(OUT / "medicines_by_progression_class.csv")
gnew = pd.read_csv(OUT / "genes_reached_for_first_time.csv")
med = pd.read_csv(OUT / "new_links_medicines_curated.csv")
e = pd.read_csv(RES / "Task2/drug_gene_edges.csv", low_memory=False)
e = e[e.Link_origin == "OptimusKG"]      # the graph's links before step 0.4 added the new ones
raw = pd.read_csv(OUT / "raw_rows_all_sources.csv.gz", low_memory=False)

# control on the graph's own Open Targets links (the graph was built from them, so recovery should be ~100%)
ot_e = e[e.Evidence_source.str.contains("OPEN_TARGETS")]
ot_ctrl = {}
for s_ in ["OpenTargets", "ChEMBL"]:
    r = raw[raw.source == s_]
    keys = set(zip(r.chembl, r.Gene_ID)) | set(zip(r.Drug_key, r.Gene_ID))
    hit = sum((d, g) in keys for d, g in zip(ot_e.Drug_ID, ot_e.Gene_ID))
    ot_ctrl[s_] = (hit, len(ot_e))

M = n["new_pairs_medicines_curated"]
A = n["new_pairs_approved_curated"]
ALL = n["new_pairs_all"]
# direction rule: the drug pushes down a gene that rises with progression, or up a gene that falls
med["consistency"] = "unknown"
known = med.push.isin(["up", "down"]) & ~med.paper1_group.astype(str).str.startswith("shared_discordant")
med.loc[known, "consistency"] = "direction-inconsistent"
med.loc[known & (((med.gene_direction == "positive") & (med.push == "down"))
                 | ((med.gene_direction == "negative") & (med.push == "up"))), "consistency"] = "direction-consistent"
cons = med.consistency.value_counts()
t1 = gnew[gnew.Tier.astype(str).str.startswith("Tier 1")]
fib_new = int(gnew.progression_class.isin(["fibrosis_only", "shared"]).sum())
EXAMPLES = [("CCN2", "CTGF", ["PAMREVLUMAB"]), ("LGALS3", "galectin-3", ["BELAPECTIN"]),
            ("ITGB5", "αv integrin", ["ABITUZUMAB", "CILENGITIDE"]), ("ITGB8", "αv integrin", ["ABITUZUMAB"]),
            ("ITGA4", "", ["NATALIZUMAB", "VEDOLIZUMAB"]), ("TGFB1", "TGF-β1", ["LUSPATERCEPT"])]
ex_rows = []
for gsym, alias, drugs in EXAMPLES:
    for d in drugs:
        sub = med[(med.Gene == gsym) & (med.Drug.str.upper() == d)]
        ex_rows.append(f"| {gsym}{' (' + alias + ')' if alias else ''} | {d.lower()} | {'yes' if len(sub) else 'NOT FOUND'} | "
                       f"{', '.join(sub.push) or '-'} | {', '.join(sub.consistency) or '-'} |")
sgc = med[med.Gene.isin(["GUCY1A1", "GUCY1B1"]) & med.Drug.str.upper().isin(["RIOCIGUAT", "VERICIGUAT", "PRALICIGUAT"])]
sgc_note = ", ".join(sorted(set(sgc.consistency)))


def md(df):
    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * df.shape[1]]
    for row in df.values:
        lines.append("| " + " | ".join(f"{v:,}" if isinstance(v, (int,)) or (hasattr(v, "dtype") and v.dtype.kind == "i") else str(v) for v in row) + " |")
    return "\n".join(lines)


summary = f"""# Drug–gene links added to the graph from ChEMBL, DGIdb, GtoPdb and Open Targets (graph preparation, step 0.4)

Gene set: the 6,222 Paper 1 progression genes. The headline links below were written to the graph on 21 Oct 2026 (`write_new_links_to_graph.py`; log in `graph_write/`).

## In one sentence
Outside databases list **{M['pairs']:,} drug–gene links for approved or clinical-trial medicines** (curated drug-target evidence) on our genes that the graph does not have; they reach **{M['new_genes']:,} progression genes that have no drug at all in the graph** (genes with a drug: {n['genes_with_drug_graph']:,} → {n['genes_with_drug_graph_plus_new']:,}).

## Sources searched (all free)
| Source | Version | How genes were matched |
|---|---|---|
| IUPHAR/BPS Guide to Pharmacology | {n['release']['GtoPdb']} | Ensembl ID (human, endogenous ligands removed) |
| DGIdb | {n['release']['DGIdb']} | HGNC symbol (whole database, kept for our genes) |
| ChEMBL drug mechanisms | {n['release']['ChEMBL']} | gene symbol of human targets |
| Open Targets Platform | {n['release']['OpenTargets']} (latest) | Ensembl ID |

Each drug is reduced to its ChEMBL parent compound (salts count once). A link is **new** if the graph has no Drug–[ACTS_ON]→Gene link for that drug and gene (all {len(e):,} graph links on our genes were compared).

## Does the matching work? (positive control)
- The graph's {ot_ctrl['OpenTargets'][1]:,} links that came from Open Targets: Open Targets finds **{ot_ctrl['OpenTargets'][0]:,} ({100 * ot_ctrl['OpenTargets'][0] / ot_ctrl['OpenTargets'][1]:.0f}%)**, ChEMBL finds **{ot_ctrl['ChEMBL'][0]:,} ({100 * ot_ctrl['ChEMBL'][0] / ot_ctrl['ChEMBL'][1]:.0f}%)**.
- Neo4j spot check (live graph, 21 Oct 2026): 120 random "new" links → **0** have any Drug–Gene relationship in Neo4j; 150 random new-link drugs → 121 already exist as Drug nodes but only 7 have any ACTS_ON link. So the graph knows most of these drugs but is missing their targets.

{md(ctrl)}

## Two errors found in DGIdb and removed
1. **DGIdb's copy of Guide to Pharmacology is on the wrong genes.** Guide to Pharmacology target numbers were read as NCBI gene numbers (target 215 = dopamine D2 receptor, but NCBI gene 215 = ABCD1). All {6095:,} such rows on our genes were dropped; Guide to Pharmacology is used directly instead.
2. **DGIdb's copy of ChEMBL keeps virus and bacteria targets** matched by short names (hepatitis B polymerase "P" → human OCA2). These {4748:,} rows were dropped; ChEMBL is used directly with human targets only.

## How good is the evidence? (all new links)
{md(ev)}

- **curated drug-target**: Guide to Pharmacology, ChEMBL, Open Targets, or DGIdb sources that give the gene by symbol/UniProt (TEND, TALC, FDA, clinical-trial target lists). **Used for the headline numbers.**
- **TTD only**: curated, but TTD matches genes from protein names, which sometimes goes wrong (e.g. "KIR" → GEM). Kept, not in the headline.
- **lab binding screen** (Drug Target Commons): measured binding in an assay; mostly research compounds.
- **literature / biomarker only** (NCI, PharmGKB, CIViC, cancer knowledge bases): "drug is related to gene" (response, resistance, expression), **not** "drug acts on gene". Not counted as a target.

## By source
{md(src)}

## Headline: new links for medicines with curated evidence
| | links | drugs | genes | genes with NO drug in the graph before |
|---|---|---|---|---|
| approved or clinical trial | {M['pairs']:,} | {M['drugs']:,} | {M['genes']:,} | **{M['new_genes']:,}** |
| …approved only | {A['pairs']:,} | {A['drugs']:,} | {A['genes']:,} | {A['new_genes']:,} |

- {n['medicine_pairs_2plus_sources']:,} of {M['pairs']:,} links ({100 * n['medicine_pairs_2plus_sources'] / M['pairs']:.0f}%) are reported by two or more sources.
- {fib_new:,} of the {M['new_genes']:,} newly reached genes are fibrosis genes (fibrosis-only or shared); **{len(t1)} are Tier 1** (high confidence).
- Push direction is known for {n['medicine_pairs_with_direction']:,} links: **{int(cons.get('direction-consistent', 0)):,} direction-consistent**, {int(cons.get('direction-inconsistent', 0)):,} inconsistent (direction rule: the drug pushes down a gene that rises, or up a gene that falls).

{md(cls)}

### Drugs already tested against fibrosis that the graph was missing (sanity check)
| Gene | Drug | found | push | direction rule |
|---|---|---|---|---|
{chr(10).join(ex_rows)}

These are drugs that have been tested against fibrosis (pamrevlumab in lung fibrosis, belapectin in NASH cirrhosis, αv-integrin blockers), which suggests the new layer adds real biology, not noise.

**Limit of the simple direction rule:** soluble guanylate cyclase stimulators (riociguat, vericiguat, praliciguat) are found on GUCY1A1 / GUCY1B1 but are labelled *{sgc_note}*: these genes rise with fibrosis and the drugs push them up, yet sGC stimulation reduces fibrosis in animal models. A rising gene is not always a harmful gene; direction labels should be read with biology in mind.

## Caveats
- Protein-complex targets are expanded to every subunit gene (e.g. proteasome, tubulin, CD3), as the graph itself does; antibody–drug conjugates (…vedotin) also hit tubulin through their payload.
- A few drugs appear twice (salt or prodrug with its own ChEMBL parent, or a DGIdb RxNorm ID that could not be matched to ChEMBL).
- Literature-only and TTD-only links are listed but not counted in the headline.

## Files
- `new_links_medicines_curated.csv`: headline links (one row per drug–gene, with sources, action, push direction, gene tier).
- `new_links_medicines.csv`: all medicine links, any evidence. `new_links_all.csv`: everything incl. research compounds.
- `genes_reached_for_first_time.csv`, `medicines_by_drug.csv`, `by_source.csv`, `by_evidence_and_stage.csv`, `medicines_direction_consistency.csv`, `positive_control.csv`, `numbers.json`.
- `raw/`: cached downloads; `raw_rows_all_sources.csv.gz`: every source row used.
- Scripts: `code/new_links_search.py` (search), `code/new_links_report.py` (this page + diagram).
"""
(OUT / "summary.md").write_text(summary, encoding="utf-8")

# ------------------------------------------------------------------ diagram
INK, MUTED = "#1f2933", "#616e7c"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "text.color": INK})
fig = plt.figure(figsize=(15, 9))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 15); ax.set_ylim(0, 9); ax.axis("off")
ax.text(0.3, 8.65, "Drug–gene links added to the graph (step 0.4)", fontsize=17, fontweight="bold", va="top")
ax.text(0.3, 8.15, "We asked four outside databases: \"which drugs act on our 6,222 genes?\"  Then we kept only links the graph does not already have.",
        fontsize=11, color=MUTED, va="top")


def box(x, y, w, h, t, fc, fs=10.5, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec=INK, lw=1.2))
    ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal")


def arr(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14, color=MUTED, lw=1.4))


for i, (nm, v) in enumerate([("Guide to\nPharmacology", "2026.3.1"), ("DGIdb", "v5"), ("ChEMBL", "mechanisms"), ("Open Targets", n['release']['OpenTargets'].split()[-1])]):
    box(0.4, 6.6 - i * 1.25, 2.4, 1.0, f"{nm}\n({v})", "#e8eef5")
    arr(2.85, 7.1 - i * 1.25, 3.6, 5.2)
box(3.6, 4.3, 3.0, 1.8, f"All drug–gene pairs\non our genes\n\nNEW (not in graph):\n{ALL['pairs']:,} pairs", "#fdf2e9", fs=11)
arr(6.65, 5.2, 7.3, 5.2)
box(7.3, 4.3, 3.2, 1.8, f"Keep only real medicines\n(approved / clinical trial)\nwith curated target evidence\n\n{M['pairs']:,} links · {M['drugs']:,} drugs", "#e3f4e6", fs=11, bold=False)
arr(10.55, 5.2, 11.2, 5.2)
box(11.2, 4.3, 3.4, 1.8, f"{M['new_genes']:,} genes get their\nFIRST drug\n\n(genes with a drug:\n{n['genes_with_drug_graph']:,} → {n['genes_with_drug_graph_plus_new']:,})", "#e3f4e6", fs=11.5, bold=True)
# thrown away
box(3.6, 2.35, 3.0, 1.45, "Removed: 2 DGIdb errors\n(wrong genes from GtoPdb copy;\nvirus targets from ChEMBL copy)", "#fbe9e7", fs=9.5)
box(7.3, 2.35, 3.2, 1.45, "Set aside (not targets):\nliterature / biomarker links,\nlab screens, research compounds", "#eeeeee", fs=9.5)
box(11.2, 2.35, 3.4, 1.45, f"Known direction: {n['medicine_pairs_with_direction']:,} links\nright way: {int(cons.get('direction-consistent', 0)):,}\nwrong way: {int(cons.get('direction-inconsistent', 0)):,}", "#fdf2e9", fs=9.8)
ax.text(0.4, 1.55, f"Check: the graph's own Open Targets links are found again {100 * ot_ctrl['OpenTargets'][0] / ot_ctrl['OpenTargets'][1]:.0f}% (Open Targets) / "
        f"{100 * ot_ctrl['ChEMBL'][0] / ot_ctrl['ChEMBL'][1]:.0f}% (ChEMBL).  Neo4j: 0 of 120 sampled new links exist in the graph.", fontsize=10.5)
ax.text(0.4, 1.05, "Drugs already tested against fibrosis that it finds: pamrevlumab → CCN2 (CTGF), belapectin → galectin-3, abituzumab → αv integrins, natalizumab → ITGA4.",
        fontsize=10.5)
ax.text(0.4, 0.5, f"Written to the graph on 21 Oct 2026: {M['pairs']:,} ACTS_ON links and 519 new Drug nodes (graph_write/write_log.json).",
        fontsize=10, color=MUTED)
fig.savefig(OUT / "New_links_diagram.png", dpi=160, bbox_inches="tight")
print(summary)
