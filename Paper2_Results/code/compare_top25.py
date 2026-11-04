"""Paper 2: compare the top 25 liver-treating drugs across the Tier 1, Tier 2 and Tier 3 runs.

Top 25 = final drugs with rank <= 25 (ties included) in each run's
<run>_BB_liver_treating_drugs_final.xlsx. Writes top25_comparison.csv (one row per drug in any
top 25, with its rank in each run) and top25_comparison.md. No Neo4j needed.

    python Paper2_Results/code/compare_top25.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "Paper2_Results" / "Metapath_runs"
SETS = {"Tier 1": "Run2_onlyconsistent_tier1", "Tier 2": "Run3_onlyconsistent_tier2",
        "Tier 3": "Run4_onlyconsistent_tier3"}

top, liver = {}, {}
for k, run in SETS.items():
    f = pd.ExcelFile(RUNS / f"{run}_BB_liver_treating_drugs_final.xlsx")
    d = f.parse(f.sheet_names[0])
    top[k] = d[d.Rank <= 25].set_index("Drug").Rank
    liver.update(d.set_index("Drug")["Liver diseases it treats"].to_dict())

t = pd.DataFrame(top)
t["in_all_three"] = t.notna().all(axis=1)
t["liver_diseases_it_treats"] = t.index.map(liver)
t = t.sort_values(["in_all_three", "Tier 3", "Tier 2", "Tier 1"], ascending=[False, True, True, True])
t.index.name = "Drug"
t.to_csv(RUNS / "top25_comparison.csv")

s = {k: set(v.index) for k, v in top.items()}
common = s["Tier 1"] & s["Tier 2"] & s["Tier 3"]
fmt = lambda x: "" if pd.isna(x) else str(int(x))
rows = "\n".join(f"| {d} | {fmt(r['Tier 1'])} | {fmt(r['Tier 2'])} | {fmt(r['Tier 3'])} | "
                 f"{'yes' if r.in_all_three else ''} | {r.liver_diseases_it_treats} |" for d, r in t.iterrows())
md = f"""# Top 25 liver-treating drugs: Tier 1 vs Tier 2 vs Tier 3

Top 25 = rank ≤ 25 by BB routes to NAFLD (ties included): Tier 1 {len(s['Tier 1'])}, Tier 2 {len(s['Tier 2'])}, Tier 3 {len(s['Tier 3'])} drugs.

| Overlap | Drugs |
|---|---|
| All three | {len(common)} |
| Tier 1 & Tier 2 | {len(s['Tier 1'] & s['Tier 2'])} |
| Tier 1 & Tier 3 | {len(s['Tier 1'] & s['Tier 3'])} |
| Tier 2 & Tier 3 | {len(s['Tier 2'] & s['Tier 3'])} |

| Drug | Rank Tier 1 | Rank Tier 2 | Rank Tier 3 | In all three | Liver diseases it treats |
|---|---|---|---|---|---|
{rows}
"""
(RUNS / "top25_comparison.md").write_text(md, encoding="utf-8")
print(md)
