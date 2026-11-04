"""Paper 2 metapath analysis on the present graph: B -> BB -> direct diseases -> ranking, for 4 gene sets.

Drugs: the direction-consistent candidates with ONLY direction-consistent links
(Paper2_Results/Direction_consistent_drugs/, direction_group = "only consistent"). Routes start only
from each drug's direction-consistent targets.

For each gene set (Gene*):
  Run1 Top 23 shortlist | Run2 Tier 1 | Run3 Tier 2 | Run4 Tier 3
  1. metapath B   Drug -ACTS_ON-> Gene -INTERACTS_WITH- Gene* -IN_PATHWAY-> Pathway -IMPLICATED_IN-> NAFLD
                  (scripts/run_six_metapaths.py --metapaths B; rules: correct direction, no hub target
                  (<= 500 partners), pathway size 5-200, drug has no TREATS/OFF_LABEL_FOR to NAFLD/NASH/cirrhosis)
  2. metapath BB  for the drugs with kept B routes:
                  Drug -> Gene1 - Gene2 - Gene* -> Pathway -> NAFLD (scripts/run_two_gene_metapath.py)
  3. follow-up    direct diseases, liver-treating drugs ranked by BB routes, stage table
                  (Paper2_Results/code/metapath_followup.py)

Outputs: Paper2_Results/Metapath_runs/. Password from NEO4J_PASSWORD (prompted if unset).
Run from the repository root:

    python Paper2_Results/code/metapath_pipeline.py            # all four runs
    python Paper2_Results/code/metapath_pipeline.py Run1 Run2  # selected runs
"""
import getpass
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "Paper2_Results" / "Metapath_runs"
OUT.mkdir(parents=True, exist_ok=True)
SRC = ROOT / "Paper2_Results/Direction_consistent_drugs/direction_consistent_drugs.csv"
DRUGS = OUT / "only_consistent_drugs.xlsx"
SHORTLIST = "Paper1_HighConfidence/11_FINAL_shortlist_25_genes.csv"   # holds the final 23 genes
ALL = "Paper1_HighConfidence/11_final_all_3845_genes.csv"
RUNS = {
    "Run1": ("Run1_onlyconsistent_top23", SHORTLIST, None),
    "Run2": ("Run2_onlyconsistent_tier1", ALL, "Tier 1"),
    "Run3": ("Run3_onlyconsistent_tier2", ALL, "Tier 2"),
    "Run4": ("Run4_onlyconsistent_tier3", ALL, "Tier 3"),
}

# drug file in the format the metapath scripts read (one row per drug, correct targets in consistent_targets)
d = pd.read_csv(SRC)
d = d[d.direction_group == "only consistent"].copy()
d["direction"] = "only consistent"
d["stable_drug"] = False     # no ranking step in this analysis
d[["Drug", "Drug_ID", "consistent_targets", "direction", "stable_drug", "n_consistent", "max_phase",
   "is_approved", "liver_contraindicated"]].to_excel(DRUGS, index=False)
print(f"drug file: {len(d):,} drugs with only direction-consistent links -> {DRUGS}")

env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
if not env.get("NEO4J_PASSWORD"):
    env["NEO4J_PASSWORD"] = getpass.getpass("Neo4j password: ")


def step(label, cmd):
    t = time.time()
    print(f"\n----- {label}", flush=True)
    r = subprocess.run([sys.executable] + cmd, cwd=ROOT, env=env)
    if r.returncode:
        raise SystemExit(f"{label} failed (exit {r.returncode})")
    print(f"----- {label} done in {time.time() - t:.0f} s", flush=True)


for key in (sys.argv[1:] or list(RUNS)):
    name, genes, tier = RUNS[key]
    tier_arg = ["--tier", tier] if tier else []
    step(f"{name}: metapath B", ["scripts/run_six_metapaths.py", "--run", name, "--drugs", str(DRUGS),
                                  "--genes", genes, *tier_arg, "--out", str(OUT), "--metapaths", "B"])
    b = pd.read_excel(OUT / f"{name}.xlsx", sheet_name="B")
    if "drug" not in b.columns or b.drug.dropna().empty:
        print(f"{name}: no kept B routes, so no BB run")
        continue
    step(f"{name}: metapath BB", ["scripts/run_two_gene_metapath.py", "--run-file", str(OUT / f"{name}.xlsx"),
                                   "--drugs", str(DRUGS), "--genes", genes, *tier_arg,
                                   "--name", f"{name}_BB", "--out", str(OUT)])
    if not (OUT / f"{name}_BB.xlsx").exists():
        print(f"{name}: no BB routes, so no follow-up")
        continue
    step(f"{name}: follow-up", ["Paper2_Results/code/metapath_followup.py", "--bb", str(OUT / f"{name}_BB"),
                                 "--genes", genes])
print("\nAll requested runs finished.")
