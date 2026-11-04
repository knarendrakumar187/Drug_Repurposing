# Paper 2 — Task 1: master progression-gene table

This table uses the corrected Paper 1 union of **6,222 genes** (`Paper1_HighConfidence/00_all_progression_genes.csv`, final combined version 10 Oct 2026). No gene was removed.

**Naming note:** columns and classes called "inflammation" hold Paper 1's **ordered disease-group** axis (control/NAFL/NASH labels; not a direct inflammation measurement). The plan's column names are kept so later tasks can read the file; the manuscript should say "disease-group".

Output: `Paper2_progression_genes.csv` (6,222 rows). Script: `Paper2_Results/code/task1_master_gene_table.py`.


## Genes per progression class

| Progression class | Genes | % of 6,222 |
|---|---|---|
| fibrosis_only | 2,023 | 32.5 |
| inflammation_only | 2,377 | 38.2 |
| shared | 1,822 | 29.3 |

Fibrosis-associated genes (fibrosis_only + shared): **3,845**.

## Genes per class x direction (Paper 1 group definitions)

Direction = fibrosis direction for fibrosis_only, inflammation direction for inflammation_only, and agreement of both for shared (as in Paper 1 Task 1).

| Group | Genes |
|---|---|
| fibrosis_only_positive | 1,040 |
| fibrosis_only_negative | 983 |
| shared_positive | 1,408 |
| shared_negative | 384 |
| shared_discordant_fib_negative | 30 |
| shared_discordant_fib_positive | 0 |
| inflammation_only_positive | 2,042 |
| inflammation_only_negative | 335 |

Shared genes split by **fibrosis** direction only: positive 1,408, negative 414.

## Paper 1 annotations available

| Annotation | Genes annotated | Source |
|---|---|---|
| Trajectory cluster (C1–C5) + membership confidence | 3,845 | Paper1_Results/Task2/task2_gene_clusters.csv |
| Largest stage transition + beta + q | 3,845 | Paper1_HighConfidence/01_master_table.csv |
| Tier | 3,845 | Paper1_HighConfidence/01_master_table.csv |
| In 23-gene shortlist | 23 | Paper1_HighConfidence/11_FINAL_shortlist_25_genes.csv |

These annotations exist only for the 3,845 fibrosis-associated genes; inflammation-only genes are blank by design (Paper 1 clustered and scored fibrosis genes only).

Trajectory clusters: C1 857, C2 1,367, C3 226, C4 397, C5 998.
Tiers: Tier 3 1,995, Tier 2 - moderate confidence 1,403, Tier 1 - high confidence 447.

## Check against Paper 1 Task 1

| Group | Direction | Paper 1 | Paper 2 | Match |
|---|---|---|---|---|
| fibrosis_only | all | 2023 | 2023 | True |
| fibrosis_only | positive | 1040 | 1040 | True |
| fibrosis_only | negative | 983 | 983 | True |
| inflammation_only | all | 2377 | 2377 | True |
| inflammation_only | positive | 2042 | 2042 | True |
| inflammation_only | negative | 335 | 335 | True |
| shared | all | 1822 | 1822 | True |
| shared | positive (both up) | 1408 | 1408 | True |
| shared | negative (both down) | 384 | 384 | True |
| shared | discordant (fib up, inf down) | 0 | 0 | True |
| shared | discordant (fib down, inf up) | 30 | 30 | True |

**All counts match Paper 1: YES.**
