# Direction-consistent repurposing candidates

Input: the Task 4 main list (2,587 mechanism-annotated repurposing candidates) and their 8,261 drug–gene links (Task 2, present graph).
Script: `Paper2_Results/code/direction_consistent_drugs.py`.

**Rule:** a link is direction-consistent when the drug switches the gene product **down** and the gene **rises** with progression, or switches it **up** and the gene **falls**. Gene direction = Paper 1 fibrosis direction; disease-group direction for disease-group-only genes.

## Drugs

| Group | Drugs |
|---|---|
| Main list | 2,587 |
| …with ≥1 link of known direction (down/up) | 2,402 |
| …only direction-consistent links | 1,131 |
| …both consistent and inconsistent links (mixed) | 501 |
| …only direction-inconsistent links | 770 |
| …no decided link | 185 |
| **Kept for the metapath analysis (≥1 direction-consistent link)** | **1,632** |

## Links of the main-list drugs

| Direction check | Links |
|---|---|
| not decided | 3,487 |
| direction-consistent | 2,960 |
| direction-inconsistent | 1,814 |

Kept drugs' direction-consistent links: 2,960 on 468 genes.

## Outputs
- `direction_consistent_drugs.csv` / `.xlsx`: the kept drugs (input for the metapath runs).
- `direction_consistent_links.csv`: their direction-consistent drug–gene links.
- `candidate_drugs_direction_summary.csv`, `candidate_links_direction_check.csv`: every main-list drug and link with its check.

A direction-consistent link is a mechanism-consistency check, not evidence that a drug would help (repurposing hypotheses only).
