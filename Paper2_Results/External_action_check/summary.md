# Drug action directions for OptimusKG links without one (graph preparation, step 0.5)

**Question:** the knowledge graph gives no push direction for 7,919 of 8,948 drug–gene links (TARGET 3,958; ENZYME / TRANSPORTER / CARRIER 3,961). Can other resources supply it?
**Answer: partly.** 297 of the 7,919 unknown links (3.8%) get a clear direction from at least one source with no disagreement. A further 3 get conflicting labels. That would raise known-direction links from **1,001 (11.2%) to 1,298 (14.5%)**.

## Positive control (does the matching work?)
The same lookup was run on the 1,001 links whose direction the graph already knows:

| source | known links found | % found | agree with graph | % agree |
|---|---|---|---|---|
| ChEMBL mechanisms | 516 | 52% | 515 | 100% |
| Guide to Pharmacology | 434 | 43% | 429 | 99% |
| DGIdb | 736 | 74% | 733 | 100% |

The sources find many known links and agree with the graph almost always, so the matching works. **The low fill rate below is real:** most unknown links have no curated direction anywhere.

## Coverage of the unknown links, by source
| index | links given a direction | % of 7,919 unknown links |
|---|---|---|
| ChEMBL mechanisms | 5 | 0.1 |
| Guide to Pharmacology | 267 | 3.4 |
| DGIdb | 226 | 2.9 |
| Any source (agreeing) | 297 | 3.8 |

Sources (all free): **ChEMBL** drug mechanisms (EBI REST API; matched by drug ChEMBL ID, including the parent compound of salts, and target gene symbol); **IUPHAR/BPS Guide to Pharmacology** 2026.3 (public interactions.csv; matched by ligand ChEMBL ID and target Ensembl ID, also gives affinity values); **DGIdb 5** (GraphQL API; combines Guide to Pharmacology, TTD, Drug Repurposing Hub, ChEMBL and others; matched by gene symbol and drug ChEMBL ID or name).

## By graph label
| graph_action_type | down | up | conflicting | still unknown |
|---|---|---|---|---|
| ENZYME/TRANSPORTER/CARRIER | 28 | 0 | 0 | 3,933 |
| TARGET (binding, no direction) | 187 | 82 | 3 | 3,686 |

ChEMBL adds almost nothing because the knowledge graph already contains the ChEMBL / Open Targets mechanisms: those links are the ones that already have INHIBITOR, AGONIST and similar labels. The graph's TARGET links are mostly DrugBank "binds to" entries for which no curated database states a direction. TARGET links (binding with no stated direction) are the ones where external sources help most. ENZYME / TRANSPORTER / CARRIER links mostly describe drug metabolism; an external "inhibitor" label there (e.g. a CYP inhibitor) is real but rarely therapeutic.

## Other resources considered
- **DrugBank** (full action fields) needs an academic licence; the graph already uses its public part.
- **CTD** (Comparative Toxicogenomics Database) gives curated chemical → gene "increases / decreases" statements, but mostly for expression, not protein activity, and is a large download.
- **TTD** is included through DGIdb.
- **ChEMBL activities / BindingDB** give affinity (Ki, IC50), not direction; useful later for potency.

## Caveats
- Directions come from GtoPdb, DGIdb and ChEMBL; they are stored on the graph links as separate external_action_* properties.
- Gene matching in ChEMBL and DGIdb uses gene symbols; Guide to Pharmacology uses Ensembl IDs.
- Labels like "modulator" or "binder" still give no direction and stay unknown.

Per-link detail: `unknown_links_external_lookup.csv`. Raw downloads: `raw/`. Script: `code/external_action_check.py`.
