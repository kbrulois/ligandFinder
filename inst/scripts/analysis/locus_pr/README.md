# Scan precision-recall with the other known ends held out (2026-10-15)

**Asked:** what are a scan's precision and recall when the other known ends are
held out instead of counted as false positives? Those are the train ends,
`ligand_list`, and UniProt peptide and propeptide ends.

**Found:** counting them as false positives understates tight candidate sets
about 2.4-fold. The trunks memorised their train ends, so those anchors sit near
the top of the ranking. `load_ends()` in `10_9z_locus_pr.py`
(`(other | train) - val`) is the reference implementation, and any scan PR
should copy it.

`10_9z_locus_pr.py`, `10_9z_locus_pr_plot.R`
