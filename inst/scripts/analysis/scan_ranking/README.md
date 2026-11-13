# How each trunk's full scan ranks the held-out ends (2026-10-14..15)

**Asked:** how do the held-out known ends rank among every scanned anchor?
- ROC and PR over the whole scan.
- ROC for exactly the series the protein pages draw.
- Recall at the top N.
- The global head's call on held-out candidate windows.

**Found:** not recorded with the scripts. The later version of this question,
with the other known ends held out, is `../locus_pr/`. The maintained metrics are
`evaluation/metrics_final.R`.

**Scripts:**
- `10_9x_auc_summary.py`, `10_9x_auc_plot.R`
- `10_9y_scan_auc_extract.R` -> `10_9y_scan_auc.py` -> `10_9y_scan_auc_plot.R`
- `10_9y_page_roc_scan.py` and `10_9y_page_roc_val.py` -> `10_9y_page_roc_plot.R`
- `10_9y_page_topn_plot.R`
