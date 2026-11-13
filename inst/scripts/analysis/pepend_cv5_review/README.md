# The 5-fold peptide-end model, reviewed (2026-10-18..21)

These 91 scripts come from the session that reviewed the 5-fold (CV5)
peptide-end model for the 2026-10-21 presentation. They were saved from the
session's scratchpad to `~/AF2_analysis/lf_session_scripts_20261021` and moved
here on 2026-11-13. That folder's HTML outputs and slide deck were not code and
stayed behind.

## Outputs
- `~/AF2_analysis/pepend_pages_cv5/`: 103 protein pages and an index.
- The talk's figures, in `~/Desktop/dt_oct20/lf_presentation/`.
- `~/AF2_analysis/pepend_cv5_pca_v4.png`, `pepend_recall_rank.*` and
  `pepend_resid_shapley_importance.*`.

## What it found
- **Never-seen sets.** The 531 annotated peptides give 37/43/45 never-seen ends.
  Length is the main filter: 383 of the 531 are longer than 40 aa (median
  120 aa). Since then "peptide" means shorter than 60 aa.
- **Contamination.** A window-overlap rule replaced the gene-level
  "contaminated" rule. It admits GDF5/GDF6 and 47 other precursors.
  - An earlier check read the `all` split (the scan pool) instead of
    train + val, and wrongly excluded PI3.
- **pep_any vs pocket.** `pep_any` (pocket and other as one class) does worse
  than the separate pocket class. pep_any: AP 0.459, 21x lift. Pocket: AP 0.377,
  96x lift, and 72% recall at 1%.
- **Metrics.** F1 moved from 0.052 to 0.205 on identical predictions just by
  changing the negative pool, so it is not reported. Use prevalence-standardised
  AP at 1% (normalised PR-AUC), reweighting negatives by
  w = n1 (1 - p') / (p' n0). The textbook (AP - p) / (1 - p) is useless at
  p ~ 1e-5.
- **Window-free Shapley** reverses the window model's ranking: AlphaMissense is
  top for pocket (49.5%) and relASA second (36.7%).
- **Recall vs rank** (anchors needed for 50% / 80% recall):

  | set | N ends | C ends |
  |---|---|---|
  | CV | 3,332 / 27,288 | 3,014 / 45,022 |
  | never-seen peptides | 36,850 / 165,498 | 49,764 / 223,492 |
  | never-seen large | 55,005 / 149,728 | 100,310 / 407,569 |

  By motif, against same-class negatives: C amidation 4 / 273 (n = 4),
  C dibasic 1,107 / 4,499, N dibasic 450 / 6,076, C monobasic 7,423 / 18,282
  (n = 3), N monobasic 10,025 / 22,335.

## Scripts, by theme
| theme | scripts |
|---|---|
| recall vs rank | `rr_curves.R`, `rr_plot.R` |
| window sets and splits | `win_stats*.R`, `win_table.py`, `split_table.R`, `splits_audit.R`, `testset_*.R`, `test_*.R` |
| CV5 pages, tracks, PCA | `cv5_*.R`, `pages_nopad.R` |
| attribution | `shap_nopad.R`, `xgbimp_nopad.R`, `xgb_*.R`, `resid_shapley.py`, `resid_shap_plot.R` |
| peaks | `peak_dist.R`, `peak_table.R`, `where_peaks.R` |
| the residue table, probed | `trace_sa*.R`, `probe*.R`, `verify_pos*.R`, `resid_table*.R`, `allres_probe*.R`, `missing_*.R` |
| ASA vs AlphaMissense | `asa_afm_plot.R`, `pep_afm_asa*.R`, `pep_chan_asa.R` |
| the gap and motifs | `gap_*.R`, `gen_gap*.R`, `motif*.R`, `gdibasic_check.R`, `amid.R` |
| figures for the talk | `cnn_arch.py` (an early `visualization/cnn_arch.py`), `aligned_third.R`, `pep_3panel.R`, `offset_*.R` |
