# analysis/

One-off questions, one folder per question. Each folder's README says what was
asked and what was found. Scripts here are records: they ran against the files
of their day and nothing else calls them. When inst/scripts was reorganised
(2026-11-13), their `inst/scripts/...` paths were updated, and the `lf_paths`
lines of the scripts that moved; nothing else in them was.

- **A new question** gets a new folder here, named for the question, not the date.
- **A script that becomes part of a pipeline** moves to `modeling/`,
  `evaluation/` or `visualization/`.

| folder | question | found |
|---|---|---|
| `pepend_cv5_review/` | the 5-fold peptide-end model, reviewed for the 2026-10-21 talk (91 scripts) | see its README: never-seen sets, recall vs rank, window-free Shapley, normalised PR |
| `head_ablation/` | attention head vs plain average pooling as the global head | indistinguishable over 10 seeds; pooling is smaller and steadier |
| `cleavage_motifs/` | how much the C-terminal cleavage motif says on its own | amidation motif AP 0.59 held out; a dibasic prior is worse than none |
| `locus_pr/` | scan precision-recall when other known ends are held out | counting them as false positives understates tight candidate sets ~2.4x |
| `scan_ranking/` | how each trunk's full scan ranks held-out ends (ROC, PR, top-N) | not recorded |
| `channel_importance/` | which channels the trunks rely on (permutation, Shapley, xgb gain) | not recorded |
| `seed_vs_ensemble/` | does one member draw the same page as the ensemble | not recorded |
| `per_index_accuracy/` | class accuracy at each window position | not recorded |
| `pepend_vs_db_model/` | peptide-end model vs the retired db-anchored model | not recorded; the peptide-end model replaced it |
| `known_peptide_audit/` | where known peptides drop out between `knowns.rds` and the page's red rings | not recorded |

"Not recorded" means that neither the scripts nor the notes of the time kept a
conclusion. The outputs may still be in `~/AF2_analysis`.

## archive/

This is retired code, kept as a record and not maintained. Delete it only with
Kevin's OK.

- **`window_model/`:** the db-anchored window model, retired in 2026-09/10 when
  the peptide-end model replaced it. It holds the model's scoring
  (`10_1dcnn_new6.R`), benchmarks, per-index profiles, importance and UMAPs.
  Two scripts from 2026-04-15 that built its v2 training sets
  (`make_training_sets_v2.R`, `split_train_val.R`) were loose in
  `~/AF2_analysis` and are here too.
- **`residue_db/`:** residue-db steps that were replaced.
  - `5_integrate_gtp_old.R`
  - `10_score_AA.R` and `10_score_AA_xgboost.R`, ported to `inst/python/lf_resid`
  - `10.1_pep_compare_metrics.R`
  - `export_parquet.R`, an older variant of `generate_residue_db/10_export_parquet.R`
