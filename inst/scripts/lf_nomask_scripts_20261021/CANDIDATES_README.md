# Candidate peptide table — how to rebuild it

One command regenerates everything, at any threshold / tolerance:

```bash
./candidate_pipeline.sh                       # the defaults: xgb_pm, T 0.5, tol 4, 8-60 aa
THRESH=0.3 TOL=2 ./candidate_pipeline.sh      # another cut
HEAD=cnn ./candidate_pipeline.sh              # the linear head instead
```

Output `lf_headcmp/candidate_peptides_<HEAD>_T<THRESH>_tol<TOL>.csv`, with a
`.manifest.txt` beside it recording parameters, input sizes and mtimes, tool
versions, the exact mmseqs command and the full column list.

## Stages

| | script | adds |
|---|---|---|
| 1 | `candidate_peptides.R` | the candidates themselves: pair N and C end segments, fan out by ±TOL, `group_id`/`group_rank`/`pair_score`, `end_score_max`, `total_offset`, motifs |
| 2 | `candidate_mmseqs_prep.R` + `mmseqs easy-search` | search against known peptides |
| 3 | `candidate_annotate.R` | `mmseqs_*` best hit and per-class bests, `length_bin`, `overlap_*` |
| 4 | `candidate_multipep.sh` | 20 `mp_*` MultiPep classes + `mp_signalling_max` |
| 5 | `candidate_logistic.R` | `lr_prob`, `lr_prob_oof` |

`WITH_PDEEPPP=1` adds a stage 4b. It is **off by default**: on our own ground
truth PDeepPP's neuro head ranks real ligands backwards (AUC 0.21), behaving
like a length cut near 29 aa. Kept only so that finding stays checkable.

## Caches — why a re-run is cheap

Both external scorers cache **by sequence**, not by candidate, so changing
THRESH or TOL only pays for sequences never seen before:

- `multipep_ext/multipep_scores.csv`
- `pdeeppp_run/pdeeppp_neuro_scores.csv`

## Environments

- `~/.virtualenvs/multipep` — numpy<2, TF, pandas, sklearn, matplotlib
- `~/.virtualenvs/pdeeppp` — torch + fair-esm (ESM-2 650M cached in ~/.cache/torch)
- **Never install into `~/.virtualenvs/r-tensorflow`** — the protein-page render
  uses it, and a numpy upgrade there breaks TensorFlow for every new render
  process.

## Two things the pipeline guards

- **Leakage.** `lr_prob`'s outcome is defined from the `mmseqs_*` columns, so
  every `mmseqs_*` / `overlap_*` column is barred from its predictors;
  `candidate_logistic.R` asserts this.
- **Correlated rows.** Up to (2·TOL+1)² near-identical candidates share one
  segment pair, so CV folds are grouped **by gene**. A random split reports a
  fantasy AUC.
