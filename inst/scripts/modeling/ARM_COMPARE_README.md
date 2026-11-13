# Two arms side by side: the comparison plots

Updated 2026-11-12. The deployed pair is **AAswap** (`curated_AAswap`, 43 ch) and
**ESM-C only** (`esmonly36`, 37 ch). Their comparison plots were first drawn on
2026-11-10, from patched copies of the scripts here. Those changes, and two new
scripts, are now in this directory, so every command below runs from a checkout.

The plots are `~/AF2_analysis/lf_headcmp/*curated_AAswap_vs_esmonly36*`.
**`~/AF2_analysis/lf_headcmp/PLOTS_INDEX_20261110.md`** has a line on every plot in
that folder, the 47 older ones too: what it compares, its takeaway and what is
wrong with it. Read it before quoting a plot.

## The scripts

| script | draws | changed 2026-11-10, and why |
|---|---|---|
| `final_metrics_plot_trunks.py` | `final_model_metrics_trunks_<A>_vs_<B>.png` (table), `final_recall_rank_trunks_<A>_vs_<B>.png` | Set / Term columns were printed only when the first arm was `flatten`, so any other pair got a blank table; they now follow the first arm of `ARMS`. An arm with no fixed colour took `hash(name)`, which Python salts per process, so its colour changed run to run and could equal the other arm's; it now takes the next validated slot that no arm of `ARMS` holds (2026-11-12: the first version of this change still gave two arms without a fixed colour, such as the two `_padg` arms, the same one). AAswap and motifind have fixed colours, AAswap a label. `HEADCMP` names the folder it reads, as for the other stage-6 plots (default `~/AF2_analysis/lf_headcmp`). |
| `segment_T_sweep_plot.R` | `<LEVEL>_T_sweep[_F1]_<ARM1>_vs_<ARM2>_tol<TOL>.png` (`curated` for an unset `ARM1`) | `ARM1`, the first curve. The default `""` is the curated model, so old commands, file names and pixels are unchanged. |
| `shap_positions_plot.py` | `shap_positions_<arms>.png`, `shap_letters_<arm>.png` | An AAswap label instead of the raw token; the N panels' legend moved off the ESM-C -1 spike. |
| `neverseen_examples_pair.R` | `curated_AAswap_vs_esmonly36_neverseen_examples.png` | New. Examples picked by rule in both directions (below). |
| `neverseen_examples_rebuilt.R` | `curated_vs_esmonly36_neverseen_examples.png`; with `ARM1`, `<arm1>_vs_esmonly36_neverseen_examples_rebuilt.png` | New. A rebuild of the lost one-off behind `esmonly36_neverseen_examples.png`; it reproduces that figure except the title. `ARM1=_curated_AAswap` draws AAswap instead of curated, under a `_rebuilt` name: the unsuffixed name is the pair script's fair figure. Its genes (`GENES`) were picked from ESM-C's wins, so it is kept to check the old figure, not to compare arms. |

## The never-seen examples rule

`neverseen_examples_pair.R` was tightened after an independent check of its first
version, whose picks had tilted toward ESM-C. An examples figure is only fair with
all of these:

- **Calls within +/-TOL**, not on the exact anchor. An arm calls an end when its
  best CNN linear-head score (all-50 mean) within 4 residues reaches its threshold.
  Exact-anchor calls made UCN N83 an "AAswap win" that ESM-C makes one residue later.
- **Each arm at its own threshold**, `THRESH_DEFAULT` in `arms/<arm>.sh`: AAswap
  0.5, ESM-C 0.65.
- **Terminus ends counted, not drawn.** An end on the first residue after the signal
  peptide or on the last residue has a window that is about half padding, so a call
  there needs no cleavage-site recognition. 6 of ESM-C's 13 one-arm calls are such
  ends.
- **Unembedded precursors left out**, with whose calls they were. ESM-C scores a
  precursor with no embedding on zeros; the subtitle gives the count (2 ENHO ends,
  both AAswap calls).
- **Picks by rule, both ways**: the `N_EACH` (3) genes with the largest interior
  wins for each arm. A win's size is the winner's margin over its own threshold
  plus the loser's distance below its own.

## Inputs

Nothing upstream needs re-running. Both arms' stage-5 and stage-8 CSVs came from
the same script content at TOL 4 (checked 2026-11-10). The one input made for this
comparison is `anchor_confusion_curated_AAswap_tol4.csv`, from `anchors_roc.R`.

The two never-seen scripts read data only, all under `~/AF2_analysis`: the tracks
rds files, `lf_pepend_residue_cache.rds`, `lf_plm/esm_c_pca36.parquet` and
`lf_run/cv5_unseen_ends.csv` ($LF_RUN). They read nothing from the
source tree, so they run from any directory and need no `LF_SRC`.

`final_metrics.csv` is one shared file, keyed by the arms of the last
`metrics_final.R` run. The two Oct-25 `final_*_trunks_flatten_vs_esmonly36.png`
were drawn while it lacked the ESM-C rows, and show none. Check that it holds both
arms before drawing a pair. The recall-vs-rank panel titles take `n` from the
curated (`flatten`) rows, so a folder without them prints `n = 0`.

## Commands

```bash
cd ~/R_projects/ligandFinder/inst/scripts                              # any checkout
H=~/AF2_analysis/lf_headcmp
PY=~/.virtualenvs/lf-plm/bin/python

# anchor-level confusion for AAswap (ESM-C's existed) -> $H/anchor_confusion_curated_AAswap_tol4.csv
LC_ALL=en_US.UTF-8 TOL=4 ARM=_curated_AAswap Rscript evaluation/anchors_roc.R

# metrics table + recall vs rank
ARMS=curated_AAswap,esmonly36 $PY visualization/final_metrics_plot_trunks.py $H

# threshold sweeps: segment auc_T, segment F1, anchor F1
LC_ALL=en_US.UTF-8 ARM1=curated_AAswap ARM2=esmonly36 Rscript visualization/segment_T_sweep_plot.R
LC_ALL=en_US.UTF-8 ARM1=curated_AAswap ARM2=esmonly36 METRIC=F1 Rscript visualization/segment_T_sweep_plot.R
LC_ALL=en_US.UTF-8 ARM1=curated_AAswap ARM2=esmonly36 METRIC=F1 LEVEL=anchor Rscript visualization/segment_T_sweep_plot.R

# attribution profile
$PY visualization/shap_positions_plot.py $H curated_AAswap esmonly36

# never-seen examples
LC_ALL=en_US.UTF-8 Rscript visualization/neverseen_examples_pair.R
```

Where the output goes:

- **Into `lf_headcmp`, by default.** The plotting R scripts take `OUT=<png>`
  instead; `final_metrics_plot_trunks.py` writes into its argument (into
  `~/AF2_analysis` when it has none).
- **`shap_positions_plot.py` writes into the folder it reads.** It also redraws
  `shap_letters_<arm>.png` (and `shap_props_<arm>.png`, for an arm with the AA
  property channels) for every arm it is given. To keep the old figures, give it a
  copy of the CSVs.

## Which venv

| job | interpreter |
|---|---|
| `final_metrics_plot_trunks.py`, `shap_positions_plot.py` | `~/.virtualenvs/lf-plm/bin/python` |
| R | `/usr/local/bin/Rscript` with `LC_ALL=en_US.UTF-8` |

## Checked 2026-11-12

Each plotting command above was run from this directory with its output sent to
a scratch folder (`anchors_roc.R` was not re-run). Each reproduced its plot in
`lf_headcmp` pixel for pixel; the metrics PNGs matched byte for byte under
`PYTHONHASHSEED` 1, 2 and 3. The letters figure `shap_positions_plot.py` redraws
differs from the Oct-26 `shap_letters_curated_AAswap.png` in its title alone,
which now names the arm "AAswap (43 ch)". With no `ARM1` or `ARMS`,
`segment_T_sweep_plot.R` and `final_metrics_plot_trunks.py` give the same output,
under the same file names, as before the changes. `neverseen_examples_rebuilt.R`
differs from `esmonly36_neverseen_examples.png` only in the title.
