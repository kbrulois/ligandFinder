# The ESM-C-only arm, and what the defaults now are

Updated 2026-10-24. Read this before passing any env var: **the defaults are now
current practice**, so the bare command is the right command.

## What changed on 2026-10-24

Every script below used to default to the settings of the October 22 curated run.
They now default to what we actually run. The old behaviour is still reachable,
and the output names encode the parameters, so nothing silently overwrites
anything.

| script | was | now | old behaviour |
|---|---|---|---|
| `segments_roc.R` | `TOL=2` | **`TOL=4`** | `TOL=2 Rscript segments_roc.R` |
| `segments_motif.R` | `TOL=2 THRESH=0.3` | **`TOL=4 THRESH=0.5`** | pass both |
| `end_calls.R` | `TOL=2` | **`TOL=4`** | `TOL=2` |
| `metrics_final.R` | `ARMS=attention,flatten` | **`ARMS=flatten,esmonly36`** | pass `ARMS` |
| `segments_roc_plot.py` | `SUF="" TOL=2` | **`SUF=_tol4 TOL=4`** | pass both |
| `final_metrics_plot_trunks.py` | flatten vs attention | **flatten vs esmonly36** | pass `ARMS` |
| the four `*_esmonly36*` / `esm_pages` shells | `PROCS=2` | **`PROCS=4`** | pass `PROCS` |

**Attention is retired.** The flatten head has been THE model since 2026-10-21.
Nothing should default to comparing against the attention trunk again.

## The arm itself

The flatten window head trained on the 36-component PCA of the ESM C per-residue
embedding and nothing else: all 25 curated channels dropped, `padding` kept
because `10_8` refuses to drop the mask channel. **37 input channels**, in the
order `padding, plm_01 … plm_36`.

```bash
cd ~/AF2_analysis/lf_nomask_scripts_20261021

# 1. the embedding (only if esm_c_pca36.parquet is missing)
~/.virtualenvs/lf-plm/bin/python -m lf_plm reduce \
    --h5 ~/AF2_analysis/lf_plm/esm_c_raw.h5 \
    --sequences ~/AF2_analysis/lf_plm/sequences.parquet \
    --k 36 --refit --out ~/AF2_analysis/lf_plm/esm_c_pca36.parquet

./esmonly36_driver.sh                                   # 2. 10 trunks
~/.virtualenvs/r-tensorflow/bin/python residues_esmonly36.py   # 3. 37-ch residue table
./fullscan_esmonly36.sh                                 # 4. 10 full scans
TAG=_nomask_flatten_esmonly36 Rscript tracks_esmonly36.R # 5. tracks rds

# evaluation (defaults are TOL 4 now)
ARM=_esmonly36 Rscript segments_roc.R                   # segment confusion + ROC
ARM=_esmonly36 THRESH=0.5 Rscript segments_motif.R      # by cleavage motif
Rscript metrics_final.R $PWD                            # ROC / normPR / recall, both arms
~/.virtualenvs/lf-plm/bin/python final_metrics_plot_trunks.py ~/AF2_analysis/lf_headcmp

# pages + gallery
./esm_pages.sh                                          # 104 known-peptide pages + index
```

`metrics_final.R` takes the scripts directory as its **first positional
argument** — it is not optional, and without it the run dies on `NA/cv5_unseen_ends.csv`.

## Operating point

ESM-C's segment optimum is **T ≈ 0.65 at tolerance 4**, stable across
cross-validated and never-seen-peptide cells (0.575–0.725). It is the more
liberal scorer — mean score ~1.8× the curated model's — so its best cut sits
above 0.5. Curated's optimum wanders between 0.175 and 0.825 depending on the
cell; T = 0.5 is its optimum in exactly **one of twelve** cells.

Never-seen *large products* optimise at T ≈ 0.025–0.125 for both arms. That is a
degenerate regime, not a finding: neither model separates those ends, so the best
threshold is the one that barely thresholds.

## Five traps, each of which cost real time

1. **The residue table is positional.** `lf_dcnn.scan.build_windows` does
   `x[ok] <- res.feat[rows[ok]]` into a `(k, seq_len, n_channels)` array. A
   wrong-width table raises; a right-width table in the wrong order scans
   silently against the wrong channels and produces plausible nonsense.
   `residues_esmonly36.py` reads the trained model's own `channel_names` and
   refuses to write unless the order it built matches name for name. Keep that.
2. **`10_8` caches by `--plm`, not by `--drop-channels`.** Two arms that differ
   only in dropped channels share a cache, and `lf_dcnn_run_isolated` reuses any
   `member_*.npz` already present on its presence alone — so one trains and
   reports as both. Always pass an explicit, distinct `--out`.
3. **Held-out triangles match a series BY NAME.** `lf_track_heldout_marks` emits
   `"C fold 1"`, which only matches the page's own series. Extra series need
   `options(lf.heldout_extra = ...)` with their own names, which is what
   `final_pages_esm.R` builds. A held-out fold draws triangles along its *whole*
   series, not just at the end residue.
4. **ggiraph builds its SVG in JavaScript at load time.** Grepping the HTML for
   markers finds nothing — the raw file has zero `<polygon>` elements. Measure in
   the live DOM, or you will conclude a working change did nothing.
5. **Never install into `~/.virtualenvs/r-tensorflow`.** Installing matplotlib
   there once upgraded numpy and broke a running render. matplotlib lives in
   `~/.virtualenvs/lf-plm` (3.8.3); use that for every plot script.

## Which venv

| job | interpreter |
|---|---|
| embedding, PCA, every `*_plot.py` | `~/.virtualenvs/lf-plm/bin/python` |
| training, scanning, `residues_esmonly36.py` | `~/.virtualenvs/r-tensorflow/bin/python` |
| R | `/usr/local/bin/Rscript` with `LC_ALL=en_US.UTF-8` |

## Results

`~/AF2_analysis/lf_headcmp/ESMONLY36_RESULTS.md` holds the write-up. Short
version: a wash on cross-validated genes by segment AUC, ahead on never-seen
genes, and **ahead on proteome-wide ranking** — roughly half the anchors to reach
50% recall on the C terminus and on never-seen peptides. It ranks better than it
calibrates, which is why the threshold-based and rank-based views disagree.
