# lf_dcnn

The 1D-CNN that scores 36-residue peptide-cleavage windows, ported from
`inst/scripts/10_1dcnn_new6.R`. Two models — terminus `N` and terminus `C` —
each carry a **global** ranking head and a supervised **per-index** class head.

Everything upstream (data prep in `9.2_add_contact_data.R`) and everything
downstream (ranks, the amidation-motif block, nearest-known retrieval, plotting)
stays in R. This package owns the modeling only.

## Install / run

No install step: the package is plain Python and lives in `inst/python`. It
targets the same virtualenv R's `keras3` already uses, which is what makes the
reticulate bridge seamless.

```bash
cd inst/python && PYTHONPATH=. /Users/kbrulois/.virtualenvs/r-tensorflow/bin/python -m lf_dcnn selftest
```

## Two entry points

**In-process, via reticulate.** `R/dcnn_bridge.R` does the whole handoff; numpy
arrays pass in memory, nothing touches disk.

```r
source("R/dcnn_bridge.R")
res <- lf_dcnn_run(nn_input, all_params3)      # after 9.2 has built nn_input
nn_input_comb$pred     <- res$pred
nn_input_comb$pred_raw <- res$pred_raw
```

**Standalone, from disk.** Runnable and testable with no R at all.

```bash
python -m lf_dcnn train --input-dir /path/to/in --output-dir /path/to/out
# 5-member ensemble, any Config field overridable as FIELD=JSON
python -m lf_dcnn train --input-dir in --output-dir out --n-seeds 5 --set position_ramp=false
# the same, one python process per member (see below), resumable
python -m lf_dcnn train --input-dir in --output-dir out --n-seeds 5 --isolated
```

**Use `--isolated` for ensembles.** The *second* model trained in one Keras/TF
process can die mid-fit in a retraced `tf.function` — an AUC-metric variable,
an Adam slot or the regularization-loss sum suddenly reads as shape `[0]`, as
if it captured resources the previous model freed. It is intermittent (maybe
one model in five) and has happened under reticulate and in plain python
alike. `--isolated` therefore runs **one model per interpreter**:
`run_member(..., terms=[T])` for every (member, terminus), each written to
`<out>/members/member_NNN_T.{npz,json}`; models already on disk are skipped
(an interrupted run resumes), then `combine()` merges them into the usual
outputs. `python -m lf_dcnn combine` redoes just that last step. `run()` itself
is `run_member()` × n + `combine()` and gives identical numbers in-process, and
a terminus trained alone reproduces the same terminus of a full member; the
selftest checks both.

`--input-dir` holds `arrays.npz`, `config.json` and an optional `meta.parquet`;
the run writes `outputs.npz`, `predictions.parquet` and `history.json`. Write
that directory from R with `lf_dcnn_export()`, or from Python with
`lf_dcnn.io.save_inputs()`.

## Array contract

Per term (`"N"`, `"C"`) and per split (`train`, `val`, `all`):

| name | shape | contents |
|---|---|---|
| `x` | `(n, 36, C)` | float32 input channels |
| `y_global` | `(n, 1)` | float32 0/1 window label |
| `y_per_index_cat` | `(n, 36, K_pi)` | cols `0:K_cat` are the `[6 real, none]` one-hot, col `K_pi-1` the padding flag |

**Do the `aperm` on the R side.** R is column-major and numpy is row-major, so
handing over the pre-`aperm` `(36, C, n)` array transposes silently and still
trains, just badly. `TermArrays.validate()` catches the shape half of that.

Returned, with rows in `term_order` (= `names(nn_input)`, so it lines up with
`bind_rows(lapply(nn_input, function(x) x$all))`):

| name | shape | contents |
|---|---|---|
| `pred_raw` | `(n,)` | uncalibrated global score |
| `pred` | `(n,)` | pooled-Platt calibrated |
| `per_index` | `(n, 36, K_pi)` | per-position softmax |
| `emb` | `(n, 16)` | L2-normalised `embed` layer output |

## Class bookkeeping

Class order is `["CT_cleavage_context", "DB", "gap", "NT_cleavage_context",
"pep_other", "pep_pocket", "padding", "none"]`. The per-index softmax reorders
to `[6 real, none, padding]`, so `K_cat = 7`, `K_pi = 8`, `none` sits at index
`6` and `padding` **last** at index `7`.

The per-index **loss** masks `none` and *keeps* `padding` — padding is a real
trained class whose label is leaked by the input padding channel, so it is
learned for free. The per-index **metric** masks both.

## Two switches where the R script diverges from its own design

`Config` defaults implement the intended architecture. `Config.r_exact()` (or
`--r-exact`) flips both back to reproduce `10_1dcnn_new6.R` bit-for-bit — use it
when comparing validation PR-AUC against an R run.

| switch | default | `r_exact` | what the R script does |
|---|---|---|---|
| `db_mask_both_termini` | `True` | `False` | The NT and CT masks are assigned to overlapping column sets in sequence, so the second assignment overwrites `DB` with the CT mask alone. `DB` ends up allowed only at positions 31–36, not at both termini. |
| `resample_each_epoch` | `True` | `False` | The `on_epoch_begin` callback rebinds `train_ds`, but `fit` already holds the original `tf.data` object, so the redraw never reaches training. The negative sample and the augmentation noise are drawn once for the whole run. |

Neither is an error the R would report; both change what the model sees.

## Layout

| file | what |
|---|---|
| `config.py` | `Config`: geometry, class bookkeeping, position masks, per-class weights, hyper-parameters |
| `model.py` | `build_model` / `compile_model` / `embed_model` and the two custom layers |
| `losses.py` | weighted BCE, focal per-index CE + smoothness, masked categorical accuracy |
| `data.py` | the array contract and the oversampled, noise-augmented sampler |
| `calibrate.py` | pooled Platt calibration (IRLS mirroring R's `glm`; scipy cross-check) |
| `pipeline.py` | `train` / `predict_all` / `embed`; `run_member` + `combine` = `run` |
| `io.py` | the `.npz` + parquet exchange format |
| `cli.py` | `python -m lf_dcnn train｜selftest` |
| `synthetic.py` | contract-shaped fake windows for the tests |
| `selftest.py` | the checks, runnable without pytest |

## Trunk architectures

`Config.trunk` selects what sits between the input and the per-residue logits.
Everything above the trunk is identical either way, so the two are directly
comparable.

| | `"unet"` (default since 2026-09-18) | `"flat"` (the R model; `Config.r_exact()`) |
|---|---|---|
| resolution | 36 → 18 → 9 → 18 → 36 | 36 positions throughout |
| pooling | `MaxPooling1D` down, `UpSampling1D` + skip concat back up | none |
| channels | 16 → 32 → 64 → 32 → 16 | 16 → 8 |
| parameters | 21,406 (21,454 with the position ramp) | 2,390 (2,438 with the ramp) |
| receptive field | most of the window | 5 positions |

```python
Config()                                                        # unet@16-32-64, no position ramp
Config(trunk="unet", unet_filters=(16, 32, 64), unet_dropout=(0.2, 0.2, 0.2))   # the same, spelled out
Config.r_exact()                                                # the pre-port R model
```

The default was set on `inst/scripts/10_5_benchmark_window_model.R` (U-Net
with vs without the position input, 5 seeds each, validation ROC/PR, top-hit
violins and rank agreement); rerun it to revisit.

`unet_filters` is one entry per level with the **last entry the bottleneck**, so
`(16, 32, 64)` means two pooling steps. Every step must divide the window
exactly — at `seq_len=36`, `pool_size=2` allows at most two levels
(36 → 18 → 9); a third would need 9/2 and `Config` rejects it up front rather
than failing at build time.

The decoder concatenates the matching encoder output at each level. Without
those skips the per-residue head only sees upsampled 9-position features and
returns blocky class boundaries — the fine positional detail that locates a
cleavage site to the residue is exactly what the pooling discards.

**The parameter count is the thing to watch.** The labelled set is tiny — 11
and 16 training positives for the N and C models — and the U-Net is ~9x the
flat trunk against the same handful, so any change to it should be judged on
held-out PR-AUC across several seeds (the benchmark script does exactly that),
not on a single run. `unet_filters=(8, 16, 32)` (~5.5k params) is the middle
ground if it ever overfits.

## The `none` class

The per-index softmax carries a `none` column for background positions. Until
2026-09-25 those positions were masked out of the loss, which had a
consequence worth knowing: a NEGATIVE window is all `none`, so it contributed
nothing at all to the per-index loss, and the `none` column never received a
positive gradient. Measured on the trained model it was effectively dead —
mean probability 0.007 at scored positions, 0.057 at background, and never the
argmax at any of 1,065,564 positions.

`Config.none_in_loss` (default `True` since 2026-09-25) trains on them
instead, with `pi_weight_none` (default 0.1) keeping the
75%-of-in-batch-positions majority from swamping the six real classes. The weight is applied after the mean-1
normalisation, so every real class keeps exactly the weight the masked model
trained with. `Config.r_exact()` pins `none_in_loss=False`.

On the C terminus this saturates validation PR AUC (1.000 vs .876 ± .066
masked), lifts held-out known peptide ends (35 vs 31 of 52 in the top 500) and
tightens the ensemble (member sd 0.11 vs 0.34 among the top 100), at a cost of
~5 points of per-residue real-class accuracy. A 0.1/0.3/0.6 sweep is within
noise on every summary metric, so the weight is not tuned; 0.1 is the cheapest
setting that gets the benefit. It stays off by default because the arms rank
CANDIDATES quite differently despite matching on every summary metric (Spearman
0.27 between w=0.1 and w=0.6; ANO8_w12-47 goes from candidate rank 282 masked to
1321 at w=0.1), and the evidence is C-terminus only.

## The position input

The windows are anchored and aligned, so absolute position within the window
is meaningful — but conv is translation-equivariant. `Config(position_ramp=True)`
appends an in-graph `[0, 1]` ramp (`PositionRamp`) as a 27th channel; the
default (`False`, since 2026-09-18) feeds the trunk the 26 raw channels only.
The first conv gains or loses exactly `conv_kernel × filters` weights
(21,406 ↔ 21,454 for the U-Net, 2,390 ↔ 2,438 flat) and nothing else moves.

Whether the ramp earns its place is an empirical question that
`inst/scripts/10_5_benchmark_window_model.R` answers on the real windows: both
arms, `n_seeds` members per terminus, validation ROC/PR quoted as mean ± sd
across seeds, plus the top-hit violin and rank agreement, written to one html.
`Result.val_scores_members` / `pred_raw_members` (`(n_seeds, m)` / `(n_seeds, n)`)
exist for that: the AUC of the averaged score is a different, usually higher,
number than the average of the per-seed AUCs, and a benchmark needs the latter.

## Testing

Three levels, cheapest first:

```bash
# 1. contract + numerics, synthetic data, ~1 min
cd inst/python && PYTHONPATH=. python -m lf_dcnn selftest

# 2. the R <-> Python boundary: axis order, one-hot layout, row order,
#    and that the standalone .npz path reproduces the in-process one
Rscript inst/python/tests/roundtrip.R

# 3. compare trunk/head variants on real windows (2 seeds, python-only)
Rscript inst/python/tests/compare_r_python.R --terms C --arms flat@39-19,unet_bneck
```

Level 3 trains variants on identical arrays with the same seed and prints
validation PR-AUC, ROC-AUC and pooled per-residue accuracy side by side.
Roughly a minute per arm/terminus/seed.

The frozen R reference is opt-in (`--arms r,flat`), not a default: parity is
established -- identical 2,438 parameters, per-residue accuracy within seed
noise -- so re-running it every sweep buys nothing. Bring it back if you change
the losses, the sampler or the array contract.

Run a sweep as one process per (arm, seed) with `--out`. A long run followed by
another model in the same process can die inside a retraced `tf.function`
summing an empty regularization-loss list; process isolation avoids it.

### Live training curves

R's keras3 viewer only draws when R drives `fit()`, which it no longer does.
The Python stand-in is TensorBoard — set `Config.tensorboard_dir` (or pass
`--tensorboard <dir>` to the comparison harness, which also logs the R
reference arm, so all arms land in one place):

```bash
Rscript inst/python/tests/compare_r_python.R --terms C --tensorboard ~/AF2_analysis/lf_dcnn_tb/run1
tensorboard --logdir ~/AF2_analysis/lf_dcnn_tb/run1 --reload_interval 15
```

Then open <http://localhost:6006> and filter tags. Start TensorBoard *after*
the first epochs have been written — pointed at an empty directory it caches
"No scalar data was found" and will not pick the run up later.

### Per-residue accuracy

`masked_cat_accuracy` is argmax accuracy over positions that are neither
`padding` nor `none`. **Do not read it off `model.evaluate()` or off the
`val_*` curve in TensorBoard.** It is a stateless metric, so Keras averages it
per batch — and validation batches mostly hold no known peptide at all, so
every position is `none`, the mask zeroes the batch, and it scores
`0/(0+eps) = 0`. Averaging those in scales the number down by the fraction of
empty batches, which is most of them. Training batches are oversampled to ~25%
positives, so they do not have this problem, and the train/val gap you see in
TensorBoard is therefore mostly an artifact.

The comparison harness reports the pooled figure instead — numerator and
denominator summed over the whole split — alongside the residue count it was
computed over (`n_res_train` / `n_res_val`).

Read the spread, not a single number: with 6–7 validation positives per
terminus, PR-AUC moves ~0.4 across seeds *within* either implementation, so a
one-seed difference is noise. Compare the between-arm gap against the
within-arm spread the summary table prints. ROC-AUC is far stabler at these
positive counts.

scikit-learn is deliberately not used — it is not installed in the target venv,
which is why calibration is a hand-rolled IRLS rather than
`LogisticRegression`.
