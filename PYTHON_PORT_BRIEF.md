# Port the 1D-CNN modeling from R/Keras3 to a native Python package

## Task
`~/R_projects/ligandFinder/inst/scripts/10_1dcnn_new6.R` (1120 lines) trains two
Keras models (terms `N` and `C`) that score 36-residue peptide-cleavage windows.
The R code drives Keras through reticulate, so it is already Python semantics in
R syntax. Port **only the modeling** to a real Python package, keep everything
else in R, and let data cross the boundary without manual conversion.

Out of scope: data prep (built by `9.2_add_contact_data.R`), the amidation-motif
block, nearest-known retrieval, and all plotting. Those stay in R.

## Environment (verified — do not create a new venv)
`/Users/kbrulois/.virtualenvs/r-tensorflow/bin/python`
  python 3.11.0 · tensorflow 2.16.2 · keras 3.6.0 · numpy 1.26.4 · pandas 2.2.3
  · scipy 1.14.1 · pyarrow 18.0.0 · **scikit-learn is NOT installed**

This is the same venv R's `keras3` already uses, which is what makes the bridge
seamless. Platt calibration must therefore use `scipy.optimize` or a hand-rolled
IRLS logistic fit, not `sklearn.linear_model.LogisticRegression`.

R side for testing: `/usr/local/bin/Rscript` (R 4.3.1). The `Rscript` first on
PATH is a pixi install with no packages — do not use it.

## What to port, by line number
| lines | what | note |
|---|---|---|
| 2–98 | six `masked_focal_loss` defs | **DEAD — only the last binding survives. Do not port.** |
| 102–155 | `per_index_loss_fn`, `masked_categorical_loss`, `masked_accuracy`, `weighted_binary_crossentropy` | only `weighted_binary_crossentropy` is live |
| 165–255 | class table, position masks, per-class weights | port |
| 257–305 | `make_oversampled_dataset` | port |
| 315–350 | `generate_keras_input` | becomes the array contract, below |
| 356–425 | architecture | port |
| 445–495 | `per_index_cat_loss`, `masked_cat_accuracy` | port |
| 520–570 | compile + fit | port |
| 577–620 | predict + pooled Platt calibration | port |
| 640–680 | `extract_per_index`, embedding extraction | port |
| 690+ | ranks, amidation, plotting | **stays in R** |

## Architecture (Keras 3 functional)
```
inputs (seq_len=36, n_channels)
  -> pos_channel: normalized ramp [0,1] appended as an extra channel
     (cumsum of ones along axis 1, divided by seq_len-1), built in-graph
  -> concat -> (36, n_channels+1)
shared trunk (L2 regularizer 1e-3 on every conv kernel):
  Conv1D(16, k=3, gelu, padding="same") -> Dropout(0.2)
  Conv1D( 8, k=3, gelu, padding="same") -> Dropout(0.3)
per_index_logits = Conv1D(K_pi, k=1)                    # K_pi = 8
per_index_cat    = softmax(per_index_logits)            # OUTPUT "per_index_cat"
masked_sum = per_index_logits[..., :K_cat]              # drop trailing padding logit
             * mask_matrix_cat                          # (36, K_cat) constant 0/1
             -> softmax
global_output = MultiHeadAttention(num_heads=2, key_dim=8)(masked_sum, masked_sum)
  -> LayerNormalization -> GlobalAveragePooling1D
  -> Dense(16, gelu, name="embed")                      # 16-d, tapped later
  -> Dropout(0.3) -> Dense(1, sigmoid, name="global")   # OUTPUT "global"
```
Class order is `["CT_cleavage_context","DB","gap","NT_cleavage_context",
"pep_other","pep_pocket","padding","none"]`; `K_cat = 7` (6 real + none),
`K_pi = 8` (6 real + none + padding, padding LAST in the softmax column order).

Position masks over 36 positions: `NT_mask = 1:5`, `CT_mask = 31:36`,
`mid_mask = 6:30` (1-based, inclusive). `mask_matrix_cat` allows `none`
everywhere; DB at both termini; pep_other/pep_pocket/gap in the middle only.

## Losses
- `global`: weighted BCE, weights currently **1.0 / 1.0** (the oversampler
  already rebalances; the old ~20x positive weight double-corrected).
  `loss_weight = 0.05`.
- `per_index_cat`: focal categorical CE over K_pi, `gamma = 2`, plus a
  smoothness term. `loss_weight = 1`.
  - mask: `keep = 1 - none_flag` — trains real **and** padding classes, masks
    only `none`. (`padding` is a real trained class; its label is leaked by the
    input padding channel, so it is learned for free.)
  - per-class weights via `pi_w_for(term)`: all 1, `DB = 2`, and 3 on
    `CT_cleavage_context` for the C model / `NT_cleavage_context` for the N
    model, then divided by the mean.
  - smoothness: mean squared difference of adjacent-position softmax vectors,
    skipping any adjacency touching a masked position, weight `0.01`.
- metric `masked_cat_accuracy`: argmax accuracy over positions that are neither
  padding nor none.

## Training
- Optimizer Adam, `lr = 3e-4`, `clipnorm = 1.0`.
- Oversampled tf.data pipeline: all positives + `3 x n_pos` negatives sampled
  with replacement, shuffled, batch 32, repeated.
- Noise augmentation, train only: jitter continuous channels
  (`cons_rs, cons_rs_n, min_afm, mean_afm, relASA`) at non-padding positions by
  `N(0, 0.1 * sd_of_that_channel_over_real_positions)`. One-hot/flag and AA_*
  property channels are left intact on purpose.
- The dataset is **rebuilt every epoch** by an `on_epoch_begin` callback so the
  negative sample and the noise are redrawn.
- `steps_per_epoch = ceil(n_pos * 4 / 32)`, `epochs = 2000`.
- EarlyStopping on `val_global_pr_auc`, `mode="max"`, `patience=300`,
  `start_from_epoch=300`, `restore_best_weights=True`.
- `metrics`: AUC and PR-AUC on the global head.

## Data contract (the R <-> Python boundary)
Per term (`"N"`, `"C"`) and per split (`train`, `val`, `all`):
- `x`               float32 `(n, 36, n_channels)`
- `y_global`        float32 `(n, 1)`
- `y_per_index_cat` float32 `(n, 36, K_pi)` — cols `0:K_cat` are the
  `[6 real, none]` one-hot, col `K_cat` (last) is the padding indicator

R builds these today in `generate_keras_input` (line 315) by unlisting
`nn_input[[term]][[split]]$data` into `(36, n_channels, n)` and `aperm`-ing to
`(n, 36, n_channels)`. **Watch the axis order** — R is column-major, numpy is
row-major; a naive `reticulate` handoff of the pre-`aperm` array will transpose
silently and still train, just badly.

Returns to R:
- `pred_raw`  `(n,)`   uncalibrated global score, concatenated in `names(nn_input)` order
- `pred`      `(n,)`   pooled-Platt calibrated
- `per_index` `(n, 36, K_pi)`
- `emb`       `(n, 16)` L2-normalised `embed` layer output

Row order must match `bind_rows(lapply(nn_input, function(x) x$all))`.

## Calibration
ONE logistic calibrator fit on the val scores/labels **pooled across both
models**, then applied to all raw scores. Per-model calibration is unstable
(few val positives). Monotone, so it leaves each model's ROC/PR-AUC unchanged
and only makes the two score scales comparable.

## Bridge
Package `lf_dcnn/` with a clean API (`build_model`, `train`, `predict_all`,
`embed`, `calibrate`). Two entry points:
1. **In-process via reticulate** — `reticulate::use_virtualenv("r-tensorflow")`
   then `source_python()`; numpy arrays pass in memory, no files.
2. **Standalone** — load `.npz` (arrays) + parquet (metadata) from disk and
   write results back the same way, so the package is runnable and testable
   without R at all. This is what makes it a real module rather than an R
   appendage.

## Cleanup, once the Python passes
1. Delete the ported R modeling code from `10_1dcnn_new6.R`, leaving the data
   prep, the bridge call, amidation, retrieval and plotting.
2. Delete the six dead `masked_focal_loss` definitions (lines 2–98) regardless.
3. `git rm` the superseded scripts (~372 KB, 34 files): `10_1dcnn.R`,
   `_new.R`, `_new2..new5`, `_new4.5`, `10_1indexOnly.R`, `10_1dcnn_feb24_restore.R`,
   `10_1dcnn_new_save.R`, `10_3b_add_params.R`, `10_3c_add_scores.R`,
   `selectable_text_test.R`, `leftover_cm_script.R`, `inst/scripts/old/`.
4. Add `.gitattributes`: `inst/scripts/putative_peptide/old/** linguist-vendored`
   — keeps the archive, drops it from the language bar.

Expected result: R ~95.6% -> ~91%, Python ~2.6% -> ~6%.

## Verification
- Python model `count_params()` matches the R model's for the same
  `seq_len`/`n_channels`.
- Same seed + same arrays -> val PR-AUC within noise of the R run.
- Round-trip test: R -> Python -> R returns arrays with the expected shapes and
  row order, and `pred` correlates ~1.0 with a rank of `pred_raw`.
- The standalone `.npz` path reproduces the reticulate path exactly.

## Gotchas
- R is 1-based and column-major; every index and `aperm` above is stated in R
  terms. Translate carefully — off-by-one on the mask rows or the padding column
  will train without error and score wrongly.
- `padding` is the LAST softmax column (`K_pi`), but `none` is at index `K_cat`.
  The loss masks `none` and keeps `padding`; the metric masks both.
- `keras3::clear_session()` is called before each term's model is built; the
  Python port needs the equivalent or the two models will share graph state.
- Don't `cd` into worktrees under `.claude/`; work in
  `/Users/kbrulois/R_projects/ligandFinder` on `master`.
