#!/bin/zsh
## The arms msa and msa_gap, end to end, unattended (Kevin, 2026-11-12: "set two
## new arms to run over night"; arms/msa.sh and arms/msa_gap.sh say what they are).
##
## Run it from a CHECKOUT, not from the ~/AF2_analysis working copy: it copies
## this directory's code into the working copy before each phase (the pipeline
## runs there, beside its data files), and a script must not be overwritten
## while zsh is still reading it.
##
##   zsh inst/scripts/lf_nomask_scripts_20261021/run_msa_arms.sh > ~/AF2_analysis/run_msa_arms.out 2>&1 &
##
## Phases, each skipping what already exists, the log in ~/AF2_analysis/run_msa_arms.log:
##   0 inputs   the v2 feature cache (10_7 --rebuild-features, the 2026-11-12
##              residue-join fix), the 10 fold inputs (_msa, _msagap), a check
##              that the two sets differ in the gap labels alone, the v2 curated
##              residue table, and the MSA channel re-masked on the pinned scale
##   1 train    both arms, one after the other (PROCS trunks at a time)
##   2 scans    stages 2, 2b, 3, 4 per arm
##   3 eval     the never-seen set of the _msa windows (cv5_unseen_sfx.R), then
##              stages 5 and 6 per arm into lf_headcmp/unseen_msa
set -u
SRC=${0:A:h}
W=${SRC:h:h:h}
[[ -f $W/DESCRIPTION && -d $W/R ]] || { print -u2 "not inside a ligandFinder checkout: $SRC"; exit 1 }
export LF_SRC=$W
A2=~/AF2_analysis
S=$A2/lf_nomask_scripts_20261021
PROCS=${PROCS:-3}
R=(env LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript)
PY_TF=~/.virtualenvs/r-tensorflow/bin/python
PY_PLM=~/.virtualenvs/lf-plm/bin/python
LOG=$A2/run_msa_arms.log
CACHE=$A2/lf_pepend_residue_cache_v2.rds
RESID_V2=$A2/lf_pepend_residues_v2.npz
PLM_V2=$A2/lf_plm/aa_onehot_msacons_v2.parquet
FLAGS=(--n-folds 5 --near-n 100 --seed 42 --pad-match --mono-match --db-match --amid-match
       --feature-cache $CACHE)

say()  { print -r -- "[$(date '+%F %H:%M:%S')] $*" | tee -a $LOG }
die()  { say "STOPPED: $*"; exit 1 }
## the code of this checkout into the working copy (code only; data untouched)
sync_code() {
  rsync -a --exclude='__pycache__/' --include='*/' --include='*.R' --include='*.sh' \
        --include='*.py' --include='*.md' --exclude='*' $SRC/ $S/ || die "rsync to $S"
  ## rsync adds and replaces but never deletes: a script this checkout removed
  ## (against origin/master) would linger in the working copy and could still be
  ## found there. Move such copies aside, beside the 2026-11-12 archive.
  local rel f gone=$A2/_archived_code_20261112/removed_from_package
  for rel in ${(f)"$(git -C $W diff --name-only --diff-filter=D origin/master -- inst/scripts/lf_nomask_scripts_20261021)"}; do
    f=$S/${rel#inst/scripts/lf_nomask_scripts_20261021/}
    [[ -e $f ]] || continue
    mkdir -p $gone/${${f#$S/}:h} && mv -n $f $gone/${f#$S/} && say "  moved aside (removed from the package): ${f#$S/}"
  done
  say "code synced: $SRC -> $S (git $(git -C $W rev-parse --short HEAD)$(git -C $W diff --quiet || echo ', with uncommitted edits'))"
}

say "=== run_msa_arms: LF_SRC=$W PROCS=$PROCS ==="

## ---- 0. inputs ----------------------------------------------------------------
sync_code
build() {   # $1 fold, $2 sfx, rest: extra flags
  local k=$1 sfx=$2; shift 2
  local out=$A2/lf_pepend_nn_input_cv5_f${k}${sfx}.rds
  [[ -s $out ]] && { say "  have $(basename $out)"; return 0 }
  $R $W/inst/scripts/10_7_pepend_windows.R $FLAGS --fold $k "$@" --out $out \
    > $A2/build${sfx}_f$k.log 2>&1 || { say "  BUILD FAILED $(basename $out), see build${sfx}_f$k.log"; return 1 }
  say "  built $(basename $out)"
}
if [[ ! -s $CACHE ]]; then
  say "phase 0: v2 feature cache (reads secretome + secretome_aa, ~10 min)"
  build 1 _msa --rebuild-features || die "the v2 cache"
fi
say "phase 0: fold inputs"
for k in 1 2 3 4 5; do      # up to four builds at a time (~3 GB each)
  build $k _msa &
  build $k _msagap --gap-class &
  (( k % 2 == 0 )) && wait
done
wait
for k in 1 2 3 4 5; do for x in _msa _msagap; do
  [[ -s $A2/lf_pepend_nn_input_cv5_f$k$x.rds ]] || die "fold input f$k$x missing"
done; done

say "phase 0: the two window sets must differ in the gap labels alone"
$R -e '
  ok <- TRUE
  for (k in 1:5) {
    a <- readRDS(sprintf("~/AF2_analysis/lf_pepend_nn_input_cv5_f%d_msa.rds", k))
    b <- readRDS(sprintf("~/AF2_analysis/lf_pepend_nn_input_cv5_f%d_msagap.rds", k))
    stopifnot(identical(a$params$dssp_index, "mapped"), !isTRUE(a$params$gap_class), isTRUE(b$params$gap_class))
    for (t in c("N", "C")) for (s in c("train", "val", "all")) {
      x <- a$nn_input[[t]][[s]]; y <- b$nn_input[[t]][[s]]
      same <- identical(x$peps, y$peps) && identical(x$known, y$known) && identical(x$data, y$data)
      la <- unlist(x$known_idx); lb <- unlist(y$known_idx); d <- la != lb
      lab_ok <- all(lb[d] == "gap") && all(la[d] %in% c("CT_cleavage_context", "NT_cleavage_context"))
      if (!same || !lab_ok) { ok <- FALSE; message(sprintf("f%d %s %s: windows %s, labels %s", k, t, s, same, lab_ok)) }
      if (s == "train") message(sprintf("f%d %s train: %d windows, %d gap residues", k, t, nrow(y), sum(lb == "gap")))
    }
  }
  if (!ok) quit(status = 1)' >> $LOG 2>&1 || die "the _msa and _msagap inputs differ beyond the gap labels"

if [[ ! -s $RESID_V2 ]]; then
  say "phase 0: v2 curated residue table"
  $R -e '
    ROOT <- Sys.getenv("LF_SRC"); source(file.path(ROOT, "R", "dcnn_bridge.R"))
    invisible(lf_dcnn_python(path = file.path(ROOT, "inst", "python")))
    fc <- readRDS(path.expand("~/AF2_analysis/lf_pepend_residue_cache_v2.rds"))
    stopifnot(identical(fc$dssp_index, "mapped"), identical(names(fc$feats), fc$prec$accession))
    lf_dcnn_pack_residues(fc$prec, fc$feats, fc$all_params3, path.expand("~/AF2_analysis/lf_pepend_residues_v2.npz"),
                          py_path = file.path(ROOT, "inst", "python"))' >> $LOG 2>&1 || die "the v2 residue table"
fi
if [[ ! -s $PLM_V2 ]]; then
  say "phase 0: MSA channel re-masked on the v2 table, scale pinned"
  $PY_PLM $W/inst/scripts/msa_cons/build_msacons_parquet.py --residues-npz $RESID_V2 \
    --scale-from $A2/lf_plm/aa_onehot_msacons.parquet.scale.json --out $PLM_V2 >> $LOG 2>&1 \
    || die "the v2 MSA parquet"
fi

## ---- 1. train ----------------------------------------------------------------
## Training is resumable, and a member can fail on a TensorFlow fault that a
## retry clears ("Incompatible shapes: [0] vs. [32,36]" in masked_cat_accuracy,
## seen in the 2026-11-10/11 runs too; lf_dcnn retries a member 3 times). On
## 2026-11-12 one member failed all three and the whole night stopped at 9 of 10
## trunks, so an arm's training is re-entered up to three times before giving up.
for ARM in msa msa_gap; do
  for attempt in 1 2 3; do
    say "phase 1: training $ARM (attempt $attempt)"
    (cd $S && ARM=$ARM STAGES=1 PROCS=$PROCS zsh ./pepend_arm_pipeline.sh) >> $A2/run_msa_arms_$ARM.log 2>&1
    n=$(ls $A2/lf_pepend_run_[NC]_cv5_f*_nomask_flatten_${ARM}.rds 2>/dev/null | wc -l | tr -d ' ')
    [[ $n -eq 10 ]] && break
    say "  $ARM: $n of 10 trunks after attempt $attempt"
  done
  [[ $n -eq 10 ]] || die "$ARM: $n of 10 trunks after 3 attempts, see run_msa_arms_$ARM.log"
  say "  $ARM: 10 trunks"
done

## ---- 2. residues, groups, scans, tracks ------------------------------------------
for ARM in msa msa_gap; do
  say "phase 2: residues, groups, scans, tracks for $ARM"
  (cd $S && ARM=$ARM STAGES=2,2b,3,4 PROCS=$PROCS zsh ./pepend_arm_pipeline.sh) >> $A2/run_msa_arms_$ARM.log 2>&1 \
    || die "stages 2-4 of $ARM, see run_msa_arms_$ARM.log"
done

## ---- 3. evaluation --------------------------------------------------------------
sync_code    # any fix to the evaluation scripts made since launch
say "phase 3: never-seen set of the _msa windows"
if [[ ! -s $S/cv5_unseen_ends_msa.csv ]]; then
  [[ -s $S/cv5_unseen_sfx.R ]] || die "no cv5_unseen_sfx.R in $S"
  (cd $S && INPUT_SFX=_msa $R cv5_unseen_sfx.R) >> $LOG 2>&1 || die "cv5_unseen_sfx.R"
fi
for ARM in msa msa_gap; do
  say "phase 3: evaluation of $ARM on cv5_unseen_ends_msa.csv"
  (cd $S && ARM=$ARM UNSEEN_SFX=_msa STAGES=5,6 PROCS=$PROCS zsh ./pepend_arm_pipeline.sh) >> $A2/run_msa_arms_$ARM.log 2>&1 \
    || say "  stages 5-6 of $ARM FAILED, see run_msa_arms_$ARM.log (trunks and scans are complete)"
done
M=$A2/lf_headcmp/unseen_msa/final_metrics.csv
if [[ -s $M ]]; then
  $R -e '
    m <- data.table::fread(commandArgs(TRUE)[1])
    print(data.table::dcast(m, term + set + head ~ trunk, value.var = "norm_pr_auc"), digits = 3)' $M >> $LOG 2>&1
  ## the morning figure: normalised PR, ROC, best F1, anchors for 50% recall --
  ## the original AAswap arm for reference, msa the comparator, msa_gap in red
  HEADCMP=$A2/lf_headcmp/unseen_msa ARMS=curated_AAswap,msa,msa_gap \
    LABELS="original windows, Aminode cons. (AAswap); MSA; MSA + gap class" \
    OUT=$A2/lf_viz_20261113 NAME=cv_metrics_msa_arms \
    TITLE="MSA arms (motif-matched negatives, mature anchors, v2 features): full-proteome scan metrics" \
    $R $S/arm_metrics_compare_plot.R >> $LOG 2>&1 && say "figure: $A2/lf_viz_20261113/cv_metrics_msa_arms.png" \
    || say "  the metrics figure FAILED (the CSVs are fine)"
  ## and every comparison figure across the reference arms (Kevin, 2026-11-13:
  ## "make this plot cv_metrics_old_vs_new_conservation. for all these plots, have
  ## it be default that it's produced"): old (Aminode) vs new (MSA) conservation,
  ## with the previous MSA arm and the original curated arm for reference; the
  ## metrics table; per-end ranks, top false positives and tracks for all five
  ## arms on the never-seen ends no arm's negatives touch (_msa AND _padg); the
  ## paired ends; the gap class. Stage 6 has drawn each arm against curated_AAswap.
  F=$A2/lf_viz_20261113
  HEADCMP=$A2/lf_headcmp/unseen_msa ARMS=flatten,curated_AAswap_msacons_padg,curated_AAswap,msa \
    LABELS="original curated arm (Aminode cons.); previous MSA arm (msacons_padg, _padg windows); curated_AAswap (Aminode cons., the site); msa (MSA cons., motif-matched negatives, v2 features)" \
    OUT=$F NAME=cv_metrics_old_vs_new_conservation \
    TITLE="Old (Aminode) vs new (MSA) conservation: full-proteome scan metrics" \
    $R $S/arm_metrics_compare_plot.R >> $LOG 2>&1 || say "  cv_metrics_old_vs_new_conservation FAILED"
  HEADCMP=$A2/lf_headcmp/unseen_msa ARMS=curated_AAswap,msa,msa_gap LABELS="AAswap (site);msa;msa_gap" OUT=$F \
    NAME=cv_metrics_msa_arms_table $R $S/arm_compare_metrics_table.R >> $LOG 2>&1 || say "  the metrics table FAILED"
  ALL5=curated,curated_AAswap,curated_AAswap_msacons_padg,msa,msa_gap
  RETICULATE_PYTHON=$PY_PLM ARMS=$ALL5 UNSEEN_SETS=_msa,_padg OUT=$F/msa_arms_data \
    $R $S/arm_compare_scan_tables.R >> $LOG 2>&1 || say "  arm_compare_scan_tables FAILED"
  DATA=$F/msa_arms_data OUT=$F TAG=msa_arms ARMS=$ALL5 TRACK_ARMS=curated_AAswap,msa,msa_gap \
    $R $S/arm_compare_scan_plots.R >> $LOG 2>&1 || say "  arm_compare_scan_plots FAILED"
  RANKS=$F/msa_arms_data/end_ranks.csv PAIRS=curated_AAswap:msa,msa:msa_gap ARMS=$ALL5 OUT=$F NAME=paired_ranks_msa_arms \
    $R $S/arm_compare_paired_ranks.R >> $LOG 2>&1 || say "  paired ranks FAILED"
  RETICULATE_PYTHON=$PY_PLM ARMS=msa_gap,msa OUT=$F NAME=gap_class_msa_arms \
    $R $S/arm_compare_gap_class.R >> $LOG 2>&1 || say "  gap class FAILED"
  say "figures: $F"
fi
say "=== done ==="
