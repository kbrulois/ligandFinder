#!/bin/zsh
## The arms msa and msa_gap, end to end, unattended (Kevin, 2026-11-12: "set two
## new arms to run over night"; arms/msa.sh and arms/msa_gap.sh say what they are).
##
## It runs every script straight from this checkout (LF_SRC), from the folders
## PACKAGE_MAP.md lists; the pipeline's input lists, blurbs and logs are in
## $LF_RUN (~/AF2_analysis/lf_run, see lf_paths.sh). Do not edit a script in
## this checkout while the run is using it: zsh and Rscript read files as they go.
##
##   zsh inst/scripts/modeling/run_msa_arms.sh > ~/AF2_analysis/run_msa_arms.out 2>&1 &
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
. ${0:A:h:h}/lf_paths.sh
A2=~/AF2_analysis
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

say "=== run_msa_arms: LF_SRC=$W PROCS=$PROCS ==="

## ---- 0. inputs ----------------------------------------------------------------
say "code: $W (git $(git -C $W rev-parse --short HEAD)$(git -C $W diff --quiet || echo ', with uncommitted edits')); run directory: $LF_RUN"
build() {   # $1 fold, $2 sfx, rest: extra flags
  local k=$1 sfx=$2; shift 2
  local out=$A2/lf_pepend_nn_input_cv5_f${k}${sfx}.rds
  [[ -s $out ]] && { say "  have $(basename $out)"; return 0 }
  $R $W/inst/scripts/modeling/10_7_pepend_windows.R $FLAGS --fold $k "$@" --out $out \
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
  $PY_PLM $W/inst/scripts/generate_residue_db/msa_cons/build_msacons_parquet.py --residues-npz $RESID_V2 \
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
    (cd $LF_RUN && ARM=$ARM STAGES=1 PROCS=$PROCS zsh $LF_MOD/pepend_arm_pipeline.sh) >> $A2/run_msa_arms_$ARM.log 2>&1
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
  (cd $LF_RUN && ARM=$ARM STAGES=2,2b,3,4 PROCS=$PROCS zsh $LF_MOD/pepend_arm_pipeline.sh) >> $A2/run_msa_arms_$ARM.log 2>&1 \
    || die "stages 2-4 of $ARM, see run_msa_arms_$ARM.log"
done

## ---- 3. evaluation --------------------------------------------------------------
say "phase 3: never-seen set of the _msa windows"
if [[ ! -s $LF_RUN/cv5_unseen_ends_msa.csv ]]; then
  [[ -s $LF_EVAL/cv5_unseen_sfx.R ]] || die "no cv5_unseen_sfx.R in $LF_EVAL"
  (cd $LF_RUN && INPUT_SFX=_msa $R $LF_EVAL/cv5_unseen_sfx.R) >> $LOG 2>&1 || die "cv5_unseen_sfx.R"
fi
for ARM in msa msa_gap; do
  say "phase 3: evaluation of $ARM on cv5_unseen_ends_msa.csv"
  (cd $LF_RUN && ARM=$ARM UNSEEN_SFX=_msa STAGES=5,6 PROCS=$PROCS zsh $LF_MOD/pepend_arm_pipeline.sh) >> $A2/run_msa_arms_$ARM.log 2>&1 \
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
    $R $LF_VIZ/arm_metrics_compare_plot.R >> $LOG 2>&1 && say "figure: $A2/lf_viz_20261113/cv_metrics_msa_arms.png" \
    || say "  the metrics figure FAILED (the CSVs are fine)"
fi
say "=== done ==="
