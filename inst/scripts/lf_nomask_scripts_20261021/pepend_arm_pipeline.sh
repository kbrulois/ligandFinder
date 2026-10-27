#!/bin/zsh
# One peptide-end "arm", end to end: extra channels -> trunks -> scans ->
# tracks -> evaluation -> plots -> pages -> per-position attribution (Kevin,
# 2026-10-24: "I want all this as a modifiable pipeline as we will be rerunning
# these trainings and visualizations repeatedly").
#
# AN ARM IS A CHANNEL SET, and it lives in ONE FILE: arms/<ARM>.sh, which says
# what the extra per-residue channels are, which curated channels to drop, how
# many channels the trunk should end up with, and how to build the two things
# that depend on the channel set (the parquet and the residue table). Everything
# else -- model names, scan npz, tracks, result CSVs, plots, page directory --
# derives from $ARM, so a new arm collides with nothing already on disk.
#
#   ./pepend_arm_pipeline.sh                            # the default arm, all stages
#   ARM=esmonly36 ./pepend_arm_pipeline.sh              # the ESM-C-only arm
#   STAGES=5,6 ARM=curated_wAAid ./pepend_arm_pipeline.sh    # re-evaluate and re-plot
#   STAGES=shap ./pepend_arm_pipeline.sh                # just the attribution
#
# TO ADD AN ARM: copy arms/curated_wAAid.sh, change PLM / DROP / N_CHAN_EXPECT
# and the two mk_ functions, and run. Nothing else in this directory changes.
#
# PARAMETERS (all have working defaults -- the bare command is the right command)
#   ARM       arm name; must have arms/$ARM.sh            curated_wAAid
#   THRESH    segment threshold for the evaluation        0.5
#   TOL       segment tolerance, in residues              4
#   PROCS     parallel jobs                               4
#   STAGES    comma list of numbers or names, or "all"    all
#   SHAP_ARMS which arms to attribute in stage 8          "$ARM curated"
#
# STAGES (number or name; both are accepted)
#   0 channels  build the arm's extra per-residue parquet (skipped if it exists)
#   1 trunks    train 5 folds x 2 termini                 (skipped per trunk)
#   2 residues  the arm's own residue table + ORDER CHECK (skipped if it exists)
#   3 scans     full-proteome scans, 10 of them           (skipped per npz)
#   4 tracks    assemble the tracks rds                   (skipped if it exists)
#   5 eval      segment confusion/ROC, motif, anchor metrics  (always re-runs)
#   6 plots     threshold scan, motif table, metrics       (always re-runs)
#   7 pages     per-gene scans, known-peptide pages, gallery (skipped per page)
#   8 shap      per-residue-position attribution + plots   (always re-runs)
#
# Stages 0-4 and 7 are resumable and skip finished work, so re-running after an
# interruption is cheap. 5, 6 and 8 are cheap and always re-run, which is what
# you want when only THRESH or TOL changed.
#
# WHY A SEPARATE RESIDUE TABLE PER ARM: lf_dcnn.scan.build_windows copies
# res.feat into a (k, seq_len, n_channels) array POSITIONALLY. A table of the
# wrong width raises; one of the right width in the wrong order scans silently
# against the wrong channels. Stage 2 runs check_resid_order.py, which compares
# the table against the trained trunk's own channel_names, name for name.
set -u
S=${0:A:h}
A2=~/AF2_analysis
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
R="LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript"
PY_PLM=~/.virtualenvs/lf-plm/bin/python
PY_TF=~/.virtualenvs/r-tensorflow/bin/python

export ARM=${ARM:-curated_wAAid}

## the arm file comes FIRST, because the segment threshold is a property of the
## arm (ESM-C's optimum is 0.65, the curated model's is nowhere near it) and
## belongs with the channel set, not in this script's defaults.
[[ -r $S/arms/$ARM.sh ]] || {
  echo "no arm definition at arms/$ARM.sh; have:"; ls $S/arms/*.sh 2>/dev/null | sed 's|.*/|  |'
  exit 1; }
source $S/arms/$ARM.sh

export THRESH=${THRESH:-${THRESH_DEFAULT:-0.5}}
export TOL=${TOL:-4}
export PROCS=${PROCS:-4}
STAGES=${STAGES:-all}
THR_CUR=${THR_CUR:-0.2}      # the curated arm's chip threshold in the gallery index
TAG=_nomask_flatten_${ARM}
PAGES=$A2/pepend_pages_${ARM}
RESID=$A2/lf_pepend_residues_${ARM}.npz
SHAP_ARMS=${SHAP_ARMS:-"$ARM curated"}

## a stage runs if STAGES is "all", or names its number, or names its name
run() {
  [[ $STAGES == all ]] && return 0
  [[ ",$STAGES," == *",$1,"* ]] && return 0
  [[ ",$STAGES," == *",$2,"* ]]
}
say() { echo "[$(date '+%H:%M:%S')] $*" }

echo "=== PEPEND ARM PIPELINE  ARM=$ARM ($ARM_LABEL) ==="
echo "    THRESH=$THRESH TOL=$TOL PROCS=$PROCS STAGES=$STAGES"
echo "    channels : $N_CHAN_EXPECT expected"
echo "    plm      : $PLM"
echo "    tag      : $TAG"
echo "    pages    : $PAGES"

# ---- 0. the arm's extra per-residue channels -------------------------------
if run 0 channels; then
  if [[ -s $PLM ]]; then say "stage 0: $(basename $PLM) exists, skipping"
  else
    say "stage 0: building $(basename $PLM)"
    mk_plm || { say "STAGE 0 FAILED"; exit 1; }
  fi
fi

# ---- 1. trunks -------------------------------------------------------------
if run 1 trunks; then
  say "stage 1: training trunks (resumable)"
  $S/train_arm.sh || { say "STAGE 1 FAILED"; exit 1; }
  n=$(ls $A2/lf_pepend_run_[NC]_cv5_f*${TAG}.rds 2>/dev/null | wc -l | tr -d ' ')
  [[ $n -eq 10 ]] || { say "STAGE 1: only $n/10 trunks, stopping"; exit 1; }
fi

# ---- 2. residue table ------------------------------------------------------
if run 2 residues; then
  if [[ -s $RESID ]]; then say "stage 2: $(basename $RESID) exists, checking its order"
  else
    say "stage 2: building the arm's $N_CHAN_EXPECT-channel residue table"
    mk_resid || { say "STAGE 2 FAILED"; exit 1; }
  fi
  $PY_PLM $S/check_resid_order.py --arm $ARM --plm $PLM --residues-npz $RESID \
    || { say "STAGE 2: residue table does not match the trunks"; exit 1; }
fi

# ---- 3. full scans ---------------------------------------------------------
if run 3 scans; then
  say "stage 3: full-proteome scans (resumable)"
  $S/fullscan_arm.sh || { say "STAGE 3 FAILED"; exit 1; }
fi

# ---- 4. tracks -------------------------------------------------------------
if run 4 tracks; then
  say "stage 4: assembling the tracks rds"
  eval TAG=$TAG $R $S/tracks_esmonly36.R || { say "STAGE 4 FAILED"; exit 1; }
fi

# ---- 5. evaluation ---------------------------------------------------------
if run 5 eval; then
  say "stage 5: segment confusion + ROC, by motif, and the anchor ranking metrics"
  eval ARM=_$ARM TOL=$TOL $R $S/segments_roc.R   > $S/pipe_segroc_$ARM.log 2>&1 \
    || { say "  segments_roc FAILED, see pipe_segroc_$ARM.log"; }
  eval ARM=_$ARM TOL=$TOL THRESH=$THRESH $R $S/segments_motif.R > $S/pipe_segmotif_$ARM.log 2>&1 \
    || { say "  segments_motif FAILED, see pipe_segmotif_$ARM.log"; }
  ## lf_headcmp/final_metrics.csv is ONE SHARED FILE keyed by the ARMS list:
  ## an arm left off a run DISAPPEARS from it. That has already happened once,
  ## when a sibling session re-ran metrics for its own arm and silently dropped
  ## another's rows. So the default is every arm that has a tracks rds on disk,
  ## not just this one. Override with METRIC_ARMS to narrow it deliberately.
  found=(flatten)
  for t in $A2/lf_pepend_cv5_tracks_C_flatten_*.rds(N); do
    a=${${t:t:r}#lf_pepend_cv5_tracks_C_flatten_}
    [[ -n $a ]] && found+=($a)
  done
  MARMS=${METRIC_ARMS:-${(j:,:)found}}
  say "  metrics arms: $MARMS"
  eval ARMS=$MARMS $R $S/metrics_final.R $S > $S/pipe_metrics_$ARM.log 2>&1 \
    || { say "  metrics_final FAILED, see pipe_metrics_$ARM.log"; }
  say "stage 5 done"
fi

# ---- 6. plots --------------------------------------------------------------
if run 6 plots; then
  say "stage 6: plots"
  O=$A2/lf_headcmp
  SUF=_${ARM}$([[ $TOL == 2 ]] && echo "" || echo "_tol$TOL")
  ## stage 6 plots whatever stage 5 wrote AT THE SAME THRESH/TOL. Check first and
  ## name the stage that produces each file, instead of letting the plot script
  ## raise FileNotFoundError on a path the reader has to decode.
  miss=0
  for f in $O/segment_confusion${SUF}.csv $O/segment_motif_T${THRESH}_${ARM}_tol${TOL}.csv; do
    [[ -s $f ]] || { say "  MISSING $(basename $f) -- run STAGES=5 at THRESH=$THRESH TOL=$TOL first"; miss=1 }
  done
  [[ $miss -eq 0 ]] || { say "stage 6: inputs missing, not plotting"; exit 1 }
  SUF=$SUF TOL=$TOL $PY_PLM $S/segments_roc_plot.py $O          >> $S/pipe_plots_$ARM.log 2>&1
  ARM=_$ARM $PY_PLM $S/segments_motif_plot.py $THRESH $TOL      >> $S/pipe_plots_$ARM.log 2>&1
  TRUNK=$ARM $PY_PLM $S/final_metrics_plot.py $O                >> $S/pipe_plots_$ARM.log 2>&1
  ARMS=flatten,$ARM $PY_PLM $S/final_metrics_plot_trunks.py $O  >> $S/pipe_plots_$ARM.log 2>&1
  say "stage 6 done; see pipe_plots_$ARM.log"
fi

# ---- 7. pages --------------------------------------------------------------
if run 7 pages; then
  say "stage 7: per-gene scans, pages and the gallery index"
  $S/cv5_post_esmonly36.sh > $S/pipe_cv5post_$ARM.log 2>&1 \
    || { say "  cv5_post FAILED, see pipe_cv5post_$ARM.log"; exit 1; }
  $S/esm_pages.sh || { say "STAGE 7 FAILED"; exit 1; }
  ## the TP/FN gallery index needs end calls for BOTH arms, each at its own
  ## threshold. They are cheap next to the pages, so make them here rather than
  ## leaving the index to a manual step that gets forgotten.
  say "stage 7: end calls for both arms, then the TP/FN index"
  eval THRESH=$THR_CUR TOL=$TOL $R $S/end_calls.R              >> $S/pipe_endcalls_$ARM.log 2>&1
  eval THRESH=$THRESH ARM=_$ARM TOL=$TOL $R $S/end_calls.R     >> $S/pipe_endcalls_$ARM.log 2>&1
  ## the index is PACKAGE code now (R/gallery_tpfn.R + inst/scripts/
  ## make_tpfn_index.R); an arm is an lf_tpfn_arm() entry in
  ## lf_tpfn_default_arms(), so adding one does not touch this pipeline.
  eval $R $W/inst/scripts/make_tpfn_index.R --dir $PAGES --tol $TOL \
       --blurb $S/final_blurb.html --footer $S/final_footer.html \
    >> $S/pipe_index_$ARM.log 2>&1 || say "  make_tpfn_index FAILED, see pipe_index_$ARM.log"
fi

# ---- 8. per-residue-position attribution -----------------------------------
if run 8 shap; then
  say "stage 8: per-position attribution for: $SHAP_ARMS"
  ## expected gradients on the trunk's logit -- the Shapley value of each
  ## (position, channel) cell with the arm's own candidate pool as the
  ## background. "curated" means the baseline arm, whose suffix is empty.
  for A in ${=SHAP_ARMS}; do
    [[ $A == curated ]] && AA="" || AA=$A
    say "  attributing arm '${AA:-curated}'"
    ## the flatten head reads the window, so it gets per-POSITION attribution;
    ## the xgb head reads 48 pooled summaries of the trunk surface, so it gets
    ## exact TreeSHAP over those, four of which are positional regions. Both,
    ## because the comparison the user asked for covers both heads.
    $PY_TF $S/shap_positions.py --arm "$AA" >> $S/pipe_shap_$ARM.log 2>&1 \
      || { say "  shap_positions FAILED for '${AA:-curated}', see pipe_shap_$ARM.log"; }
    $PY_TF $S/shap_xgb_head.py  --arm "$AA" >> $S/pipe_shap_$ARM.log 2>&1 \
      || { say "  shap_xgb_head FAILED for '${AA:-curated}', see pipe_shap_$ARM.log"; }
  done
  $PY_PLM $S/shap_positions_plot.py $A2/lf_headcmp ${=SHAP_ARMS} \
    >> $S/pipe_shap_$ARM.log 2>&1 || say "  shap plot FAILED, see pipe_shap_$ARM.log"
  say "stage 8 done; see pipe_shap_$ARM.log"
fi

say "=== PIPELINE DONE  ARM=$ARM ==="
echo "    trunks : $(ls $A2/lf_pepend_run_[NC]_cv5_f*${TAG}.rds 2>/dev/null | wc -l | tr -d ' ')/10"
echo "    scans  : $(ls $A2/lf_winxgb_fullscan_[NC]_cv5_f*${TAG}.npz 2>/dev/null | wc -l | tr -d ' ')/10"
echo "    pages  : $(ls $PAGES/*.html 2>/dev/null | grep -vc index.html)"
echo "    results: $A2/lf_headcmp/segment_*_${ARM}_tol${TOL}.csv"
echo "    shap   : $A2/lf_headcmp/shap_position_${ARM}.csv"
