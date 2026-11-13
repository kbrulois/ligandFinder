#!/bin/zsh
# The visualization pipeline: every comparison figure for a set of peptide-end
# arms, from whatever the arm pipeline has produced so far (Kevin, 2026-11-13:
# "all these scripts should be part of a visualization pipeline", after "for all
# these plots, have it be default that it's produced").
#
#   ARMS=curated_AAswap,msa ./viz_pipeline.sh            # msa against the site's model
#   ARMS=curated,curated_AAswap,curated_AAswap_msacons_padg,msa,msa_gap \
#     PAIRS=curated_AAswap:msa,msa:msa_gap UNSEEN_SETS=_msa,_padg TAG=msa_arms ./viz_pipeline.sh
#   STAGES=shap SHAP_RUN=0 ARMS=curated,curated_AAswap ./viz_pipeline.sh   # plots from saved attributions
#
# pepend_arm_pipeline.sh stage 6 runs it for every arm against CMP_ARMS, and
# run_msa_arms.sh at the end of phase 3 for the msa arms. Every stage checks its
# own inputs and says what is missing instead of failing, so it can run at any
# point of an arm's life: right after training (trunks, gap), after the scans
# (scan), after stage 5 (metrics).
#
# STAGES (comma list of names, or "all")
#   trunks   held-out ranks from the trunk summaries alone, paired by end
#            (trunk_val_compare.R, arm_compare_paired_ranks.R)
#   gap      per-residue class probabilities at the held-out ends, for arms with
#            a gap class (arm_compare_gap_class.R)
#   metrics  normalised PR / ROC / best F1 / n50, figure and table
#            (arm_metrics_compare_plot.R, arm_compare_metrics_table.R); needs
#            stage 5's final_metrics.csv with every arm in HEADCMP
#   scan     per-end ranks (standard and held-out negatives), what the top non-end
#            anchors are, tracks along example precursors, the paired ends
#            (arm_compare_scan_tables.R, _scan_plots.R, _paired_ranks.R); needs
#            each arm's tracks rds and window groups
#   shap     attribution (shap_positions.py, TensorFlow) for any arm that has none
#            yet, then per-position by family, the AA identity heatmap and global
#            SHAP per model (arm_compare_shap_plots.R)
#   motif    amino-acid enrichment of the known ends against the mature secretome,
#            as a heatmap and as Kullback-Leibler logos (also zoomed past the
#            cleavage motif), drawn under the shap stage's AA-identity heatmap
#            (arm_compare_motif_enrichment.R)
#   index    index.html over every figure in OUT, with OUT/takeaways.tsv's line
#            (file <TAB> takeaway) under each one that has it
#
# PARAMETERS (all have working defaults)
#   ARMS         comma list; the LAST is the arm under test, the one before it
#                its comparator (dark grey), any earlier ones references
#                                                       curated_AAswap,msa
#   PAIRS        x:y[,x:y] for the paired-end figures   <second-last>:<last>
#   LABELS       ";"-separated, one per arm, for the metrics figure    the names
#   TAG          in every file name                     <last arm>
#   OUT          ~/AF2_analysis/lf_viz_<today>
#   UNSEEN_SETS  the never-seen set(s), intersected     the last arm's INPUT_SFX
#   HEADCMP      stage 5's metrics directory            ~/AF2_analysis/lf_headcmp[/unseen<UNSEEN_SFX>]
#   SHAP_DIR     where new attribution tables go        ~/AF2_analysis/lf_shap
#   SHAP_RUN     0 = never compute attributions, plot what exists   1
#   STAGES       comma list or "all"                    all
#
# "curated" and "flatten" both name the original curated trunk (no suffix).
# Logs: OUT/viz_pipeline_<TAG>.log.
set -u
S=${0:A:h}
A2=~/AF2_analysis
. $S/lf_paths.sh                     # W = the ligandFinder checkout
export LF_SRC=${LF_SRC:-$W}
R=(env LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript)
PY_PLM=~/.virtualenvs/lf-plm/bin/python
PY_TF=~/.virtualenvs/r-tensorflow/bin/python

ARMS=${ARMS:-curated_AAswap,msa}
arms=(${(s:,:)ARMS})
(( ${#arms} >= 2 )) || { echo "viz_pipeline: ARMS needs at least two arms (got '$ARMS')"; exit 1 }
LAST=${arms[-1]}
PAIRS=${PAIRS:-${arms[-2]}:$LAST}
TAG=${TAG:-$LAST}
OUT=${OUT:-$A2/lf_viz_$(date +%Y%m%d)}
STAGES=${STAGES:-all}
SHAP_DIR=${SHAP_DIR:-$A2/lf_shap}
SHAP_RUN=${SHAP_RUN:-1}
## an arm's training-window suffix, from its arm file (INPUT_SFX); "" without one
in_sfx() { [[ -r $S/arms/$1.sh ]] && sed -n 's/^INPUT_SFX=//p' $S/arms/$1.sh | tail -1 }
UNSEEN_SETS=${UNSEEN_SETS-$(in_sfx $LAST)}
U1=${UNSEEN_SETS%%,*}
HEADCMP=${HEADCMP:-$A2/lf_headcmp${U1:+/unseen$U1}}
LABELS=${LABELS:-${ARMS//,/;}}
## the metrics figure takes one label per arm; the other figures take arm=label
## pairs (and must not see the plain list, or an arm can come out unlabelled)
labs=(${(s:;:)LABELS}); ARM_LABELS=""
for i in {1..${#arms}}; do ARM_LABELS+="${ARM_LABELS:+;}${arms[$i]}=${${labs[$i]:-${arms[$i]}}## }"; done
GROUPS=""; for a in $arms; do s=$(in_sfx $a); [[ -n $s ]] && GROUPS+="${GROUPS:+,}$a:$s"; done
DATA=$OUT/${TAG}_data
mkdir -p $OUT $DATA
LOG=$OUT/viz_pipeline_$TAG.log
say() { print -r -- "[$(date '+%H:%M:%S')] viz $TAG: $*" | tee -a $LOG }
run() { [[ $STAGES == all || ",$STAGES," == *",$1,"* ]] }
## metrics_final.R and the trunk tables call the curated trunk "flatten" / "curated"
mname() { [[ $1 == curated ]] && print flatten || print $1 }
say "ARMS=$ARMS PAIRS=$PAIRS OUT=$OUT never-seen=${UNSEEN_SETS:-<original>} HEADCMP=$HEADCMP"

# ---- trunks: held-out ranks from the trunk summaries -----------------------------
if run trunks; then
  ## trunk_val_compare.R's short names for three arms of 2026-11-10
  tshort() { case $1 in curated_AAswap) print AAswap_old;; curated_AAswap_padg) print AAswap_padg;;
                       curated_AAswap_msacons_padg) print msacons_padg;; *) print $1;; esac }
  tarms=(); for a in $arms; do tarms+=$(tshort $a); done
  tpairs=(); for p in ${(s:,:)PAIRS}; do tpairs+="$(tshort ${p%%:*}):$(tshort ${p##*:})"; done
  (cd $S && ARMS=$ARMS PAIR= HEADCMP=$DATA $R trunk_val_compare.R) >> $LOG 2>&1 \
    && RANKS=$DATA/trunk_val_ranks.csv PAIRS=${(j:,:)tpairs} ARMS=${(j:,:)tarms} LABELS= OUT=$OUT NAME=paired_ranks_${TAG}_trunks \
         RANK_DESC="from the trunk summaries: rank among the ~40k candidate windows (superseded by the scan stage)" \
         $R $S/arm_compare_paired_ranks.R >> $LOG 2>&1 \
    && say "trunks: paired_ranks_${TAG}_trunks{,_by_motif}.png" || say "trunks FAILED, see $LOG"
fi

# ---- gap: class probabilities at the held-out ends ---------------------------------
if run gap; then
  for a in $arms; do
    grep -q '"gap"' $A2/lf_pepend_run_C_cv5_f1_nomask_flatten_${a}_isolated/in/config.json 2>/dev/null || continue
    ## drawn next to its pair's first arm if it has one, else the arm before it
    cmp=""; for p in ${(s:,:)PAIRS}; do [[ ${p##*:} == $a ]] && cmp=${p%%:*}; done
    [[ -z $cmp ]] && cmp=${arms[$(( ${arms[(i)$a]} - 1 ))]}
    RETICULATE_PYTHON=$PY_PLM ARMS=$a,$cmp OUT=$OUT NAME=gap_class_$a $R $S/arm_compare_gap_class.R >> $LOG 2>&1 \
      && say "gap: gap_class_$a.png" || say "gap FAILED for $a, see $LOG"
  done
fi

# ---- metrics: the headline figure and table ----------------------------------------
if run metrics; then
  marms=(); for a in $arms; do marms+=$(mname $a); done
  M=$HEADCMP/final_metrics.csv
  miss=(); for a in $marms; do [[ -s $M ]] && cut -d, -f3 $M | grep -qx $a || miss+=$a; done
  if (( ${#miss} )); then
    say "metrics: skipped, $M lacks ${(j:, :)miss} (stage 5 of pepend_arm_pipeline.sh writes it)"
  else
    HEADCMP=$HEADCMP ARMS=${(j:,:)marms} LABELS=$LABELS OUT=$OUT NAME=cv_metrics_$TAG \
      TITLE="${TITLE:-${LABELS//;/ / }: full-proteome scan metrics}" $R $S/arm_metrics_compare_plot.R >> $LOG 2>&1 \
      && say "metrics: cv_metrics_$TAG.png" || say "metrics figure FAILED, see $LOG"
    HEADCMP=$HEADCMP ARMS=${(j:,:)marms} LABELS=$LABELS OUT=$OUT NAME=cv_metrics_${TAG}_table \
      $R $S/arm_compare_metrics_table.R >> $LOG 2>&1 \
      && say "metrics: cv_metrics_${TAG}_table.png" || say "metrics table FAILED, see $LOG"
  fi
fi

# ---- scan: per-end ranks, false positives, tracks, paired ends -----------------------
if run scan; then
  miss=()
  for a in $arms; do
    t=$A2/lf_pepend_cv5_tracks_C_flatten${${a:#(curated|flatten)}:+_$a}.rds
    [[ -s $t ]] || miss+=$a
  done
  if (( ${#miss} )); then
    say "scan: skipped, no tracks for ${(j:, :)miss} (stages 3-4 of pepend_arm_pipeline.sh)"
  else
    RETICULATE_PYTHON=$PY_PLM ARMS=$ARMS UNSEEN_SETS=$UNSEEN_SETS GROUPS=$GROUPS OUT=$DATA \
      $R $S/arm_compare_scan_tables.R >> $LOG 2>&1 || say "scan tables FAILED, see $LOG"
    if [[ -s $DATA/end_ranks.csv ]]; then
      tr=(${(u)${(s:,:)${PAIRS//:/,}}}); (( ${#tr} > 3 )) && tr=(${tr[-3,-1]})
      DATA=$DATA OUT=$OUT TAG=$TAG ARMS=$ARMS TRACK_ARMS=${(j:,:)tr} LABELS=$ARM_LABELS $R $S/arm_compare_scan_plots.R >> $LOG 2>&1 \
        || say "scan plots FAILED, see $LOG"
      RANKS=$DATA/end_ranks.csv PAIRS=$PAIRS ARMS=$ARMS LABELS=$ARM_LABELS OUT=$OUT NAME=paired_ranks_$TAG \
        $R $S/arm_compare_paired_ranks.R >> $LOG 2>&1 || say "paired ranks FAILED, see $LOG"
      RANKS=$DATA/end_ranks.csv RANK_COL=rank_ho PAIRS=$PAIRS ARMS=$ARMS LABELS=$ARM_LABELS OUT=$OUT NAME=paired_ranks_${TAG}_heldout \
        RANK_DESC="rank against held-out negatives only (anchors no training window of any compared arm reaches)" \
        $R $S/arm_compare_paired_ranks.R >> $LOG 2>&1 || say "paired ranks (held-out) FAILED, see $LOG"
      say "scan: recall_rank_$TAG, fp_motif_$TAG, tracks_examples_$TAG, paired_ranks_$TAG{,_heldout}"
    fi
  fi
fi

# ---- shap: attribution where missing, then the figures --------------------------------
if run shap; then
  mkdir -p $SHAP_DIR
  for a in $arms; do
    sfx=${${a:#(curated|flatten)}:+_$a}
    f=""; for d in $SHAP_DIR $A2/lf_headcmp; do [[ -s $d/shap_channel$sfx.csv ]] && { f=$d/shap_completeness$sfx.csv; break } done
    ## all 10 trunks attributed? (one completeness row per terminus and fold)
    [[ -n $f && $(( $(wc -l < $f 2>/dev/null || echo 0) - 1 )) -ge 10 ]] && continue
    if [[ $SHAP_RUN == 0 ]]; then say "shap: no complete attribution for $a (SHAP_RUN=0)"; continue; fi
    n=$(ls $A2/lf_pepend_run_[NC]_cv5_f*_nomask_flatten$sfx.rds 2>/dev/null | wc -l | tr -d ' ')
    say "shap: attributing $a ($n of 10 trunks; TensorFlow, ~10-20 min)"
    $PY_TF $S/shap_positions.py --arm "${sfx#_}" --out-dir $SHAP_DIR >> $LOG 2>&1 || say "shap_positions FAILED for $a, see $LOG"
  done
  SHAP_DIRS=$SHAP_DIR:$A2/lf_headcmp ARMS=$ARMS OUT=$OUT TAG=$TAG $R $S/arm_compare_shap_plots.R >> $LOG 2>&1 \
    && say "shap: shap_global_$TAG, shap_positions_by_family_$TAG{,_signed}, shap_aa_heatmap_$TAG" \
    || say "shap plots FAILED, see $LOG"
fi

# ---- motif: what is enriched at real ends, under what the models use ------------------
if run motif; then
  SEQS=$DATA/prec_seqs.rds OUT=$OUT TAG=$TAG $R $S/arm_compare_motif_enrichment.R >> $LOG 2>&1 \
    && say "motif: motif_logo_$TAG{,_zoom}, motif_enrichment_$TAG, motif_logo_known_ends{,_zoom}" \
    || say "motif FAILED, see $LOG"
fi

# ---- index: one page over every figure in OUT ----------------------------------------
if run index; then
  T=$OUT/takeaways.tsv
  {
    print '<!doctype html><meta charset="utf-8"><title>'${OUT:t}'</title>'
    print '<style>body{font:14px/1.45 -apple-system,Helvetica,sans-serif;max-width:1500px;margin:24px auto;padding:0 16px;color:#222}'
    print 'figure{margin:0 0 36px}img{max-width:100%;border:1px solid #ddd}figcaption{margin:6px 0}code{color:#555}</style>'
    print "<h1>${OUT:t}</h1><p>Figures from <code>viz_pipeline.sh</code>; newest first.</p>"
    for p in $OUT/*.png(Nom); do
      b=${p:t}; tk=""
      [[ -r $T ]] && tk=$(awk -F'\t' -v f=$b '$1 == f {print $2}' $T | sed 's/&/\&amp;/g; s/</\&lt;/g; s/>/\&gt;/g')
      print "<figure><a href=\"$b\"><img src=\"$b\" loading=\"lazy\"></a><figcaption><code>$b</code>${tk:+ &mdash; $tk}</figcaption></figure>"
    done
  } > $OUT/index.html
  say "index: $OUT/index.html"
fi
say "done"
