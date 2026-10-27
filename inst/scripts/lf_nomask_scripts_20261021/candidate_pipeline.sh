#!/bin/zsh
# Candidate peptides, end to end, with a provenance record (Kevin, 2026-10-23:
# "keep a record of how this table is constructed so the analysis can be redone
# with different threshold and tolerance").
#
# Every stage reads the SAME five parameters from the environment, and the
# output file name carries three of them, so two parameter sets never collide:
#
#   HEAD     which window head's profile to segment   (xgb_pm | cnn)
#   THRESH   score at or above which an anchor is in an end segment
#   TOL      half-width of the fan-out around each es centre
#   LEN_MIN  shortest candidate peptide to keep
#   LEN_MAX  longest
#   ARM      which arm's tracks ("" = curated, "_esmonly36", "_curated_wAAid");
#            suffixes every file name. Pass the arm's own THRESH with it.
#   ARMS_MAX segment the per-anchor MAX of several arms instead ("curated,
#            esmonly36,curated_AAswap"); each es and candidate end then says
#            how many arms support it. Needs its own ARM suffix.
#   KNOWN_LEN_MIN / KNOWN_LEN_MAX  the length range of the known peptides
#            MMseqs searches against (8 / 60)
#
#   ./candidate_pipeline.sh                              # the defaults below
#   THRESH=0.3 TOL=2 ./candidate_pipeline.sh             # redo at another cut
#   HEAD=cnn THRESH=0.5 TOL=4 ./candidate_pipeline.sh    # the linear head
#   ARM=_esmonly36 THRESH=0.65 ./candidate_pipeline.sh   # the ESM-C-only arm
#   ARMS_MAX=curated,esmonly36,curated_AAswap ARM=_max3 ./candidate_pipeline.sh
#
# Stages (each is a script you can run on its own with the same env):
#   1 candidate_peptides.R      pair N and C end segments -> candidates
#   2 candidate_mmseqs_prep.R   build known.fasta + cand.fasta
#     mmseqs easy-search        search candidates against known peptides
#   3 candidate_annotate.R      merge the best hit + derived columns
#   4 candidate_multipep.sh     MultiPep bioactivity classes (cached by sequence)
#   5 candidate_logistic.R      the 5-variable ranker -> lr_prob
#
# WITH_PDEEPPP=1 also runs the PDeepPP neuro head. Off by default: measured on
# our own ground truth it ranks real ligands BACKWARDS (AUC 0.21) because it
# behaves like a length cut near 29 aa. Kept only so the finding is checkable.
#
# Writes <output>.manifest.txt beside the CSV: parameters, inputs with their
# size and mtime, tool versions, the exact mmseqs command, and row counts.
set -u
S=${0:A:h}
A2=~/AF2_analysis
R="LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript"

export HEAD=${HEAD:-xgb_pm}
export THRESH=${THRESH:-0.5}
export TOL=${TOL:-4}
export LEN_MIN=${LEN_MIN:-8}
export LEN_MAX=${LEN_MAX:-60}
export ARM=${ARM:-}
export ARMS_MAX=${ARMS_MAX:-}
export KNOWN_LEN_MIN=${KNOWN_LEN_MIN:-8}
export KNOWN_LEN_MAX=${KNOWN_LEN_MAX:-60}
[[ -n $ARMS_MAX && -z $ARM ]] && { echo "ARMS_MAX needs its own ARM suffix (e.g. ARM=_max3)"; exit 1 }

OUT=$A2/lf_headcmp/candidate_peptides_${HEAD}_T${THRESH}_tol${TOL}${ARM}.csv
MAN=$A2/lf_headcmp/candidate_peptides_${HEAD}_T${THRESH}_tol${TOL}${ARM}.manifest.txt
O=$A2/mmseq_cand

# the exact search settings, and why: candidate_mmseqs_args.sh
. $S/candidate_mmseqs_args.sh

echo "=== CANDIDATE PIPELINE  HEAD=$HEAD THRESH=$THRESH TOL=$TOL LEN=[$LEN_MIN,$LEN_MAX] ARM=${ARM:-(curated)} ==="
echo "--- stage 1: assemble candidates from paired end segments"
eval $R $S/candidate_peptides.R || { echo "STAGE 1 FAILED"; exit 1; }

echo "--- stage 2: MMseqs2 against known peptides"
eval $R $S/candidate_mmseqs_prep.R || { echo "STAGE 2a FAILED"; exit 1; }
rm -rf $O/tmp $O/res.m8 $O/res_all.m8
mmseqs easy-search $O/cand.fasta $O/known.fasta $O/res_all.m8 $O/tmp "${MMSEQS_ARGS[@]}" > $O/mmseqs.log 2>&1 \
  || { echo "mmseqs FAILED, see $O/mmseqs.log"; exit 1; }
echo "    $(wc -l < $O/res_all.m8 | tr -d ' ') hits, $(cut -f1 $O/res_all.m8 | sort -u | wc -l | tr -d ' ') queries with one"
# Stage 3 uses only each candidate's best hit overall and its best within each
# target type, so keep exactly those -- every row tied on top bits for its
# (candidate, type), ties left for stage 3 to break as before. ~80 hits per
# candidate become ~3, which keeps tens of millions of rows out of R's memory.
awk -F'\t' 'FNR == 1 { p++ }
  p == 1 { ty[$1] = $2; next }
  p == 2 { k = $1 SUBSEP ty[$2]; b = $5 + 0; if (!(k in mx) || b > mx[k]) mx[k] = b; next }
  { k = $1 SUBSEP ty[$2]; if ($5 + 0 == mx[k]) print }' \
  $O/known_types.tsv $O/res_all.m8 $O/res_all.m8 > $O/res.m8 \
  || { echo "reducing the hits FAILED"; exit 1; }
echo "    kept $(wc -l < $O/res.m8 | tr -d ' ') best-per-type hits for stage 3"
rm -f $O/res_all.m8

echo "--- stage 3: merge best hit + derived columns"
eval $R $S/candidate_annotate.R || { echo "STAGE 3 FAILED"; exit 1; }

echo "--- stage 4: MultiPep bioactivity classes"
$S/candidate_multipep.sh || { echo "STAGE 4 FAILED"; exit 1; }

if [[ ${WITH_PDEEPPP:-0} == 1 ]]; then
  echo "--- stage 4b: PDeepPP neuro head (off by default; see the header)"
  TASK=neuro ~/.virtualenvs/pdeeppp/bin/python -I $A2/pdeeppp_run/score_candidates.py \
    && eval TASK=neuro $R $S/candidate_pdeeppp_merge.R || echo "  stage 4b failed, continuing"
fi

echo "--- stage 5: the logistic ranker"
eval $R $S/candidate_logistic.R || { echo "STAGE 5 FAILED"; exit 1; }

# ---- provenance (its own script, so it can be refreshed alone) --------------
HEAD=$HEAD THRESH=$THRESH TOL=$TOL LEN_MIN=$LEN_MIN LEN_MAX=$LEN_MAX ARM=$ARM $S/candidate_manifest.sh

echo
echo "=== DONE ==="
echo "  table    : $OUT"
echo "  segments : $A2/lf_headcmp/end_segments_${HEAD}_T${THRESH}${ARM}.csv"
echo "  manifest : $MAN"
