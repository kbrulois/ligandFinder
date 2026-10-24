#!/bin/zsh
# Writes the provenance record for a candidate table. Called by
# candidate_pipeline.sh, but safe to run on its own after re-running only a
# later stage (e.g. candidate_annotate.R added columns):
#   HEAD=xgb_pm THRESH=0.5 TOL=4 ./candidate_manifest.sh
set -u
S=${0:A:h}
. $S/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
A2=~/AF2_analysis
HEAD=${HEAD:-xgb_pm}; THRESH=${THRESH:-0.5}; TOL=${TOL:-4}
LEN_MIN=${LEN_MIN:-8}; LEN_MAX=${LEN_MAX:-60}
OUT=$A2/lf_headcmp/candidate_peptides_${HEAD}_T${THRESH}_tol${TOL}.csv
MAN=$A2/lf_headcmp/candidate_peptides_${HEAD}_T${THRESH}_tol${TOL}.manifest.txt
MMSEQS_ARGS=(-s 7.5 --min-seq-id 0 -e 10000 -k 5 --comp-bias-corr 0 --max-seqs 300
             --alignment-mode 3 --format-output "query,target,pident,alnlen,bits,evalue,qcov,tcov")
stat_line() { [[ -e $1 ]] && echo "    $(basename $1)  $(stat -f '%z bytes, modified %Sm' -t '%F %T' $1)" || echo "    $(basename $1)  MISSING" }
{
  echo "candidate peptide table -- how it was built"
  echo "built            : $(date '+%F %T %Z') by $(whoami) on $(hostname -s)"
  echo
  echo "PARAMETERS"
  echo "    HEAD    = $HEAD      (window head whose all-50 profile is segmented)"
  echo "    THRESH  = $THRESH     (anchor score at or above which it joins an end segment)"
  echo "    TOL     = $TOL        (fan-out: every anchor within +/-TOL of an es centre is offered)"
  echo "    LEN_MIN = $LEN_MIN, LEN_MAX = $LEN_MAX"
  echo
  echo "RULE"
  echo "    An end segment (es) is a run of consecutive anchors scoring >= THRESH, centred"
  echo "    at the run's maximum. An N anchor is a peptide's FIRST residue and a C anchor"
  echo "    its LAST, so one (N es, C es) pair on a precursor spans n_pos..c_pos. TOL is a"
  echo "    fan-out, not a filter: each es offers every anchor within +/-TOL of its centre,"
  echo "    so a pair yields up to (2*TOL+1)^2 overlapping candidates, grouped by group_id."
  echo "    Scores are <HEAD>__all50 taken from each terminus' OWN profile."
  echo
  echo "INPUTS"
  for f in $A2/lf_pepend_cv5_tracks_N_flatten.rds $A2/lf_pepend_cv5_tracks_C_flatten.rds \
           $A2/lf_pepend_residue_cache.rds $A2/pepend_locus/uniprot_peptides.csv \
           $W/inst/extdata/ligand_list.rds; do
    stat_line $f
  done
  echo
  echo "STAGES"
  echo "    1  candidate_peptides.R"
  echo "    2  candidate_mmseqs_prep.R  +  mmseqs easy-search"
  echo "    3  candidate_annotate.R"
  echo
  echo "MMSEQS COMMAND"
  echo "    mmseqs easy-search cand.fasta known.fasta res.m8 tmp \\"
  echo "        ${MMSEQS_ARGS[@]}"
  echo "    (defaults find nothing at these lengths; the prefilter is tuned for proteins)"
  echo "    best hit per candidate chosen by BIT SCORE, not percent identity"
  echo
  echo "VERSIONS"
  echo "    $(mmseqs version 2>/dev/null | head -1 | sed 's/^/mmseqs /')"
  echo "    $(/usr/local/bin/Rscript -e 'cat(R.version.string, "| data.table", as.character(packageVersion("data.table")))' 2>/dev/null)"
  echo
  echo "OUTPUT"
  echo "    $(basename $OUT)"
  echo "    rows     : $(( $(wc -l < $OUT) - 1 ))"
  echo "    groups   : $(tail -n +2 $OUT | cut -d, -f2 | sort -u | wc -l | tr -d ' ')  (one per N es x C es pair)"
  echo "    columns  : $(head -1 $OUT | tr ',' '\n' | wc -l | tr -d ' ')"
  head -1 $OUT | tr ',' '\n' | sed 's/^/        /'
} > $MAN


echo "manifest -> $MAN"
