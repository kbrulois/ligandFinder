#!/bin/zsh
# Deferred work for the two _padg arms (Kevin, 2026-11-10: "train both arms now,
# scan after the API finishes").
#
#   1. wait for the ColabFold API fetch to finish
#   2. rescore every MSA and rebuild the conservation channel, with the raw-JSD
#      min/max PINNED to the build the trunks were trained on -- recomputing it
#      over a larger MSA set would put the scan's channel on a different scale
#      from the one the model learned
#   3. verify that only the previously-uncovered precursors changed
#   4. stages 2b/2/3 (window groups, residue table, full-proteome scans) for both
#      arms. NOT pages: those are a separate, much larger job.
set -u
S=${0:A:h}
A2=~/AF2_analysis
. ${0:A:h}/lf_paths.sh          # W = the ligandFinder checkout (LF_SRC, or the one this sits in)
PY_PLM=~/.virtualenvs/lf-plm/bin/python

## Wait for the MSAs to BE THERE, not merely for the fetcher to be gone: the
## first run of this chain started rescoring at 380/388 because the client had
## died on a network timeout and "process absent" read as "finished".
## MIN_A3M lets a deliberate short run through (e.g. a precursor the API refuses).
MIN_A3M=${MIN_A3M:-388}
while true; do
  n=$(ls $A2/lf_msa/api_colabfold/*.a3m.gz 2>/dev/null | wc -l | tr -d ' ')
  [[ $n -ge $MIN_A3M ]] && break
  if ! pgrep -qf "fetch_colabfold_api.py"; then
    echo "API fetcher is not running and only $n/$MIN_A3M alignments are present."
    echo "Re-run it (inst/scripts/msa_cons/fetch_colabfold_api.py) or set MIN_A3M=$n to accept this."
    exit 1
  fi
  sleep 300
done
echo "API alignments complete: $(ls $A2/lf_msa/api_colabfold/*.a3m.gz | wc -l | tr -d ' ')/388, $(date '+%F %H:%M:%S')"

cp $A2/lf_plm/aa_onehot_msacons.parquet $A2/lf_plm/aa_onehot_msacons_trainbuild.parquet
$PY_PLM $W/inst/scripts/msa_cons/msa_conservation.py --workers 4 --allow-missing \
  > $A2/lf_msa/score_run_full.log 2>&1 || { echo "RESCORE FAILED"; exit 1; }
tail -3 $A2/lf_msa/score_run_full.log
$PY_PLM $W/inst/scripts/msa_cons/build_msacons_parquet.py --allow-missing \
  --scale-from $A2/lf_plm/aa_onehot_msacons.parquet.scale.json || { echo "BUILD FAILED"; exit 1; }

## only the precursors that had no MSA before may have changed
$PY_PLM $S/check_msacons_rebuild.py || { echo "CHANNEL CHANGED WHERE IT SHOULD NOT -- not scanning"; exit 1; }

for ARM in curated_AAswap_padg curated_AAswap_msacons_padg; do
  echo "=== $ARM: groups, residue table, scans  $(date '+%F %H:%M:%S')"
  (cd $S && ARM=$ARM STAGES=2b,2,3 PROCS=3 ./pepend_arm_pipeline.sh) \
    > $A2/pipe_${ARM}_scan.log 2>&1 || echo "  $ARM STAGE FAILED, see pipe_${ARM}_scan.log"
  ls $A2/lf_winxgb_fullscan_[NC]_cv5_f*_nomask_flatten_${ARM}.npz 2>/dev/null | wc -l | tr -d ' '
done
echo "DONE $(date '+%F %H:%M:%S')"
