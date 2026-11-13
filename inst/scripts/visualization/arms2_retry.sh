#!/bin/zsh
# Re-runs arms2_allsec.sh until every page exists. A chunk can die on an
# intermittent TensorFlow segfault (first seen 2026-11-03, c_009); its finished
# pages are kept and its other genes simply stay to do, so a later pass picks
# them up. Stage 4 and 5 run again, finding little to do.
S=~/AF2_analysis/lf_src_arms2_20261103/inst/scripts/lf_nomask_scripts_20261021
O=~/AF2_analysis/pepend_pages_allsec_arms2
L=~/AF2_analysis/arms2_allsec.log
while kill -0 21847 2>/dev/null || kill -0 21849 2>/dev/null; do sleep 60; done
for pass in 1 2 3 4 5; do
  (( $(df -g ~ | tail -1 | awk '{print $4}') < 60 )) && { echo "RETRY STOPPED: under 60 GB free" >> $L; exit 0; }
  missing=0
  for g in $(cat $O/genes.txt $O/giant_genes.txt); do [[ -s $O/$g.html ]] || (( missing++ )); done
  (( missing )) || { echo "RETRY: every page exists" >> $L; break; }
  echo "[$(date '+%F %H:%M:%S')] RETRY pass $pass: $missing pages to do" >> $L
  ( cd $S && caffeinate -is ./arms2_allsec.sh ) >> $L 2>&1
done
echo "RETRY DONE $(date '+%F %H:%M:%S')" >> $L
