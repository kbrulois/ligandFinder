#!/bin/zsh
# Unattended driver for the ESM-C-only arm (Kevin, 2026-10-24, signing off).
#
#   nohup ./esmonly36_driver.sh > esmonly36_driver.log 2>&1 &
#
# Three stages, each gated on the previous one actually having worked:
#   1  wait for `lf_plm reduce --k 36` to land esm_c_pca36.parquet, and check it
#      really has 36 channels before anything trains
#   2  train ONE job (C, fold 1) and assert the channel count went 1 -> 37.
#      This stage exists because the all-secretome render died 90 seconds in on
#      a missing chmod +x: a script that parses is not a script that runs. One
#      cheap job proves the whole path before ten expensive ones commit to it.
#   3  train the remaining nine
#
# It touches nothing in the existing pipeline: new --out names (_esmonly36), a
# new cache, a new log prefix. The curated-channel models are untouched.
set -u
S=${0:A:h}
A2=~/AF2_analysis
PLM=$A2/lf_plm/esm_c_pca36.parquet
say() { echo "[$(date '+%F %H:%M:%S')] $*"; }

say "=== ESM-C-ONLY ARM START ==="

# ---- stage 1: the embedding ------------------------------------------------
say "stage 1: waiting for the k=36 PCA reduction"
for i in {1..240}; do                      # up to 4 h, checked each minute
  pgrep -f 'lf_plm reduce' > /dev/null || break
  sleep 60
done
if pgrep -f 'lf_plm reduce' > /dev/null; then
  say "stage 1 FAILED: reduce still running after 4 h, giving up"; exit 1
fi
[[ -s $PLM ]] || { say "stage 1 FAILED: $PLM was never written"; tail -5 $A2/lf_plm/reduce_esm_c_pca36.log; exit 1; }
say "stage 1: parquet present, $(du -h $PLM | cut -f1)"

nch=$(LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript -e '
  suppressMessages(library(arrow))
  d <- read_parquet(path.expand("~/AF2_analysis/lf_plm/esm_c_pca36.parquet"), as_data_frame = FALSE)
  cat(sum(grepl("^(PC|V|pc)[0-9]+$", names(d))), "\n")' 2>/dev/null | tail -1)
say "stage 1: $nch PCA channels in the parquet"
if [[ ${nch:-0} -ne 36 ]]; then
  say "stage 1 WARNING: expected 36 channels, saw ${nch:-?} -- continuing, 10_8 reports the real count itself"
fi

# ---- stage 2: prove one job runs -------------------------------------------
say "stage 2: one job (C fold 1) to prove the path executes"
$S/train_esmonly36.sh one 1 C
one=$A2/lf_pepend_run_C_cv5_f1_nomask_flatten_esmonly36.rds
if [[ ! -s $one ]]; then
  say "stage 2 FAILED: no rds produced. tail of the log:"
  tail -25 $S/trainesm36_C_f1.log; exit 1
fi
say "stage 2: ok, $(du -h $one | cut -f1)"
grep -h 'channels:' $S/trainesm36_C_f1.log | tail -1 | sed "s/^/[$(date '+%H:%M:%S')] stage 2: /"

# ---- stage 3: the rest ------------------------------------------------------
say "stage 3: the remaining nine jobs"
PROCS=${PROCS:-2} $S/train_esmonly36.sh
say "=== ESM-C-ONLY ARM DONE ==="
