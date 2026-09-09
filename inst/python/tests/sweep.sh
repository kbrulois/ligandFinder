#!/bin/zsh
# Run compare_r_python.R as ONE PROCESS PER (arm, seed), then combine the CSVs.
#
#   inst/python/tests/sweep.sh                                  # defaults below
#   inst/python/tests/sweep.sh "flat@39-19 unet_bneck" "1 2"
#   ARMS="flat unet" SEEDS="1 2 3 42" TERMS=C inst/python/tests/sweep.sh
#
# Why one process each: a long run followed by another model in the same
# process can die inside a retraced tf.function summing an empty
# regularization-loss list. Isolation avoids it, and one crash costs one cell.
set -u
ARMS=${1:-${ARMS:-"flat unet"}}
SEEDS=${2:-${SEEDS:-"1 2"}}
TERMS=${TERMS:-C}
EPOCHS=${EPOCHS:-2000}
ROOT=${ROOT:-$(cd "$(dirname "$0")/../../.." && pwd)}
OUTDIR=${OUTDIR:-~/AF2_analysis/lf_dcnn_sweep_$(date +%Y%m%d_%H%M%S)}
RSCRIPT=${RSCRIPT:-/usr/local/bin/Rscript}   # NOT the pixi Rscript on PATH

mkdir -p "${OUTDIR}"
cd "$ROOT"
echo "arms=[$ARMS] seeds=[$SEEDS] terms=$TERMS epochs=$EPOCHS -> $OUTDIR"

for arm in ${=ARMS}; do
  for seed in ${=SEEDS}; do
    tag=$(echo "$arm" | tr '@-' '__')_s${seed}
    $RSCRIPT inst/python/tests/compare_r_python.R \
      --terms "$TERMS" --seeds "$seed" --arms "$arm" --epochs "$EPOCHS" \
      --out "${OUTDIR}/${tag}.csv" > "${OUTDIR}/${tag}.log" 2>&1
    echo "  $arm seed=$seed rc=$?"
  done
done

$RSCRIPT -e '
  d <- commandArgs(TRUE)[1]
  fs <- list.files(d, "\\.csv$", full.names = TRUE)
  fs <- fs[basename(fs) != "combined.csv"]
  if (!length(fs)) { cat("no results\n"); quit(status = 1) }
  x <- do.call(rbind, lapply(fs, read.csv))
  write.csv(x, file.path(d, "combined.csv"), row.names = FALSE)
  s <- aggregate(cbind(acc_val, acc_train, roc_auc, pr_auc) ~ arm + params, x, mean)
  sd_ <- aggregate(acc_val ~ arm, x, sd); names(sd_)[2] <- "acc_val_sd"
  n_  <- aggregate(acc_val ~ arm, x, length); names(n_)[2] <- "n"
  s <- merge(merge(s, sd_, by = "arm"), n_, by = "arm")
  s <- s[order(s$params), c("arm","params","n","acc_val","acc_val_sd","acc_train","roc_auc","pr_auc")]
  print(format(s, digits = 3), row.names = FALSE)
  cat("\ncombined ->", file.path(d, "combined.csv"), "\n")
' "${OUTDIR}"
