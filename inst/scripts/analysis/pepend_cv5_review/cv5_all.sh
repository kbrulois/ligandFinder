#!/bin/zsh
set -u
S2=/private/tmp/claude-501/-Users-kbrulois-R-projects-ligandFinder--claude-worktrees-modest-satoshi-87b2c8/92c0469f-235a-4310-b3a9-55060cf86b72/scratchpad
# wait for the lf_resid arms to clear so the folds get the cores
while pgrep -f "lf_resid (train|member)" >/dev/null 2>&1; do sleep 60; done
echo "=== lf_resid clear, starting CV $(date +%H:%M:%S)"
# two folds at a time: each training run uses ~2 cores, 8 performance cores here
for pair in "1 2" "3 4" "5"; do
  for k in ${=pair}; do $S2/cv5_fold.sh $k >> $S2/cv5_driver.log 2>&1 & done
  wait
  echo "=== pair [$pair] complete $(date +%H:%M:%S)"
done
echo "=== ALL FOLDS DONE $(date +%H:%M:%S)"
