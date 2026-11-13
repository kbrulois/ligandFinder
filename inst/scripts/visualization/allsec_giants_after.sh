#!/bin/zsh
# Wait for the main all-secretome render to finish, then render the 19 giants.
#
# Why this exists: page rendering is MEMORY-bound on this 32 GB machine, not
# CPU-bound (see the note at the top of allsec_giants.sh), so the giants cannot
# share the machine with the main run's five workers -- two of them beside it
# took swap to 97% and halved its throughput. They have to go afterwards, and
# the main run has hours left, so this waits instead of someone watching for it.
#
#   nohup ./allsec_giants_after.sh > giants_after.log 2>&1 &
#
# Safe to start, stop and restart: it only ever waits and then calls
# allsec_giants.sh, which skips any gene whose .html already exists.
set -u
S=${0:A:h}
A2=~/AF2_analysis
POLL=${POLL:-120}
GIANT_PROCS=${GIANT_PROCS:-1}

say() { echo "[$(date '+%F %H:%M:%S')] $*" }

# "finished" = no final_pages.R worker left AND the driver has written its
# closing line. Both, because there is a gap between the last worker exiting
# and the driver finishing its index.
main_running() {
  ps aux | grep -q "[-]-file=$A2/lf_nomask_scripts_20261021/final_pages.R"
}

say "waiting for the main render to finish (polling every ${POLL}s)"
say "  at start: $(ls $A2/pepend_pages_allsec/*.html 2>/dev/null | grep -vc index.html) pages, main running: $(main_running && echo yes || echo no)"

while main_running; do sleep $POLL; done

# let the driver finish its own index pass before adding to the directory
say "no render workers left; settling for 3 min in case the driver is indexing"
sleep 180
while main_running; do sleep $POLL; done     # a late chunk would restart one

say "main render done: $(ls $A2/pepend_pages_allsec/*.html | grep -vc index.html) pages"
say "free swap before the giants: $(sysctl -n vm.swapusage | sed -E 's/.*free = ([0-9.]+M).*/\1/')"
say "starting allsec_giants.sh PROCS=$GIANT_PROCS"

PROCS=$GIANT_PROCS $S/allsec_giants.sh
say "allsec_giants.sh exited rc=$?"
