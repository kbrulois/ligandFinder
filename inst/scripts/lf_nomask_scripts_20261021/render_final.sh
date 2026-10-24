#!/bin/zsh
# Render the final (flatten, 5-fold) protein pages: 104 genes, four processes.
# The five largest precursors go at popup max 600 (as the cv5 gallery did: at
# 2000 they produce 80-120 MB pages); everything else at the default 2000.
set -u
S=${0:A:h}
BIG=(ANO8 C3 C5 EGF COPA)
REST=($(grep -v -x -F -e ANO8 -e C3 -e C5 -e EGF -e COPA $S/page_genes.txt))
n=${#REST}; per=$(( (n + 2) / 3 ))
cd ~
for i in 1 2 3; do
  chunk=(${REST[$(( (i-1)*per + 1 )),$(( i*per ))]})
  ( MODE=cv5 PROV_N=2 LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $S/final_pages.R ${chunk[@]} \
      > $S/render_final_$i.log 2>&1; echo "chunk $i (${#chunk} genes) rc=$?" ) &
done
( MODE=cv5 PROV_N=1 LF_POPUP_MAX=600 LC_ALL=en_US.UTF-8 /usr/local/bin/Rscript $S/final_pages.R ${BIG[@]} \
    > $S/render_final_big.log 2>&1; echo "big (${#BIG} genes) rc=$?" ) &
wait
echo "RENDER DONE $(date +%H:%M:%S)"
