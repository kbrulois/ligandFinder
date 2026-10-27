# ARM: motifind -- "curated + motif indicators" (2026-10-26)
#
# The usual 26 hand-built channels PLUS THREE binary indicators of the residue's
# own identity at each of the 36 window positions -- is it G, is it basic (K/R),
# is it P -- for 26 + 3 = 29 channels. Nothing is dropped, so like
# curated_wAAid this arm is a strict superset of the curated model.
#
# WHAT IT TESTS. curated_wAAid (26 + a 21-letter one-hot, 47 channels) roughly
# doubled C-terminal AP (4 of 4 non-degenerate cells, mean 2.07x) and halved the
# anchors needed for 50% recall, but made anchors-for-80%-recall WORSE in 5 of 6
# C cells. Its attribution keyed on G at offset +1 and R at +2/+3 -- the
# amidation glycine then the dibasic site -- so the signal it found is a
# three-letter fact carried by 21 channels. This arm asks whether ~3 indicator
# channels recover the AP gain without the cost in the tail.
#
# Coverage is 100%: every residue has a letter, so like curated_wAAid and unlike
# the ESM-C arm there is no zero-embedding subset to qualify the result with.
# 77% of residues are all-zero across the three indicators, which is the
# ordinary case rather than missingness -- the curated `padding` channel is what
# marks "outside the sequence".

ARM_LABEL="curated+motif-ind"
PLM=${PLM:-$A2/lf_plm/motif_ind.parquet}
DROP=""                                   # keep every curated channel
N_CHAN_EXPECT=29

# how to build $PLM if it is missing (stage 0). The R step exports the sequences
# of all 5,616 precursors from the residue cache; the python step assigns the
# indicators and refuses to write unless the letters reproduce the residue
# table's own four AA property channels, all 3M rows.
mk_plm() {
  local tsv=$A2/lf_plm/precursor_seqs.tsv
  if [[ ! -s $tsv ]]; then
    eval $R -e \''p <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")$prec
      write.table(data.frame(accession=p$accession, gene=p$gene, seq=p$seq),
                  path.expand("~/AF2_analysis/lf_plm/precursor_seqs.tsv"),
                  sep="\t", quote=FALSE, row.names=FALSE)'\' || return 1
  fi
  $PY_PLM $S/motif_ind_build.py --out $PLM
}

# how to build the arm's residue table (stage 2). This arm keeps all 26 curated
# channels, so lf_winxgb.plm_residues -- which appends the plm_* columns to the
# full curated table in its own row order -- gives exactly the trained order
# (26 curated, then the 3 indicators). No custom builder needed.
mk_resid() {
  (cd $W/inst/python && $PY_TF -m lf_winxgb.plm_residues \
      --plm $PLM --out $A2/lf_pepend_residues_${ARM}.npz)
}

# No sweep has been run for this arm; this is the curated model's conventional
# cut, matching curated_wAAid so the two are compared at the same T rather than
# a measured optimum. Sweep it the way esmonly36 was swept and change this line.
THRESH_DEFAULT=0.5
