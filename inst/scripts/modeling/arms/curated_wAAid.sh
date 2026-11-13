# ARM: curated_wAAid -- "curated-wAAid" (Kevin, 2026-10-25)
#
# The usual 26 hand-built channels PLUS a one-hot of the residue's own identity
# at each of the 36 window positions: 26 + 21 = 47 channels. Nothing is dropped,
# so this arm is a strict superset of the curated model -- the comparison is
# "does the letter itself add anything over the four AA properties".
#
# The one-hot rides the `--plm` mechanism, which is why no feature builder, no
# window builder and no scan code needed an edit. Its channels are named
# plm_AA_A .. plm_AA_Y (21 letters, U included), sorted by letter.
#
# Coverage is 100%: every residue has a letter, so unlike the ESM-C arm there is
# no zero-embedding subset to qualify the result with.

ARM_LABEL="curated-wAAid"
PLM=${PLM:-$A2/lf_plm/aa_onehot.parquet}
DROP=""                                   # keep every curated channel
N_CHAN_EXPECT=47

# how to build $PLM if it is missing (stage 0). The R step exports the sequences
# of all 5,616 precursors from the residue cache; the python step one-hots them
# and refuses to write unless the letters reproduce the residue table's own four
# AA property channels, all 3M rows.
mk_plm() {
  local tsv=$A2/lf_plm/precursor_seqs.tsv
  if [[ ! -s $tsv ]]; then
    eval $R -e \''p <- readRDS("~/AF2_analysis/lf_pepend_residue_cache.rds")$prec
      write.table(data.frame(accession=p$accession, gene=p$gene, seq=p$seq),
                  path.expand("~/AF2_analysis/lf_plm/precursor_seqs.tsv"),
                  sep="\t", quote=FALSE, row.names=FALSE)'\' || return 1
  fi
  $PY_PLM $S/aa_onehot_build.py --out $PLM
}

# how to build the arm's residue table (stage 2). This arm keeps all 26 curated
# channels, so lf_winxgb.plm_residues -- which appends the plm_* columns to the
# full curated table in its own row order -- gives exactly the trained order
# (26 curated, then the 21 one-hots). No custom builder needed.
mk_resid() {
  (cd $W/inst/python && $PY_TF -m lf_winxgb.plm_residues \
      --plm $PLM --out $A2/lf_pepend_residues_${ARM}.npz)
}

# No sweep has been run for this arm yet, so this is the curated model's
# conventional cut rather than a measured optimum. Sweep it the way esmonly36
# was swept and change this line.
THRESH_DEFAULT=0.5
