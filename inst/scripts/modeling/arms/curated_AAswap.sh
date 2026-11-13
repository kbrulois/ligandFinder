# ARM: curated_AAswap -- "curated, AA properties -> one-hot" (Kevin, 2026-10-26)
#
# The curated channels with the FOUR AA PROPERTY channels REMOVED and the
# 21-letter one-hot of residue identity put in their place:
#
#   26 curated - 4 properties = 22, + 21 one-hot = 43 channels
#
# WHAT IT TESTS, and why it is a sharper question than curated_wAAid's. That arm
# ADDED identity beside the properties and asked "does the letter add anything".
# This one asks "does the letter REPLACE them" -- and the attribution says it
# might do better than replace. On the curated model the four properties sum to
# a NEGATIVE contribution at offset +1 (-0.051, with AA_pI at -0.044 the single
# strongest property cell on the C panel), which is precisely where the one-hot
# puts its largest POSITIVE cell (G, +0.446): glycine is small, neutral and
# middling-pI, so in property space the amidation donor looks like weak negative
# evidence. Carrying the properties alongside identity may therefore be worse
# than not carrying them. In curated_wAAid their wrong-signed +1 contribution
# nearly vanished once identity was available (-0.051 -> -0.003), which is the
# model routing around them; dropping them tests whether it should not have to.
#
# The one-hot also strictly dominates the properties in information: the four
# property values are a BIJECTION over the 21 letters (aa_onehot_build.py checks
# exactly this, at max |diff| 0 over all 3,056,424 residues), so nothing is lost
# by the swap -- any function of the properties is a function of the letter.
# What changes is how easily the trunk can express it: a property is a smooth
# scalar the convolution can interpolate along, a one-hot is not.
#
# Coverage is 100%: every residue has a letter.

ARM_LABEL="curated, AA props -> one-hot"
PLM=${PLM:-$A2/lf_plm/aa_onehot.parquet}   # the same parquet curated_wAAid uses
DROP='AA_hydro,AA_charge,AA_mw,AA_pI'      # the four properties, and only those
N_CHAN_EXPECT=43

# the parquet already exists (curated_wAAid built it); rebuild only if missing
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

# This arm CANNOT use lf_winxgb.plm_residues: that appends the parquet to ALL 26
# curated channels and gives 47, while these trunks expect 43 with the four
# properties gone. residues_arm.py resolves each channel by NAME against the
# trained trunk's own channel_names, so it handles a drop and an add together.
mk_resid() { $PY_TF $S/residues_arm.py --arm $ARM --plm $PLM }

# No sweep yet; the curated model's conventional cut, matching curated_wAAid so
# the two identity arms are compared at the same T rather than a measured
# optimum. Sweep it the way esmonly36 was swept and change this line.
THRESH_DEFAULT=0.5
