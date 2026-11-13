# The C-terminal cleavage motif on its own (2026-10-16)

**Asked:** how well does the cleavage motif at an anchor rank C ends by itself
(dibasic, G + dibasic for amidation, monobasic)? Does a solvent-exposure floor add
anything?

**Found:**
- The amidation motif alone reaches AP 0.59 on held-out C ends.
- A dibasic prior does worse than no prior, because dibasic sites are exactly the
  hard negatives.
- Motif classes over the 2,980,852 C anchors: dibasic 33,970; G + dibasic 2,444;
  monobasic 269,303; terminus 5,616; other 2,669,519.

**Outputs:** figures `~/AF2_analysis/pepend_motif_prior_pr.*` and
`pepend_motif_relasa_pr.*`; tables in `~/AF2_analysis/pepend_motif/`.

**Scripts:**
- `10_9z_motif_codes.R`: the motif class at every scanned anchor.
- `10_9z_motif_pr.py` and `_plot.R`: the motif as a ranking.
- `10_9z_motif_ppv.py` and `_plot.R`: precision of each motif.
- `10_9z_relasa_window.R` and `10_9z_motif_asa_plot.R`: the solvent-exposure floor.
