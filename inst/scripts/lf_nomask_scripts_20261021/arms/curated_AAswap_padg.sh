# ARM: curated_AAswap_padg -- curated_AAswap's channels on the _padg windows
# (Kevin, 2026-11-10). Same 43 channels, same PLM parquet, same DROP as
# curated_AAswap; ONLY the training windows differ, built by 10_7 with
#
#   --amid-g-anchor  C ends followed by G + dibasic are anchored ON the G (the G
#                    labelled peptide); db_amid leaves the training negatives
#   --pad-match      negatives topped up until every padding profile is as
#                    common among negatives as among positives (the AAswap trunks
#                    scored 17% of unlabelled extracellular precursor-terminus
#                    anchors > 0.5, vs ~0.1% of interior ones)
#   --mono-match     monobasic-anchored negatives topped up to the positives' share
#
# Against curated_AAswap this measures the window fixes; against
# curated_AAswap_msacons_padg it is the baseline for the conservation swap.
# The C folds are pinned to the old ones (10_7 split_by_mature), so all three
# arms share one cross-validation partition.
ARM_LABEL="curated, AA props -> one-hot (padg windows)"
PLM=${PLM:-$A2/lf_plm/aa_onehot.parquet}
DROP='AA_hydro,AA_charge,AA_mw,AA_pI'
N_CHAN_EXPECT=43
mk_plm() {
  [[ -s $PLM ]] || { echo "missing $PLM: build it with ARM=curated_wAAid stage 0"; return 1; }
}
mk_resid() { $PY_TF $S/residues_arm.py --arm $ARM --plm $PLM }
THRESH_DEFAULT=0.5

## the windows need lf_dcnn / 10_7 from a source tree that has them (10_7
## --amid-g-anchor, 2026-11-10): $W, i.e. LF_SRC or the checkout this sits in
grep -q 'amid-g-anchor' ${W:-/nonexistent}/inst/scripts/10_7_pepend_windows.R 2>/dev/null || {
  echo "arm curated_AAswap_padg: W=${W:-<unset>} predates the _padg windows (its 10_7 has no --amid-g-anchor)"
  echo "  run from a checkout that has them, or point LF_SRC at one"
  exit 1; }
INPUT_SFX=_padg
INPUT_FLAGS="amid_g=TRUE pad_match=TRUE mono_match=TRUE"
