# ARM: curated_AAswap_msacons_padg -- "AAswap, Aminode conservation -> MSA conservation",
# trained on the _padg windows (see INPUT_SFX below)
# (Kevin, 2026-11-10: "replace the existing conservation channel with
# conservation scores computed from alphafold db msa")
#
# curated_AAswap's channels with ONE change: `cons_rs_n` (Aminode evolutionary
# rate) is dropped and `plm_cons_msa` takes its place, so the comparison
# isolates the conservation source -- against curated_AAswap_padg, which has the
# same _padg windows (curated_AAswap itself has the old ones):
#
#   26 curated - 4 AA properties - cons_rs_n = 21, + 21 one-hot + 1 = 43 channels
#
# THE MSA. AFDB's own MSAs (msaUrl, v6) have returned 403 to anonymous callers
# since Aug 2026 (google-deepmind/alphafold#1137), so the alignments are the
# ColabFold/MMseqs2 ones (UniRef30 2023-02 + ColabFold env DB) from the LiteFold
# human-proteome set for 5,228 precursors whose sequence it holds exactly, and
# from the ColabFold API (same search, mode=env) for the other 388. The scorer
# takes any .a3m, so AFDB MSAs can be swapped in by listing their dir first.
#
# THE SCORE. Capra & Singh JSD, Henikoff weights, BLOSUM62 background, gap
# penalised, no smoothing, shrunk toward the background by a prior of 10
# pseudo-sequences; global min-max to [0,1]; 0 on the residues excluded
# by design, as cons_rs_n is. Code: $W/inst/scripts/msa_cons/.
#
# NOISE. plm_cons_msa is in lf_dcnn's DEFAULT_CONT_CHANNELS (since the
# 2026-11-10 conservation work), so it is augmented as cons_rs_n was. That needs
# W -- LF_SRC, or the checkout these scripts sit in -- to be a source tree that
# has it:
#   LF_SRC=<checkout> ARM=curated_AAswap_msacons_padg ./pepend_arm_pipeline.sh

## GUARD: in a source tree from before that, plm_cons_msa is NOT a
## noise-augmented channel; training there would silently differ from the
## baseline in one more way than intended. Refuse rather than warn. The test is
## the capability, not the checkout's name: plm_cons_msa inside the
## DEFAULT_CONT_CHANNELS tuple of $W's lf_dcnn/config.py.
sed -n '/^DEFAULT_CONT_CHANNELS/,/^)/p' ${W:-/nonexistent}/inst/python/lf_dcnn/config.py 2>/dev/null |
  grep -q '"plm_cons_msa"' || {
  echo "arm curated_AAswap_msacons_padg: W=${W:-<unset>} lacks plm_cons_msa in DEFAULT_CONT_CHANNELS;"
  echo "  run from a checkout that has it, or point LF_SRC at one"
  exit 1; }

ARM_LABEL="AAswap, MSA conservation (padg windows)"
PLM=${PLM:-$A2/lf_plm/aa_onehot_msacons.parquet}
DROP='cons_rs_n,AA_hydro,AA_charge,AA_mw,AA_pI'
N_CHAN_EXPECT=43

mk_plm() {
  local M=$W/inst/scripts/msa_cons
  local raw=$A2/lf_msa/msa_cons_raw.parquet
  ## rescore when any MSA is newer than the scores (a late download, a re-fetch)
  if [[ ! -s $raw || -n $(find $A2/lf_msa -name '*.a3m.gz' -newer $raw | head -1) ]]; then
    $PY_PLM $M/msa_conservation.py --workers 6 || return 1
  fi
  $PY_PLM $M/build_msacons_parquet.py --out $PLM
}

mk_resid() { $PY_TF $S/residues_arm.py --arm $ARM --plm $PLM }

# AAswap's cut, so the two are compared at the same T.
THRESH_DEFAULT=0.5

# trunks must noise-augment the replacement channel, as the baseline did cons_rs_n
CONT_EXPECT=plm_cons_msa

## the training windows: lf_pepend_nn_input_cv5_f<k>${INPUT_SFX}.rds, built by 10_7 with
##   --amid-g-anchor --pad-match --mono-match   (Kevin, 2026-11-10)
## INPUT_FLAGS is asserted against each input's own params before training.
INPUT_SFX=_padg
INPUT_FLAGS="amid_g=TRUE pad_match=TRUE mono_match=TRUE"
