# ARM: msa_gap -- "MSA + gap": arm msa exactly, plus a `gap` class. AAswap's
# channels with MSA conservation in place of Aminode,
# trained on motif-matched negatives, every end anchored on its own last / first
# residue (Kevin, 2026-11-12: "i want both arms to use plm_cons_msa instead of
# (not in addition to) cons_rs and cons_rs_n. both arms should have equal
# representation of padded, monobasic, dibasic, amidation motif negative
# windows. both arms should have the peptide end be the anchor point").
#
#   26 curated - 4 AA properties - cons_rs_n = 21, + 21 one-hot + plm_cons_msa = 43
#
# cons_rs is no curated channel (cons_rs_n is the one the trunks see), so
# dropping cons_rs_n leaves plm_cons_msa the only conservation channel.
#
# WINDOWS (INPUT_SFX=_msa). 10_7 with --pad-match --mono-match --db-match
# --amid-match: negatives topped up until padded windows and windows at a
# monobasic, dibasic or amidation (G + dibasic) site are at least as common
# among negatives as among positives. Mature anchors throughout -- the G-anchored
# C end of the _padg arms was removed from the code on 2026-11-12.
#
# FEATURES (v2). Built on lf_pepend_residue_cache_v2.rds, the 2026-11-12 fix of
# the residue-feature join: each secretome_aa row is placed at its UniProt
# position instead of its AlphaFold-model index, which differed for ~160
# precursors whose model was built on another sequence (AGT, VEGFA, ...). So the
# curated residue table is lf_pepend_residues_v2.npz (RESID_BASE) and the
# conservation channel is re-masked against it, on the SCALE the _padg msacons
# trunks used (aa_onehot_msacons.parquet.scale.json, pinned).
#
# GAP (Kevin, 2026-11-12: "add an additional class prediction called 'gap' and
# define it as any residue between a dibasic site and the peptide (usually
# just a single G residue)"). 10_7 --gap-class labels those residues `gap`
# instead of cleavage context (lf_pepend_gap(): 31 of 86 C windows -- 30 single
# G, one GFL -- and 3 of 92 N windows) and saves the 7-class vocabulary
# LF_PEPEND_CLASSES_GAP. Windows, splits and negatives are arm msa's exactly
# (same seed and flags), so the two arms differ in the labels alone.
#
# THE SCORE and THE MSA: as inst/scripts/generate_residue_db/msa_cons/ documents (ColabFold MSAs,
# Capra-Singh JSD, Henikoff weights, BLOSUM62 background, a 10-pseudo-sequence
# prior; 0 on residues excluded by design).

## GUARD: plm_cons_msa must be a noise-augmented channel of this source tree
sed -n '/^DEFAULT_CONT_CHANNELS/,/^)/p' ${W:-/nonexistent}/inst/python/lf_dcnn/config.py 2>/dev/null |
  grep -q '"plm_cons_msa"' || {
  echo "arm msa_gap: W=${W:-<unset>} lacks plm_cons_msa in DEFAULT_CONT_CHANNELS;"
  echo "  run from a checkout that has it, or point LF_SRC at one"
  exit 1; }

ARM_LABEL="MSA + gap"
PLM=${PLM:-$A2/lf_plm/aa_onehot_msacons_v2.parquet}
DROP='cons_rs_n,AA_hydro,AA_charge,AA_mw,AA_pI'
N_CHAN_EXPECT=43
## the curated residue table this arm's channels come from (pepend_arm_pipeline.sh
## stage 2 checks the arm's table against it)
RESID_BASE=$A2/lf_pepend_residues_v2.npz

mk_plm() {
  $PY_PLM $W/inst/scripts/generate_residue_db/msa_cons/build_msacons_parquet.py --residues-npz $RESID_BASE \
    --scale-from $A2/lf_plm/aa_onehot_msacons.parquet.scale.json --out $PLM
}

mk_resid() { $PY_TF $LF_RDB/residues_arm.py --arm $ARM --plm $PLM --residues-npz $RESID_BASE }

THRESH_DEFAULT=0.5

# trunks must noise-augment the conservation channel
CONT_EXPECT=plm_cons_msa

## the training windows: lf_pepend_nn_input_cv5_f<k>${INPUT_SFX}.rds, built by 10_7 with
##   --n-folds 5 --fold k --near-n 100 --seed 42 --pad-match --mono-match --db-match
##   --amid-match --gap-class --feature-cache ~/AF2_analysis/lf_pepend_residue_cache_v2.rds
## INPUT_FLAGS is asserted against each input's own params before training.
INPUT_SFX=_msagap
INPUT_FLAGS="pad_match=TRUE mono_match=TRUE db_match=TRUE amid_match=TRUE gap_class=TRUE dssp_index=mapped"
