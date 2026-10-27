# ARM: esmonly36 -- "ESM-C only" (Kevin, 2026-10-24)
#
# The flatten window head on the 36-component PCA of the ESM C per-residue
# embedding and NOTHING ELSE: all 25 curated channels dropped, `padding` kept
# because 10_8 refuses to drop the mask channel. 37 channels, 36 informative.
#
# This file states what train_esmonly36.sh / fullscan_esmonly36.sh already did,
# so the generalised runner reproduces that arm exactly. Those scripts remain on
# disk and still work; this is the same thing as a parameter set.

ARM_LABEL="ESM-C only"
PLM=${PLM:-$A2/lf_plm/esm_c_pca36.parquet}
# all of all_params3 except `padding`
DROP='cons_rs_n,min_afm,mean_afm,relASA,SS_P,SS_S,SS_E,SS_-,SS_T,SS_G,SS_B,SS_H,SS_I,Phi_cos,Psi_cos,Phi_sin,Psi_sin,NH->O_1_energy,O->NH_1_energy,NH->O_2_energy,O->NH_2_energy,AA_hydro,AA_charge,AA_mw,AA_pI'
N_CHAN_EXPECT=37

mk_plm() {
  (cd $W/inst/python && $PY_PLM -m lf_plm reduce \
      --h5 ${H5:-$A2/lf_plm/esm_c_raw.h5} --sequences $A2/lf_plm/sequences.parquet \
      --k ${K:-36} --refit --out $PLM)
}

# this arm CANNOT use lf_winxgb.plm_residues: that appends to all 26 curated
# channels and gives 62, while these trunks expect 37. residues_esmonly36.py
# builds `padding` + the 36 and asserts the order against the model.
mk_resid() { $PY_TF $S/residues_esmonly36.py }

# ESM-C's segment optimum, found by the T sweep at tolerance 4 (stable across
# cells over 0.575-0.725). It is the more liberal scorer -- mean score ~1.8x the
# curated model's -- so its best cut sits above 0.5.
THRESH_DEFAULT=0.65
