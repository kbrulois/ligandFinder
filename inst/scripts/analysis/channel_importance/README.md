# Which channels the trunks rely on (2026-10-15..26)

**Asked:** how much does each input channel, or group of channels, matter to the
peptide-end trunks? Measured three ways: grouped permutation importance, Shapley
values, and the xgboost head's gain.

**Found:** not recorded with the scripts. The figures of the time were
`~/AF2_analysis/pepend_shapley_importance.svg` and `pepend_xgb_importance.svg`,
redrawn without the padded arms. The maintained attribution is
`evaluation/shap_positions.py` and `shap_channel_table.py`.

**Scripts:** `10_9z_grouped_importance.py`, `10_9z_shapley_importance.py` and
`10_9z_xgb_importance.py`, each with its `_plot.R`.
