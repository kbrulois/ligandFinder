"""The two residue classifiers, ported from 10_score_AA.R and its xgboost twin.

Both are per-residue binary classifiers. Neither sees a window: the MLP is a
plain stack of dense layers over one residue's feature row, and xgboost is
trees over the same row.

Kept faithful to the R, including the odd bits:
  * MLP 256 -> 64 -> 4 -> 1, leaky_relu, L2 1e-2 on every layer, BatchNorm
    before each and Dropout 0.5 after the first three, Adam at 1e-4.
  * xgboost eta 0.08, subsample 0.5, max_depth 6, lambda 0, alpha 1,
    min_child_weight 0, tree_method exact, 500 rounds, 100-round early stop.
  * class weights balanced from the training counts, as the R does.

ONE DELIBERATE DEPARTURE, in `MONITOR_DEFAULT`. The R early-stops the MLP on
``val_binary_accuracy``. These labels run at ~4% positive, where predicting "no"
for every residue scores 96% -- so that monitor rewards the degenerate model and
restores its weights. The default here is ``val_auc``; pass
``--monitor val_binary_accuracy`` to reproduce the R exactly.
"""
from __future__ import annotations

import numpy as np

MONITOR_DEFAULT = "val_auc"


def class_weights(y: np.ndarray) -> dict[int, float]:
    """Balanced weights, as the R computes them from the training counts."""
    n0 = float((y == 0).sum())
    n1 = float((y == 1).sum())
    if n1 == 0:
        return {0: 1.0, 1: 1.0}
    return {0: 1.0, 1: n0 / n1}


def build_mlp(n_feat: int, l2: float = 0.01, dropout: float = 0.5,
              units=(256, 64, 4), lr: float = 1e-4):
    import keras
    from keras import layers, regularizers

    reg = lambda: regularizers.l2(l2)
    m = keras.Sequential(name="lf_resid_mlp")
    m.add(layers.Input(shape=(n_feat,)))
    m.add(layers.BatchNormalization())
    for u in units:
        m.add(layers.Dense(u, activation="leaky_relu", kernel_regularizer=reg()))
        m.add(layers.BatchNormalization())
        m.add(layers.Dropout(dropout))
    m.add(layers.Dense(1, activation="sigmoid", kernel_regularizer=reg()))
    m.compile(
        optimizer=keras.optimizers.Adam(learning_rate=lr),
        loss=keras.losses.BinaryCrossentropy(),
        metrics=[keras.metrics.BinaryAccuracy(name="binary_accuracy"),
                 keras.metrics.SpecificityAtSensitivity(0.8, name="spec_at_sens"),
                 keras.metrics.SensitivityAtSpecificity(0.8, name="sens_at_spec"),
                 keras.metrics.AUC(name="auc")],
    )
    return m


def fit_mlp(x_tr, y_tr, x_val, y_val, *, epochs=100, patience=40,
            start_from_epoch=20, monitor=MONITOR_DEFAULT, verbose=0, seed=0):
    import keras
    # An ensemble builds many models in one process, and keras state accumulates
    # across them: the 13th MLP of a 20-seed run died inside the PR-curve AUC
    # metric with "Incompatible shapes: [0] vs. [199]" while the 10th of a
    # 5-seed run was fine. lf_dcnn hit the same class of bug and answered it with
    # one process per model; here clearing the session between models is enough,
    # and the PR metric that actually blew up is gone from compile() -- it was
    # redundant anyway, since PR AUC is computed from the predictions with
    # sklearn.
    keras.backend.clear_session()
    keras.utils.set_random_seed(int(seed))
    m = build_mlp(x_tr.shape[1])
    mode = "min" if monitor.endswith("loss") else "max"
    cb = keras.callbacks.EarlyStopping(monitor=monitor, mode=mode, patience=patience,
                                       start_from_epoch=start_from_epoch,
                                       restore_best_weights=True)
    h = m.fit(x_tr, y_tr, epochs=epochs, class_weight=class_weights(y_tr),
              validation_data=(x_val, y_val), callbacks=[cb], verbose=verbose,
              batch_size=256)
    return m, {k: [float(v) for v in vs] for k, vs in h.history.items()}


def fit_xgb(x_tr, y_tr, x_val, y_val, *, nrounds=500, early_stopping_rounds=100,
            verbose=False, seed=0, feat_names=None):
    import xgboost as xgb

    dtr = xgb.DMatrix(x_tr, label=y_tr, feature_names=feat_names)
    dva = xgb.DMatrix(x_val, label=y_val, feature_names=feat_names)
    params = {
        "objective": "binary:logistic", "eta": 0.08, "subsample": 0.5,
        "min_child_weight": 0, "gamma": 0, "lambda": 0, "alpha": 1,
        "tree_method": "exact", "max_depth": 6, "eval_metric": ["auc", "aucpr"],
        "seed": int(seed),
        # the R leaves this at 1 and leans on class_weight only for the MLP;
        # without it the trees see a 4% positive rate and barely split
        "scale_pos_weight": float((y_tr == 0).sum()) / max(1.0, float((y_tr == 1).sum())),
    }
    ev = {}
    bst = xgb.train(params, dtr, num_boost_round=nrounds,
                    evals=[(dtr, "train"), (dva, "eval")],
                    early_stopping_rounds=early_stopping_rounds,
                    evals_result=ev, verbose_eval=verbose)
    return bst, ev


def predict_mlp(model, x, batch_size: int = 8192) -> np.ndarray:
    return np.asarray(model.predict(x, verbose=0, batch_size=batch_size)).reshape(-1)


def predict_xgb(bst, x, feat_names=None) -> np.ndarray:
    import xgboost as xgb
    d = xgb.DMatrix(x, feature_names=feat_names)
    best = getattr(bst, "best_iteration", None)
    kw = {"iteration_range": (0, best + 1)} if best is not None else {}
    return np.asarray(bst.predict(d, **kw)).reshape(-1)
