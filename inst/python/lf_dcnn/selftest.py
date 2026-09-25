"""Built-in checks -- ``python -m lf_dcnn selftest``.

pytest is not installed in the target venv, so these run standalone.  The same
functions are collected by pytest if it is ever available (see
``inst/python/tests/test_lf_dcnn.py``).
"""

from __future__ import annotations

import shutil
import tempfile
import traceback
from pathlib import Path

import numpy as np

from .calibrate import PlattCalibrator
from .config import Config
from .data import OversampledWindows, TermArrays
from .losses import PerIndexCatLoss, WeightedBinaryCrossentropy, make_masked_cat_accuracy
from keras import ops

from .model import build_model
from .synthetic import make_data

#: parameter count of the default architecture (U-Net 16-32-64, no position ramp)
EXPECTED_PARAMS = 21406
#: the pre-port R model (flat trunk + position ramp), what Config.r_exact() builds
FLAT_PARAMS = 2438


def check_class_indices():
    cfg = Config()
    assert cfg.K_cat == 7 and cfg.K_pi == 8, (cfg.K_cat, cfg.K_pi)
    assert cfg.pi_names == (
        "CT_cleavage_context", "DB", "gap", "NT_cleavage_context",
        "pep_other", "pep_pocket", "none", "padding",
    ), cfg.pi_names
    # `none` sits at K_cat - 1, `padding` last -- the loss masks none and keeps
    # padding, the metric masks both.
    assert cfg.none_index == cfg.K_cat - 1 == 6
    assert cfg.padding_index == cfg.K_pi - 1 == 7
    assert cfg.pad_channel_id == cfg.n_channels - 1 == 25
    assert cfg.cont_channel_ids == (0, 1, 2, 3)


def check_position_masks():
    cfg = Config()
    m = cfg.mask_matrix_cat
    assert m.shape == (36, 7), m.shape
    names = [cfg.class_names[c] for c in cfg.cat_cols]
    col = {n: m[:, i] for i, n in enumerate(names)}
    assert np.array_equal(np.flatnonzero(col["NT_cleavage_context"]), np.arange(0, 5))
    assert np.array_equal(np.flatnonzero(col["CT_cleavage_context"]), np.arange(30, 36))
    assert np.array_equal(np.flatnonzero(col["gap"]), np.arange(5, 30))
    assert col["none"].sum() == 36
    # documented intent: DB at both termini
    assert col["DB"].sum() == 11, col["DB"].sum()
    # r_exact reproduces the R script's overwrite: DB gets the CT mask alone
    r_col = Config.r_exact().mask_matrix_cat[:, names.index("DB")]
    assert np.array_equal(np.flatnonzero(r_col), np.arange(30, 36))


def check_pi_weights():
    masked = Config()                        # the default: `none` masked out
    trained = Config(none_in_loss=True)      # opt-in: `none` trained, weighted down
    for term, anchor in (("C", "CT_cleavage_context"), ("N", "NT_cleavage_context")):
        w = masked.pi_weights(term)
        # with `none` at its nominal 1, the vector is normalised to mean 1
        assert abs(w.mean() - 1.0) < 1e-6, w.mean()
        d = dict(zip(masked.pi_names, w))
        assert d[anchor] > d["DB"] > d["gap"]
        assert abs(d[anchor] / d["gap"] - 3.0) < 1e-5
        assert abs(d["DB"] / d["gap"] - 2.0) < 1e-5

        # pi_weight_none is applied AFTER that normalisation, so the mean is no
        # longer 1 -- deliberately: what must hold is that every other class
        # keeps EXACTLY the weight the masked model trains with, so the two
        # differ by the `none` positions alone and not by a shifted loss scale
        t = trained.pi_weights(term)
        keep = [i for i, n in enumerate(trained.pi_names) if n != "none"]
        assert np.allclose(t[keep], w[keep]), (t, w)
        i_none = trained.pi_names.index("none")
        assert np.isclose(t[i_none], w[i_none] * trained.pi_weight_none)
        assert t.mean() < 1.0


def check_losses():
    cfg = Config()
    rng = np.random.default_rng(0)
    n, seq, K = 4, cfg.seq_len, cfg.K_pi
    cls = rng.integers(0, K, size=(n, seq))
    y_true = np.eye(K, dtype="float32")[cls]
    y_pred = rng.random((n, seq, K)).astype("float32")
    y_pred /= y_pred.sum(-1, keepdims=True)

    loss = PerIndexCatLoss(cfg.pi_weights("C"), cfg.none_index, cfg.gamma, cfg.smoothness_weight)
    base = float(loss(y_true, y_pred))
    assert np.isfinite(base) and base > 0

    # `none` positions are masked out: perturbing predictions only there must
    # not move the main term (it can still move smoothness, so compare a config
    # with the smoothness term off).
    hard = PerIndexCatLoss(cfg.pi_weights("C"), cfg.none_index, cfg.gamma, 0.0)
    before = float(hard(y_true, y_pred))
    y2 = y_pred.copy()
    none_pos = cls == cfg.none_index
    y2[none_pos] = np.eye(K, dtype="float32")[(cls[none_pos] + 1) % K]
    assert abs(float(hard(y_true, y2)) - before) < 1e-5

    # ...whereas padding positions ARE trained, so perturbing them must move it
    y3 = y_pred.copy()
    pad_pos = cls == cfg.padding_index
    y3[pad_pos] = np.eye(K, dtype="float32")[(cls[pad_pos] + 1) % K]
    assert abs(float(hard(y_true, y3)) - before) > 1e-4

    # metric masks BOTH padding and none
    acc = make_masked_cat_accuracy(cfg.none_index, cfg.padding_index)
    perfect = np.eye(K, dtype="float32")[cls]
    assert abs(float(acc(y_true, perfect)) - 1.0) < 1e-5
    wrong = np.eye(K, dtype="float32")[(cls + 1) % K]
    assert float(acc(y_true, wrong)) < 1e-5
    # only-real-classes-wrong drives it to 0; only-none/padding-wrong leaves it 1
    mixed = perfect.copy()
    mask_np = none_pos | pad_pos
    mixed[mask_np] = np.eye(K, dtype="float32")[(cls[mask_np] + 1) % K]
    assert abs(float(acc(y_true, mixed)) - 1.0) < 1e-5

    bce = WeightedBinaryCrossentropy(1.0, 1.0)
    yt = np.array([[1.0], [0.0]], dtype="float32")
    yp = np.array([[0.9], [0.2]], dtype="float32")
    want = float(np.mean([-np.log(0.9), -np.log(0.8)]))
    assert abs(float(bce(yt, yp)) - want) < 1e-6
    up = WeightedBinaryCrossentropy(1.0, 5.0)
    assert float(up(yt[:1], yp[:1])) > float(bce(yt[:1], yp[:1]))


def check_model_shapes():
    cfg = Config()
    model = build_model(cfg)
    assert model.count_params() == EXPECTED_PARAMS, model.count_params()
    x = np.zeros((3, cfg.seq_len, cfg.n_channels), dtype="float32")
    out = model.predict(x, verbose=0)
    assert set(out) == {"global", "per_index_cat"}, set(out)
    assert out["global"].shape == (3, 1), out["global"].shape
    assert out["per_index_cat"].shape == (3, cfg.seq_len, cfg.K_pi)
    # per_index_cat is a softmax
    assert np.allclose(out["per_index_cat"].sum(-1), 1.0, atol=1e-5)
    assert model.get_layer("embed").output.shape[-1] == cfg.embed_units


def check_include_none():
    """Dropping the `none` column must change only the softmax, not the masking.

    Uses ``none_in_loss = False`` throughout: with `none` masked out of the
    loss it is never a training target either way, so the set of scored
    positions and the loss on them must be identical and only the column count
    moves.
    """
    from .losses import PerIndexCatLoss, make_masked_cat_accuracy

    on, off = Config(none_in_loss=False), Config(include_none=False, none_in_loss=False)
    assert (on.K_cat, on.K_pi, on.none_index) == (7, 8, 6)
    assert (off.K_cat, off.K_pi, off.none_index) == (6, 7, -1)
    assert "none" not in off.pi_names and off.pi_names[-1] == "padding"
    assert off.mask_matrix_cat.shape == (off.seq_len, off.K_cat)
    # the first conv loses the dropped column's weights and nothing else moves
    assert build_model(off).count_params() < build_model(on).count_params()
    assert Config.r_exact().none_in_loss is False        # the R model masked `none`

    # same windows, the two label layouts: with-none one-hot vs an all-zero row
    rng = np.random.default_rng(3)
    n, seq = 4, on.seq_len
    cls = rng.integers(0, on.K_pi, size=(n, seq))
    y_on = np.eye(on.K_pi, dtype="float32")[cls]
    y_off = np.delete(y_on, on.none_index, axis=-1)          # none rows -> all zero
    assert np.allclose(y_off.sum(-1), (cls != on.none_index).astype("float32"))

    # A prediction that never puts its argmax on `none`, so the only thing left
    # that could differ between the two layouts is WHICH POSITIONS are scored.
    # (With mass on the `none` column the two genuinely differ, and in the
    # dropped-column model's favour: a scored position can no longer lose its
    # argmax to a class it is never trained to emit. That difference is the
    # point of the ablation, not an invariant.)
    p_on = np.asarray(rng.random((n, seq, on.K_pi)), dtype="float32")
    p_on[:, :, on.none_index] = 0.0
    p_on /= p_on.sum(-1, keepdims=True)
    p_off = np.delete(p_on, on.none_index, axis=-1)
    p_off /= p_off.sum(-1, keepdims=True)

    acc_on = make_masked_cat_accuracy(on.none_index, on.padding_index)
    acc_off = make_masked_cat_accuracy(off.none_index, off.padding_index)
    keep = (cls != on.none_index) & (cls != on.padding_index)
    assert keep.sum() > 0
    assert abs(float(acc_on(y_on, p_on)) - float(acc_off(y_off, p_off))) < 1e-5, (
        float(acc_on(y_on, p_on)), float(acc_off(y_off, p_off)))

    # and the loss masks the same positions: a loss with every `none` position
    # removed from the mask is unchanged when those rows go all-zero
    l_off = PerIndexCatLoss(off.pi_weights("C"), off.none_index, off.gamma, 0.0)
    v_all = float(ops.mean(l_off(y_off, p_off)))
    y_zero = y_off.copy(); y_zero[cls == on.none_index] = 0.0
    assert abs(v_all - float(ops.mean(l_off(y_zero, p_off)))) < 1e-6


def check_none_in_loss():
    """Training on `none` must add those positions and change nothing else."""
    from .losses import PerIndexCatLoss

    base, on = Config(), Config(none_in_loss=True)       # `none` masked by default
    # the real classes keep exactly the weights the default trains with; only
    # `none` moves, so the two arms differ by the `none` positions alone
    wb, wo = base.pi_weights("C"), on.pi_weights("C")
    keep = [i for i, n in enumerate(base.pi_names) if n != "none"]
    assert np.allclose(wb[keep], wo[keep]), (wb, wo)
    assert np.isclose(wo[on.none_index], wb[on.none_index] * on.pi_weight_none)
    assert on.pi_weight_none < 1.0, "a weight-1 `none` would swamp the real classes"

    rng = np.random.default_rng(11)
    n, seq = 6, base.seq_len
    cls = rng.integers(0, base.K_pi, size=(n, seq))
    y = np.eye(base.K_pi, dtype="float32")[cls]
    p_ = np.asarray(rng.random((n, seq, base.K_pi)), dtype="float32")
    p_ /= p_.sum(-1, keepdims=True)

    masked = PerIndexCatLoss(wb, base.none_index, base.gamma, 0.0, mask_none=True)
    trained = PerIndexCatLoss(wo, on.none_index, on.gamma, 0.0, mask_none=False)
    # with no `none` position at all the two must agree; introduce some and they
    # must not (the whole point), and the trained one must stay finite
    no_none = cls.copy(); no_none[no_none == base.none_index] = 0
    y2 = np.eye(base.K_pi, dtype="float32")[no_none]
    assert np.allclose(float(ops.mean(masked(y2, p_))), float(ops.mean(trained(y2, p_))), atol=1e-5)
    assert not np.isclose(float(ops.mean(masked(y, p_))), float(ops.mean(trained(y, p_))))
    assert np.isfinite(float(ops.mean(trained(y, p_))))


def check_position_ramp():
    import keras

    cfg = Config(position_ramp=True)
    model = build_model(cfg)
    ramp = keras.Model(model.input, model.get_layer("pos_ramp").output)
    got = np.asarray(ramp.predict(np.zeros((2, cfg.seq_len, cfg.n_channels), "float32"), verbose=0))
    want = np.linspace(0.0, 1.0, cfg.seq_len)
    assert np.allclose(got[0, :, 0], want, atol=1e-6), got[0, :3, 0]

    # position_ramp=False (the default): no ramp layer, and the first conv loses
    # exactly the ramp's kernel_size x filters weights -- nothing else moves
    off = build_model(cfg.evolve(position_ramp=False))
    assert "pos_ramp" not in [l.name for l in off.layers]
    first = cfg.unet_filters[0] if cfg.trunk == "unet" else cfg.conv_filters[0]
    assert model.count_params() - off.count_params() == cfg.conv_kernel * first, (
        model.count_params(), off.count_params())
    assert off.count_params() == EXPECTED_PARAMS
    assert off.get_layer("per_index_cat").output.shape[1:] == (cfg.seq_len, cfg.K_pi)
    # the R model: flat trunk, ramp on
    assert build_model(Config.r_exact()).count_params() == FLAT_PARAMS


def check_resample_reaches_training():
    """The per-epoch redraw must actually reach fit(), not just exist.

    This is the failure the port was written to fix: the R script's
    `on_epoch_begin` callback rebound its own `train_ds` variable, which `fit`
    no longer read, so the redraw never happened and every epoch trained on one
    fixed draw. A test that calls `_resample()` by hand would pass even if
    Keras never invoked it, so drive a real `fit()` and fingerprint the draw
    ORDER-INDEPENDENTLY -- `resample_each_epoch = False` still reshuffles, so a
    per-batch fingerprint cannot tell the two apart.
    """
    import keras

    from .data import OversampledWindows, TermArrays
    from .model import compile_model

    for flag, want in ((True, 4), (False, 1)):
        cfg = Config(epochs=4, patience=3, start_from_epoch=1, seed=42,
                     resample_each_epoch=flag)
        arrays = TermArrays.from_mapping(make_data(cfg, seed=5)["C"]["train"]).validate(cfg)
        ds = OversampledWindows(arrays, cfg, rng=np.random.default_rng(0))
        seen = []

        class Probe(keras.callbacks.Callback):
            def on_epoch_begin(self, epoch, logs=None):
                seen.append((round(float(ds._x.sum()), 2), float(ds._yg.sum())))

        compile_model(build_model(cfg), cfg, "C").fit(ds, epochs=4, verbose=0, callbacks=[Probe()])
        draws = {x for x, _ in seen}
        assert len(draws) == want, (flag, len(draws), want)
        # every positive is in every draw either way; only the negatives move
        assert len({p for _, p in seen}) == 1, seen


def check_oversampler():
    cfg = Config()
    data = make_data(cfg, seed=3)
    arrays = TermArrays.from_mapping(data["N"]["train"]).validate(cfg)
    ds = OversampledWindows(arrays, cfg, rng=np.random.default_rng(1))

    n_pos = int(arrays.y_global.sum())
    assert ds.n_selected == n_pos * 4
    assert ds.steps_per_epoch == int(np.ceil(n_pos * 4 / cfg.batch_size))
    assert len(ds) == ds.steps_per_epoch

    ys = np.concatenate([ds[i][1]["global"][:, 0] for i in range(len(ds))])
    assert ys.size == ds.n_selected
    assert abs(ys.sum() - n_pos) < 1e-6, ys.sum()   # every positive, exactly once

    xs = np.concatenate([ds[i][0] for i in range(len(ds))])
    # the padding channel itself is never jittered, so it stays a clean flag
    assert np.all(np.isin(xs[:, :, cfg.pad_channel_id], (0.0, 1.0)))

    first = xs.copy()
    ds.on_epoch_end()
    second = np.concatenate([ds[i][0] for i in range(len(ds))])
    assert not np.allclose(first, second), "resample_each_epoch did not redraw"

    fixed = OversampledWindows(arrays, cfg.evolve(resample_each_epoch=False), rng=np.random.default_rng(1))
    a = np.concatenate([fixed[i][0] for i in range(len(fixed))])
    fixed.on_epoch_end()
    b = np.concatenate([fixed[i][0] for i in range(len(fixed))])
    # same windows, re-ordered: the draw is fixed but batch composition still varies
    assert not np.allclose(a, b), "fixed draw did not reshuffle between epochs"
    assert np.allclose(np.sort(a, axis=0), np.sort(b, axis=0)), \
        "resample_each_epoch=False redrew instead of merely reshuffling"


def check_noise_only_on_continuous_channels():
    cfg = Config()
    data = make_data(cfg, seed=5)
    arrays = TermArrays.from_mapping(data["C"]["train"]).validate(cfg)
    plain = OversampledWindows(arrays, cfg.evolve(noise_frac=0.0), rng=np.random.default_rng(7))
    noisy = OversampledWindows(arrays, cfg, rng=np.random.default_rng(7))
    # same rng seed => the index draw is identical, so the two differ only by noise
    a = np.concatenate([plain[i][0] for i in range(len(plain))])
    b = np.concatenate([noisy[i][0] for i in range(len(noisy))])
    changed = np.flatnonzero(np.abs(a - b).sum(axis=(0, 1)) > 0)
    assert set(changed.tolist()) <= set(cfg.cont_channel_ids), changed
    assert changed.size > 0, "noise augmentation changed nothing"

    # ...and never at padding positions, on any channel
    pad = a[:, :, cfg.pad_channel_id] == 1
    assert pad.any(), "synthetic data has no padding positions to check"
    assert np.array_equal(a[pad], b[pad]), "noise leaked onto padding positions"


def check_validate_catches_transpose():
    cfg = Config()
    bad = {
        "x": np.zeros((cfg.seq_len, cfg.n_channels, 5), "float32"),
        "y_global": np.zeros((5, 1), "float32"),
        "y_per_index_cat": np.zeros((5, cfg.seq_len, cfg.K_pi), "float32"),
    }
    try:
        TermArrays.from_mapping(bad).validate(cfg)
    except ValueError as e:
        assert "aperm" in str(e)
    else:
        raise AssertionError("validate() accepted a pre-aperm array")


def check_calibrator():
    rng = np.random.default_rng(11)
    s = rng.normal(size=4000)
    a, b = -0.7, 2.3
    y = (rng.random(4000) < 1 / (1 + np.exp(-(a + b * s)))).astype("float64")

    irls = PlattCalibrator.fit(s, y)
    assert irls.converged, irls
    assert abs(irls.intercept - a) < 0.15, irls.intercept
    assert abs(irls.slope - b) < 0.2, irls.slope

    lbfgs = PlattCalibrator.fit(s, y, method="lbfgs")
    assert abs(lbfgs.intercept - irls.intercept) < 1e-4
    assert abs(lbfgs.slope - irls.slope) < 1e-4

    # monotone: calibration cannot reorder, so AUC is untouched
    raw = rng.normal(size=500)
    cal = irls.predict(raw)
    assert np.array_equal(np.argsort(raw), np.argsort(cal))

    # perfectly separable data diverges; R stops at 25 iterations and so do we
    s2 = np.array([-3.0, -2.0, -1.0, 1.0, 2.0, 3.0])
    y2 = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    sep = PlattCalibrator.fit(s2, y2)
    assert sep.iterations <= 25
    assert np.all(np.diff(sep.predict(np.sort(s2))) >= 0)

    try:
        PlattCalibrator.fit([0.1, 0.2], [1.0, 1.0])
    except ValueError as e:
        assert "one class" in str(e)
    else:
        raise AssertionError("calibrator accepted single-class labels")


def _tiny_cfg() -> Config:
    return Config(epochs=3, patience=2, start_from_epoch=1, seed=42)


def check_end_to_end():
    from .pipeline import run

    cfg = _tiny_cfg()
    data = make_data(cfg, seed=17)
    res = run(data, cfg, verbose=0)

    n = sum(res.n_by_term[t] for t in res.term_order)
    assert res.term_order == ("N", "C"), res.term_order
    assert res.pred.shape == (n,)
    assert res.pred_raw.shape == (n,)
    assert res.per_index.shape == (n, cfg.seq_len, cfg.K_pi)
    assert res.emb.shape == (n, cfg.embed_units)
    assert np.allclose((res.emb**2).sum(1), 1.0, atol=1e-6), "embeddings not L2-normalised"
    assert np.all((res.pred >= 0) & (res.pred <= 1))
    # Calibration is MONOTONE in the raw score, so it cannot reorder. The sign
    # is not guaranteed: on a degenerate fit (few val positives, or a short run
    # where val scores anti-correlate with labels) the logistic slope comes out
    # negative, which reverses the ranking without breaking monotonicity.
    order_raw = np.argsort(res.pred_raw, kind="stable")
    order_cal = np.argsort(res.pred, kind="stable")
    assert np.array_equal(order_raw, order_cal) or np.array_equal(
        order_raw, np.argsort(-res.pred, kind="stable")
    ), f"calibration reordered (slope {res.calibrator.slope:.3f})"
    assert abs(abs(np.corrcoef(
        np.argsort(order_raw).astype(float), np.argsort(order_cal).astype(float)
    )[0, 1]) - 1.0) < 1e-9
    assert res.term_index("N").size == res.n_by_term["N"]
    for term in res.term_order:
        assert "val_global_pr_auc" in res.histories[term], list(res.histories[term])
    # single member: the member rows ARE the means
    assert res.val_scores_members.shape == (1, res.val_scores.size)
    assert np.allclose(res.val_scores_members[0], res.val_scores)
    assert np.allclose(res.pred_raw_members[0], res.pred_raw)

    # ensemble: mean/sd over the member rows must reproduce pred_raw / pred_sd,
    # and the members must actually differ (different seeds)
    res3 = run(data, cfg, verbose=0, n_seeds=3)
    assert res3.pred_raw_members.shape == (3, n)
    assert res3.val_scores_members.shape == (3, res3.val_scores.size)
    assert np.allclose(res3.pred_raw_members.mean(0), res3.pred_raw, atol=1e-6)
    assert np.allclose(res3.pred_raw_members.std(0, ddof=1), res3.pred_sd, atol=1e-5)
    assert np.allclose(res3.val_scores_members.mean(0), res3.val_scores, atol=1e-6)
    assert not np.allclose(res3.pred_raw_members[0], res3.pred_raw_members[1])
    assert res3.params == {t: res3.models[t].count_params() for t in res3.term_order}


def check_members_roundtrip():
    """run() == run_member() x n + combine(), also through the on-disk member files."""
    import tempfile

    from .io import load_members, save_member
    from .pipeline import combine, run, run_member, set_seed

    cfg = _tiny_cfg()
    data = make_data(cfg, seed=17)
    ref = run(data, cfg, verbose=0, n_seeds=2)

    set_seed(cfg.seed)
    members = [run_member(data, cfg, k=k, n_seeds=2, verbose=0, keep_models=False) for k in range(2)]
    assert [m.seed for m in members] == [cfg.seed, cfg.seed + 1]
    with tempfile.TemporaryDirectory() as d:
        for m in members:
            assert len(save_member(d, m)) == len(m.terms)     # one file pair per terminus
        back = load_members(d, term_order=cfg.term_order)
    assert [m.k for m in back] == [0, 1]
    assert back[0].terms == ref.term_order

    # a terminus trained on its own (one model per process, the CLI's
    # --isolated unit) reproduces the same terminus of the full member: the
    # per-term rng offset keys on the FULL term order, not on what is trained
    set_seed(cfg.seed)
    solo = run_member(data, cfg, k=1, n_seeds=2, verbose=0, keep_models=False, terms=["C"])
    assert solo.terms == ("C",)
    assert np.allclose(solo.global_all["C"], members[1].global_all["C"], atol=1e-6)
    with tempfile.TemporaryDirectory() as d:
        save_member(d, run_member(data, cfg, k=0, n_seeds=2, verbose=0, keep_models=False, terms=["N"]))
        save_member(d, run_member(data, cfg, k=0, n_seeds=2, verbose=0, keep_models=False, terms=["C"]))
        merged = load_members(d, term_order=cfg.term_order)
    assert len(merged) == 1 and merged[0].terms == ref.term_order
    assert not back[0].models                                   # nothing live survives disk
    res = combine(back, cfg)
    # numerics through combine must match run()'s -- same seeds, same reduction
    assert res.n_seeds == 2 and res.term_order == ref.term_order
    assert np.allclose(res.pred_raw, ref.pred_raw, atol=1e-6)
    assert np.allclose(res.pred_sd, ref.pred_sd, atol=1e-6)
    assert np.allclose(res.per_index, ref.per_index, atol=1e-6)
    assert np.allclose(res.val_scores_members, ref.val_scores_members, atol=1e-6)
    assert np.allclose(res.emb, ref.emb, atol=1e-6)
    assert res.params == ref.params
    assert res.histories.keys() == ref.histories.keys()
    # members in any order combine the same
    assert np.allclose(combine(back[::-1], cfg).pred_raw, res.pred_raw)
    return res


def check_roundtrip_npz(reference=None):
    from .io import load_inputs, load_outputs, save_inputs, save_outputs
    from .pipeline import run
    import pandas as pd

    cfg = _tiny_cfg()
    data = make_data(cfg, seed=17)
    n = sum(len(data[t]["all"]["y_global"]) for t in cfg.term_order)
    meta = pd.DataFrame({"peps": [f"w{i}" for i in range(n)]})

    tmp = Path(tempfile.mkdtemp(prefix="lf_dcnn_"))
    try:
        save_inputs(tmp / "in", data, cfg, meta=meta)
        data2, cfg2, meta2 = load_inputs(tmp / "in")
        assert cfg2 == cfg, "config did not survive the round trip"
        assert meta2 is not None and list(meta2["peps"]) == list(meta["peps"])
        for t in cfg.term_order:
            for s in ("train", "val", "all"):
                for k in ("x", "y_global", "y_per_index_cat"):
                    assert np.array_equal(data[t][s][k], data2[t][s][k]), (t, s, k)

        res2 = run(data2, cfg2, verbose=0)
        out = save_outputs(tmp / "out", res2, meta=meta2)
        back = load_outputs(out)
        assert np.allclose(back["pred"], res2.pred)
        assert back["predictions"].shape[0] == n
        assert list(back["predictions"].columns[:1]) == ["peps"]
        assert back["calibrator"]["method"] == "irls"

        if reference is not None:
            # the standalone path must reproduce the in-process path exactly
            assert np.array_equal(reference.pred_raw, res2.pred_raw), (
                "npz path diverged from the in-process path"
            )
            assert np.array_equal(reference.pred, res2.pred)
            assert np.array_equal(reference.emb, res2.emb)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def check_unet_trunk():
    """The pooling trunk must still deliver per-residue output at seq_len."""
    cfg = Config(trunk="unet")
    model = build_model(cfg)
    by_name = {l.name: l for l in model.layers}

    # down the encoder, then back up: 36 -> 18 -> 9 -> 18 -> 36
    assert by_name["enc1"].output.shape[1] == 36
    assert by_name["pool1"].output.shape[1] == 18
    assert by_name["enc2"].output.shape[1] == 18
    assert by_name["pool2"].output.shape[1] == 9
    assert by_name["bottleneck"].output.shape[1] == 9
    assert by_name["bottleneck"].output.shape[2] == cfg.unet_filters[-1]
    assert by_name["up2"].output.shape[1] == 18
    assert by_name["up1"].output.shape[1] == 36
    # skips concatenate decoder + encoder channels
    assert by_name["skip2"].output.shape[2] == cfg.unet_filters[2] + cfg.unet_filters[1]
    assert by_name["skip1"].output.shape[2] == cfg.unet_filters[1] + cfg.unet_filters[0]

    x = np.zeros((3, cfg.seq_len, cfg.n_channels), dtype="float32")
    out = model.predict(x, verbose=0)
    assert out["per_index_cat"].shape == (3, cfg.seq_len, cfg.K_pi), out["per_index_cat"].shape
    assert out["global"].shape == (3, 1)
    assert np.allclose(out["per_index_cat"].sum(-1), 1.0, atol=1e-5)
    assert model.get_layer("embed").output.shape[-1] == cfg.embed_units

    # the flat trunk must be untouched by any of this
    assert build_model(Config.r_exact()).count_params() == FLAT_PARAMS

    # a depth that does not divide the window is rejected up front, not at build
    try:
        Config(trunk="unet", unet_filters=(8, 16, 32, 64), unet_dropout=(0.2,) * 4)
    except ValueError as e:
        assert "divisible" in str(e), e
    else:
        raise AssertionError("accepted a U-Net depth that does not divide seq_len")


def check_global_head_variants():
    """The window score can read the class softmax, the bottleneck, or both."""
    flat = build_model(Config.r_exact())
    assert flat.count_params() == FLAT_PARAMS          # the R model is untouched
    assert "gap_bneck" not in {l.name for l in flat.layers}

    seen = {}
    for head, has_attn, embed_in in [
        ("attn", True, 7),          # masked class softmax -> attention -> pool
        ("bottleneck", False, 64),  # straight off the U-Net bottleneck
        ("both", True, 71),         # concatenation of the two
    ]:
        cfg = Config(trunk="unet", global_head=head)
        m = build_model(cfg)
        names = {l.name for l in m.layers}
        assert ("attn" in names) is has_attn, (head, names & {"attn"})
        assert ("gap_bneck" in names) is (head != "attn"), head
        assert m.get_layer("embed").input.shape[-1] == embed_in, (head, m.get_layer("embed").input.shape)
        # the per-residue head is unaffected by where the score reads from
        out = m.predict(np.zeros((2, cfg.seq_len, cfg.n_channels), "float32"), verbose=0)
        assert out["per_index_cat"].shape == (2, cfg.seq_len, cfg.K_pi)
        assert out["global"].shape == (2, 1)
        seen[head] = m.count_params()
    assert seen["both"] > seen["bottleneck"] > seen["attn"], seen

    # a bottleneck head needs a trunk that has one
    try:
        Config(trunk="flat", global_head="bottleneck")
    except ValueError as e:
        assert "trunk='unet'" in str(e), e
    else:
        raise AssertionError("accepted a bottleneck head on the flat trunk")


def check_unet_trains():
    """A short end-to-end run on the pooling trunk."""
    from .pipeline import run

    cfg = _tiny_cfg().evolve(trunk="unet")
    res = run(make_data(cfg, seed=17), cfg, verbose=0)
    n = sum(res.n_by_term[t] for t in res.term_order)
    assert res.per_index.shape == (n, cfg.seq_len, cfg.K_pi)
    assert res.emb.shape == (n, cfg.embed_units)
    assert np.all(np.isfinite(res.pred)) and np.all(np.isfinite(res.pred_raw))


CHECKS = [
    check_class_indices,
    check_position_masks,
    check_pi_weights,
    check_losses,
    check_model_shapes,
    check_unet_trunk,
    check_global_head_variants,
    check_unet_trains,
    check_position_ramp,
    check_include_none,
    check_none_in_loss,
    check_oversampler,
    check_resample_reaches_training,
    check_noise_only_on_continuous_channels,
    check_validate_catches_transpose,
    check_calibrator,
    check_members_roundtrip,
]


def main(verbose: int = 1) -> int:
    failures = 0
    reference = None
    for fn in CHECKS:
        try:
            fn()
        except Exception:
            failures += 1
            print(f"FAIL {fn.__name__}")
            traceback.print_exc()
        else:
            if verbose:
                print(f"ok   {fn.__name__}")

    # these two are ordered: the round trip compares against the in-process run
    try:
        reference = check_end_to_end()
    except Exception:
        failures += 1
        print("FAIL check_end_to_end")
        traceback.print_exc()
    else:
        if verbose:
            print("ok   check_end_to_end")

    try:
        check_roundtrip_npz(reference)
    except Exception:
        failures += 1
        print("FAIL check_roundtrip_npz")
        traceback.print_exc()
    else:
        if verbose:
            print("ok   check_roundtrip_npz")

    total = len(CHECKS) + 2
    print(f"\n{total - failures}/{total} checks passed")
    return 1 if failures else 0
