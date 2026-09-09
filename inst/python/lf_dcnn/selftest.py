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
from .model import build_model
from .synthetic import make_data

#: parameter count of the default architecture; must equal the R model's
EXPECTED_PARAMS = 2438


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
    cfg = Config()
    for term, anchor in (("C", "CT_cleavage_context"), ("N", "NT_cleavage_context")):
        w = cfg.pi_weights(term)
        assert abs(w.mean() - 1.0) < 1e-6, w.mean()
        d = dict(zip(cfg.pi_names, w))
        assert d[anchor] > d["DB"] > d["gap"]
        assert abs(d[anchor] / d["gap"] - 3.0) < 1e-5
        assert abs(d["DB"] / d["gap"] - 2.0) < 1e-5


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


def check_position_ramp():
    import keras

    cfg = Config()
    model = build_model(cfg)
    ramp = keras.Model(model.input, model.get_layer("pos_ramp").output)
    got = np.asarray(ramp.predict(np.zeros((2, cfg.seq_len, cfg.n_channels), "float32"), verbose=0))
    want = np.linspace(0.0, 1.0, cfg.seq_len)
    assert np.allclose(got[0, :, 0], want, atol=1e-6), got[0, :3, 0]


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
    # calibration is monotone in the raw score
    order_raw = np.argsort(res.pred_raw, kind="stable")
    order_cal = np.argsort(res.pred, kind="stable")
    assert np.array_equal(order_raw, order_cal)
    assert abs(np.corrcoef(
        np.argsort(order_raw).astype(float), np.argsort(order_cal).astype(float)
    )[0, 1] - 1.0) < 1e-9
    assert res.term_index("N").size == res.n_by_term["N"]
    for term in res.term_order:
        assert "val_global_pr_auc" in res.histories[term], list(res.histories[term])
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
    assert build_model(Config(trunk="flat")).count_params() == EXPECTED_PARAMS

    # a depth that does not divide the window is rejected up front, not at build
    try:
        Config(trunk="unet", unet_filters=(8, 16, 32, 64), unet_dropout=(0.2,) * 4)
    except ValueError as e:
        assert "divisible" in str(e), e
    else:
        raise AssertionError("accepted a U-Net depth that does not divide seq_len")


def check_global_head_variants():
    """The window score can read the class softmax, the bottleneck, or both."""
    flat = build_model(Config(trunk="flat"))
    assert flat.count_params() == EXPECTED_PARAMS      # default is untouched
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
    check_oversampler,
    check_noise_only_on_continuous_channels,
    check_validate_catches_transpose,
    check_calibrator,
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
