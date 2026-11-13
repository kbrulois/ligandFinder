#!/usr/bin/env python3
"""Read the deployed peptide-end CNN's layers out of its saved members.

Builds each arm/terminus from the run's own config.json with the code snapshot
that rendered the public gallery, loads real member weights, and writes
cnn_layers.json for cnn_arch.py. Also tests, on real windows, the claims the
figure makes in words. PYTHONPATH must point at the snapshot's inst/python.
"""
import json, pathlib, sys, os
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")
import h5py, keras, numpy as np
import lf_dcnn
from lf_dcnn.config import Config
from lf_dcnn.model import build_model
A = pathlib.Path.home()/"AF2_analysis"; HERE = pathlib.Path(__file__).parent
print("lf_dcnn from", lf_dcnn.__file__)
ARMS, TERMS, FOLDS = ("curated_AAswap", "esmonly36"), ("C", "N"), (1, 2, 3, 4, 5)
res = {}

def h5_shapes(p):
    out = {}
    with h5py.File(p, "r") as f:
        f.visititems(lambda n, o: out.__setitem__(n, o.shape) if isinstance(o, h5py.Dataset) else None)
    return out

for arm in ARMS:
    for T in TERMS:
        run = lambda k: A/f"lf_pepend_run_{T}_cv5_f{k}_nomask_flatten_{arm}_isolated"
        cfg_raw = json.load(open(run(1)/"in/config.json"))
        cfg = Config.from_dict(cfg_raw)
        m = build_model(cfg)
        layers = [dict(name=l.name, cls=type(l).__name__,
                       shape=str(tuple(l.output.shape[1:])) if hasattr(l, "output") else "?",
                       params=int(l.count_params())) for l in m.layers]
        # every deployed member must have the same weight layout, on every fold
        mem = sorted(p for k in FOLDS for p in (run(k)/"out/members").glob("*.weights.h5"))
        ref = h5_shapes(mem[0])
        assert all(h5_shapes(p) == ref for p in mem), "member layouts differ"
        cfgs = [json.load(open(run(k)/"in/config.json")) for k in FOLDS]
        drift = sorted({c for c in cfgs[0] for x in cfgs if x.get(c) != cfgs[0].get(c)})
        m.load_weights(mem[0])
        # end to end: this rebuild must reproduce the member's own saved predictions
        saved = np.load(str(mem[0]).replace(".weights.h5", ".npz"), allow_pickle=True)
        Xall = np.load(run(1)/"in/arrays.npz", allow_pickle=True)[f"{T}__all__x"][:2000].astype("float32")
        got = m.predict(Xall, verbose=0, batch_size=256)
        pred_err = float(np.abs(got["global"].ravel() - saved[f"{T}__global__x"][:2000]).max())
        pi_err = float(np.abs(got["per_index_cat"] - saved[f"{T}__per_index__x"][:2000]).max())
        # real windows, and the member's own saved predictions as a reference
        arrays = np.load(run(1)/"in/arrays.npz", allow_pickle=True)
        xkey = [k for k in arrays.files if arrays[k].ndim == 3 and arrays[k].shape[1:] == (36, cfg.n_channels)]
        X = arrays[xkey[0]][:256].astype("float32")
        sub = keras.Model(m.input, {n: m.get_layer(n).output for n in
              ("per_index_logits", "per_index_cat", "class_position_mask", "masked_sum", "flat_cls")})
        o = sub.predict(X, verbose=0, batch_size=64)
        L, ms = o["per_index_logits"], o["masked_sum"]
        sm = lambda a: np.exp(a - a.max(-1, keepdims=True)) / np.exp(a - a.max(-1, keepdims=True)).sum(-1, keepdims=True)
        names = list(cfg.pi_names)
        ip = names.index("padding"); keep = [i for i in range(len(names)) if i != ip]
        inone = names.index("none"); keep2 = [i for i in range(len(names)) if i != inone]
        check = dict(
            pi_names=names,
            head_classes=[names[i] for i in keep],
            members_checked=len(mem), cfg_keys_varying_across_folds=drift,
            mask_all_ones=bool(np.all(m.get_layer("class_position_mask").mask_matrix == 1)),
            mask_shape=list(m.get_layer("class_position_mask").mask_matrix.shape),
            sum_over_classes_is_1=float(np.abs(ms.sum(-1) - 1).max()),
            sum_over_positions_minus_1_max=float(np.abs(ms.sum(1) - 1).max()),
            drops_padding_err=float(np.abs(ms - sm(L[..., keep])).max()),
            drops_none_err=float(np.abs(ms - sm(L[..., keep2])).max()),
            per_index_softmax_over_classes=float(np.abs(o["per_index_cat"].sum(-1) - 1).max()),
            flat_dim=int(o["flat_cls"].shape[1]),
            flat_equals_ms=float(np.abs(o["flat_cls"] - ms.reshape(len(ms), -1)).max()),
            head_reads_trunk_directly=any(n.get("name") == "gap_bneck" for n in layers),
            total_params=int(m.count_params()),
            saved_global_max_err=pred_err, saved_per_index_max_err=pi_err,
            xkey=xkey[0],
        )
        res[f"{arm}|{T}"] = dict(layers=layers, check=check, config={k: cfg_raw[k] for k in cfg_raw
                                 if k not in ("mask_matrix_cat",)})
        print(f"\n{arm} {T}: total {check['total_params']:,}, {len(mem)} members")
        print({k: v for k, v in check.items() if k not in ("pi_names",)})
json.dump(res, open(HERE/"cnn_layers.json", "w"), indent=1)
print("\nLAYERS (curated_AAswap C):")
for l in res["curated_AAswap|C"]["layers"]:
    print(f"  {l['name']:22}{l['cls']:22}{l['shape']:14}{l['params']:>7,}")
print("\nparams differing between arms:", [(a['name'], a['params'], b['params']) for a, b in
      zip(res['curated_AAswap|C']['layers'], res['esmonly36|C']['layers']) if a['params'] != b['params']])
print("C vs N layer tables equal:", all(res[f"{a}|C"]["layers"] == res[f"{a}|N"]["layers"] for a in ARMS))
