"""PCA-reduce the raw embeddings to k channels and write the long residue table.

Why reduce: 1024 channels beside the 26 hand-built ones would put ~50k
parameters in the first conv of a model that has ~21k in total, against a few
dozen positives.  PCA fit ONCE over a residue sample of the whole set (the
same global-not-per-protein rule 9.2 applies to its min-max ranges), cached,
and applied identically to every protein.

Scaling: each component is z-scored by its own sd from the fit, clipped to
+/- ``clip`` sd, then mapped to [0, 1] -- the range every other channel is on.
Without the clip, min-max would let a handful of outlier residues squeeze the
bulk into a sliver.  Padding positions get 0 on the R side, like every other
channel, and the ``padding`` channel flags them.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def fit_pca(h5_path, k: int = 32, n_residues: int = 500_000, seed: int = 1,
            batch_size: int = 20_000, verbose: int = 1):
    """IncrementalPCA over a random ``n_residues`` sample of all residues.

    Returns ``(pca, sd)`` where ``sd`` is the per-component sd of the sampled
    scores, used for the z-scoring in :func:`transform_to_table`.
    """
    from sklearn.decomposition import IncrementalPCA

    rng = np.random.default_rng(seed)
    with h5py.File(Path(h5_path).expanduser(), "r") as h5:
        keys = list(h5.keys())
        total = sum(h5[a].shape[0] for a in keys)
        frac = min(1.0, n_residues / max(total, 1))
        if verbose:
            print(f"PCA: {len(keys)} proteins, {total} residues, sampling {frac:.1%} -> k={k}", flush=True)
        # two passes: fit on the sample, then score the sample for its sd
        pca = IncrementalPCA(n_components=k)
        buf = []
        n_buf = 0

        def flush():
            nonlocal buf, n_buf
            if n_buf >= k:                       # IncrementalPCA needs >= k rows per batch
                pca.partial_fit(np.concatenate(buf).astype("float32"))
            buf, n_buf = [], 0

        sample_rows = []
        for a in keys:
            x = h5[a][...]
            keep = rng.random(x.shape[0]) < frac
            if keep.any():
                s = x[keep].astype("float32")
                buf.append(s); n_buf += s.shape[0]
                sample_rows.append(s)
                if n_buf >= batch_size:
                    flush()
        flush()
    sample = np.concatenate(sample_rows)
    sd = pca.transform(sample).std(axis=0, ddof=1)
    if verbose:
        ev = pca.explained_variance_ratio_
        print(f"PCA: explained variance {ev[:5].round(3)} ... total {ev.sum():.3f} over {k} components",
              flush=True)
    return pca, sd


def save_pca(path, pca, sd, meta: dict | None = None) -> Path:
    import joblib

    p = Path(path).expanduser()
    p.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({"pca": pca, "sd": np.asarray(sd), "meta": meta or {}}, p)
    return p


def load_pca(path):
    import joblib

    d = joblib.load(Path(path).expanduser())
    return d["pca"], d["sd"], d.get("meta", {})


def transform_to_table(h5_path, sequences: pd.DataFrame, pca, sd, out_parquet,
                       clip: float = 3.0, prefix: str = "plm_", verbose: int = 1) -> Path:
    """``(accession, index, AA, plm_01..plm_k)`` for every residue, scaled to [0, 1].

    ``sequences`` supplies the letters (``accession``, ``sequence_uni``) so the
    table carries the same ``(index, AA)`` join keys ``secretome_aa`` has --
    a sequence-version mismatch then surfaces as an unmatched row, not as
    silently misaligned channels.
    """
    k = pca.n_components_
    cols = [f"{prefix}{i + 1:02d}" for i in range(k)]
    seq_of = dict(zip(sequences["accession"], sequences["sequence_uni"]))
    frames = []
    with h5py.File(Path(h5_path).expanduser(), "r") as h5:
        keys = [a for a in h5.keys() if a in seq_of]
        for i, a in enumerate(keys, 1):
            x = h5[a][...].astype("float32")
            seq = seq_of[a]
            if x.shape[0] != len(seq):
                raise ValueError(f"{a}: {x.shape[0]} embedded residues vs {len(seq)}-aa sequence")
            z = pca.transform(x) / sd
            z = np.clip(z, -clip, clip)
            u = (z + clip) / (2 * clip)                       # [0, 1]
            f = pd.DataFrame(u.astype("float32"), columns=cols)
            f.insert(0, "AA", list(seq))
            f.insert(0, "index", np.arange(1, len(seq) + 1, dtype="int32"))
            f.insert(0, "accession", a)
            frames.append(f)
            if verbose and i % 1000 == 0:
                print(f"  transformed {i}/{len(keys)}", flush=True)
    tbl = pd.concat(frames, ignore_index=True)
    out = Path(out_parquet).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    tbl.to_parquet(out, index=False)
    with open(out.with_suffix(".json"), "w") as fh:
        json.dump({"k": int(k), "clip": clip, "prefix": prefix, "columns": cols,
                   "n_proteins": len(keys), "n_residues": int(len(tbl)),
                   "explained_variance_ratio": [float(v) for v in pca.explained_variance_ratio_]},
                  fh, indent=2)
    if verbose:
        print(f"wrote {out}  ({len(tbl)} residues x {k} channels)", flush=True)
    return out
