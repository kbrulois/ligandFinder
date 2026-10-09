"""Per-residue embeddings for a table of sequences, one (L, D) array per accession.

Whole precursors are embedded (signal peptide included): the models were
trained on full UniProt entries, and a dibasic-site residue needs the
surrounding domain context that a 36-residue window would strip.  Slicing to
windows happens later, on the R side, by (accession, index).

Sequences longer than ``max_len`` are embedded in overlapping chunks and the
overlap is blended linearly, so no residue's vector comes from a chunk edge.
"""

from __future__ import annotations

import re
import time
from pathlib import Path
from typing import Callable, Iterable

import h5py
import numpy as np
import pandas as pd

PROT_T5 = "Rostlab/prot_t5_xl_half_uniref50-enc"   # encoder only, fp16 weights, 1024-d
ESM_C = "biohub/ESMC-600M"                          # 1152-d; ungated on HF, `pip install esm`


def pick_device(device: str | None = None) -> str:
    import torch

    if device and device != "auto":
        return device
    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def _prot_t5(device: str) -> Callable[[str], np.ndarray]:
    import torch
    from transformers import T5EncoderModel, T5Tokenizer

    tok = T5Tokenizer.from_pretrained(PROT_T5, do_lower_case=False)
    # fp16 on a GPU/MPS, fp32 on CPU (fp16 matmuls are not implemented there)
    dtype = torch.float32 if device == "cpu" else torch.float16
    mdl = T5EncoderModel.from_pretrained(PROT_T5, torch_dtype=dtype).to(device).eval()

    def run(seq: str) -> np.ndarray:
        # ProtT5 wants space-separated residues; rare letters map to X
        spaced = " ".join(re.sub(r"[UZOB]", "X", seq))
        ids = tok(spaced, return_tensors="pt", add_special_tokens=True).to(device)
        with torch.no_grad():
            h = mdl(**ids).last_hidden_state[0]
        return h[: len(seq)].float().cpu().numpy()         # drop the trailing </s>

    return run


def _esm_c(device: str) -> Callable[[str], np.ndarray]:
    """ESM C 600M through the esm package's HF-style API (esm >= 3.4)."""
    import torch
    from esm.models.esmc import EsmcForMaskedLM
    from esm.tokenization import EsmSequenceTokenizer

    tok = EsmSequenceTokenizer()
    mdl = EsmcForMaskedLM.from_pretrained(ESM_C).to(device).eval()
    if device != "cpu":
        mdl = mdl.to(torch.bfloat16)      # bf16 halves memory; ESM C was trained in bf16

    def run(seq: str) -> np.ndarray:
        ids = torch.tensor([tok.encode(seq)], device=device)       # BOS + residues + EOS
        with torch.no_grad():
            out = mdl(input_ids=ids, return_dict=True)
        h = out.last_hidden_state[0, 1:-1]                          # drop BOS / EOS
        if h.shape[0] != len(seq):
            raise RuntimeError(f"ESM C returned {h.shape[0]} tokens for {len(seq)} residues")
        return h.float().cpu().numpy()

    return run


BACKENDS = {"prot_t5": _prot_t5, "esm_c": _esm_c}


def embed_chunked(run: Callable[[str], np.ndarray], seq: str,
                  max_len: int = 2000, overlap: int = 200) -> np.ndarray:
    """Embed ``seq`` in windows of ``max_len`` with ``overlap``, blending the
    overlaps with a linear ramp so every residue sees context on both sides."""
    n = len(seq)
    if n <= max_len:
        return run(seq)
    step = max_len - overlap
    starts = list(range(0, max(n - max_len, 0) + 1, step))
    if starts[-1] + max_len < n:
        starts.append(n - max_len)
    out = None
    weight = None
    for s in starts:
        e = min(s + max_len, n)
        h = run(seq[s:e]).astype("float32")
        if out is None:
            out = np.zeros((n, h.shape[1]), dtype="float32")
            weight = np.zeros(n, dtype="float32")
        w = np.ones(e - s, dtype="float32")
        if s > 0:
            w[:overlap] = np.linspace(0.0, 1.0, overlap, endpoint=False)
        if e < n:
            w[-overlap:] = np.linspace(1.0, 0.0, overlap, endpoint=False)
        out[s:e] += h * w[:, None]
        weight[s:e] += w
    return out / np.maximum(weight, 1e-6)[:, None]


def read_sequences(path) -> pd.DataFrame:
    """``accession`` + ``sequence_uni`` (parquet or csv), one row per accession."""
    p = Path(path).expanduser()
    df = pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p)
    need = {"accession", "sequence_uni"}
    if not need <= set(df.columns):
        raise ValueError(f"{p} needs columns {sorted(need)}, has {list(df.columns)}")
    df = df.dropna(subset=["sequence_uni"]).drop_duplicates("accession")
    return df.reset_index(drop=True)


def embed_sequences(sequences, out_h5, backend: str = "prot_t5", device: str | None = None,
                    max_len: int = 2000, overlap: int = 200, verbose: int = 1) -> Path:
    """Embed every accession in ``sequences`` into ``out_h5`` (resumable).

    ``sequences`` is a path or a DataFrame with ``accession`` and
    ``sequence_uni``.  Accessions already in the file are skipped.  Sequences are
    processed longest first so a memory problem shows up in the first minute,
    not the last.
    """
    df = read_sequences(sequences) if not isinstance(sequences, pd.DataFrame) else sequences
    df = df.assign(_len=df["sequence_uni"].str.len()).sort_values("_len", ascending=False)
    device = pick_device(device)
    out = Path(out_h5).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(out, "a") as h5:
        done = set(h5.keys())
        todo = df[~df["accession"].isin(done)]
        if verbose:
            print(f"[{backend}] {len(todo)} to embed, {len(done)} already in {out.name}; "
                  f"device={device}, longest {int(todo['_len'].max()) if len(todo) else 0} aa",
                  flush=True)
        if not len(todo):
            return out
        run = BACKENDS[backend](device)
        t0 = time.time()
        n_res = 0
        for i, (acc, seq) in enumerate(zip(todo["accession"], todo["sequence_uni"]), 1):
            emb = embed_chunked(run, seq, max_len=max_len, overlap=overlap)
            if emb.shape[0] != len(seq):
                raise RuntimeError(f"{acc}: embedded {emb.shape[0]} residues for a {len(seq)}-aa sequence")
            ds = h5.create_dataset(acc, data=emb.astype("float16"), compression="gzip", compression_opts=4)
            ds.attrs["backend"] = backend
            n_res += len(seq)
            if verbose and (i % 50 == 0 or i == len(todo)):
                dt = time.time() - t0
                print(f"  {i}/{len(todo)}  {n_res} residues  {dt / 60:.1f} min  "
                      f"({n_res / max(dt, 1e-9):.0f} res/s)", flush=True)
            if i % 50 == 0:
                h5.flush()
    return out


def iter_embeddings(h5_path) -> Iterable[tuple[str, np.ndarray]]:
    with h5py.File(Path(h5_path).expanduser(), "r") as h5:
        for acc in h5.keys():
            yield acc, h5[acc][...]
