"""Per-residue conservation from each precursor's MSA (Kevin, 2026-11-10).

Replaces the Aminode evolutionary-rate channel (`cons_rs_n`). The score is the
Jensen-Shannon divergence of Capra & Singh (2007, Bioinformatics 23:1875):

  * sequences are weighted by Henikoff & Henikoff position-based weights, so a
    family that is over-represented in UniRef/env does not dominate;
  * a column's weighted amino-acid distribution is compared with the BLOSUM62
    background distribution (JSD, log base 2, so the score is in [0, 1]);
  * the score is multiplied by the column's weighted fraction of NON-gap
    sequences, so a column only a handful of homologues cover cannot look
    conserved. This is what keeps signal peptides and disordered linkers, which
    many hits do not reach, from scoring as conserved on thin evidence.

No window smoothing (Capra & Singh's default lambda=0.5, window 3 is off): the
trunk already convolves over neighbours, and smoothing would leak a residue's
neighbours into the one channel meant to describe that residue.

The MSA is read as A3M: lower-case letters are insertions relative to the query
and are deleted, so every row has exactly one column per query residue. The
query row MUST equal the precursor's sequence in the residue cache; a mismatch
is a hard error, because the columns would then describe other residues.

Input: one `<accession>.a3m[.gz]` per precursor in any of --msa-dirs (first hit
wins, so a better source can be listed first and shadow the others; AFDB's own
MSAs, if they come back, would go first).

Output (--out, parquet), one row per residue of every precursor that has an MSA:
  accession, index (1-based), AA, jsd, gap_frac (weighted), n_cov (raw count of
  sequences with a residue in the column)
and --per-protein (tsv): accession, length, depth, neff (exp of the entropy of
the Henikoff weights, which sum to 1), source.

Evidence shrinkage (--prior, default 10): each column's distribution is mixed
with `prior` pseudo-sequences of BLOSUM62 background before the JSD, against an
effective count of neff x (1 - gap fraction). Without it a query-only MSA scores
~0.85 everywhere; with it a protein with no homologues reads as unconserved.

    ~/.virtualenvs/lf-plm/bin/python msa_conservation.py --workers 6
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

A2 = os.path.expanduser("~/AF2_analysis")

AA20 = "ARNDCQEGHILKMFPSTWYV"
GAP = 20                                   # code for '-' and anything not in AA20
_LUT = np.full(256, GAP, dtype=np.int8)
for _i, _c in enumerate(AA20):
    _LUT[ord(_c)] = _i

# BLOSUM62 background amino-acid frequencies, in AA20 order (the distribution
# Capra & Singh's score_conservation.py uses).
BLOSUM62_BG = np.array([
    0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062,
    0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072])
BLOSUM62_BG = BLOSUM62_BG / BLOSUM62_BG.sum()
PSEUDO = 1e-6


def read_a3m(path: str) -> list[str]:
    """Rows of an A3M, insertions deleted. The query is row 0 and appears ONCE:
    ColabFold's combined files (the LiteFold set) repeat the query record at the
    head of the env-DB block; that copy is dropped so both sources weight the
    query the same (fetch_colabfold_api.py drops it when concatenating)."""
    op = gzip.open if path.endswith(".gz") else open
    seqs, cur, hdr, qhdr = [], [], None, None

    def flush():
        if hdr is not None and cur:
            if qhdr is None or hdr != qhdr:
                seqs.append("".join(cur))

    with op(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n").replace("\x00", "")
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                flush()
                if qhdr is None and hdr is not None:
                    qhdr = hdr
                hdr, cur = line.split("\t")[0], []
            else:
                cur.append(line)
    flush()
    # delete insertions (lower case) and the a3m '.' padding
    return ["".join(ch for ch in s if not ch.islower() and ch != ".") for s in seqs]


def encode(seqs: list[str], L: int) -> np.ndarray:
    keep = [s for s in seqs if len(s) == L]
    buf = np.frombuffer("".join(keep).encode("ascii"), dtype=np.uint8)
    return _LUT[buf].reshape(len(keep), L), len(seqs) - len(keep)


def henikoff_weights(m: np.ndarray) -> np.ndarray:
    """Position-based weights (gap as a 21st symbol), normalised to sum 1."""
    n, L = m.shape
    w = np.zeros(n)
    for c in range(L):
        col = m[:, c]
        counts = np.bincount(col, minlength=21)
        k = np.count_nonzero(counts)
        w += 1.0 / (k * counts[col])
    w /= L
    return w / w.sum()


def jsd_columns(m: np.ndarray, w: np.ndarray, neff: float, prior: float):
    n, L = m.shape
    freq = np.zeros((L, 21))
    for c in range(L):             # per column: an (n x L) float temporary is ~1 GB on deep MSAs
        freq[c] = np.bincount(m[:, c], weights=w, minlength=21)
    gap_frac = freq[:, GAP]
    aa = freq[:, :20] + PSEUDO
    p = aa / aa.sum(axis=1, keepdims=True)
    q = BLOSUM62_BG[None, :]
    # evidence shrinkage: a Dirichlet prior of `prior` pseudo-sequences at the
    # background, against the column's effective count neff * (1 - gap_frac).
    # A query-only MSA otherwise scores ~0.85 (its own one-hot vs BLOSUM62) at
    # EVERY residue -- "no homologues" read as "maximally conserved"; with the
    # prior it scores near 0. Deep MSAs (neff in the hundreds) are unchanged.
    n_c = (neff * (1.0 - gap_frac))[:, None]
    p = (n_c * p + prior * q) / (n_c + prior)
    r = 0.5 * (p + q)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = 0.5 * np.sum(p * np.log2(p / r), axis=1) + 0.5 * np.sum(q * np.log2(q / r), axis=1)
    n_cov = (m != GAP).sum(axis=0)
    return d * (1.0 - gap_frac), gap_frac, n_cov


def score_one(acc: str, seq: str, path: str, max_depth: int, prior: float):
    seqs = read_a3m(path)
    # ColabFold writes selenocysteine as X; the cache keeps U (SELENOP, GPX3, GPX6,
    # DIO3). Both encode to the gap code, so the scores are the same either way.
    if not seqs or seqs[0] != seq.replace("U", "X"):
        q = seqs[0] if seqs else ""
        raise ValueError(f"{acc}: MSA query ({len(q)} aa) != cache sequence ({len(seq)} aa)")
    m, n_bad = encode(seqs, len(seq))
    if n_bad:
        raise ValueError(f"{acc}: {n_bad} MSA rows are not {len(seq)} columns after deleting insertions")
    depth = m.shape[0]
    if depth > max_depth:          # keep the query + the first hits (sorted by e-value per DB)
        m = m[:max_depth]
    w = henikoff_weights(m)
    neff = float(np.exp(-(w * np.log(w)).sum()))
    jsd, gap, ncov = jsd_columns(m, w, neff, prior)
    return acc, jsd.astype("float32"), gap.astype("float32"), ncov.astype("int32"), depth, neff


def main(argv=None):
    ap = argparse.ArgumentParser(prog="msa_conservation")
    ap.add_argument("--seqs", default=f"{A2}/lf_plm/precursor_seqs.tsv",
                    help="accession/gene/seq of every precursor in lf_pepend_residue_cache.rds")
    ap.add_argument("--msa-dirs", nargs="+",
                    default=[f"{A2}/lf_msa/hf_colabfold", f"{A2}/lf_msa/api_colabfold"])
    ap.add_argument("--out", default=f"{A2}/lf_msa/msa_cons_raw.parquet")
    ap.add_argument("--per-protein", default=f"{A2}/lf_msa/msa_cons_per_protein.tsv")
    ap.add_argument("--max-depth", type=int, default=50000)
    ap.add_argument("--alias", default=f"{A2}/lf_msa/hf_needed_alias.tsv",
                    help="tsv acc<TAB>sequence for MSA files named by non-cache accessions")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--prior", type=float, default=10.0,
                    help="pseudo-sequences of BLOSUM62 background mixed into every column")
    ap.add_argument("--allow-missing", action="store_true",
                    help="write even if some precursors have no MSA (they become 0 downstream)")
    a = ap.parse_args(argv)

    import pyarrow as pa
    import pyarrow.parquet as pq

    with open(a.seqs) as f:
        prec = [(r["accession"], r["seq"]) for r in csv.DictReader(f, delimiter="\t")]
    # an MSA is found under the precursor's own accession or, failing that, under
    # any accession whose sequence is IDENTICAL (the fetch lists keep one file per
    # unique sequence, e.g. AMY1A/AMY1C share P0DTE7's). --alias maps a file's
    # accession to its sequence, for files named by accessions outside the cache
    # (e.g. P0DY56, the HF set's IFNA1, which the cache holds as P01562).
    def find(acc):
        for d in a.msa_dirs:
            for ext in (".a3m.gz", ".a3m"):
                p = os.path.join(d, acc + ext)
                if os.path.exists(p) and os.path.getsize(p) > 0:
                    return p, os.path.basename(d)
        return None
    accs_of_seq = {}
    for acc, seq in prec:
        accs_of_seq.setdefault(seq, []).append(acc)
    if a.alias and os.path.exists(a.alias):
        with open(a.alias) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                accs_of_seq.setdefault(r["sequence"], []).append(r["acc"])
    jobs, missing = [], []
    for acc, seq in prec:
        hit = find(acc)
        if hit is None:
            for other in accs_of_seq.get(seq, []):
                if other != acc and (hit := find(other)):
                    hit = (hit[0], hit[1] + f" (via {other})"); break
        (jobs.append((acc, seq, *hit)) if hit else missing.append(acc))
    print(f"{len(prec):,} precursors; MSA for {len(jobs):,}; none for {len(missing):,}", flush=True)
    if missing and not a.allow_missing:
        raise SystemExit(f"no MSA for {len(missing)} precursors, e.g. {missing[:10]}; "
                         "finish the downloads or pass --allow-missing")

    seq_of = dict(prec)
    src_of = {j[0]: j[3] for j in jobs}
    cols = {k: [] for k in ("accession", "index", "AA", "jsd", "gap_frac", "n_cov")}
    per = []
    errors = []
    with ProcessPoolExecutor(a.workers) as ex:
        futs = {ex.submit(score_one, acc, seq, p, a.max_depth, a.prior): acc for acc, seq, p, _ in jobs}
        for i, fu in enumerate(as_completed(futs), 1):
            try:
                acc, jsd, gap, ncov, depth, neff = fu.result()
            except Exception as e:      # collected, then fatal: never write a partial table silently
                errors.append(str(e)); continue
            L = jsd.size
            cols["accession"] += [acc] * L
            cols["index"].append(np.arange(1, L + 1, dtype="int32"))
            cols["AA"] += list(seq_of[acc])
            cols["jsd"].append(jsd); cols["gap_frac"].append(gap); cols["n_cov"].append(ncov)
            per.append((acc, L, depth, neff, src_of[acc]))
            if i % 250 == 0:
                print(f"  {i:,}/{len(jobs):,}", flush=True)
    if errors:
        print("\n".join(errors[:20]))
        raise SystemExit(f"{len(errors)} precursors failed; nothing written")

    t = pa.table({
        "accession": cols["accession"],
        "index": np.concatenate(cols["index"]),
        "AA": cols["AA"],
        "jsd": np.concatenate(cols["jsd"]),
        "gap_frac": np.concatenate(cols["gap_frac"]),
        "n_cov": np.concatenate(cols["n_cov"]),
    })
    pq.write_table(t, a.out, compression="zstd")
    with open(a.per_protein, "w") as f:
        f.write("accession\tlength\tdepth\tneff\tsource\n")
        for r in sorted(per):
            f.write("\t".join(map(str, r)) + "\n")
        for acc in missing:
            f.write(f"{acc}\t{len(seq_of[acc])}\t0\t0\tnone\n")
    j = t.column("jsd").to_numpy()
    print(f"-> {a.out}: {t.num_rows:,} residues; jsd min {j.min():.3f} "
          f"median {np.median(j):.3f} max {j.max():.3f}", flush=True)
    if missing:
        print(f"   no MSA: {', '.join(missing[:20])}{' ...' if len(missing) > 20 else ''}")


if __name__ == "__main__":
    main()
