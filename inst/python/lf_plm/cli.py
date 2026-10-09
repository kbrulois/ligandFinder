"""Command line entry point.

    python -m lf_plm embed  --sequences seqs.parquet --out plm_raw.h5 [--backend prot_t5] [--device auto]
    python -m lf_plm reduce --h5 plm_raw.h5 --sequences seqs.parquet --k 32 --out plm_residues.parquet
"""

from __future__ import annotations

import argparse
import sys


def cmd_embed(args) -> int:
    from .embed import embed_sequences

    out = embed_sequences(args.sequences, args.out, backend=args.backend, device=args.device,
                          max_len=args.max_len, overlap=args.overlap, verbose=args.verbose)
    print(f"wrote {out}")
    return 0


def cmd_reduce(args) -> int:
    from pathlib import Path

    from .embed import read_sequences
    from .reduce import fit_pca, load_pca, save_pca, transform_to_table

    pca_path = Path(args.pca or (Path(args.out).expanduser().with_suffix("") .as_posix() + f"_pca{args.k}.joblib"))
    if pca_path.exists() and not args.refit:
        pca, sd, _ = load_pca(pca_path)
        print(f"reusing {pca_path}")
    else:
        pca, sd = fit_pca(args.h5, k=args.k, n_residues=args.n_residues, seed=args.seed,
                          verbose=args.verbose)
        save_pca(pca_path, pca, sd, meta={"h5": str(args.h5), "k": args.k,
                                          "n_residues": args.n_residues, "seed": args.seed})
        print(f"wrote {pca_path}")
    transform_to_table(args.h5, read_sequences(args.sequences), pca, sd, args.out,
                       clip=args.clip, verbose=args.verbose)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="lf_plm", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)

    e = sub.add_parser("embed", help="per-residue embeddings -> h5, one dataset per accession")
    e.add_argument("--sequences", required=True, help="parquet/csv with accession, sequence_uni")
    e.add_argument("--out", required=True, help="output .h5 (resumable)")
    e.add_argument("--backend", default="prot_t5", choices=["prot_t5", "esm_c"])
    e.add_argument("--device", default="auto", help="auto | cuda | mps | cpu")
    e.add_argument("--max-len", type=int, default=2000, help="chunk longer sequences")
    e.add_argument("--overlap", type=int, default=200)
    e.add_argument("--verbose", type=int, default=1)
    e.set_defaults(func=cmd_embed)

    r = sub.add_parser("reduce", help="PCA to k channels -> long residue table")
    r.add_argument("--h5", required=True)
    r.add_argument("--sequences", required=True, help="the same table embed ran on")
    r.add_argument("--k", type=int, default=32)
    r.add_argument("--out", required=True, help="output .parquet")
    r.add_argument("--pca", default=None, help="fitted PCA (.joblib); default beside --out")
    r.add_argument("--refit", action="store_true")
    r.add_argument("--n-residues", type=int, default=500_000, help="residue sample for the fit")
    r.add_argument("--clip", type=float, default=3.0, help="+/- sd clip before [0,1] scaling")
    r.add_argument("--seed", type=int, default=1)
    r.add_argument("--verbose", type=int, default=1)
    r.set_defaults(func=cmd_reduce)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
