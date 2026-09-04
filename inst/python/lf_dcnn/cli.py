"""Standalone command line entry point -- runnable with no R in the loop.

    python -m lf_dcnn train --input-dir DIR --output-dir DIR
    python -m lf_dcnn selftest
"""

from __future__ import annotations

import argparse
import sys

from .config import Config


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--epochs", type=int, help="override Config.epochs")
    p.add_argument("--seed", type=int, help="override Config.seed")
    p.add_argument(
        "--r-exact",
        action="store_true",
        help="reproduce 10_1dcnn_new6.R exactly, including its DB-mask overwrite "
        "and its non-functional per-epoch resample",
    )
    p.add_argument("--verbose", type=int, default=1)


def _apply(cfg: Config, args) -> Config:
    over = {}
    if args.epochs is not None:
        over["epochs"] = args.epochs
    if args.seed is not None:
        over["seed"] = args.seed
    if getattr(args, "r_exact", False):
        over["db_mask_both_termini"] = False
        over["resample_each_epoch"] = False
    return cfg.evolve(**over) if over else cfg


def cmd_train(args) -> int:
    from .io import load_inputs, save_outputs
    from .pipeline import run

    data, cfg, meta = load_inputs(args.input_dir)
    cfg = _apply(cfg, args)
    result = run(data, cfg, verbose=args.verbose)
    out = save_outputs(args.output_dir, result, meta=meta)
    print(f"wrote {out}")
    print(f"  windows scored : {len(result.pred)}")
    print(f"  calibrator     : {result.calibrator}")
    return 0


def cmd_selftest(args) -> int:
    from .selftest import main as selftest_main

    return selftest_main(verbose=args.verbose)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="lf_dcnn", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)

    t = sub.add_parser("train", help="train both models from a .npz/parquet directory")
    t.add_argument("--input-dir", required=True)
    t.add_argument("--output-dir", required=True)
    _add_common(t)
    t.set_defaults(func=cmd_train)

    s = sub.add_parser("selftest", help="run the built-in checks on synthetic data")
    s.add_argument("--verbose", type=int, default=1)
    s.set_defaults(func=cmd_selftest)

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
