"""Standalone command line entry point -- runnable with no R in the loop.

    python -m lf_dcnn train --input-dir DIR --output-dir DIR
    python -m lf_dcnn train --input-dir DIR --output-dir DIR --n-seeds 5 --set position_ramp=false
    python -m lf_dcnn train --input-dir DIR --output-dir DIR --n-seeds 5 --isolated
    python -m lf_dcnn combine --input-dir DIR --output-dir DIR
    python -m lf_dcnn selftest

``--isolated`` trains every (member, terminus) model in its own python process
(they land under ``<output-dir>/members/``) and then combines them; a model
already on disk is not retrained, so an interrupted run resumes. The second
model trained in ONE process can die mid-fit in a retraced tf.function, which
is what this avoids.
"""

from __future__ import annotations

import argparse
import subprocess
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
    p.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="FIELD=JSON",
        help="override any Config field, value parsed as JSON: "
        "--set position_ramp=false --set conv_filters=[16,8]",
    )
    p.add_argument("--verbose", type=int, default=1)


def _apply(cfg: Config, args) -> Config:
    import json

    over = {}
    if args.epochs is not None:
        over["epochs"] = args.epochs
    if args.seed is not None:
        over["seed"] = args.seed
    if getattr(args, "r_exact", False):
        over["db_mask_both_termini"] = False
        over["resample_each_epoch"] = False
    for item in getattr(args, "set", []):
        field, _, raw = item.partition("=")
        if not _ or not field:
            raise SystemExit(f"--set expects FIELD=JSON, got {item!r}")
        try:
            over[field] = json.loads(raw)
        except json.JSONDecodeError:
            over[field] = raw          # bare strings: --set trunk=unet
    return cfg.evolve(**over) if over else cfg


def _report(out, result) -> None:
    print(f"wrote {out}")
    print(f"  windows scored : {len(result.pred)}")
    print(f"  members        : {result.n_seeds}")
    print(f"  calibrator     : {result.calibrator}")


def cmd_train(args) -> int:
    from .io import load_inputs, member_paths, save_member, save_outputs
    from .pipeline import run, run_member

    if args.member is not None:
        # one member (one terminus with --term), this process: write it and stop
        data, cfg, _ = load_inputs(args.input_dir)
        cfg = _apply(cfg, args)
        m = run_member(data, cfg, k=args.member, n_seeds=args.n_seeds,
                       verbose=args.verbose, keep_models=(args.member == 0),
                       terms=[args.term] if args.term else None)
        for f in save_member(args.output_dir, m):
            print(f"wrote {f}")
        # member 0 is the one combine() keeps for embeddings/models; its
        # weights let a caller rebuild live models (R's lf_dcnn_run does)
        for term, model in m.models.items():
            w = member_paths(args.output_dir, m.k, term)[0].with_suffix(".weights.h5")
            model.save_weights(w)
            print(f"wrote {w}")
        return 0

    if args.isolated:
        # one subprocess per (member, terminus) -- one model per process --
        # then combine. The subprocess argv is this call's argv with
        # --isolated swapped for --member K --term T, so every override
        # (--set, --epochs, --seed, ...) reaches it unchanged.
        _, cfg, _ = load_inputs(args.input_dir)
        cfg = _apply(cfg, args)
        terms = [t for t in cfg.term_order]
        base = [a for a in sys.argv[1:] if a != "--isolated"]
        n = args.n_seeds * len(terms)
        i = 0
        for k in range(args.n_seeds):
            for term in terms:
                i += 1
                npz, _ = member_paths(args.output_dir, k, term)
                if npz.exists():
                    print(f"model {i}/{n} (member {k + 1}, {term}): reusing {npz}")
                    continue
                print(f"model {i}/{n} (member {k + 1}, {term}): training in a subprocess", flush=True)
                rc = subprocess.call([sys.executable, "-u", "-m", "lf_dcnn", *base,
                                      "--member", str(k), "--term", term])
                if rc != 0 or not npz.exists():
                    raise SystemExit(f"member {k} {term} failed (rc={rc}); "
                                     f"rerun to resume from {npz.parent}")
        return cmd_combine(args)

    data, cfg, meta = load_inputs(args.input_dir)
    cfg = _apply(cfg, args)
    result = run(data, cfg, verbose=args.verbose, n_seeds=args.n_seeds)
    _report(save_outputs(args.output_dir, result, meta=meta), result)
    return 0


def cmd_combine(args) -> int:
    """Combine the members under ``<output-dir>/members/`` into the outputs."""
    from .io import load_inputs, load_members, save_outputs
    from .pipeline import combine

    _, cfg, meta = load_inputs(args.input_dir)
    cfg = _apply(cfg, args)
    members = load_members(args.output_dir, term_order=cfg.term_order)
    if not members:
        raise SystemExit(f"no members under {args.output_dir}/members")
    want = getattr(args, "n_seeds", None)
    if want and len(members) != want:
        raise SystemExit(f"{len(members)} members on disk but --n-seeds {want}")
    result = combine(members, cfg)
    _report(save_outputs(args.output_dir, result, meta=meta), result)
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
    t.add_argument("--n-seeds", type=int, default=1,
                   help="ensemble members per terminus (seeds seed..seed+n-1)")
    t.add_argument("--isolated", action="store_true",
                   help="train each member in its own process, then combine")
    t.add_argument("--member", type=int, default=None, metavar="K",
                   help="train only member K (0-based) and write it under "
                        "<output-dir>/members/; used by --isolated")
    t.add_argument("--term", default=None, metavar="T",
                   help="with --member: train only this terminus (N or C)")
    _add_common(t)
    t.set_defaults(func=cmd_train)

    c = sub.add_parser("combine", help="combine on-disk members into the ensemble outputs")
    c.add_argument("--input-dir", required=True, help="for config.json and meta.parquet")
    c.add_argument("--output-dir", required=True, help="holds members/, receives the outputs")
    c.add_argument("--n-seeds", type=int, default=None, help="expected member count (optional check)")
    _add_common(c)
    c.set_defaults(func=cmd_combine)

    s = sub.add_parser("selftest", help="run the built-in checks on synthetic data")
    s.add_argument("--verbose", type=int, default=1)
    s.set_defaults(func=cmd_selftest)

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
