"""Where the ligandFinder package's own python, R code and data live --
resolved ONCE, here, instead of once per script.

The python twin of inst/scripts/lf_paths.R, and the same reasoning: these
scripts import ``lf_winxgb`` and ``lf_dcnn`` from a ligandFinder SOURCE TREE,
not from the installed package. That is deliberate. The install lags badly
(built 2026-09-10) and its ``inst/python`` carries neither ``lf_winxgb`` nor
``lf_resid`` nor ``lf_dcnn/scan.py`` at all, so there is nothing there to
import; what used to come from a throwaway git worktree now comes from here.

The source tree is found as:

  1. ``$LF_SRC``, if set (validated, so a typo fails loudly)
  2. the checkout this file sits in -- walking up from ``__file__`` for a
     DESCRIPTION naming ligandFinder beside an ``R/`` directory. These scripts
     are committed under ``inst/scripts/``, so a copy that travels with the
     repository needs no configuration at all.

A copy kept OUTSIDE a checkout (the working copy in ~/AF2_analysis, say)
therefore needs ``LF_SRC``, and says so when it does not have it::

    export LF_SRC=~/R_projects/ligandFinder

``__file__`` makes this strictly more reliable than the R side, which has to
dig the script's path out of ``commandArgs()``.

Using it is a two-line preamble, the same in every script, and it has to come
BEFORE the first ``lf_winxgb`` / ``lf_dcnn`` import::

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import lf_paths; lf_paths.add_python()

The first line is what makes ``import lf_paths`` work. ``python some_script.py``
already puts the script's own directory on ``sys.path``, so that line is
usually redundant -- but ``python -m``, ``-I`` and ``-P`` do not, and being
explicit costs one line and removes the question.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

__all__ = ["lf_src_root", "lf_src", "add_python"]

_root: Path | None = None


def _is_src(d: Path) -> bool:
    f = d / "DESCRIPTION"
    if not f.is_file() or not (d / "R").is_dir():
        return False
    try:
        with open(f, encoding="utf-8", errors="replace") as fh:
            return any(ln.strip() == "Package: ligandFinder" or
                       ln.rstrip().replace(" ", "") == "Package:ligandFinder"
                       for ln in fh)
    except OSError:
        return False


def lf_src_root() -> Path:
    """The top of the ligandFinder source tree, resolved once and cached."""
    global _root
    if _root is not None:
        return _root

    env = os.environ.get("LF_SRC", "")
    if env:
        d = Path(os.path.abspath(os.path.expanduser(env)))
        if not _is_src(d):
            raise RuntimeError(
                f"LF_SRC is not a ligandFinder source tree (no DESCRIPTION / R/): {env}")
        _root = d
        return _root

    here = Path(os.path.abspath(__file__)).parent
    for d in (here, *here.parents):
        if _is_src(d):
            _root = d
            return _root
    raise RuntimeError(
        f"no ligandFinder source tree at or above {here}\n"
        "  these scripts import lf_winxgb and lf_dcnn from a checkout, not from the\n"
        "  installed package, which has neither. Run them from inst/scripts/ in a\n"
        "  checkout, or set LF_SRC:\n"
        "      export LF_SRC=~/R_projects/ligandFinder")


def lf_src(*parts: str) -> Path:
    """One file or directory inside the source tree, checked -- a missing path
    here means the tree is not the one the script expects, worth stopping for."""
    p = lf_src_root().joinpath(*parts)
    if not os.path.exists(p):
        raise FileNotFoundError(
            f"not in the ligandFinder source tree {lf_src_root()}: {os.path.join(*parts)}")
    return p


def add_python() -> str:
    """Put the source tree's ``inst/python`` first on ``sys.path``; idempotent."""
    p = str(lf_src("inst", "python"))
    while p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)
    return p
