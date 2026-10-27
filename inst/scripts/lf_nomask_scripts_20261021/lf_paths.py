"""A forwarder. The real thing is one directory up, at inst/scripts/lf_paths.py,
and this is its python counterpart of lf_paths.R beside it: the scripts here
use the same two-line preamble as the ones up there,

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import lf_paths; lf_paths.add_python()

and land on this file, which hands straight over -- ``import lf_paths`` returns
the real module, so ``lf_src_root`` / ``lf_src`` / ``add_python`` are its own.

Two ways to find the real file, matching lf_paths.R's forwarder: this directory
kept inside a checkout, or copied out of one (the working copy in
~/AF2_analysis, say) with LF_SRC set. The second is why this cannot just put
``..`` on sys.path: a copied directory has no parent inst/scripts.
"""
import importlib.util
import os
import sys


def _real() -> str:
    here = os.path.dirname(os.path.abspath(__file__))
    cand = os.path.join(here, "..", "lf_paths.py")
    if not os.path.isfile(cand):
        e = os.environ.get("LF_SRC", "")
        cand = os.path.join(os.path.expanduser(e), "inst", "scripts", "lf_paths.py") if e else ""
    if not cand or not os.path.isfile(cand):
        raise RuntimeError(
            f"cannot find inst/scripts/lf_paths.py from {here}\n"
            "  this file only forwards to it. Run these scripts from a checkout,\n"
            "  or point LF_SRC at one:    export LF_SRC=~/R_projects/ligandFinder")
    return os.path.realpath(cand)


_spec = importlib.util.spec_from_file_location(__name__, _real())
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
sys.modules[__name__] = _mod
