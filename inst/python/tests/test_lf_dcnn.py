"""pytest shim over the built-in checks.

The target venv has no pytest, so the checks themselves live in
``lf_dcnn.selftest`` and run via ``python -m lf_dcnn selftest``.  This module
exposes them to pytest wherever it happens to be available.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from lf_dcnn import selftest  # noqa: E402


@pytest.mark.parametrize("check", selftest.CHECKS, ids=lambda f: f.__name__)
def test_check(check):
    check()


def test_end_to_end_and_roundtrip():
    reference = selftest.check_end_to_end()
    selftest.check_roundtrip_npz(reference)
