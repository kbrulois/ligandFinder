#!/usr/bin/env python3
"""Teach AlphaPulldown's run_multimer_jobs.py to forward initial-guess flags.

run_multimer_jobs.py is a wrapper: it reads --protein_lists, turns each line
into an --input string, assembles a fixed argument list, and shells out to
run_structure_prediction.py.  The list is hardcoded, so --initial_guess_dir /
--initial_guess_map cannot reach the inner script and a screen cannot be
seeded through the wrapper.

This adds the two flags to the wrapper and appends them to the command it
builds, when set.  run_structure_prediction.py must already be patched to
understand them (apply_initial_guess_patches_v2.py).

    python patch_run_multimer_jobs_initial_guess.py \
        --wrapper $ENV/lib/python3.11/site-packages/alphapulldown/scripts/run_multimer_jobs.py

Idempotent.  Writes <wrapper>.pre_initial_guess once; --revert restores it.
Nothing is written unless every anchor is found, so a failed match leaves the
file untouched.

Note the wrapper is often installed twice -- under scripts/ and again under
bin/ -- and `run_multimer_jobs.py` on PATH is usually the bin/ copy.  Patch
whichever one you actually invoke, or both.  --wrapper accepts several paths.
"""
import argparse
import os
import re
import shutil
import sys

FLAG_BLOCK = '''
# Initial guess settings (added by patcher).  Guarded because absl flags are
# process-global: if anything in this process has already imported
# run_structure_prediction, these names exist and redefining them raises
# DuplicateFlagError.
if 'initial_guess_dir' not in flags.FLAGS:
    flags.DEFINE_string(
        'initial_guess_dir', None,
        'Directory of PDB files for AlphaFold initial guess. Each PDB '
        'filename stem must match the job description AlphaPulldown builds '
        '(e.g. P21452_and_P20366_72-107.pdb).')
if 'initial_guess_map' not in flags.FLAGS:
    flags.DEFINE_string(
        'initial_guess_map', None,
        'TSV mapping job descriptions to initial-guess PDB paths. '
        'Two columns: description<TAB>pdb_path.')

'''

# Injected as the first entries of `constant_args`.  That dict is the single
# source for `command_args`, which every subprocess path consumes -- the
# wrapper has at least two (a batched one when --desired_num_res and
# --desired_num_msa are both set, and a per-fold loop otherwise), and patching
# the call sites individually silently misses whichever one is not taken.
# Entries whose value is None are dropped when command_args is built, so an
# unset flag adds nothing to the command.
ARGS_TMPL = '''{i}# Added by patcher: forwarded to run_structure_prediction via command_args.
{i}"--initial_guess_dir": FLAGS.initial_guess_dir,
{i}"--initial_guess_map": FLAGS.initial_guess_map,
'''

CONST_ARGS_RE = re.compile(r'^(?P<indent>[ \t]*)constant_args\s*=\s*\{\s*$')


def patch_one(path):
    print(f"\n=== {path}")
    if not os.path.isfile(path):
        print("  [ERROR] not a file")
        return False

    with open(path) as f:
        src = f.read()

    if 'initial_guess_dir' in src:
        print("  [skip] already patched")
        return True

    lines = src.splitlines(keepends=True)

    # --- anchor 1: the constant_args dict that feeds every subprocess path
    hit = None
    for idx, line in enumerate(lines):
        m = CONST_ARGS_RE.match(line)
        if m:
            hit = (idx, m.group('indent'))
            break
    if hit is None:
        print("  [ERROR] no `constant_args = {` line found.")
        print("          Lines mentioning constant_args or subprocess:")
        for idx, line in enumerate(lines):
            if 'constant_args' in line or 'subprocess' in line:
                print(f"            {idx + 1}: {line.rstrip()}")
        print("          Not modified. Send those lines and the patcher can be adjusted.")
        return False
    const_idx, indent = hit
    entry_indent = indent + '    '
    n_calls = sum(1 for l in lines if 'subprocess.run' in l)
    print(f"  constant_args : line {const_idx + 1}  (covers {n_calls} subprocess call site(s))")
    insert_idx = const_idx + 1

    # --- anchor 2: module-level insertion point for the flag definitions.
    # Before `def main(` keeps them at column 0 and after every existing
    # DEFINE, without having to parse multi-line DEFINE statements.
    main_idx = None
    for idx, line in enumerate(lines):
        if re.match(r'^def main\s*\(', line):
            main_idx = idx
            break
    if main_idx is None:
        print("  [ERROR] no module-level `def main(` found. Not modified.")
        return False
    if main_idx > const_idx:
        print("  [ERROR] `def main(` appears after constant_args; unexpected layout.")
        return False
    print(f"  flags inserted: before `def main(` (line {main_idx + 1})")

    if 'import flags' not in src and 'from absl import' not in src:
        print("  [ERROR] absl flags do not appear to be imported. Not modified.")
        return False

    # Apply bottom-up so the earlier index stays valid.
    lines.insert(insert_idx, ARGS_TMPL.format(i=entry_indent))
    lines.insert(main_idx, FLAG_BLOCK)

    backup = path + '.pre_initial_guess'
    if not os.path.exists(backup):
        shutil.copy2(path, backup)
        print(f"  backup        : {backup}")

    with open(path, 'w') as f:
        f.write(''.join(lines))

    # Syntax-check the result; restore the backup if the edit broke the file.
    try:
        with open(path) as f:
            compile(f.read(), path, 'exec')
    except SyntaxError as e:
        shutil.copy2(backup, path)
        print(f"  [ERROR] patched file does not compile, reverted: {e}")
        return False

    print("  [ok] patched and compiles")
    return True


def revert_one(path):
    backup = path + '.pre_initial_guess'
    if os.path.exists(backup):
        shutil.copy2(backup, path)
        print(f"reverted {path} from {backup}")
        return True
    print(f"no backup for {path}")
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--wrapper', required=True, nargs='+',
                    help='Path(s) to run_multimer_jobs.py')
    ap.add_argument('--revert', action='store_true',
                    help='Restore from the .pre_initial_guess backup')
    args = ap.parse_args()

    fn = revert_one if args.revert else patch_one
    ok = all([fn(p) for p in args.wrapper])

    if not args.revert:
        print("\n" + "=" * 60)
        if ok:
            print("Done. run_multimer_jobs.py now accepts:")
            print("  --initial_guess_dir DIR")
            print("  --initial_guess_map FILE")
            print("\nA job whose description matches no seed runs UNGUIDED and logs")
            print("nothing. Check 'Using initial guess from' line count against the")
            print("number of folds before trusting a screen.")
        else:
            print("Some files were NOT patched -- see errors above.")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
