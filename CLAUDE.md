# ligandFinder: rules for Claude

`PACKAGE_MAP.md` says where each kind of code lives, what the pipelines are, and
where their data comes from and goes. Read it before adding, moving or running
code.

## Where code goes
1. **New code goes where the map says.** If no folder fits, add one and update
   `PACKAGE_MAP.md` in the same commit. **Update the map whenever you add a
   folder** or change what a folder holds.
2. **No code outside the package.** This covers:
   - scripts left in `~/AF2_analysis`, on the Desktop or in a scratchpad;
   - edited copies of package scripts;
   - a working copy that runs instead of the checkout.

   If a script produced a number, a figure or a file that anyone may use or
   cite, commit it to the package. Run pipelines from a checkout
   (`LF_SRC=<checkout>`), never from a copy.
3. **No one-off scripts outside `inst/scripts/analysis/`.**
   - A one-off goes in `inst/scripts/analysis/<question>/`, one folder per
     question.
   - Its README says what was asked and what was found, or "not recorded".
   - Prefer extending a pipeline script to writing a one-off.
   - When a one-off becomes part of a pipeline, move it to `modeling/`,
     `evaluation/` or `visualization/`.
4. **No dated folder names.** Name a folder for what it holds. Dates go in
   READMEs, commit messages and output directories, not in code paths.
5. **No data in code folders.**
   - Inputs and outputs live under `~/AF2_analysis`. The arm pipeline's lists,
     blurbs and logs go in `$LF_RUN`.
   - Only small, versioned reference tables go in `inst/extdata/` or `data/`.
   - A script reads its data paths from `lf_paths` / `LF_RUN` or from its
     arguments; never from its own directory.
6. **AlphaFold complex-model code goes in fap** (`~/R_projects/fap`).
   ligandFinder must never depend on fap. What fap produces (for example the
   pocket-label CSV) reaches ligandFinder as data, listed in the map.

## How scripts find things
- **Shell:** `. ${0:A:h:h}/lf_paths.sh`, then `$LF_RDB`, `$LF_MOD`, `$LF_EVAL`,
  `$LF_VIZ`, `$LF_RUN`.
- **R:** `source(file.path(.lf_here, "..", "lf_paths.R"))`, then
  `lf_src("inst/scripts/<folder>/<file>")` or `lf_src("R/<file>")`.
- **Python:** put `inst/scripts/` on `sys.path`, then
  `import lf_paths; lf_paths.add_python()`.
- **Never** a path to a sibling script relative to the caller's directory, and
  never a hard-coded checkout or working-copy path.

## Before running anything
- **Do not edit a script while a run is using it.** zsh and Rscript read files
  as they go, and an edit mid-run kills the run at the end. Check with
  `pgrep -f <script>` first.
- **`evaluation/end_calls.R`, `metrics_final.R` and `segments_*.R` write into
  the shared `~/AF2_analysis/lf_headcmp`.** Do not run them for a comparison or
  a check unless the output is pointed elsewhere (`HEADCMP=`).
- **Ask Kevin before any multi-agent workflow or fan-out of agents.** He is on
  the Pro plan.
