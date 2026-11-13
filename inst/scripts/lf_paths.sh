# W = the ligandFinder source tree these scripts read python, R code and data
# from: $LF_SRC if set, else the checkout this file sits in (at inst/scripts/).
# See lf_paths.R for why this is a checkout and not the installed package -- in
# short, the install's inst/python has no lf_winxgb and no lf_resid. A script
# in inst/scripts/<folder>/ sources it from one level up:
#   . ${0:A:h:h}/lf_paths.sh
# zsh gives a sourced file its own $0, so ${0:A:h} below is THIS file's
# directory, whichever script sourced it.
lf_is_src() { [[ -f $1/DESCRIPTION && -d $1/R ]] &&
  grep -q '^Package:[[:space:]]*ligandFinder[[:space:]]*$' $1/DESCRIPTION }
if [[ -n ${LF_SRC:-} ]]; then
  W=${~LF_SRC}; W=${W:A}
  lf_is_src $W || { print -u2 "LF_SRC is not a ligandFinder source tree (no DESCRIPTION / R/): $LF_SRC"; exit 1 }
else
  W=${0:A:h}
  while [[ $W != / ]] && ! lf_is_src $W; do W=${W:h}; done
  lf_is_src $W || { print -u2 \
    "no ligandFinder source tree at or above ${0:A:h} -- run these from a checkout, or set LF_SRC=~/R_projects/ligandFinder"; exit 1 }
fi

# Where each kind of script lives (PACKAGE_MAP.md). A script calls another
# through these, so it does not matter which folder the caller sits in.
LF_SCRIPTS=$W/inst/scripts
LF_RDB=$LF_SCRIPTS/generate_residue_db
LF_MOD=$LF_SCRIPTS/modeling
LF_EVAL=$LF_SCRIPTS/evaluation
LF_VIZ=$LF_SCRIPTS/visualization

# LF_RUN = the arm pipeline's run directory: the lists, blurbs and tables it
# reads (page_genes.txt, cv5_unseen_ends*.csv, final_blurb*.html, index_cards.rds,
# ...) and the logs it writes. Until 2026-11 these sat beside the code in the
# working copy ~/AF2_analysis/lf_nomask_scripts_20261021. R scripts read the
# same variable (S <- Sys.getenv("LF_RUN", ...)).
export LF_RUN=${LF_RUN:-$HOME/AF2_analysis/lf_run}
