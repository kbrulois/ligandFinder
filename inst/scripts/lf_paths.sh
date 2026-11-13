# W = the ligandFinder source tree these scripts read python, R code and data
# from: $LF_SRC if set, else the checkout this script directory sits in (the
# scripts are committed at inst/scripts/lf_nomask_scripts_20261021/). See
# lf_paths.R for why this is a checkout and not the installed package -- in
# short, the install's inst/python has no lf_winxgb and no lf_resid.
#   . ${0:A:h}/lf_paths.sh
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
