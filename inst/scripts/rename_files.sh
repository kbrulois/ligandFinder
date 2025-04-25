#!/bin/bash

ml R/4.1.2

export R_LIBS=/home/groups/ebutcher/programs/pipeline/R_libs4.1

conda activate spoc_venv

exec_path=$(Rscript -e 'cat(system.file("exec", package = "ligandFinder"))')

python "$exec_path/generate_spoc_json.py" $1

Rscript "$exec_path/afpd_rename_files_cli.R" \
--input=$1 \
--name=deepX14 \
--person=KB \
--location=SU
