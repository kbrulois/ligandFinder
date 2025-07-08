#!/bin/bash

ml R/4.1.2

export R_LIBS=/home/groups/ebutcher/programs/pipeline/R_libs4.1

exec_path=$(Rscript -e 'cat(system.file("exec", package = "ligandFinder"))')

Rscript "$exec_path/af3_rename_files_cli.R" \
--input=$1 \
--name=$2 \
--person=SE \
--location=SD
