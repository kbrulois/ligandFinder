#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --partition=normal
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=test.out

ml R/4.1.2

Rscript /home/groups/ebutcher/programs/pipeline/R_libs4.1/ligandFinder/exec/afpd_rename_files_cli.R \
--input=/oak/stanford/groups/ebutcher/deorphan-AI-ze/models/deeperCXCL142 \
--name=test \
--person=KB \
--location=SU
