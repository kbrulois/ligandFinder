#!/bin/bash
#SBATCH --job-name=bm_jobs
#SBATCH --partition=normal
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=test.out
#SBATCH --qos long

ml R/4.1.2  

Rscript /home/groups/ebutcher/programs/pipeline/R_libs4.1/ligandFinder/exec/afpd_submit_jobs.R \
--first=722
--last=999