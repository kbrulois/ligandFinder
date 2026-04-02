#!/bin/bash
#SBATCH --job-name=extract_metrics
#SBATCH --partition=normal
#SBATCH --time=2-00:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

ml R/4.1.2

JOB_DIR=$1

Rscript "$(dirname "$0")/afpd_extract_metrics_worker.R" "${JOB_DIR}" "${SLURM_ARRAY_TASK_ID}"
