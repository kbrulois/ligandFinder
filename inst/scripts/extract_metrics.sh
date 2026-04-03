#!/bin/bash
#SBATCH --job-name=extract_metrics
#SBATCH --partition=normal
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=extract_metrics_%A_%a.out
#SBATCH --error=extract_metrics_%A_%a.err

JOB_DIR=$1

ml R/4.1.2

WORKER_SCRIPT=$(Rscript -e ".libPaths('/home/groups/ebutcher/programs/pipeline/R_libs4.1'); cat(system.file('scripts/afpd_extract_metrics_worker.R', package='ligandFinder'))")

Rscript "${WORKER_SCRIPT}" "${JOB_DIR}" "${SLURM_ARRAY_TASK_ID}"
