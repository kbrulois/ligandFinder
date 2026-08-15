#!/bin/bash
## ---------------------------------------------------------------------------
## AlphaPulldown screening run WITH AlphaFold initial guess.
##
## Drop-in replacement for alphapulldown.sh.  Same arguments:
##     sbatch alphapulldown_initial_guess.sh <protein_list.txt> <output_path> [seed_dir]
##
## WHY THIS DOESN'T USE run_multimer_jobs.py
## -----------------------------------------
## run_multimer_jobs.py is a thin wrapper that reads --protein_lists, converts
## each line to an --input string, and shells out to run_structure_prediction.py
## with a fixed set of flags.  --initial_guess_dir is not among them, and the
## wrapper offers no way to add it.  So this script does the wrapper's job
## directly: same conversion, same flags (copied from the command the wrapper
## actually emits), plus the initial-guess flag.
##
## JOB FILE FORMAT
## ---------------
## AlphaPulldown custom mode, one fold per line:
##     NK2R;TKN1,72-107
##     GPR15;GP15L,24-81
## ";" separates proteins, "," separates a protein from its residue range.
## run_structure_prediction.py wants "+" and ":" instead, and takes all folds
## as one comma-separated --input list -- the conversion is done below.
##
## SEEDS
## -----
## Generate them first, from the same job file, so every description has a
## match:
##     Rscript $R_LIBS/ligandFinder/scripts/ap_initial_guess_from_job_file.R \
##         --jobs <protein_list.txt> --out <seed_dir>
##
## A job whose description is not a key in the seed directory runs UNGUIDED
## and says nothing about it.  Check the run's "Initial guess PDB for" lines
## against the number of folds before trusting a whole screen -- the count
## printed at the end of this script is there for that.
## ---------------------------------------------------------------------------

#SBATCH --job-name=ap_iguess
#SBATCH --time=16:00:00
#SBATCH -e run_dir/ap_iguess_%A_%a_err.txt
#SBATCH -o run_dir/ap_iguess_%A_%a_out.txt
#SBATCH -p gpu
#SBATCH --gpu_cmode=shared
#SBATCH --gpus=1
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH --mem=16000

set -euo pipefail

source $GROUP_HOME/programs/anaconda3/bin/activate $GROUP_HOME/programs/anaconda3/envs/AlphaPulldown

protein_list="$1"
output_path="$2"
seed_dir="${3:-}"

mkdir -p "$output_path"
umask 0002

MAXRAM=$(echo `ulimit -m` '/ 1024.0'|bc)
GPUMEM=`nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits|tail -1`
export XLA_PYTHON_CLIENT_MEM_FRACTION=`echo "scale=3;$MAXRAM / $GPUMEM"|bc`
export TF_FORCE_UNIFIED_MEMORY='1'

AP_ENV=$GROUP_HOME/programs/anaconda3/envs/AlphaPulldown
RSP=$AP_ENV/lib/python3.11/site-packages/alphapulldown/scripts/run_structure_prediction.py

# Job-file format -> run_structure_prediction format, all folds in one list.
#   NK2R;TKN1,72-107   ->   NK2R+TKN1:72-107
# Blank lines and # comments are dropped; lines are then joined with ",".
INPUTS=$(grep -vE '^[[:space:]]*(#|$)' "$protein_list" \
         | tr -d '\r' \
         | sed 's/;/+/g; s/,/:/g' \
         | paste -sd, -)

N_FOLDS=$(grep -cvE '^[[:space:]]*(#|$)' "$protein_list")
echo "protein_list : $protein_list  ($N_FOLDS folds)"
echo "output_path  : $output_path"
echo "seed_dir     : ${seed_dir:-<none - running UNGUIDED>}"

# Flags below mirror exactly what run_multimer_jobs.py emits, so results stay
# comparable to previous unguided screens.
ARGS=(
    --input "$INPUTS"
    --output_directory "$output_path"
    --num_cycle 3
    --num_predictions_per_model 1
    --data_directory /oak/stanford/groups/ebutcher/catherine/alphafold_db
    --features_directory /oak/stanford/groups/ebutcher/deorphan-AI-ze/alphapulldown/input_features/Homo_sapiens
    --pair_msa
    --nomsa_depth_scan
    --nomultimeric_template
    --fold_backend alphafold
    --nocompress_result_pickles
    --noremove_result_pickles
    --remove_keys_from_pickles
    --use_ap_style
    --use_gpu_relax
    --protein_delimiter +
    --models_to_relax None
)

if [[ -n "$seed_dir" ]]; then
    if [[ ! -d "$seed_dir" ]]; then
        echo "ERROR: seed_dir '$seed_dir' does not exist." >&2
        exit 1
    fi
    N_SEEDS=$(find "$seed_dir" -maxdepth 1 -name '*.pdb' | wc -l | tr -d ' ')
    echo "seeds found  : $N_SEEDS"
    if [[ "$N_SEEDS" -eq 0 ]]; then
        echo "ERROR: no .pdb files in '$seed_dir'." >&2
        exit 1
    fi
    ARGS+=(--initial_guess_dir "$seed_dir")
fi

python "$RSP" "${ARGS[@]}" 2>&1 | tee "$output_path/run.log"

# A silent map miss is the failure mode that matters here: the job runs
# unguided and the output looks entirely normal.  Compare counts.
if [[ -n "$seed_dir" ]]; then
    N_USED=$(grep -c 'Using initial guess from' "$output_path/run.log" || true)
    echo
    echo "folds requested        : $N_FOLDS"
    echo "folds given a seed     : $N_USED"
    if [[ "$N_USED" -lt "$N_FOLDS" ]]; then
        echo "WARNING: $(( N_FOLDS - N_USED )) fold(s) ran UNGUIDED -- no seed matched"
        echo "         their description. Descriptions AP looked for:"
        grep -oP "Now running predictions on \K[^ ]+" "$output_path/run.log" | sort -u | sed 's/^/           /'
        echo "         Seeds present:"
        find "$seed_dir" -maxdepth 1 -name '*.pdb' -printf '           %f\n' | sort
    fi
fi
