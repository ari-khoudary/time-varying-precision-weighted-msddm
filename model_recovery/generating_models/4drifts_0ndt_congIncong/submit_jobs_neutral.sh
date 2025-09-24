#!/bin/bash
#SBATCH --job-name=neutral
#SBATCH -A bornstea_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 16
#SBATCH --array=1-11250  # 1 * 6 * 6 * 6 * 6 * 8 = 11,250 jobs
#SBATCH --output=slurm_messages/slurm-%A_%a.out
#SBATCH --error=slurm_messages/slurm-%A_%a.err
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=makhouda@uci.edu

echo "Job started at: $(date)"
mkdir -p slurm_messages
source ~/.bashrc
module load python/3.10.2

# Define parameter ranges (same as generator)
CUE_LEVELS=(0.5)
N1=($(seq 0 2 10))
S1=($(seq 0 2 10))
N2=($(seq 0 2 10))
S2=($(seq 0 2 10))
B_LEVELS=($(seq 1 2 15))

# Calculate array dimensions
NUM_CUE=${#CUE_LEVELS[@]}
NUM_N1=${#N1[@]}
NUM_S1=${#S1[@]}
NUM_N2=${#N2[@]}
NUM_S2=${#S2[@]}
NUM_B=${#B_LEVELS[@]}

# Convert SLURM_ARRAY_TASK_ID to parameter combination
TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))

# Calculate indices (same order as generator)
B_INDEX=$((TASK_ID % NUM_B))
TASK_ID=$((TASK_ID / NUM_B))

S2_INDEX=$((TASK_ID % NUM_S2))
TASK_ID=$((TASK_ID / NUM_S2))

N2_INDEX=$((TASK_ID % NUM_N2))
TASK_ID=$((TASK_ID / NUM_N2))

S1_INDEX=$((TASK_ID % NUM_S1))
TASK_ID=$((TASK_ID / NUM_S1))

N1_INDEX=$((TASK_ID % NUM_N1))
CUE_INDEX=$((TASK_ID / NUM_N1))

# Extract parameter values
CUE=${CUE_LEVELS[$CUE_INDEX]}
N1_VAL=${N1[$N1_INDEX]}
S1_VAL=${S1[$S1_INDEX]}
N2_VAL=${N2[$N2_INDEX]}
S2_VAL=${S2[$S2_INDEX]}
B=${B_LEVELS[$B_INDEX]}

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Parameters: CUE=$CUE, N1=$N1_VAL, S1=$S1_VAL, N2=$N2_VAL, S2=$S2_VAL, B=$B"

# Setup checkpoint directory
CHECKPOINT_DIR="checkpoints/task_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$CHECKPOINT_DIR"

# Run the fitting code
python -u runModel_neutral.py --cue $CUE --n1 $N1_VAL --s1 $S1_VAL --n2 $N2_VAL --s2 $S2_VAL --thresh $B --checkpoint-dir "$CHECKPOINT_DIR"

echo "Completed: CUE=$CUE, N1=$N1_VAL, S1=$S1_VAL, N2=$N2_VAL, S2=$S2_VAL, B=$B"
echo "Job ended at: $(date)"