#!/bin/bash
#SBATCH --job-name=biased
#SBATCH -A bornstea_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 16
#SBATCH --array=1-1492992 # Total jobs with 2 cue levels
#SBATCH --output=slurm_messages/slurm-%A_%a.out
#SBATCH --error=slurm_messages/slurm-%A_%a.err
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=makhouda@uci.edu

echo "Job started at: $(date)"
mkdir -p slurm_messages
source ~/.bashrc
module load python/3.10.2

## define ranges of parameter values (same as job list generator)
CUE_LEVELS=(0.65 0.8)
# drift levels
N1=($(seq 0 2 10))
S1=($(seq 0 2 10))
S1_INCONG=($(seq 0 2 10))
N2=($(seq 0 2 10))
N2_INCONG=($(seq 0 2 10))
S2=($(seq 0 2 10))
S2_INCONG=($(seq 0 2 10))
# threshold levels
B_LEVELS=($(seq 1 2 15))

# Calculate array dimensions
NUM_CUE=${#CUE_LEVELS[@]}
NUM_N1=${#N1[@]}
NUM_S1=${#S1[@]}
NUM_S1_INCONG=${#S1_INCONG[@]}
NUM_N2=${#N2[@]}
NUM_N2_INCONG=${#N2_INCONG[@]}
NUM_S2=${#S2[@]}
NUM_S2_INCONG=${#S2_INCONG[@]}
NUM_B=${#B_LEVELS[@]}

# Convert SLURM_ARRAY_TASK_ID to parameter combination
TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))  # Convert to 0-based indexing

# Calculate indices for each parameter (same order as job list generator)
B_INDEX=$((TASK_ID % NUM_B))
TASK_ID=$((TASK_ID / NUM_B))

S2_INCONG_INDEX=$((TASK_ID % NUM_S2_INCONG))
TASK_ID=$((TASK_ID / NUM_S2_INCONG))

S2_INDEX=$((TASK_ID % NUM_S2))
TASK_ID=$((TASK_ID / NUM_S2))

N2_INCONG_INDEX=$((TASK_ID % NUM_N2_INCONG))
TASK_ID=$((TASK_ID / NUM_N2_INCONG))

N2_INDEX=$((TASK_ID % NUM_N2))
TASK_ID=$((TASK_ID / NUM_N2))

S1_INCONG_INDEX=$((TASK_ID % NUM_S1_INCONG))
TASK_ID=$((TASK_ID / NUM_S1_INCONG))

S1_INDEX=$((TASK_ID % NUM_S1))
TASK_ID=$((TASK_ID / NUM_S1))

N1_INDEX=$((TASK_ID % NUM_N1))
CUE_INDEX=$((TASK_ID / NUM_N1))

# Extract parameter values
CUE=${CUE_LEVELS[$CUE_INDEX]}
N1_VAL=${N1[$N1_INDEX]}
S1_VAL=${S1[$S1_INDEX]}
S1_INCONG_VAL=${S1_INCONG[$S1_INCONG_INDEX]}
N2_VAL=${N2[$N2_INDEX]}
N2_INCONG_VAL=${N2_INCONG[$N2_INCONG_INDEX]}
S2_VAL=${S2[$S2_INDEX]}
S2_INCONG_VAL=${S2_INCONG[$S2_INCONG_INDEX]}
B=${B_LEVELS[$B_INDEX]}

echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Parameters: CUE=$CUE, N1=$N1_VAL, S1=$S1_VAL, S1_INCONG=$S1_INCONG_VAL, N2=$N2_VAL, N2_INCONG=$N2_INCONG_VAL, S2=$S2_VAL, S2_INCONG=$S2_INCONG_VAL, B=$B"

# Setup checkpoint directory
CHECKPOINT_DIR="checkpoints/task_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$CHECKPOINT_DIR"

# Run the fitting code with all parameters
python -u runModel_biased.py --cue $CUE --n1 $N1_VAL --s1 $S1_VAL --s1_incong $S1_INCONG_VAL --n2 $N2_VAL --n2_incong $N2_INCONG_VAL --s2 $S2_VAL --s2_incong $S2_INCONG_VAL --b $B --checkpoint-dir "$CHECKPOINT_DIR"

# timestamp completion
echo "Completed: CUE=$CUE, N1=$N1_VAL, S1=$S1_VAL, S1_INCONG=$S1_INCONG_VAL, N2=$N2_VAL, N2_INCONG=$N2_INCONG_VAL, S2=$S2_VAL, S2_INCONG=$S2_INCONG_VAL, B=$B"
echo "Job ended at: $(date)"