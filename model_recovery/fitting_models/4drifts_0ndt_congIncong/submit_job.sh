#!/bin/bash
#SBATCH --job-name=fit_recovery
#SBATCH -A bornstea_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 16
#SBATCH --array=1-7  # 1-N, where N = nSubs * nCoherences
#SBATCH --output=slurm_messages/slurm-%A_%a.out
#SBATCH --error=slurm_messages/slurm-%A_%a.err
#SBATCH -t 1-12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=makhouda@uci.edu


# print job start time
echo "Job started at: $(date)"
# setup
mkdir -p slurm_messages
source ~/.bashrc
module load python/3.10.2

# Define coherence levels in the simulated data
COHERENCE_LEVELS=(0.61 0.56 0.54 0.53 0.52 0.515 0.505)
NUM_COHERENCES=${#COHERENCE_LEVELS[@]}

# Calculate subject index and coherence index from SLURM_ARRAY_TASK_ID
SUBJECT_INDEX=$(((SLURM_ARRAY_TASK_ID - 1) / NUM_COHERENCES + 1))
COHERENCE_INDEX=$(((SLURM_ARRAY_TASK_ID - 1) % NUM_COHERENCES))

# Get the specific coherence value
COHERENCE=${COHERENCE_LEVELS[$COHERENCE_INDEX]}

# Get the subject ID 
SUBJECT="sub${SUBJECT_INDEX}"

# Skip if no subject found
if [ -z "$SUBJECT" ]; then
    echo "No subject found for task ID $SLURM_ARRAY_TASK_ID (subject index: $SUBJECT_INDEX). Exiting."
    exit 0
fi

echo "Processing $SUBJECT, Coherence: $COHERENCE"
echo "Task ID: $SLURM_ARRAY_TASK_ID, Subject Index: $SUBJECT_INDEX, Coherence Index: $COHERENCE_INDEX"

# Run the fitting code with subject and coherence parameters
python -u fit_model.py ${SUBJECT} --coherence ${COHERENCE}

# timestamp completion
echo "Completed Subject: $SUBJECT, Coherence: $COHERENCE"
echo "Job ended at: $(date)"
