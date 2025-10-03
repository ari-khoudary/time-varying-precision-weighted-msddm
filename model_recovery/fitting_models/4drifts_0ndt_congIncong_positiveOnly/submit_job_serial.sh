#!/bin/bash
#SBATCH --job-name=serial+
#SBATCH -A bornstea_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 16
#SBATCH --array=1-90 # 1-N where N = nSubs(5) * nCoherence (6) * nCue (3)
#SBATCH --output=slurm_messages/slurm-%A_%a.out
#SBATCH --error=slurm_messages/slurm-%A_%a.err
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=makhouda@uci.edu

echo "Job started at: $(date)"
mkdir -p slurm_messages
source ~/.bashrc
module load python/3.10.2

# Define coherence levels from the simulation
COHERENCE_LEVELS=(0.515 0.52 0.53 0.54 0.56 0.61)
CUE_LEVELS=(0.5 0.65 0.8)
NUM_COHERENCES=${#COHERENCE_LEVELS[@]}
NUM_CUES=${#CUE_LEVELS[@]}

# Calculate subject index and coherence index from SLURM_ARRAY_TASK_ID
TOTAL_COMBINATIONS=$((NUM_COHERENCES * NUM_CUES))
SUBJECT_INDEX=$(((SLURM_ARRAY_TASK_ID - 1) / TOTAL_COMBINATIONS + 1))
COMBINATION_INDEX=$(((SLURM_ARRAY_TASK_ID - 1) % TOTAL_COMBINATIONS))
COHERENCE_INDEX=$((COMBINATION_INDEX / NUM_CUES))
CUE_INDEX=$((COMBINATION_INDEX % NUM_CUES))

# Get the specific coherence and cue values, as well as subID
COHERENCE=${COHERENCE_LEVELS[$COHERENCE_INDEX]}
CUE=${CUE_LEVELS[$CUE_INDEX]}
SUBJECT="sub${SUBJECT_INDEX}"

# Skip if no subject found
if [ -z "$SUBJECT" ]; then
    echo "No subject found for task ID $SLURM_ARRAY_TASK_ID (subject index: $SUBJECT_INDEX). Exiting."
    exit 0
fi

echo "Processing Subject: $SUBJECT, Coherence: $COHERENCE, Cue: $CUE"
echo "Task ID: $SLURM_ARRAY_TASK_ID, Subject Index: $SUBJECT_INDEX, Combination Index: $COMBINATION_INDEX"

# Run the fitting code with subject and coherence parameters
python -u fit_model_serial.py ${SUBJECT} --coherence ${COHERENCE} --cue ${CUE}

# Move SLURM files to subject directory after job completes
mv "slurm_messages/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" "${SUBJECT_RESULTS_DIR}/${COHERENCE}coh_${CUE}cue_slurm.out" 2>/dev/null || true
mv "slurm_messages/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" "${SUBJECT_RESULTS_DIR}/${COHERENCE}coh_${CUE}cue_slurm.err" 2>/dev/null || true

# timestamp completion
echo "Completed Subject: $SUBJECT, Coherence: $COHERENCE, Cue: $CUE"
echo "Job ended at: $(date)"
