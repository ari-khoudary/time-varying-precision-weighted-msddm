#!/bin/bash

#SBATCH --job-name=negAll
#SBATCH -A bornstea_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 8
#SBATCH --output=slurm_messages/slurm-%A_%a.out 
#SBATCH --error=slurm_messages/slurm-%A_%a.err
#SBATCH -t 3-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=makhouda@uci.edu
#SBATCH --array=0-89  # Adjust based on total number of subject-cue combinations

mkdir -p slurm_messages

source ~/.bashrc
module load python/3.10.2

# Get the subject ID and cue value for this task
SUBJECT_CUE=$(python -c "
import pandas as pd
df = pd.read_csv('inference_test.csv')
# Get all unique subject-cue combinations
subject_cue_pairs = []
for subject in df['subID'].unique():
    subject_df = df[df['subID'] == subject]
    for cue in subject_df['trueCue'].unique():
        subject_cue_pairs.append((subject, cue))

if $SLURM_ARRAY_TASK_ID < len(subject_cue_pairs):
    subject, cue = subject_cue_pairs[$SLURM_ARRAY_TASK_ID]
    print(f'{subject} {cue}')
")

# Skip if no subject-cue pair found
if [ -z "$SUBJECT_CUE" ]; then
    echo "No subject-cue pair found for task ID $SLURM_ARRAY_TASK_ID. Exiting."
    exit 0
fi

# Parse subject and cue
SUBJECT=$(echo $SUBJECT_CUE | cut -d' ' -f1)
CUE=$(echo $SUBJECT_CUE | cut -d' ' -f2)

# Run the fitting code
python -u fit_model_negAll.py ${SUBJECT} ${CUE}
