#!/bin/bash

#SBATCH --job-name=neutral
#SBATCH -A bornstea_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 8
#SBATCH --output=slurm_messages/slurm-%A_%a.out 
#SBATCH --error=slurm_messages/slurm-%A_%a.err
#SBATCH -t 1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=makhouda@uci.edu
#SBATCH --array=0-45  # Adjust based on total number of subject-cue combinations

mkdir -p slurm_messages

source ~/.bashrc
module load python/3.10.2

SUBJECT=$(python -c "
import pandas as pd
df = pd.read_csv('inference_test.csv')
# Get subjects that have trueCue == 0.5 data
subjects_with_neutral = df[df['trueCue'] == 0.5]['subID'].unique()

if $SLURM_ARRAY_TASK_ID < len(subjects_with_neutral):
    print(subjects_with_neutral[$SLURM_ARRAY_TASK_ID])
")

# Skip if no subject found
if [ -z "$SUBJECT" ]; then
    echo "No subject found for task ID $SLURM_ARRAY_TASK_ID. Exiting."
    exit 0
fi

# Set cue to 0.5 since we're only fitting neutral cue data
CUE=0.5

# Run the fitting code
python -u fit_model.py ${SUBJECT} ${CUE}
