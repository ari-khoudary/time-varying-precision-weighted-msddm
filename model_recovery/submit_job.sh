#!/bin/bash
#SBATCH --job-name=recovery
#SBATCH -A bornstea_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 16
#SBATCH --array=1-50%10
#SBATCH --output=slurm_messages/slurm-%A_%a.out
#SBATCH --error=slurm_messages/slurm-%A_%a.err
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=makhouda@uci.edu

# Create necessary directories
mkdir -p slurm_messages
module load MATLAB/R2023b

# Print job information
echo "Starting SLURM array job"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"

# Run MATLAB simulation
matlab -nodisplay -nosplash -nodesktop -r "doSampling_wrapper; exit(0);"

echo "Array task $SLURM_ARRAY_TASK_ID completed at $(date)"
