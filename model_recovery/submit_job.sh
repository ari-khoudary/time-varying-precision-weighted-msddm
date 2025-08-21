#!/bin/bash
#SBATCH --job-name=recovery
#SBATCH -A bornstea_lab
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -c 8
#SBATCH --output=slurm_messages/slurm-%A_%a.out 
#SBATCH --error=slurm_messages/slurm-%A_%a.err
#SBATCH -t 2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=makhouda@uci.edu

mkdir -p slurm_messages
module load MATLAB/R2023b

matlab -nodisplay -nosplash -nodesktop -r "doSampling_wrapper; exit(0);"

echo "Job completed."

