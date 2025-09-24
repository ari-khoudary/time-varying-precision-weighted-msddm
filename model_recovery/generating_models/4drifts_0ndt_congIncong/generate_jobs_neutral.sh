#!/bin/bash

# Script to generate a complete list of all job combinations
# Run this before submitting SLURM array job

echo "Generating complete job list..."

# Define the same parameter ranges as in submit_jobs.sh
CUE_LEVELS=(0.5)
N1=($(seq 0 2 10))
S1=($(seq 0 2 10))
N2=($(seq 0 2 10))
S2=($(seq 0 2 10))
B_LEVELS=($(seq 1 2 15))

# Calculate dimensions
NUM_CUE=${#CUE_LEVELS[@]}
NUM_N1=${#N1[@]}
NUM_S1=${#S1[@]} 
NUM_N2=${#N2[@]}
NUM_S2=${#S2[@]}
NUM_B=${#B_LEVELS[@]}

TOTAL_JOBS=$((NUM_CUE * NUM_N1 * NUM_S1 * NUM_N2 * NUM_S2 * NUM_B))

# Create output files
JOB_LIST_FILE="neutral_job_list_complete.txt"
STATUS_LOG_FILE="neutral_job_status_log.txt"

# Create header for complete job list
echo "# Complete Job List - Generated on $(date)" > "$JOB_LIST_FILE"
echo "# Total Jobs: $TOTAL_JOBS" >> "$JOB_LIST_FILE"
echo "# Format: TASK_ID,CUE,N1,S1,N2,S2,B" >> "$JOB_LIST_FILE"
echo "TASK_ID,CUE,N1,S1,N2,S2,B" >> "$JOB_LIST_FILE"

# Create header for status log (this will track completed jobs)
echo "# Job Status Log - Started on $(date)" > "$STATUS_LOG_FILE"
echo "# This file tracks completed jobs. Jobs not listed here are still running or pending." >> "$STATUS_LOG_FILE"
echo "# Format: TASK_ID,CUE,N1,S1,N2,S2,B,STATUS,COMPLETION_TIME" >> "$STATUS_LOG_FILE"
echo "TASK_ID,CUE,N1,S1,N2,S2,B,STATUS,COMPLETION_TIME" >> "$STATUS_LOG_FILE"

echo "Generating $TOTAL_JOBS job combinations..."

# Generate all combinations
for ((task_id=1; task_id<=TOTAL_JOBS; task_id++)); do
    # Convert to 0-based indexing for calculations
    calc_id=$((task_id - 1))
    
    # Calculate indices for each parameter (same logic as submit_jobs.sh)
    
    b_idx=$((calc_id % NUM_B))
    calc_id=$((calc_id / NUM_B))
    
    s2_idx=$((calc_id % NUM_S2))
    calc_id=$((calc_id / NUM_S2))
    
    n2_idx=$((calc_id % NUM_N2))
    calc_id=$((calc_id / NUM_N2))
    
    s1_idx=$((calc_id % NUM_S1))
    calc_id=$((calc_id / NUM_S1))
    
    n1_idx=$((calc_id % NUM_N1))
    cue_idx=$((calc_id / NUM_N1))
    
    # Extract parameter values
    cue=${CUE_LEVELS[$cue_idx]}
    n1_val=${N1[$n1_idx]}
    s1_val=${S1[$s1_idx]}
    n2_val=${N2[$n2_idx]}
    s2_val=${S2[$s2_idx]}
    b=${B_LEVELS[$b_idx]}
    
    # Write to job list file
    echo "$task_id,$cue,$n1_val,$s1_val,$n2_val,$s2_val,$b" >> "$JOB_LIST_FILE"
    
    # Progress indicator
    if (( task_id % 5000 == 0 )); then
        echo "Generated $task_id/$TOTAL_JOBS jobs..."
    fi
done

echo "Job list generation complete!"
echo "Files created:"
echo "  - $JOB_LIST_FILE: Complete list of all $TOTAL_JOBS job combinations"
echo "  - $STATUS_LOG_FILE: Status log file (will be updated as jobs complete)"
echo ""
echo "You can now submit your SLURM job array and monitor progress by checking $STATUS_LOG_FILE"