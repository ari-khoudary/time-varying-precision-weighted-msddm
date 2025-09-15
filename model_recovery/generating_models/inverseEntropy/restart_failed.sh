#!/bin/bash

# Script to restart failed or incomplete array jobs
# Usage: ./restart_failed.sh [results_directory]

RESULTS_DIR=${1:-"results/$(date +%Y-%m-%d)"}
LOG_DIR="$RESULTS_DIR/logs"
TOTAL_SUBJECTS=500  # Adjust based on your nSub setting
SUBJECTS_PER_JOB=10  # Must match subsPerJob in MATLAB script

if [ ! -d "$LOG_DIR" ]; then
    echo "Log directory not found: $LOG_DIR"
    echo "Cannot determine failed jobs without logs."
    exit 1
fi

echo "=== Checking for failed/incomplete jobs ==="
echo "Results directory: $RESULTS_DIR"
echo "Total subjects expected: $TOTAL_SUBJECTS"
echo

# Find completed subjects
if [ -f "$LOG_DIR/completed_subjects.txt" ]; then
    COMPLETED_SUBJECTS=($(grep -o "subject [0-9]\+" "$LOG_DIR/completed_subjects.txt" | grep -o "[0-9]\+" | sort -n | uniq))
    echo "Found ${#COMPLETED_SUBJECTS[@]} completed subjects"
else
    COMPLETED_SUBJECTS=()
    echo "No completed subjects log found"
fi

# Find missing subjects
MISSING_SUBJECTS=()
for ((i=1; i<=TOTAL_SUBJECTS; i++)); do
    if [[ ! " ${COMPLETED_SUBJECTS[@]} " =~ " ${i} " ]]; then
        MISSING_SUBJECTS+=($i)
    fi
done

echo "Missing subjects: ${#MISSING_SUBJECTS[@]}"

if [ ${#MISSING_SUBJECTS[@]} -eq 0 ]; then
    echo "All subjects completed successfully!"
    exit 0
fi

# Determine which array jobs need to be restarted
RESTART_JOBS=()
for subject in "${MISSING_SUBJECTS[@]}"; do
    job_id=$(( (subject - 1) / SUBJECTS_PER_JOB + 1 ))
    if [[ ! " ${RESTART_JOBS[@]} " =~ " ${job_id} " ]]; then
        RESTART_JOBS+=($job_id)
    fi
done

echo "Array jobs that need restarting: ${RESTART_JOBS[@]}"
echo

# Ask for confirmation
read -p "Do you want to restart these array jobs? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Restart cancelled."
    exit 0
fi

# Create array specification
ARRAY_SPEC=""
for job in "${RESTART_JOBS[@]}"; do
    if [ -z "$ARRAY_SPEC" ]; then
        ARRAY_SPEC="$job"
    else
        ARRAY_SPEC="$ARRAY_SPEC,$job"
    fi
done

echo "Submitting restart job with array specification: $ARRAY_SPEC"

# Submit the restart job
sbatch --array="$ARRAY_SPEC" submit_job.sh

echo "Restart jobs submitted. Monitor progress with: ./monitor_progress.sh"
