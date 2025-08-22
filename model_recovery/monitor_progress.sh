#!/bin/bash

# Progress monitoring script
# Usage: ./monitor_progress.sh [results_directory]

RESULTS_DIR=${1:-"results/$(date +%Y-%m-%d)"}
LOG_DIR="$RESULTS_DIR/logs"

if [ ! -d "$LOG_DIR" ]; then
    echo "Log directory not found: $LOG_DIR"
    echo "Make sure the simulation has started."
    exit 1
fi

echo "=== Simulation Progress Monitor ==="
echo "Results directory: $RESULTS_DIR"
echo "Log directory: $LOG_DIR"
echo "Current time: $(date)"
echo

# Count completed CSV files
if [ -d "$RESULTS_DIR" ]; then
    COMPLETED_FILES=$(ls -1 "$RESULTS_DIR"/sub*.csv 2>/dev/null | wc -l)
    echo "Completed subject files: $COMPLETED_FILES"
else
    COMPLETED_FILES=0
    echo "No results directory found yet."
fi

# Show recent log entries if available
if [ -f "$LOG_DIR/completed_subjects.txt" ]; then
    echo
    echo "=== Recent completions ==="
    tail -10 "$LOG_DIR/completed_subjects.txt"
fi

# Show active SLURM jobs
echo
echo "=== Active SLURM jobs ==="
squeue -u $(whoami) -n recovery_array --format="%.10i %.12j %.8T %.10M %.6D %.20S"

# Show array job status if available
JOB_IDS=$(squeue -u $(whoami) -n recovery_array -h -o "%A" | sort -u)
for JOB_ID in $JOB_IDS; do
    if [ ! -z "$JOB_ID" ]; then
        echo
        echo "Array job $JOB_ID status:"
        squeue -j $JOB_ID --format="%.15i %.8T %.10M %.6D" | head -20
    fi
done

# Show disk usage
echo
echo "=== Disk usage ==="
if [ -d "$RESULTS_DIR" ]; then
    du -sh "$RESULTS_DIR"
else
    echo "Results directory not created yet."
fi

# Estimate completion
if [ -f "$LOG_DIR/completed_subjects.txt" ]; then
    echo
    echo "=== Progress estimation ==="
    TOTAL_SUBJECTS=500  # Adjust based on your nSub setting
    COMPLETION_RATE=$(echo "scale=2; $COMPLETED_FILES / $TOTAL_SUBJECTS * 100" | bc -l)
    echo "Progress: $COMPLETED_FILES/$TOTAL_SUBJECTS subjects ($COMPLETION_RATE%)"
    
    if [ $COMPLETED_FILES -gt 0 ]; then
        # Estimate completion time based on recent completions
        RECENT_COMPLETIONS=$(tail -10 "$LOG_DIR/completed_subjects.txt" | wc -l)
        if [ $RECENT_COMPLETIONS -gt 1 ]; then
            FIRST_TIME=$(tail -10 "$LOG_DIR/completed_subjects.txt" | head -1 | cut -d: -f1-2)
            LAST_TIME=$(tail -1 "$LOG_DIR/completed_subjects.txt" | cut -d: -f1-2)
            echo "Recent completion rate: $RECENT_COMPLETIONS subjects in recent timeframe"
        fi
    fi
fi
