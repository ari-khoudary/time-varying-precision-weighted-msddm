#!/bin/bash
# Find failed jobs and resubmit them
python3 -c "
import json
import os
from pathlib import Path

failed_jobs = []
checkpoint_dir = Path('checkpoints')
for status_file in checkpoint_dir.glob('*/status_task_*.json'):
    try:
        with open(status_file) as f:
            status = json.load(f)
        if status.get('status') == 'failed':
            task_id = status.get('task_id')
            failed_jobs.append(task_id)
    except:
        pass

if failed_jobs:
    print(f'Found {len(failed_jobs)} failed jobs')
    print('Task IDs:', ','.join(map(str, failed_jobs[:10])))  # Show first 10
    # You can modify this to generate a new SBATCH array for just these jobs
else:
    print('No failed jobs found')
"