#!/usr/bin/env python3
"""
Job monitoring and progress tracking system
Run this script periodically to monitor job progress
"""

import json
import time
import os
from pathlib import Path
import pandas as pd
from collections import defaultdict, Counter
import argparse

class JobMonitor:
    def __init__(self, checkpoint_base_dir="checkpoints", results_base_dir="results"):
        self.checkpoint_dir = Path(checkpoint_base_dir)
        self.results_dir = Path(results_base_dir)
        self.total_jobs = 1492992  # Update this to match your array size
        
    def scan_job_status(self):
        """Scan all job status files and return summary"""
        status_files = list(self.checkpoint_dir.glob("*/status_task_*.json"))
        
        job_statuses = []
        stages = Counter()
        statuses = Counter()
        
        for status_file in status_files:
            try:
                with open(status_file, 'r') as f:
                    status_data = json.load(f)
                job_statuses.append(status_data)
                stages[status_data.get('current_stage', 'unknown')] += 1
                statuses[status_data.get('status', 'unknown')] += 1
            except Exception as e:
                print(f"Error reading {status_file}: {e}")
        
        return job_statuses, stages, statuses
    
    def get_slurm_job_info(self, job_id=None):
        """Get SLURM job information"""
        try:
            if job_id:
                cmd = f"squeue -j {job_id} --format='%.10i %.8T %.10M %.6D %.20S'"
            else:
                cmd = "squeue -u $USER --format='%.10i %.8T %.10M %.6D %.20S'"
            
            import subprocess
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            return result.stdout
        except Exception as e:
            return f"Error getting SLURM info: {e}"
    
    def generate_report(self, detailed=False):
        """Generate comprehensive monitoring report"""
        print("=" * 60)
        print(f"DDM PARAMETER RECOVERY JOB MONITOR - {time.strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 60)
        
        # Scan job statuses
        job_statuses, stages, statuses = self.scan_job_status()
        
        # Overall progress
        jobs_started = len(job_statuses)
        jobs_completed = statuses.get('completed', 0)
        jobs_failed = statuses.get('failed', 0)
        jobs_running = statuses.get('running', 0)
        
        print(f"\nOVERALL PROGRESS:")
        print(f"  Total jobs in array: {self.total_jobs:,}")
        print(f"  Jobs started: {jobs_started:,} ({jobs_started/self.total_jobs*100:.1f}%)")
        print(f"  Jobs completed: {jobs_completed:,} ({jobs_completed/self.total_jobs*100:.1f}%)")
        print(f"  Jobs running: {jobs_running:,} ({jobs_running/self.total_jobs*100:.1f}%)")
        print(f"  Jobs failed: {jobs_failed:,} ({jobs_failed/self.total_jobs*100:.1f}%)")
        print(f"  Jobs not started: {self.total_jobs - jobs_started:,}")
        
        # Stage breakdown
        print(f"\nCURRENT STAGES:")
        for stage, count in stages.most_common():
            print(f"  {stage}: {count:,} jobs")
        
        # Recent activity
        if job_statuses:
            recent_updates = sorted(job_statuses, key=lambda x: x.get('last_update', 0), reverse=True)[:10]
            print(f"\nRECENT ACTIVITY (last 10 updates):")
            for job in recent_updates:
                task_id = job.get('task_id', 'unknown')
                stage = job.get('current_stage', 'unknown')
                status = job.get('status', 'unknown')
                update_time = time.strftime('%H:%M:%S', time.localtime(job.get('last_update', 0)))
                print(f"  Task {task_id}: {stage} ({status}) at {update_time}")
        
        # Failed jobs
        failed_jobs = [job for job in job_statuses if job.get('status') == 'failed']
        if failed_jobs:
            print(f"\nFAILED JOBS ({len(failed_jobs)} total):")
            for job in failed_jobs[:5]:  # Show first 5
                task_id = job.get('task_id', 'unknown')
                error = job.get('error_message', 'No error message')
                print(f"  Task {task_id}: {error[:100]}...")
            if len(failed_jobs) > 5:
                print(f"  ... and {len(failed_jobs) - 5} more failed jobs")
        
        # Performance metrics
        if jobs_completed > 0:
            completed_jobs = [job for job in job_statuses if job.get('status') == 'completed']
            if completed_jobs:
                # Estimate completion times (you'll need to track start times)
                print(f"\nPERFORMANCE METRICS:")
                print(f"  Completed jobs: {len(completed_jobs)}")
                # Add more metrics as needed
        
        # SLURM queue status
        print(f"\nSLURM QUEUE STATUS:")
        slurm_info = self.get_slurm_job_info()
        print(slurm_info)
        
        # Detailed report
        if detailed and job_statuses:
            print(f"\nDETAILED STATUS (showing all {len(job_statuses)} active jobs):")
            df = pd.DataFrame(job_statuses)
            if not df.empty:
                summary = df.groupby(['current_stage', 'status']).size().unstack(fill_value=0)
                print(summary)
        
        print("=" * 60)
    
    def find_stuck_jobs(self, hours_threshold=6):
        """Find jobs that have been stuck in the same stage for too long"""
        job_statuses, _, _ = self.scan_job_status()
        current_time = time.time()
        stuck_jobs = []
        
        for job in job_statuses:
            if job.get('status') == 'running':
                last_update = job.get('last_update', 0)
                hours_since_update = (current_time - last_update) / 3600
                
                if hours_since_update > hours_threshold:
                    stuck_jobs.append({
                        'task_id': job.get('task_id'),
                        'stage': job.get('current_stage'),
                        'hours_stuck': hours_since_update,
                        'parameters': job.get('parameters', {})
                    })
        
        return stuck_jobs
    
    def save_progress_log(self, filename=None):
        """Save progress to a CSV file for later analysis"""
        if filename is None:
            filename = f"progress_log_{time.strftime('%Y%m%d_%H%M%S')}.csv"
        
        job_statuses, _, _ = self.scan_job_status()
        
        if job_statuses:
            df = pd.DataFrame(job_statuses)
            df.to_csv(filename, index=False)
            print(f"Progress saved to {filename}")
        else:
            print("No job data to save")

def main():
    parser = argparse.ArgumentParser(description='Monitor DDM parameter recovery jobs')
    parser.add_argument('--detailed', action='store_true', help='Show detailed report')
    parser.add_argument('--stuck-threshold', type=int, default=6, help='Hours threshold for stuck jobs')
    parser.add_argument('--save-log', action='store_true', help='Save progress to CSV')
    parser.add_argument('--watch', type=int, help='Watch mode: refresh every N seconds')
    
    args = parser.parse_args()
    
    monitor = JobMonitor()
    
    def run_monitor():
        monitor.generate_report(detailed=args.detailed)
        
        # Check for stuck jobs
        stuck_jobs = monitor.find_stuck_jobs(hours_threshold=args.stuck_threshold)
        if stuck_jobs:
            print(f"\nWARNING: {len(stuck_jobs)} jobs may be stuck:")
            for job in stuck_jobs[:5]:  # Show first 5
                print(f"  Task {job['task_id']}: {job['stage']} for {job['hours_stuck']:.1f} hours")
        
        if args.save_log:
            monitor.save_progress_log()
    
    if args.watch:
        try:
            while True:
                os.system('clear')  # Clear screen
                run_monitor()
                print(f"\nRefreshing in {args.watch} seconds... (Ctrl+C to exit)")
                time.sleep(args.watch)
        except KeyboardInterrupt:
            print("\nMonitoring stopped.")
    else:
        run_monitor()

if __name__ == "__main__":
    main()