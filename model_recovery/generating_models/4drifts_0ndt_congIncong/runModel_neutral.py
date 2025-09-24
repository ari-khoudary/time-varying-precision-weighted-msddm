# initialize
from asyncio import Condition
import pyddm
import pandas as pd
import numpy as np
from pyddm import Sample
import os
import sys
import matplotlib.pyplot as plt
import gc
import argparse
import pickle
import time
import json

# Simple checkpoint functions (minimal addition)
def save_checkpoint(stage, data, task_id, checkpoint_dir):
    """Save checkpoint data"""
    checkpoint_file = os.path.join(checkpoint_dir, f"checkpoint_task_{task_id}.pkl")
    status_file = os.path.join(checkpoint_dir, f"status_task_{task_id}.json")
    
    # Save checkpoint data
    checkpoint_data = {
        'stage': stage,
        'data': data,
        'timestamp': time.time(),
        'readable_time': time.strftime('%Y-%m-%d %H:%M:%S')
    }
    with open(checkpoint_file, 'wb') as f:
        pickle.dump(checkpoint_data, f)
    
    # Update status
    status_data = {
        'task_id': task_id,
        'current_stage': stage,
        'status': 'running',
        'last_update': time.time(),
        'readable_time': time.strftime('%Y-%m-%d %H:%M:%S')
    }
    with open(status_file, 'w') as f:
        json.dump(status_data, f)
    
    print(f"Checkpoint saved: {stage}")

def load_checkpoint(task_id, checkpoint_dir):
    """Load checkpoint if exists"""
    checkpoint_file = os.path.join(checkpoint_dir, f"checkpoint_task_{task_id}.pkl")
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'rb') as f:
            return pickle.load(f)
    return None

def mark_completed(task_id, checkpoint_dir):
    """Mark job as completed and clean up"""
    status_file = os.path.join(checkpoint_dir, f"status_task_{task_id}.json")
    checkpoint_file = os.path.join(checkpoint_dir, f"checkpoint_task_{task_id}.pkl")
    
    # Update status to completed
    status_data = {
        'task_id': task_id,
        'current_stage': 'completed',
        'status': 'completed',
        'last_update': time.time(),
        'readable_time': time.strftime('%Y-%m-%d %H:%M:%S')
    }
    with open(status_file, 'w') as f:
        json.dump(status_data, f)
    
    # Remove checkpoint file
    if os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)

# Get number of CPUs from SLURM environment variable
n_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
pyddm.set_N_cpus(n_cpus)

# Get the cue, drifts, and threshold from command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--cue', type=float, required=True)
parser.add_argument('--n1', type=float, required=True)
parser.add_argument('--s1', type=float, required=True)
parser.add_argument('--n2', type=float, required=True)
parser.add_argument('--s2', type=float, required=True)
parser.add_argument('--thresh', type=float, required=True)
parser.add_argument('--checkpoint-dir', type=str, default='checkpoints')
args = parser.parse_args()

cue = args.cue
n1_drift = args.n1
s1_drift = args.s1
n2_drift = args.n2
s2_drift = args.s2
thresh = args.thresh

# Checkpoint setup
task_id = os.environ.get('SLURM_ARRAY_TASK_ID', 'test')
checkpoint_dir = os.path.join(args.checkpoint_dir, f"task_{task_id}")
os.makedirs(checkpoint_dir, exist_ok=True)

# Create results directory if it doesn't exist
results_dir = 'results/'
os.makedirs(results_dir, exist_ok=True)

# define one drift function for neutral and one for biased cues
def drift_neutral(t, signal1_onset, noise2_onset, signal2_onset,
                  n1_neut, s1_neut, n2_neut, s2_neut):
    """Drift function for neutral cue (0.5)"""
    if t < signal1_onset:
        return n1_neut
    elif t >= signal1_onset and t < noise2_onset:
        return s1_neut
    elif t >= noise2_onset and t < signal2_onset:
        return n2_neut
    elif t >= signal2_onset:
        return s2_neut

# define arrays of onsets
onsets_5q = {
    'signal1_onset': 0.65164947953539,
    'noise2_onset': 1.12322767774387,
    'signal2_onset': 1.82007829376656}

onsets_25q = {
    'signal1_onset': 0.680087375592259,
    'noise2_onset': 1.19414194805919,
    'signal2_onset': 1.92522084540362
}

onsets_50q = {
    'signal1_onset': 0.733298731308342,
    'noise2_onset': 1.27795054776867,
    'signal2_onset': 2.03740795226059
}

onsets_75q = {
    'signal1_onset': 0.817374716889882,
    'noise2_onset': 1.39719884239008,
    'signal2_onset': 2.19037702319693
}

onsets_95q = {
    'signal1_onset': 1.01164563696949,
    'noise2_onset': 1.64670498464725,
    'signal2_onset': 2.47961739538071
}

onsets_100q = {
    'signal1_onset': 1.51837789329293,
    'noise2_onset': 2.63630447406904,
    'signal2_onset': 3.52063445587701
}

# total number of trials per cue
nTrial = 500
quantiles = np.array([0, 0.05, 0.25, 0.50, .75, 0.95, 1])
bin_proportion = np.diff(quantiles)
trial_bins = np.round(bin_proportion * nTrial).astype(int)

try:
    # Check for existing checkpoint
    checkpoint = load_checkpoint(task_id, checkpoint_dir)
    
    if checkpoint:
        print(f"Resuming from stage: {checkpoint['stage']}")
        stage = checkpoint['stage']
        data = checkpoint['data']
    else:
        print("Starting fresh job")
        stage = 'initialization'
        data = {}
    
    # Stage 1: Model initialization
    if stage == 'initialization':
        print("Stage 1: Initializing simulation model...")
        # initialize generating model
        sim_model = pyddm.gddm(
                drift = drift_neutral,
                starting_position= 0,
                bound='B',
                T_dur = 4.3,
                nondecision=0,
                parameters={'B': thresh, 
                            'n1_neut': n1_drift, 's1_neut': s1_drift, 
                            'n2_neut': n2_drift, 's2_neut': s2_drift},
                conditions = ['signal1_onset', 'noise2_onset', 'signal2_onset'])
        
        data['sim_model_params'] = {
            'B': thresh, 'n1_neut': n1_drift, 's1_neut': s1_drift, 
            'n2_neut': n2_drift, 's2_neut': s2_drift
        }
        save_checkpoint('model_initialized', data, task_id, checkpoint_dir)
        stage = 'model_initialized'
    
    # Stage 2: Data simulation
    if stage == 'model_initialized':
        print("Stage 2: Running simulation...")
        # Recreate model if resuming
        sim_model = pyddm.gddm(
                drift = drift_neutral,
                starting_position= 0,
                bound='B',
                T_dur = 4.3,
                nondecision=0,
                parameters=data['sim_model_params'],
                conditions = ['signal1_onset', 'noise2_onset', 'signal2_onset'])
        
        # generate data
        onsets_list = [onsets_5q, onsets_25q, onsets_50q, onsets_75q, onsets_95q, onsets_100q]
        samples = []
        for i, onsets in enumerate(onsets_list):
            samp = sim_model.solve(conditions = onsets).sample(trial_bins[i])
            samples.append(samp)

        # concatenate into one sample
        sample = sum(samples)
        data['sample'] = sample
        save_checkpoint('simulation_complete', data, task_id, checkpoint_dir)
        stage = 'simulation_complete'
    
    # Stage 3: Model fitting
    if stage == 'simulation_complete':
        print("Stage 3: Fitting model...")
        sample = data['sample']
        
        # initialize fitting models 
        fitting_model = pyddm.gddm(
            drift = drift_neutral,
            starting_position = 0,
            bound="B",
            T_dur = 4.3,
            nondecision=0,
            parameters={'B': (1,15), 
                        'n1_neut': (0, 10), 's1_neut': (0, 10), 
                        'n2_neut': (0, 10), 's2_neut': (0, 10)},
            conditions = ['signal1_onset', 'noise2_onset', 'signal2_onset']
        )
        
        # fit
        fitting_model.fit(sample, verbose=True)
        
        # Gather results
        loss = fitting_model.get_fit_result().value()
        params = fitting_model.parameters()
        
        data['loss'] = loss
        data['fitted_params'] = params
        data['fitting_model_params'] = fitting_model.parameters()
        save_checkpoint('fitting_complete', data, task_id, checkpoint_dir)
        stage = 'fitting_complete'
    
    # Stage 4: Save results
    if stage == 'fitting_complete':
        print("Stage 4: Saving results...")
        loss = data['loss']
        params = data['fitted_params']
        
        # Recreate models for display if needed
        sim_model = pyddm.gddm(
                drift = drift_neutral,
                starting_position= 0,
                bound='B',
                T_dur = 4.3,
                nondecision=0,
                parameters=data['sim_model_params'],
                conditions = ['signal1_onset', 'noise2_onset', 'signal2_onset'])
        
        fitting_model = pyddm.gddm(
            drift = drift_neutral,
            starting_position = 0,
            bound="B",
            T_dur = 4.3,
            nondecision=0,
            parameters=data['fitting_model_params'],
            conditions = ['signal1_onset', 'noise2_onset', 'signal2_onset']
        )
        
        # Save text results (basic summary) - updated filename convention
        with open(os.path.join(results_dir, f'{cue}cue_{thresh}thresh_{n1_drift}n1_{s1_drift}s1_{n2_drift}n2_{s2_drift}s2_fits.txt'), 'w') as f:
            f.write(f'\nLoss: {loss}\nParameters:\n')
            for param, value in params.items():
                f.write(f'{param}: {value}\n')
        
        original_stdout = sys.stdout
        with open(os.path.join(results_dir, f'{cue}cue_{thresh}thresh_{n1_drift}n1_{s1_drift}s1_{n2_drift}n2_{s2_drift}s2_summary.txt'), 'w') as f:
            sys.stdout = f
            f.write(f'\nGROUND TRUTH: \n')
            sim_model.show()
            f.write(f'\n\n FITTED VALUES: \n')
            fitting_model.show()
        sys.stdout = original_stdout
        
        mark_completed(task_id, checkpoint_dir)
        print("Job completed successfully!")
    
except Exception as e:
    error_msg = f"Error processing job: {str(e)}"
    print(error_msg)
    
    # Save error to status file
    status_file = os.path.join(checkpoint_dir, f"status_task_{task_id}.json")
    status_data = {
        'task_id': task_id,
        'parameters': {  # Add this
            'cue': cue,
            'n1_drift': n1_drift,
            's1_drift': s1_drift,
            'n2_drift': n2_drift,
            's2_drift': s2_drift,
            'thresh': thresh
        },
        'current_stage': stage if 'stage' in locals() else 'unknown',
        'status': 'failed',
        'error_message': error_msg,
        'last_update': time.time(),
        'readable_time': time.strftime('%Y-%m-%d %H:%M:%S')
    }
    with open(status_file, 'w') as f:
        json.dump(status_data, f)
    
    # Save error information - updated filename convention
    with open(os.path.join(results_dir, f'{cue}cue_{thresh}thresh_{n1_drift}n1_{s1_drift}s1_{n2_drift}n2_{s2_drift}s2_error.txt'), 'w') as f:
        f.write(error_msg)
    
    raise