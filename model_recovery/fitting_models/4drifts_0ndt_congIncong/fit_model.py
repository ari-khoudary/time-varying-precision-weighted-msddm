# initialize
import pyddm
import pandas as pd
import numpy as np
from pyddm import Sample
import os
import sys
import matplotlib.pyplot as plt
import gc
import argparse

# Get number of CPUs from SLURM environment variable
n_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
pyddm.set_N_cpus(n_cpus)

# Get the subject ID and coherence level from command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('subject_id')
parser.add_argument('--coherence', type=float, required=True)
args = parser.parse_args()

subject_id = args.subject_id
coherence = args.coherence

# Create results directory if it doesn't exist
results_dir = f'results/{coherence}coh'
os.makedirs(results_dir, exist_ok=True)

# make df filtering robust against type mismatches in bash and python
subject_id_str = str(subject_id)
subject_id_int = int(subject_id) 

# define drift function
def drift(t, congruent, cue, signal1_onset, noise2_onset, signal2_onset,
                n1_weak, n1_strong, n1_neut, 
                s1_cong_weak, s1_cong_strong, s1_incong_weak, s1_incong_strong, s1_neut,
                n2_cong_weak, n2_cong_strong, n2_incong_weak, n2_incong_strong, n2_neut, 
                s2_cong_weak, s2_cong_strong, s2_incong_weak, s2_incong_strong, s2_neut):
    
    # drift rate during first noise period
    if t < signal1_onset:
        if congruent == 'congruent': 
            if cue == 0.65:
                return n1_weak
            else:
                return n1_strong
        elif congruent == 'incongruent':
            if cue == 0.65:
                return -n1_weak
            else:
                return -n1_strong
        else:
            return n1_neut

    # drift rates during first signal period
    if t >= signal1_onset and t < noise2_onset:
        if congruent == 'congruent':
            if cue == 0.65:
                return s1_cong_weak
            else:
                return s1_cong_strong
        elif congruent == 'incongruent':
            if cue == 0.65:
                return -s1_incong_weak
            else:
                return -s1_incong_strong
        else:
            return s1_neut

    # drift rates during the second noise period
    if t >= noise2_onset and t < signal2_onset:
        if congruent == 'congruent':
            if cue == 0.65:
                return n2_cong_weak
            else:
                return n2_cong_strong
        elif congruent == 'incongruent':
            if cue == 0.65:
                return -n2_incong_weak
            else:
                return -n2_incong_strong
        else:
            return n2_neut

    # drift rates during the second signal period
    if t >= signal2_onset:
        if congruent == 'congruent':
            if cue == 0.65:
                return s2_cong_weak
            else:
                return s2_cong_strong
        elif congruent == 'incongruent':
            if cue == 0.65:
                return -s2_incong_weak
            else:
                return -s2_incong_strong
        else:
            return s2_neut

try:
    # Load and filter data
    df = pd.read_csv(f'../../../{args.subject_id}.csv') 
    df = df.dropna(subset=['freeChoice'])
    df[['signal1_onset', 'noise2_onset', 'signal2_onset']] = df[['signal1_onset', 'noise2_onset', 'signal2_onset']].fillna(0)
    subject_df = df[(df['subID'] == subject_id_str) | 
                (df['subID'] == subject_id_int)].copy()
    
    if len(subject_df) == 0:
        error_msg = f"Error: No data found for subject {subject_id}"
        with open(os.path.join(results_dir, f's{subject_id}_error.txt'), 'w') as f:
            f.write(error_msg)
        sys.exit(1)
    
    # Create sample
    sample = pyddm.Sample.from_pandas_dataframe(
        subject_df, rt_column_name='RTs', choice_column_name='freeChoice'
    )
    
    # initialize model
    model = pyddm.gddm(
        drift = drift,
        starting_position = 0,
        bound="B",
        T_dur = 4.1,
        nondecision=0,
        parameters={'B': (0.5, 12), 
                    'noise1_cong': (0, 10), 'noise1_neut': (0,10),
                    'signal1_cong': (0, 10), 'signal1_incong': (0, 10), 'signal1_neut': (0, 10),
                    'noise2_cong': (0, 10), 'noise2_incong': (0,10), 'noise2_neut': (0,10),
                    'signal2_cong': (0, 10), 'signal2_incong': (0, 10), 'signal2_neut': (0,10)},
        conditions = ['trueCongruence', 'signal1_onset', 'noise2_onset', 'signal2_onset']
    )

    # fit
    model.fit(sample, verbose=True)
    
    # Gather results
    loss = model.get_fit_result().value()
    params = model.parameters()
    
    # Save text results (basic summary)
    with open(os.path.join(results_dir, f's{subject_id}_results.txt'), 'w') as f:
        f.write(f'Subject: {subject_id}\nLoss: {loss}\nParameters:\n')
        for param, value in params.items():
            f.write(f'{param}: {value}\n')

    original_stdout = sys.stdout
    with open(os.path.join(results_dir, f's{subject_id}_summary.txt'), "w") as f:
        sys.stdout = f
        model.show()
    sys.stdout = original_stdout
    
except Exception as e:
    error_msg = f"Error processing subject {subject_id}: {str(e)}"
    # Save error information
    with open(os.path.join(results_dir, f's{subject_id}_error.txt'), 'w') as f:
        f.write(error_msg)