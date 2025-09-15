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
    df = pd.read_csv(f'../../simulated_data/{args.subject_id}.csv') 
    # remove trials with no free choice
    df = df.dropna(subset=['freeChoice'])
    # update variable names to match existing code
    df.rename(columns={'RTs': 'RT', 'signal1Onsets': 'signal1_onset', 'noise2Onsets': 'noise2_onset', 'signal2Onsets': 'signal2_onset'}, inplace=True)
    # convert RTs & changepoints to seconds
    df[['RT', 'signal1_onset', 'noise2_onset', 'signal2_onset']] = df[['RT', 'signal1_onset', 'noise2_onset', 'signal2_onset']] / 60
    # change any unobserved changepoints to 0
    df[['signal1_onset', 'noise2_onset', 'signal2_onset']] = df[['signal1_onset', 'noise2_onset', 'signal2_onset']].fillna(0)
    # convert congruent to a string
    df['congruent'] = df.apply(lambda row: 'neutral' if row['cue'] == 0.5 
                          else ('incongruent' if row['congruent'] == 0 
                               else 'congruent'), axis=1)

    # write out tidied data for sanity checking
    df.to_csv(f'../../simulated_data/{args.subject_id}_tidy.csv', index=False)
    
    # Create sample
    sample = pyddm.Sample.from_pandas_dataframe(
        df, rt_column_name='RT', choice_column_name='freeChoice'
    )
    
    # initialize model
    model = pyddm.gddm(
        drift = drift,
        starting_position = 0,
        bound="B",
        T_dur = 4.1,
        nondecision=0,
        parameters={'B': (0.5, 35), 
                    'n1_weak': (0, 20), 'n1_strong': (0,20), 'n1_neut': (0,20),
                    's1_cong_weak': (0, 20), 's1_cong_strong': (0, 20), 's1_incong_weak': (0, 20), 's1_incong_strong': (0, 20), 's1_neut': (0, 20),
                    'n2_cong_weak': (0, 20), 'n2_cong_strong': (0,20), 'n2_incong_weak': (0,20), 'n2_incong_strong': (0,20), 'n2_neut': (0,20),
                    's2_cong_weak': (0, 20), 's2_cong_strong': (0, 20), 's2_incong_weak': (0, 20), 's2_incong_strong': (0, 20), 's2_neut': (0,20)},
        conditions = ['congruent', 'cue', 'signal1_onset', 'noise2_onset', 'signal2_onset']
    )

    # fit
    model.fit(sample, verbose=True)
    
    # Gather results
    loss = model.get_fit_result().value()
    params = model.parameters()
    
    # Save text results (basic summary)
    with open(os.path.join(results_dir, f's{subject_id}_{coherence}coh_results.txt'), 'w') as f:
        f.write(f'Subject: {subject_id}\n Coherence: {coherence}\n Loss: {loss}\nParameters:\n')
        for param, value in params.items():
            f.write(f'{param}: {value}\n')

    original_stdout = sys.stdout
    with open(os.path.join(results_dir, f's{subject_id}_{coherence}coh_summary.txt'), "w") as f:
        sys.stdout = f
        model.show()
    sys.stdout = original_stdout
    
except Exception as e:
    error_msg = f"Error processing subject {subject_id}: {str(e)}"
    # Save error information
    with open(os.path.join(results_dir, f's{subject_id}_error.txt'), 'w') as f:
        f.write(error_msg)