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
parser.add_argument('--cue', type=float, required=True)
args = parser.parse_args()

subject_id = args.subject_id
coherence = args.coherence
cue = args.cue

# Create results directory if it doesn't exist
results_dir = f'results/{cue}cue/'
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

def drift_biased(t, congruent, signal1_onset, noise2_onset, signal2_onset,
                 n1_biased, s1_cong, s1_incong, 
                 n2_cong, n2_incong, s2_cong, s2_incong):
    # drift rate during first noise period
    if t < signal1_onset:
        if congruent == 'congruent':
            return n1_biased
        else:
            return -n1_biased

    # drift rates during first signal period
    elif t >= signal1_onset and t < noise2_onset:
        if congruent == 'congruent':
            return s1_cong
        else:  # incongruent
            return s1_incong
    # drift rates during the second noise period
    elif t >= noise2_onset and t < signal2_onset:
        if congruent == 'congruent':
            return n2_cong
        else:  # incongruent
            return n2_incong
    # drift rates during the second signal period
    elif t >= signal2_onset:
        if congruent == 'congruent':
            return s2_cong
        else:  # incongruent
            return s2_incong

try:
    # Check if tidy version of dataframe already exists
    tidy_file_path = f'../../simulated_data/parallel_data/{subject_id}_tidy.csv'
    
    if os.path.exists(tidy_file_path):
        # Load the existing tidy dataframe
        df = pd.read_csv(tidy_file_path)
    else:
        # Load and tidy the raw data
        df = pd.read_csv(f'../../simulated_data/{subject_id}.csv') 
        # remove trials with no free choice
        df = df.dropna(subset=['freeChoice'])
        # update variable names to match existing code
        df.rename(columns={'RTs': 'RT', 'signal1Onsets': 'signal1_onset', 'noise2Onsets': 'noise2_onset', 'signal2Onsets': 'signal2_onset'}, inplace=True)
        # convert RTs & changepoints to seconds
        df[['RT', 'signal1_onset', 'noise2_onset', 'signal2_onset']] = df[['RT', 'signal1_onset', 'noise2_onset', 'signal2_onset']] / 60
        # change any unobserved changepoints to 0
        for prefix in ['signal1', 'noise2', 'signal2']:
            df[f'{prefix}_onset'] = df[f'{prefix}Frames_obs'].fillna(0)
        # add congruence as a string        
        df['congruent'] = df.apply(lambda row: 'neutral' if row['cue'] == 0.5 
                              else ('incongruent' if row['congruent'] == 0 
                                   else 'congruent'), axis=1)
        # filter to include only memoryThinning values in the (broad) theta range
        # df = df[df['memoryThinning'] < 25]
        # save tidied data
        df.to_csv(tidy_file_path, index=False)
    
    # filter to trials with specified coherence
    #df['coherence'] = df['coherence'].round(2)
    df = df[(df['coherence'] == coherence) & (df['cue'] == cue)]
    
    # Create sample
    sample = pyddm.Sample.from_pandas_dataframe(
        df, rt_column_name='RT', choice_column_name='freeChoice')
    
   # initialize model with parameters depending on cue level
    if cue == 0.5: 
        model = pyddm.gddm(
            drift = drift_neutral,
            starting_position = 0,
            bound="B",
            T_dur = 4.3,
            nondecision=0,
            parameters={'B': (1, 15), 
                        'n1_neut': (-10, 10), 's1_neut': (-10, 10), 
                        'n2_neut': (-10, 10), 's2_neut': (-10, 10)},
            conditions = ['signal1_onset', 'noise2_onset', 'signal2_onset']
        )
    else:  # biased cues
        model = pyddm.gddm(
            drift = drift_biased,
            starting_position = 0,
            bound="B",
            T_dur = 4.3,
            nondecision=0,
            parameters={'B': (1, 15), 
                        'n1_biased': (-10, 10), 
                        's1_cong': (-1, 10), 's1_incong': (-10, 1),
                        'n2_cong': (-1, 10), 'n2_incong': (-10, 1),
                        's2_cong': (-1, 10), 's2_incong': (-10, 1)},
            conditions = ['congruent', 'signal1_onset', 'noise2_onset', 'signal2_onset']
        )

    # fit
    model.fit(sample, verbose=True)
    
    # Gather results
    loss = model.get_fit_result().value()
    params = model.parameters()
    
    # Save text results (basic summary)
    with open(os.path.join(results_dir, f'{subject_id}_{coherence}coh_results.txt'), 'w') as f:
        f.write(f'Subject: {subject_id}\nCoherence: {coherence}\nCue: {cue}\nLoss: {loss}\nParameters:\n')
        for param, value in params.items():
            f.write(f'{param}: {value}\n')
    original_stdout = sys.stdout
    
    with open(os.path.join(results_dir, f'{subject_id}_{coherence}coh_summary.txt'), "w") as f:
        sys.stdout = f
        print(f'Subject: {subject_id}, Coherence: {coherence}, Cue: {cue}')
        model.show()
    sys.stdout = original_stdout
    
except Exception as e:
    error_msg = f"Error processing subject {subject_id}: {str(e)}"
    # Save error information with cue and coherence in filename
    with open(os.path.join(results_dir, f'{subject_id}_{coherence}coh_error.txt'), 'w') as f:
        f.write(error_msg)
