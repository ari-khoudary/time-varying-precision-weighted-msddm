# initialize
import pyddm
import pandas as pd
import numpy as np
from pyddm import Sample
import os
import sys
import matplotlib.pyplot as plt
import gc

# Get number of CPUs from SLURM environment variable
n_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
pyddm.set_N_cpus(n_cpus)

# Create results directory if it doesn't exist
results_dir = f'results/'
os.makedirs(results_dir, exist_ok=True)

# Process command line arguments
if len(sys.argv) < 3:
    print("Error: Subject ID and cue value not provided. Use as: python fit_model.py [subject_id] [cue_value]")
    sys.exit(1)

subject_id = sys.argv[1]
cue_value = float(sys.argv[2])
subject_id_str = str(subject_id)
subject_id_int = int(subject_id) 

print(f"Fitting model for subject {subject_id}, cue {cue_value}")

# Define drift functions
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
                 n1_biased, s1_incong, s1_cong, 
                 n2_cong, n2_incong, s2_cong, s2_incong):
    """Drift function for biased cue (>0.5)"""
    # drift rate during first noise period
    if t < signal1_onset:
        if congruent == 'congruent': 
            return n1_biased
        else:  # incongruent
            return -n1_biased

    # drift rates during first signal period
    elif t >= signal1_onset and t < noise2_onset:
        if congruent == 'congruent':
            return s1_cong
        else:  # incongruent
            return -s1_incong

    # drift rates during the second noise period
    elif t >= noise2_onset and t < signal2_onset:
        if congruent == 'congruent':
            return n2_cong
        else:  # incongruent
            return -n2_incong

    # drift rates during the second signal period
    elif t >= signal2_onset:
        if congruent == 'congruent':
            return s2_cong
        else:  # incongruent
            return -s2_incong

try:
    # Load and prepare data
    df = pd.read_csv('inference_test.csv')
    df = df.dropna(subset=['RT'])
    df = df.rename(columns={'trueCongruence': 'congruent'})
    df[['signal1_onset', 'noise2_onset', 'signal2_onset']] = df[['signal1_onset', 'noise2_onset', 'signal2_onset']].fillna(0)
    
    # Filter for specific subject and cue
    subject_df = df[(df['subID'] == subject_id_str) | 
                    (df['subID'] == subject_id_int)].copy()
    
    if len(subject_df) == 0:
        error_msg = f"Error: No data found for subject {subject_id}"
        cue_str = f"{cue_value:.2f}".replace('.', 'p')
        with open(os.path.join(results_dir, f's{subject_id}_cue{cue_str}_error.txt'), 'w') as f:
            f.write(error_msg)
        sys.exit(1)
    
    # Filter for specific cue value
    cue_df = subject_df[subject_df['trueCue'] == cue_value].copy()
    
    if len(cue_df) == 0:
        error_msg = f"Error: No data found for subject {subject_id} with cue {cue_value}"
        cue_str = f"{cue_value:.2f}".replace('.', 'p')
        with open(os.path.join(results_dir, f's{subject_id}_cue{cue_str}_error.txt'), 'w') as f:
            f.write(error_msg)
        sys.exit(1)
    
    print(f"Found {len(cue_df)} trials for subject {subject_id}, cue {cue_value}")
    
    # Create sample
    sample = pyddm.Sample.from_pandas_dataframe(
        cue_df, rt_column_name='RT', choice_column_name='accuracy'
    )
    
    # Initialize model based on cue type
    if cue_value == 0.5:  # neutral cue
        model = pyddm.gddm(
            drift=drift_neutral,
            starting_position=0,
            bound="B",
            T_dur=4.3,
            nondecision=0,
            parameters={'B': (1, 15), 
                        'n1_neut': (0, 10), 's1_neut': (0, 10), 
                        'n2_neut': (0, 10), 's2_neut': (0, 10)},
            conditions=['signal1_onset', 'noise2_onset', 'signal2_onset']
        )
    else:  # biased cues (0.65 or 0.8)
        model = pyddm.gddm(
            drift=drift_biased,
            starting_position=0,
            bound="B",
            T_dur=4.3,
            nondecision=0,
            parameters={'B': (1, 15), 
                        'n1_biased': (0, 10), 
                        's1_cong': (0, 10), 's1_incong': (0, 10),
                        'n2_cong': (0, 10), 'n2_incong': (0, 10),
                        's2_cong': (0, 10), 's2_incong': (0, 10)},
            conditions=['congruent', 'signal1_onset', 'noise2_onset', 'signal2_onset']
        )
    
    # Fit the model
    print("Starting model fitting...")
    model.fit(sample, verbose=True)
    
    # Gather results
    loss = model.get_fit_result().value()
    params = model.parameters()
    
    # Create filename-safe cue string
    cue_str = f"{cue_value:.2f}".replace('.', 'p')
    
    # Save parameter results
    with open(os.path.join(results_dir, f's{subject_id}_cue{cue_str}_results.txt'), 'w') as f:
        f.write(f'Subject: {subject_id}\nCue: {cue_value}\nTrials: {len(cue_df)}\nLoss: {loss}\nParameters:\n')
        for param, value in params.items():
            f.write(f'{param}: {value}\n')
    
    # Save model summary
    original_stdout = sys.stdout
    with open(os.path.join(results_dir, f's{subject_id}_cue{cue_str}_summary.txt'), "w") as f:
        sys.stdout = f
        print(f'Subject: {subject_id}, Cue: {cue_value}, Trials: {len(cue_df)}')
        model.show()
    sys.stdout = original_stdout
    
    print(f"Successfully fitted model for subject {subject_id}, cue {cue_value}")
    print(f"Loss: {loss:.4f}")
    
    # Clean up memory
    gc.collect()

except Exception as e:
    error_msg = f"Error processing subject {subject_id}, cue {cue_value}: {str(e)}"
    print(error_msg)
    cue_str = f"{cue_value:.2f}".replace('.', 'p')
    with open(os.path.join(results_dir, f's{subject_id}_cue{cue_str}_error.txt'), 'w') as f:
        f.write(error_msg)
