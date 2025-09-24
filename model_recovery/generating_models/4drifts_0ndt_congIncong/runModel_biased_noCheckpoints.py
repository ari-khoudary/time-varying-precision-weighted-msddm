# initialize
from asyncio import Condition
import enum
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

# Get the cue, drifts, and threshold from command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--cue', type=float, required=True)
parser.add_argument('--n1', type=float, required=True)
parser.add_argument('--s1', type=float, required=True)
parser.add_argument('--s1_incong', type=float, required=True)
parser.add_argument('--n2', type=float, required=True)
parser.add_argument('--n2_incong', type=float, required=True)
parser.add_argument('--s2', type=float, required=True)
parser.add_argument('--s2_incong', type=float, required=True)
parser.add_argument('--b', type=float, required=True)
args = parser.parse_args()

cue = args.cue
n1 = args.n1
s1_cong = args.s1
s1_incong = args.s1_incong
n2_cong = args.n2
n2_incong = args.n2_incong
s2_cong = args.s2
s2_incong = args.s2_incong
thresh = args.b

# Create results directory if it doesn't exist
results_dir = 'results/'
os.makedirs(results_dir, exist_ok=True)

def drift_biased(t, congruent, signal1_onset, noise2_onset, signal2_onset,
                 n1_cong, s1_cong, s1_incong, 
                 n2_cong, n2_incong, s2_cong, s2_incong):
    # drift rate during first noise period
    if t < signal1_onset:
        if congruent == 'congruent': 
            return n1_cong
        else:  # incongruent
            return -n1_cong

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
nTrial = 1000
nCongTrials = round(cue * nTrial)
nIncongTrials = nTrial - nCongTrials

quantiles = np.array([0, 0.05, 0.25, 0.50, .75, 0.95, 1])
bin_proportion = np.diff(quantiles)
trial_bins_congruent = np.round(bin_proportion * nCongTrials).astype(int)
trial_bins_incongruent = np.round(bin_proportion * nIncongTrials).astype(int)


try:    
    # initialize generating model
    sim_model = pyddm.gddm(
            drift = drift_biased,
            starting_position = 0,
            bound="B",
            T_dur = 4.3,
            nondecision=0,
            parameters={'B': thresh, 
                        'n1': n1,
                        's1_cong': s1_cong, 's1_incong': s1_incong,
                        'n2_cong': n2_cong, 'n2_incong': n2_incong,
                        's2_cong': s2_cong, 's2_incong': s2_incong},
            conditions = ['congruent', 'signal1_onset', 'noise2_onset', 'signal2_onset'])

    # simulate data 
    onsets_list = [onsets_5q, onsets_25q, onsets_50q, onsets_75q, onsets_95q, onsets_100q]
    congruence_list = ['congruent', 'incongruent']
    samples = []
    # loop over all conditions
    for congruent in congruence_list:
        trial_bins = trial_bins_congruent if congruent == "congruent" else trial_bins_incongruent
        
        for i, onsets in enumerate(onsets_list):
            onsets_with_congruent = onsets.copy()
            onsets_with_congruent['congruent'] = congruent
            print(onsets_with_congruent)
            samp = sim_model.solve(conditions=onsets_with_congruent).sample(trial_bins[i])
            samples.append(samp)

    # concatenate into one sample
    sample = sum(samples)

    # specify fitting model
    fitting_model = pyddm.gddm(
        drift = drift_biased,
        starting_position = 0,
        bound="B",
        T_dur = 4.3,
        nondecision=0,
        parameters={'B': (1, 15), 
                    'n1_cong': (0, 10), 'n1_incong': (0, 10),
                    's1_cong': (0, 10), 's1_incong': (0, 10),
                    'n2_cong': (0, 10), 'n2_incong': (0, 10),
                    's2_cong': (0, 10), 's2_incong': (0, 10)},
        conditions = ['congruent', 'signal1_onset', 'noise2_onset', 'signal2_onset'])
        

    # fit
    fitting_model.fit(sample, verbose=True)
    
    # Gather results
    loss = fitting_model.get_fit_result().value()
    params = fitting_model.parameters()
    
    # Save text results (basic summary)
    with open(os.path.join(results_dir, f'{cue}cue_{thresh}thresh_{n1}n1_{s1_cong}s1_{s1_incong}s1Incong_{n2_cong}n2_{n2_incong}n2Incong_{s2_cong}s2_{s2_incong}s2Incong_fits.txt'), 'w') as f:
        f.write(f'\nLoss: {loss}\nParameters:\n')
        for param, value in params.items():
            f.write(f'{param}: {value}\n')
    original_stdout = sys.stdout
    
    with open(os.path.join(results_dir, f'{cue}cue_{thresh}thresh_{n1}n1_{s1_cong}s1_{s1_incong}s1Incong_{n2_cong}n2_{n2_incong}n2Incong_{s2_cong}s2_{s2_incong}s2Incong_summary.txt'), 'w') as f:
        sys.stdout = f
        f.write(f'\nGROUND TRUTH: \n')
        sim_model.show()
        f.write(f'\n\n FITTED VALUES: \n')
        fitting_model.show()
    sys.stdout = original_stdout
    
except Exception as e:
    error_msg = f"Error processing job"
    # Save error information with cue and coherence in filename
    with open(os.path.join(results_dir, f'{cue}cue_{thresh}thresh_{n1}n1_{s1_cong}s1_{s1_incong}s1Incong_{n2_cong}n2_{n2_incong}n2Incong_{s2_cong}s2_{s2_incong}s2Incong_error.txt'), 'w') as f:
        f.write(error_msg)