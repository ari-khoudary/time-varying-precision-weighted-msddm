import pyddm
import pyddm.plot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import re
import argparse
import glob

def load_subject_data(data_dir, subID):
    """Load and preprocess subject data from individual CSV file."""
    csv_file = os.path.join(data_dir, f'sub{subID}_tidy.csv')
    if not os.path.exists(csv_file):
        return None
    df = pd.read_csv(csv_file)
    return df

def extract_parameters(summary_file):
    """Extract fitted parameters and metadata from summary file."""
    with open(summary_file, 'r') as f:
        content = f.read()
    
    loss = float(re.search(r'Loss function value: ([\d\.]+)', content).group(1))
    cue = float(re.search(r'Cue: ([\d\.]+)', content).group(1))
    coherence = float(re.search(r'Coherence: ([\d\.]+)', content).group(1))
    sample_size = int(re.search(r'- samplesize: (\d+)', content).group(1))
    n_params = int(re.search(r'- nparams: (\d+)', content).group(1))
    
    params = {}
    for param_name, param_value in re.findall(r'^\s*- ([^:]+): ([\d\.\-]+)', content, re.MULTILINE):
        if param_name not in ['samplesize', 'nparams', 'noise', 'x0', 'nondectime', 'umixturecoef']:
            params[param_name] = float(param_value)
    
    return params, loss, cue, coherence, sample_size, n_params

def get_model(cue, params):
    """Create pyDDM model with fitted parameters."""
    
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
                    n1_cong, n1_incong, s1_cong, s1_incong, 
                    n2_cong, n2_incong, s2_cong, s2_incong):
        # drift rate during first noise period
        if t < signal1_onset:
            if congruent == 'congruent': 
                return n1_cong
            else:  # incongruent
                return -n1_incong

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

    drift_fn = drift_neutral if cue == 0.5 else drift_biased
    conditions = ['signal1_onset', 'noise2_onset', 'signal2_onset']
    if cue != 0.5:
        conditions = ['congruent'] + conditions
    
    return pyddm.gddm(drift=drift_fn, starting_position=0, bound="B", 
                      T_dur=4.3, nondecision=0, parameters=params, conditions=conditions)

def plot_fit(model, sample, subID, cue, coherence, params, loss, sample_size, n_params, 
             subject_df, plots_dir):
    """Create and save diagnostic plot."""
    
    # Build parameter text 
    param_text = f"Loss: {loss:.4f}\n"
    param_text += '\n'.join([f"{k}: {v:.3f}" for k, v in params.items()])
    param_text += f"\nSample size: {sample_size}\nn_params: {n_params}"
    
    # Get timing info
    df = subject_df.copy()
    df[['noise2_onset', 'signal2_onset']] = df[['noise2_onset', 'signal2_onset']].replace(0, np.nan)
    s1_onset = df['signal1_onset'].median()
    n2_onset = df['noise2_onset'].median()
    s2_onset = df['signal2_onset'].median()
    median_rt = df.loc[df['RT'] > s2_onset, 'RT'].median()
    max_rt = df.loc[df['RT'] > s2_onset, 'RT'].max()
    
    # Create plot
    fig = plt.figure(figsize=(12, 10))
    pyddm.plot.plot_fit_diagnostics(model, sample, fig, data_dt=0.05)
    
    # Add timing regions to subplots
    for i, ax in enumerate(fig.get_axes()):
        ylim = ax.get_ylim()
        ax.axvspan(s1_onset, n2_onset, alpha=0.25, color='lightblue', zorder=0)
        ax.axvspan(s2_onset, median_rt, alpha=0.25, color='lightblue', zorder=0)
        if max_rt > median_rt:
            ax.axvspan(s2_onset, max_rt, alpha=0.25, color='lightpink', zorder=0)
        
        if i > 0:  # Add labels to all but first subplot
            y_pos = ylim[0] + (ylim[1] - ylim[0]) * 0.07
            ax.text(s1_onset/2, y_pos, 'median noise 1', ha='center', va='top', 
                   backgroundcolor='white', alpha=0.7, zorder=5)
            ax.text((s1_onset+n2_onset)/2, y_pos, 'median signal 1', ha='center', va='top',
                   backgroundcolor='white', alpha=0.7, zorder=5)
            ax.text((n2_onset+s2_onset)/2, y_pos, 'median noise2', ha='center', va='top',
                   backgroundcolor='white', alpha=0.7, zorder=5)
            ax.text((s2_onset+median_rt)/2, y_pos, 'median signal 2', ha='center', va='top',
                   backgroundcolor='white', alpha=0.7, zorder=5)
            if max_rt > median_rt:
                ax.text((s2_onset+max_rt)/2, y_pos, 'max signal 2', ha='center', va='top',
                       backgroundcolor='white', alpha=0.7, zorder=5)
    
    # Add parameter text
    ax = fig.get_axes()[-1]
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    ax.text(xlim[0], ylim[0] - (ylim[1] - ylim[0]) * 0.15, param_text,
            ha='left', va='top', fontsize=14, fontfamily='monospace', transform=ax.transData)
    
    plt.tight_layout(rect=[0, 0.25, 1, 0.95])
    plt.suptitle(f"Subject {subID}: Cue = {cue}, Coherence = {coherence}", 
                 fontsize=16, y=1.025)
    plt.savefig(os.path.join(plots_dir, f's{subID}_{cue}cue_{coherence}coh.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()

def process_subjects(results_dir, data_dir, plot_results=False, subjects_to_process=None):
    """Process all subjects and generate fits.csv files."""

    # Find all summary files across cue subdirectories
    cue_dirs = ['0.5cue', '0.65cue', '0.8cue']
    all_summary_files = []
    
    for cue_dir in cue_dirs:
        cue_path = os.path.join(results_dir, cue_dir)
        if os.path.exists(cue_path):
            summary_files = glob.glob(os.path.join(cue_path, 'ssub*_*coh_summary.txt'))
            all_summary_files.extend([(cue_dir, os.path.basename(f)) for f in summary_files])
    
    # Extract unique subjects
    all_subjects = sorted(set(int(re.match(r'ssub(\d+)_', f[1]).group(1)) 
                              for f in all_summary_files if re.match(r'ssub(\d+)_', f[1])))
    
    subjects = [s for s in all_subjects if s in subjects_to_process] if subjects_to_process else all_subjects
    print(f"Processing {len(subjects)} subjects across {len(cue_dirs)} cue directories")
    
    if plot_results:
        plots_dir = os.path.join(results_dir, "plots")
        os.makedirs(plots_dir, exist_ok=True)
    
    # Process each subject
    all_param_data = []
    for subID in subjects:
        print(f"Processing subject {subID}...")
        subject_df = load_subject_data(data_dir, subID)
        if subject_df is None:
            print(f"  No data file found, skipping")
            continue
        
        # Process all summary files for this subject
        subject_files = [(cue_dir, fname) for cue_dir, fname in all_summary_files 
                        if fname.startswith(f'ssub{subID}_')]
        
        for cue_dir, summary_filename in subject_files:
            try:
                summary_path = os.path.join(results_dir, cue_dir, summary_filename)
                params, loss, cue, coherence, sample_size, n_params = extract_parameters(summary_path)
                
                # Filter data for this cue and coherence
                cue_df = subject_df[(subject_df['cue'] == cue) & (subject_df['coherence'] == coherence)]
                if cue_df.empty:
                    continue
                
                # Store parameters
                param_entry = {'subID': subID, 'cue': cue, 'coherence': coherence, 
                              'model': 'diff noise1, min drift = -1', 'loss': loss, 
                              'sample_size': sample_size, 'n_params': n_params, **params}
                all_param_data.append(param_entry)
                
                # Generate plot
                if plot_results:
                    model = get_model(cue, params)
                    sample = pyddm.Sample.from_pandas_dataframe(
                        cue_df, rt_column_name='RT', choice_column_name='freeChoice')
                    plot_fit(model, sample, subID, cue, coherence, params, loss, 
                            sample_size, n_params, cue_df, plots_dir)
                    
            except Exception as e:
                print(f"  Error processing {summary_filename}: {e}")
    
    if not all_param_data:
        return None, None
    
    # Save wide format
    param_df = pd.DataFrame(all_param_data)
    param_df.to_csv(os.path.join(results_dir, 'fits.csv'), index=False)
    print(f"Wide-format data saved to {os.path.join(results_dir, 'fits.csv')}")
    
    # Save long format
    long_data = []
    for row in all_param_data:
        for key, value in row.items():
            if key not in ['subID', 'cue', 'coherence', 'model', 'loss', 'sample_size', 'n_params'] \
               and isinstance(value, (int, float)):
                long_data.append({'subID': row['subID'], 'cue': row['cue'], 
                                 'coherence': row['coherence'], 'model': row['model'], 
                                 'parameter': key, 'fitted_value': value})
    
    long_df = pd.DataFrame(long_data)
    long_df.to_csv(os.path.join(results_dir, 'fits_long.csv'), index=False)
    print(f"Long-format data saved to {os.path.join(results_dir, 'fits_long.csv')}")
    
    return param_df, long_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process DDM results and generate fits.csv')
    parser.add_argument('--results_dir', type=str, default='./',
                       help='Directory containing results and summary files')
    parser.add_argument('--data_dir', type=str, default='../../../simulated_data/parallel_data/',
                       help='Directory containing subject CSV files (sub*_tidy.csv)')
    parser.add_argument('--plot_results', action='store_true',
                       help='Generate diagnostic plots for each subject')
    parser.add_argument('--subjects', type=int, nargs='+',
                       help='Specific subject IDs to process (default: all)')
    
    args = parser.parse_args()
    
    print(f"Data directory: {args.data_dir}")
    print(f"Results directory: {args.results_dir}")
    print(f"Generate plots: {args.plot_results}")
    if args.subjects:
        print(f"Processing subjects: {args.subjects}")
    
    process_subjects(args.results_dir, args.data_dir, args.plot_results, args.subjects)
    print("Processing complete!")