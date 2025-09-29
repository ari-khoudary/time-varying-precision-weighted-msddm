import pyddm
from pyddm import Sample
import pyddm.plot
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re

def extract_parameters_from_summary(summary_file):
    """Extract parameters from a summary file."""
    if not os.path.exists(summary_file):
        return None, None, None, None, None
        
    with open(summary_file, 'r') as f:
        content = f.read()
    
    # Extract loss
    loss_match = re.search(r'Loss function value: ([\d\.]+)', content)
    if not loss_match:
        return None, None, None, None, None
    loss = float(loss_match.group(1))
    
    # Extract cue value from the first line
    cue_match = re.search(r'Cue: ([\d\.]+)', content)
    cue = float(cue_match.group(1)) if cue_match else None
    
    # Extract sample size
    sample_match = re.search(r'- samplesize: (\d+)', content)
    sample_size = int(sample_match.group(1)) if sample_match else None
        
    # Extract number of parameters
    nparams_match = re.search(r'- nparams: (\d+)', content)
    n_params = int(nparams_match.group(1)) if nparams_match else None
    
    # Extract fitted parameters - looking for lines with "- paramname: value"
    params = {}
    fitted_param_lines = re.findall(r'^\s*- ([^:]+): ([\d\.\-]+)', content, re.MULTILINE)
    
    for param_name, param_value in fitted_param_lines:
        # Skip non-parameter lines like samplesize, nparams, etc.
        if param_name not in ['samplesize', 'nparams', 'noise', 'x0', 'nondectime', 'umixturecoef']:
            params[param_name] = float(param_value)
    
    return params, loss, cue, sample_size, n_params

def create_long_format_df(all_param_data):
    """Create long-format dataframe where each row is a subject-parameter combination."""
    long_data = []
    
    for subject_data in all_param_data:
        subID = subject_data['subID']
        cue = subject_data['cue']
        model = subject_data['model']
        condition = subject_data.get('condition', '')
        
        for key, value in subject_data.items():
            if key not in ['subID', 'cue', 'model', 'loss', 'condition', 'sample_size', 'n_params']:
                if isinstance(value, (int, float)):  # Only include numeric parameters
                    long_data.append({
                        'subID': subID,
                        'cue': cue,
                        'model': model,
                        'condition': condition,
                        'parameter': key,
                        'fitted_value': value
                    })
    
    return pd.DataFrame(long_data)

def get_model_with_params(cue, params):
    """Create a pyDDM model with fitted parameters."""
    
    # Define drift functions (same as in fit_model.py)
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

    # Create model based on cue type
    if cue == 0.5:  # neutral cue
        model = pyddm.gddm(
            drift=drift_neutral,
            starting_position=0,
            bound="B",
            T_dur=4.3,
            nondecision=0,
            parameters=params,
            conditions=['signal1_onset', 'noise2_onset', 'signal2_onset']
        )
    else:  # biased cues (0.65 or 0.8)
        model = pyddm.gddm(
            drift=drift_biased,
            starting_position=0,
            bound="B",
            T_dur=4.3,
            nondecision=0,
            parameters=params,
            conditions=['congruent', 'signal1_onset', 'noise2_onset', 'signal2_onset']
        )
    
    return model

def create_diagnostic_plot(model, sample, subID, cue, params, loss, sample_size, n_params, subject_df, plots_dir, model_name, condition):
    """Create diagnostic plot for a specific subject and cue."""
    
    # Create parameter text
    param_text = f"Loss: {loss:.4f}\n"
    for key, value in params.items():
        param_text += f"{key}: {value:.3f}\n"
    
    if sample_size is not None:
        param_text += f"Sample size: {sample_size}\n"
        
    if n_params is not None:
        param_text += f"n_params: {n_params}"

    # Extract timing information
    subject_df = subject_df.copy()
    subject_df[['noise2_onset', 'signal2_onset']] = subject_df[['noise2_onset', 'signal2_onset']].replace(0, np.nan)
    signal1_onset = subject_df['signal1_onset'].median()
    noise2_onset = subject_df['noise2_onset'].median()
    signal2_onset = subject_df['signal2_onset'].median()
    median_rt = subject_df.loc[subject_df['RT'] > subject_df['signal2_onset'], 'RT'].median()
    max_rt = subject_df.loc[subject_df['RT'] > subject_df['signal2_onset'], 'RT'].max()

    # Create figure
    fig = plt.figure(figsize=(12, 10))
    pyddm.plot.plot_fit_diagnostics(model, sample, fig, data_dt=0.05)

    # Add noise regions to each subplot
    second_ax = False
    for ax in fig.get_axes():
        try:
            ylim = ax.get_ylim()
            xlim = ax.get_xlim()
            # Add boxes
            ax.axvspan(signal1_onset, noise2_onset, alpha=0.25, color='lightblue', zorder=0)
            ax.axvspan(signal2_onset, median_rt, alpha=0.25, color='lightblue', zorder=0)
            if max_rt > median_rt:
                ax.axvspan(signal2_onset, max_rt, alpha=0.25, color='lightpink', zorder=0)
            # Add labels
            if second_ax:
                label_y_pos = ylim[0] + (ylim[1] - ylim[0]) * 0.07
                # Label signal 1
                ax.text((signal1_onset+noise2_onset)/2, label_y_pos, 'median signal 1', 
                        ha='center', va='top', backgroundcolor='white', alpha=0.7, zorder=5)
                # Label median signal 2
                ax.text((signal2_onset+median_rt)/2, label_y_pos, 'median signal 2', 
                        ha='center', va='top', backgroundcolor='white', alpha=0.7, zorder=5)
                # Label max signal 2
                if max_rt > median_rt:
                    ax.text((signal2_onset+max_rt)/2, label_y_pos, 'max signal 2', 
                            ha='center', va='top', backgroundcolor='white', alpha=0.7, zorder=5)
                # Label noise 1 
                ax.text((signal1_onset+0)/2, label_y_pos, 'median noise 1', 
                        ha='center', va='top', backgroundcolor='white', alpha=0.7, zorder=5)
                # Label noise 2 
                ax.text((noise2_onset+signal2_onset)/2, label_y_pos, 'median noise2', 
                        ha='center', va='top', backgroundcolor='white', alpha=0.7, zorder=5)
            # Force redraw
            ax.figure.canvas.draw()
            second_ax = True
        except Exception as e:
            print(f"Warning: Could not add noise regions to a subplot: {e}")

    # Add the parameter text to the bottom right corner
    ax = fig.get_axes()[-1]
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    # Add the text below the x-axis
    ax.text(xlim[0],
            ylim[0] - (ylim[1] - ylim[0]) * 0.15,
            param_text,
            ha='left',
            va='top',
            fontsize=14, fontfamily='monospace',
            transform=ax.transData)

    # Ensure there's enough room at the bottom
    plt.tight_layout(rect=[0, 0.25, 1, 0.95])
    
    # Create title
    cue_str = f"{cue:.2f}".replace('.', 'p')
    plt.suptitle(f"Subject {subID} Fit: {condition} - Cue {cue} - {model_name}", fontsize=16, y=1.025)
    
    # Save the figure
    plt.savefig(os.path.join(plots_dir, f's{subID}_cue{cue_str}_fit.png'), dpi=300, bbox_inches='tight')
    plt.close()

def process_subjects(results_dir, data_csv, plot_results=False, subjects_to_process=None):
    """Process subjects data and generate fits.csv and optional plots.
    
    Parameters:
        results_dir (str): Directory containing results files
        data_csv (str): Path to the CSV file with subject data
        plot_results (bool): Whether to generate diagnostic plots for each subject
        subjects_to_process (list or None): List of subject IDs to process; if None, process all subjects
    """
    # Extract model name from results directory path
    model_name_match = re.search(r'cluster_fitting/([^/]+)/results', results_dir)
    model_name = model_name_match.group(1) if model_name_match else "unknown_model"
    
    # Read the data file
    print(f"Reading data from {data_csv}...")
    df = pd.read_csv(data_csv)
    df = df.dropna(subset=['RT'])
    df = df.rename(columns={'trueCongruence': 'congruent'})
    df[['signal1_onset', 'noise2_onset', 'signal2_onset']] = df[['signal1_onset', 'noise2_onset', 'signal2_onset']].fillna(0)
    
    # Get all subjects from the results directory
    all_summary_files = [f for f in os.listdir(results_dir) if f.endswith('_summary.txt')]
    all_subjects = set()
    for filename in all_summary_files:
        # Extract subject ID from filename like "s78_cue0p65_summary.txt"
        match = re.match(r's(\d+)_cue', filename)
        if match:
            all_subjects.add(int(match.group(1)))
    all_subjects = sorted(all_subjects)
    
    # Determine which subjects to process
    if subjects_to_process is not None:
        subjects_to_process_set = set(subjects_to_process)
        subjects = [s for s in all_subjects if s in subjects_to_process_set]
        if len(subjects) < len(subjects_to_process):
            missing = set(subjects_to_process) - set(subjects)
            print(f"Warning: Some requested subjects were not found in the data: {missing}")
    else:
        subjects = all_subjects
    
    print(f"Will process {len(subjects)} subjects from results directory")
    
    # Create output directory for plots if needed
    if plot_results:
        plots_dir = os.path.join(results_dir, "plots/")
        os.makedirs(plots_dir, exist_ok=True)
    
    # Initialize list to store parameter data
    all_param_data = []
    
    # Process each subject
    for subID in subjects:
        print(f"Processing subject {subID}...")
        
        # Get subject-specific data
        subject_df = df[df['subID'] == subID]
        if subject_df.empty:
            print(f"No data found for subject {subID}, skipping.")
            continue

        # Find all summary files for this subject
        subject_summary_files = [f for f in all_summary_files if f.startswith(f's{subID}_cue')]
        
        if not subject_summary_files:
            print(f"No summary files found for subject {subID}, skipping.")
            continue
        
        # Process each cue for this subject
        for summary_filename in subject_summary_files:
            summary_file = os.path.join(results_dir, summary_filename)
            params, loss, cue, sample_size, n_params = extract_parameters_from_summary(summary_file)
            
            if params is None or cue is None:
                print(f"Could not extract parameters or cue from {summary_filename}, skipping.")
                continue
            
            # Filter data for this specific cue
            cue_df = subject_df[subject_df['trueCue'] == cue]
            if cue_df.empty:
                print(f"No data found for subject {subID}, cue {cue}, skipping.")
                continue
            
            # Store parameter data
            param_entry = {
                'subID': subID, 
                'cue': cue,
                'model': model_name, 
                'loss': loss
            }
            
            # Add condition information
            if 'condition' in subject_df.columns:
                param_entry['condition'] = subject_df['condition'].iloc[0]
            else:
                param_entry['condition'] = 'condition_80cue' if subID < 55 else 'condition_65cue'
            
            if sample_size is not None:
                param_entry['sample_size'] = sample_size
            if n_params is not None:
                param_entry['n_params'] = n_params
                
            # Add fitted parameters
            param_entry.update(params)
            all_param_data.append(param_entry)
            
            # Generate diagnostic plot if requested
            if plot_results:
                try:
                    # Create model with fitted parameters
                    model = get_model_with_params(cue, params)
                    
                    # Create pyDDM sample object
                    sample = pyddm.Sample.from_pandas_dataframe(
                        cue_df, rt_column_name='RT', choice_column_name='accuracy'
                    )
                    
                    # Create diagnostic plot
                    create_diagnostic_plot(
                        model, sample, subID, cue, params, loss, 
                        sample_size, n_params, cue_df, plots_dir, 
                        model_name, param_entry['condition']
                    )
                    
                except Exception as e:
                    print(f"Error creating plot for subject {subID}, cue {cue}: {e}")

    if all_param_data:
        # Create and save parameter dataframe (wide format)
        param_df = pd.DataFrame(all_param_data)
        csv_path = os.path.join(results_dir, 'fits.csv')
        param_df.to_csv(csv_path, index=False)
        print(f"Wide-format parameter data saved to {csv_path}")
        
        # Create and save long-format dataframe
        long_param_df = create_long_format_df(all_param_data)
        long_csv_path = os.path.join(results_dir, 'fits_long.csv')
        long_param_df.to_csv(long_csv_path, index=False)
        print(f"Long-format parameter data saved to {long_csv_path}")
        
        return param_df, long_param_df
    
    return None, None

# Run the analysis when script is executed
if __name__ == "__main__":
    import argparse
    
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Process DDM results and generate fits.csv')
    parser.add_argument('--results_dir', type=str, default='results/',
                        help='Directory containing results and summary files (default: results/)')
    parser.add_argument('--data_csv', type=str, default='inference_test.csv',
                        help='Path to the CSV file with subject data (default: inference_test.csv)')
    parser.add_argument('--plot_results', action='store_true',
                        help='Whether to generate diagnostic plots for each subject')
    parser.add_argument('--subjects', type=int, nargs='+',
                        help='Optional list of specific subject IDs to process (default: process all subjects)')
    
    # Parse arguments
    args = parser.parse_args()
    
    print(f"Processing data from {args.data_csv}...")
    print(f"Using results from {args.results_dir}...")
    print(f"Generating plots: {args.plot_results}")
    
    if args.subjects:
        print(f"Processing only the following subjects: {args.subjects}")
    else:
        print("Processing all subjects")
    
    param_df, long_param_df = process_subjects(
        args.results_dir,
        args.data_csv, 
        args.plot_results,
        subjects_to_process=args.subjects
    )
    
    print("Processing complete!")