import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import os
import glob

from scipy.stats import ttest_rel
from scipy.stats import gaussian_kde

from Bio.PDB import *


plt.rcParams.update({'font.size': 30})

def get_kde_iqr(kde, x_values):
    # Calculate the cumulative density function (CDF)
    cdf = np.cumsum(kde) / np.sum(kde)
    
    # Find the indices for the 25th and 75th percentiles
    q1_index = np.searchsorted(cdf, 0.25)
    q3_index = np.searchsorted(cdf, 0.75)
    
    # Get the corresponding x-values for these percentiles
    q1 = x_values[q1_index]
    q3 = x_values[q3_index]
    
    return q1, q3


# Use glob to find all *_pdb_stats.txt files in the specified directory. This has 
pdb_stats_files = glob.glob('*_pdb_stats.txt')

# Initialize a list to store the DataFrames
pdb_stats_dfs = []

# Loop through each file and read it into a DataFrame
for file in pdb_stats_files:
    try:
        # Read the TXT file into a DataFrame
        df = pd.read_csv(file, sep='\t')
    except:
        continue
    
    # Append the DataFrame to the list
    pdb_stats_dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
all_pdb_stats_data = pd.concat(pdb_stats_dfs, ignore_index=True)

#Clean up files
all_pdb_stats_data.rename(columns={'PDB ID': 'PDB'}, inplace=True)
all_pdb_stats_data['RFree'] = pd.to_numeric(all_pdb_stats_data['RFree'], errors='coerce')

close_resi = pd.read_csv('close_resi.csv', sep=',', dtype={'resi': 'int32', 'chain': 'str', 'PDB': 'str'})

# Use glob to find all *_residue_metrics.txt files in the specified directory
###_____________________________________________###
# These files include information on a residue basis for RSCCS, OPIA, EDIAm, RSR, SRSR

b = glob.glob('*_residue_metrics.txt')
dfs = []

# Loop through each file and read it into a DataFrame
for file in b:
    try:
        # Read the CSV file into a DataFrame
        df = pd.read_csv(file, sep='\t', dtype={'resi': 'int32', 'chain': 'str', 'PDB': 'str'})
    except:
        continue
    
    # Extract the first four digits of the file name for the PDB column
    pdb_id = os.path.basename(file)[:4]  # Get the first four characters from the filename
    
    # Add a new column 'PDB' with the extracted PDB ID
    df['PDB'] = pdb_id
    
    # Append the DataFrame to the list
    dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
residue_metrics = pd.concat(dfs, ignore_index=True)

# Clean up column names
residue_metrics = residue_metrics.rename(columns={'seqNum': 'resi', 'AsymID': 'chain'})
residue_metrics['PDB'] = residue_metrics['PDB'].str.upper()
residue_metrics = residue_metrics.rename(columns={'seqNum': 'resi', 'AsymID': 'chain'})

results = []
for pdb, group in close_resi_metrics.groupby('PDB'):
    # Store the result for the current PDB
    results.append({
        'PDB': pdb,
        'RSCCS_med': (group['RSCCS'].median()),
        'RSCCS_mean': (group['RSCCS'].mean()), 
        'OPIA_med': (group['OPIA'].median()),
        'OPIA_mean': (group['OPIA'].mean()),
        'EDIA_med': (group['EDIAm'].median()),
        'EDIA_mean': (group['EDIAm'].mean()),
        'RSR_med': (group['RSR'].median()),
        'RSR_mean': (group['RSR'].mean()),
        'SRSR_med': (group['SRSR'].median()),
        'SRSR_mean': (group['SRSR'].mean()),
    })

metrics_close = pd.DataFrame(results)


# Merge residue_metrics with close_resi to get all residues
not_close_resi_metrics = dd.merge(
    residue_metrics[['resi', 'chain', 'PDB', 'RSCCS', 'OPIA', 'EDIAm', 'RSR', 'SRSR']], 
    close_resi[['resi', 'chain', 'PDB']], 
    on=['resi', 'chain', 'PDB'], 
    how='left',  # Left join to keep all residues from residue_metrics
    indicator=True  # Add an indicator column to track the join status
)

# Filter out rows where the residue is in close_resi (i.e., where '_merge' column is 'both')
not_close_resi_metrics = not_close_resi_metrics[not_close_resi_metrics['_merge'] == 'left_only']

# Drop the '_merge' column, as it's no longer needed
not_close_resi_metrics = not_close_resi_metrics.drop(columns=['_merge'])

# Drop any rows with NaN values in the metrics columns
not_close_resi_metrics = not_close_resi_metrics.dropna(subset=['RSCCS', 'OPIA', 'EDIAm', 'RSR', 'SRSR'])

# Remove any duplicate rows
not_close_resi_metrics = not_close_resi_metrics.drop_duplicates()

# Initialize results list to store aggregated data
results = []

# Group the data by PDB
for pdb, group in not_close_resi_metrics.groupby('PDB'):
    # Store the result for the current PDB
    results.append({
        'PDB': pdb,
        'RSCCS_med': group['RSCCS'].median(),
        'RSCCS_mean': group['RSCCS'].mean(),
        'OPIA_med': group['OPIA'].median(),
        'OPIA_mean': group['OPIA'].mean(),
        'EDIA_med': group['EDIAm'].median(),
        'EDIA_mean': group['EDIAm'].mean(),
        'RSR_med': group['RSR'].median(),
        'RSR_mean': group['RSR'].mean(),
        'SRSR_med': group['SRSR'].median(),
        'SRSR_mean': group['SRSR'].mean(),
    })

# Convert the results list to a DataFrame
metrics_not_close = pd.DataFrame(results)
metrics_not_close.to_csv('metrics_not_close.csv')                               
# Filter residue_metrics by PDBs present in close_resi_metrics, this reduces the size of residue_metrics
residue_metrics = metrics_not_close[metrics_not_close['PDB'].isin(close_resi_metrics['PDB'])]

# Perform one merge between close_resi_metrics and filtered residue_metrics, including all necessary columns
combined_metrics = pd.merge(metrics_close, metrics_not_close, on=['PDB'], suffixes=('_close', '_all'))
combined_metrics.to_csv('combined_metrics.csv')

# Merge combined_metrics with all_pdb_stats_data (which contains resolution and Rfree) in one step
combined_metrics2 = pd.merge(combined_metrics, all_pdb_stats_data[['PDB', 'Resolution', 'RFree']], on='PDB')


# Plot and compute statistics for each metric
for metric in ['RSCCS_med', 'OPIA_med', 'EDIA_med', 'RSR_med', 'SRSR_med', 'RSCCS_mean', 'OPIA_mean', 'EDIA_mean', 'RSR_mean', 'SRSR_mean']:

    # Compute KDE for the 'All Residues' data
    data_all = combined_metrics[f'{metric}_all'].dropna()
    kde_all = gaussian_kde(data_all)
    x_values_all = np.linspace(min(data_all), max(data_all), 1000)
    y_values_all = kde_all(x_values_all)

    # Compute KDE for the 'Close Residues' data
    data_close = combined_metrics[f'{metric}_close'].dropna()
    kde_close = gaussian_kde(data_close)
    x_values_close = np.linspace(min(data_close), max(data_close), 1000)
    y_values_close = kde_close(x_values_close)

    # Calculate IQR for both distributions
    iqr_all = get_kde_iqr(y_values_all, x_values_all)
    iqr_close = get_kde_iqr(y_values_close, x_values_close)

    # Plotting the KDEs with lineplot
    fig = plt.figure(figsize=(10, 6))
    sns.lineplot(x=x_values_all, y=y_values_all, color='#006C6C', linewidth=2, label='All Residues')
    sns.lineplot(x=x_values_close, y=y_values_close, color='#FF8C00', linewidth=2, label='Close Residues')

    # Shade the IQR regions
    plt.fill_between(x_values_all, 0, y_values_all, where=(x_values_all >= iqr_all[0]) & (x_values_all <= iqr_all[1]),
                 color='#006C6C', alpha=0.3, label='IQR (All Residues)')
    plt.fill_between(x_values_close, 0, y_values_close, where=(x_values_close >= iqr_close[0]) & (x_values_close <= iqr_close[1]),
                 color='#FF8C00', alpha=0.3, label='IQR (Close Residues)')

    # Set labels and legend
    plt.xlabel(f'{metric}', fontsize=16)
    plt.ylabel('Density', fontsize=16)
    plt.legend()
    plt.savefig(f"{metric}_binding_nonbinding_distribution.png")
    plt.close(fig)
    

    # Calculate and print statistical information
    mean = combined_metrics[f'{metric}_all'].mean()
    median = combined_metrics[f'{metric}_all'].median()
    std_dev = combined_metrics[f'{metric}_all'].std()
    iqr = combined_metrics[f'{metric}_all'].quantile(0.75) - combined_metrics[f'{metric}_all'].quantile(0.25)
    print(metric)
    print('ALL:')
    print(f"Mean: {mean}")
    print(f"Median: {median}")
    print(f"Standard Deviation: {std_dev}")
    print(f"IQR: {iqr}")

    mean = combined_metrics[f'{metric}_close'].mean()
    median = combined_metrics[f'{metric}_close'].median()
    std_dev = combined_metrics[f'{metric}_close'].std()
    iqr = combined_metrics[f'{metric}_close'].quantile(0.75) - combined_metrics[f'{metric}_close'].quantile(0.25)

    print('CLOSE:')
    print(f"Mean: {mean}")
    print(f"Median: {median}")
    print(f"Standard Deviation: {std_dev}")
    print(f"IQR: {iqr}")

    t_stat, p_value = ttest_rel(combined_metrics[f'{metric}_all'], combined_metrics[f'{metric}_close'])

    # Print the results of the paired t-test
    print(f'Paired t-test results on {metric}:')
    print(f'T-statistic: {t_stat}')
    print(f'P-value: {p_value}')


# Create a grouped boxplot for each metric
for metric in ['RSCCS_med', 'OPIA_med', 'EDIA_med', 'RSR_med', 'SRSR_med']:
    # Perform paired t-tests

    combined_metrics['resolution_group'] = (combined_metrics['Resolution'] // 0.5) * 0.5

    plt.figure(figsize=(12, 8))
    sns.boxplot(x='resolution_group', y=f'{metric}_close', data=combined_metrics, color='#FF8C00')
    sns.boxplot(x='resolution_group', y=f'{metric}_all', data=combined_metrics, color='#006C6C')
    plt.xlabel('Resolution (grouped by 0.5 increments)', fontsize=16)
    plt.ylabel(f'{metric}', fontsize=16)
    plt.legend(['Close Residues', 'All Residues'])
    plt.savefig(f"{metric}_resolution_grouped_boxplot.png")
    plt.close()



