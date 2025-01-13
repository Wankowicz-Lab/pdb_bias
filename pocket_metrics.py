import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
from scipy import stats
import Bio.PDB
from Bio.PDB import *
from scipy.stats import ttest_rel
import gc
import dask.dataframe as dd

print('######################################################################################')
print('_________________________________POCKET_______________________________________________')

plt.rcParams.update({'font.size': 30})

# Use glob to find all *_pdb_stats.txt files in the specified directory
pdb_stats_files = glob.glob('/dors/wankowicz_lab/PDBRedo/pdb-redo3/output/*_pdb_stats.txt')

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
all_pdb_stats_data.rename(columns={'PDB ID': 'PDB'}, inplace=True)
all_pdb_stats_data.to_csv('all_pdb_stats_data.csv')

# Convert RFree to numeric, replacing 'None' with NaN
all_pdb_stats_data['RFree'] = pd.to_numeric(all_pdb_stats_data['RFree'], errors='coerce')

close_resi = pd.read_csv('pocket_resi_close.csv', sep=',', dtype={'resi': 'int32', 'chain': 'str', 'PDB': 'str'})
close_resi = close_resi.rename(columns={'Resi': 'resi'})
close_resi['chain'] = 'A'
close_resi['PDB'] = close_resi['PDB'].str.upper()
# Remove spaces before numbers in 'resi' column
#close_resi['resi'] = close_resi['resi'].astype(str).str.strip()

# Use glob to find all *_residue_metrics.txt files in the specified directory
b = glob.glob('/dors/wankowicz_lab/PDBRedo/pdb-redo3/output/*_residue_metrics.txt')

# Initialize a list to store the DataFrames
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

# Display the first few rows of the combined DataFrame
residue_metrics = residue_metrics.rename(columns={'seqNum': 'resi', 'AsymID': 'chain'})
residue_metrics['PDB'] = residue_metrics['PDB'].str.upper()


# Clean up column names
residue_metrics = residue_metrics.rename(columns={'seqNum': 'resi', 'AsymID': 'chain'})


gc.collect()
# Now perform the merge with consistent data types
#close_resi_metrics = pd.merge(close_resi[['resi', 'chain', 'PDB']], residue_metrics[['resi', 'chain', 'PDB', 'RSCCS', 'OPIA', 'EDIAm', 'RSR', 'SRSR']], 
#                            on=['resi', 'chain', 'PDB'], how='inner')

print(close_resi.head())

close_resi_metrics = dd.merge(close_resi[['resi', 'chain', 'PDB']], 
                                 residue_metrics[['resi', 'chain', 'PDB', 'RSCCS', 'OPIA', 'EDIAm', 'RSR', 'SRSR']],
                                 on=['resi', 'chain', 'PDB'], how='inner')

print(residue_metrics.head())

close_resi_metrics = close_resi_metrics.drop_duplicates()
print(close_resi_metrics.head())

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

metrics_close.to_csv('metrics_pocket_close.csv')
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
metrics_not_close.to_csv('metrics_not_pocket.csv')                               
# Filter residue_metrics by PDBs present in close_resi_metrics, this reduces the size of residue_metrics
residue_metrics = metrics_not_close[metrics_not_close['PDB'].isin(close_resi_metrics['PDB'])]

# Perform one merge between close_resi_metrics and filtered residue_metrics, including all necessary columns
combined_metrics = pd.merge(metrics_close, metrics_not_close, on=['PDB'], suffixes=('_close', '_all'))
combined_metrics.to_csv('combined_metrics_pocket.csv')



# Merge combined_metrics with all_pdb_stats_data (which contains resolution and Rfree) in one step
combined_metrics2 = pd.merge(combined_metrics, all_pdb_stats_data[['PDB', 'Resolution', 'RFree']], on='PDB')


combined_metrics = combined_metrics.dropna()
#residue_metrics = residue_metrics[residue_metrics['PDB'].isin(close_resi_metrics['PDB'])]
#combined_metrics = pd.merge(close_resi_metrics, residue_metrics, on='PDB', suffixes=('_close', '_all'))
#combined_metrics = pd.merge(combined_metrics, all_pdb_stats_data[['PDB', 'resolution', 'Rfree']], on='PDB')
del residue_metrics
del close_resi_metrics

# Plot and compute statistics for each metric
for metric in ['RSCCS_med', 'OPIA_med', 'EDIA_med', 'RSR_med', 'SRSR_med', 'RSCCS_mean', 'OPIA_mean', 'EDIA_mean', 'RSR_mean', 'SRSR_mean']:
    
    fig = plt.figure(figsize=(10, 6))
    sns.kdeplot(combined_metrics[f'{metric}_all'].dropna(), color='#006C6C', linewidth=8)
    sns.kdeplot(combined_metrics[f'{metric}_close'].dropna(), color='#FF8C00', linewidth=8)
    plt.xlabel(f'{metric}', fontsize=16)
    plt.ylabel('Density', fontsize=16)
    plt.legend(['All Residues', 'Pocket Residues'],fontsize=12)
    plt.savefig(f"{metric}_pocket_distribution.png")
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


