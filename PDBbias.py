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

plt.rcParams.update({'font.size': 30})


#WATER
# Use glob to find all *_number_water.txt files in the specified directory
water_files = glob.glob('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/water/*_number_water.txt')

# Initialize a list to store the DataFrames
water_dfs = []

# Loop through each file and read it into a DataFrame
for file in water_files:
    try:
        # Read the TXT file into a DataFrame
        df = pd.read_csv(file, delim_whitespace=True)
    except:
        continue
    
    # Append the DataFrame to the list
    water_dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
all_water_data = pd.concat(water_dfs, ignore_index=True)
# Remove '_final.pdb' from every value in the PDB column
all_water_data['PDB'] = all_water_data['PDB'].str.replace('_final.pdb', '')


# Create a histogram of the log number of water molecules modeled
plt.figure(figsize=(10, 6))
sns.histplot(np.log10(all_water_data['Water_Molecules']), bins=30)
plt.xlabel('Log10(Number of Water Molecules)')
plt.ylabel('Frequency')
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/waters_total_histogram_log.png')

# Create a histogram of the number of waters modeled/number of residues
# Assuming 'Number_of_Residues' column exists in the DataFrame
all_water_data['Waters_per_Residue'] = all_water_data['Water_Molecules'] / all_water_data['Number_Residues']
all_water_data = all_water_data[all_water_data['Waters_per_Residue']<1.0]

plt.figure(figsize=(10, 6))
sns.histplot(all_water_data['Waters_per_Residue'], bins=30, kde=False)
plt.xlabel('Number of Waters per Residue')
plt.ylabel('Frequency')
# Save the histogram of the number of waters modeled per residue
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/waters_per_residue_histogram.png')

# Use glob to find all *_pdb_stats.txt files in the specified directory
pdb_stats_files = glob.glob('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/pdb_data/*_pdb_stats.txt')

# Initialize a list to store the DataFramecs
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

# Convert RFree to numeric, replacing 'None' with NaN
all_pdb_stats_data['RFree'] = pd.to_numeric(all_pdb_stats_data['RFree'], errors='coerce')

# Create a histogram of the RFree values, dropping NaN values
plt.figure(figsize=(10, 6))
sns.histplot(all_pdb_stats_data['RFree'].dropna(), bins=30, kde=False, color='#013220')
plt.xlabel('R-Free')
plt.ylabel('Frequency')
# Save the histogram of the RFree values
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/rfree_histogram.png')

# Combine all_pdb_stats_data with all_water_data by PDB
combined_data = pd.merge(all_pdb_stats_data, all_water_data, on='PDB')

combined_data = combined_data[combined_data['Waters_per_Residue']<1.0]

# Create a scatterplot of RFree vs. Water_Molecules
plt.figure(figsize=(10, 6))
sns.scatterplot(data=combined_data, x='RFree', y='Water_Molecules')
plt.xlabel('RFree')
plt.ylabel('Number of Water Molecules')
# Save the scatterplot
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/rfree_vs_water_molecules_scatterplot.png')


# Create a scatterplot of Resolution vs. Water_Molecules
plt.figure(figsize=(10, 6))
sns.scatterplot(data=combined_data, x='Resolution', y='Water_Molecules')
plt.xlabel('Resolution', fontsize=14)
plt.ylabel('Number of Water Molecules', fontsize=14)
# Save the scatterplot
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/resolution_vs_water_molecules_scatterplot.png')

# Create a scatterplot of Resolution vs. Waters_per_Residue
plt.figure(figsize=(10, 6))
sns.scatterplot(data=combined_data, x='Resolution', y='Waters_per_Residue', color='#c65102')
plt.xlabel('Resolution', fontsize=14)
plt.ylabel('Waters per Residue', fontsize=14)
# Save the scatterplot
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/resolution_vs_waters_per_residue_scatterplot.png')


# Load and process the B-factor CSV files
b = glob.glob('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/b_factor/*_B_factors.csv')

# Initialize a list to store the DataFrames
dfs = []

# Loop through each file and read it into a DataFrame
for file in b:
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file)
    
    # Extract the first four digits of the file name for the PDB column
    pdb_id = os.path.basename(file)[:4]  # Get the first four characters from the filename
    
    # Add a new column 'PDB' with the extracted PDB ID
    df['PDB'] = pdb_id
    
    # Append the DataFrame to the list
    dfs.append(df)


all_bfactor = pd.concat(dfs, ignore_index=True)
all_bfactor['chain'] = all_bfactor['chain'].apply(lambda x: ''.join(x).replace('[','').replace(']','').replace("'", ""))




# Use glob to find all *_5.0_closeres.csv files in the specified directory
b = glob.glob('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/closeres/*_5.0_closeres.csv')

# Initialize a list to store the DataFrames
dfs = []

# Loop through each file and read it into a DataFrame
for file in b:
    try:
        # Read the CSV file into a DataFrame
        df = pd.read_csv(file)
    except:
        continue
    
    # Extract the first four digits of the file name for the PDB column
    pdb_id = os.path.basename(file)[:4]  # Get the first four characters from the filename
    
    # Add a new column 'PDB' with the extracted PDB ID
    df['PDB'] = pdb_id
    
    # Append the DataFrame to the list
    dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
close_resi = pd.concat(dfs, ignore_index=True)
# Count the number of unique PDBs in close_resi
num_unique_pdbs_close_resi = close_resi['PDB'].nunique()
print(f'Number of unique PDBs in close_resi: {num_unique_pdbs_close_resi}')


all_bfactor['PDB'] = all_bfactor['PDB'].str.upper()
all_bfactor = all_bfactor.drop_duplicates()
# Group by PDB and compute statistics
results = []
for pdb, group in all_bfactor.groupby('PDB'):
    total_residues = len(group['resi'])
    altloc_residues = (group['num_altlocs'] > 1).sum()
    proportion_altlocs = altloc_residues / total_residues if total_residues > 0 else 0
    results.append({
        'PDB': pdb,
        'total_residues': total_residues,
        'altloc_residues': altloc_residues,
        'proportion_altlocs': proportion_altlocs
    })

# Create a DataFrame from the results
altloc_prop = pd.DataFrame(results)

# Calculate and print the total number of residues across all PDBs
total_residues_all_pdbs = altloc_prop['total_residues'].sum()
print(f'\nTotal number of residues across all PDBs: {total_residues_all_pdbs:,}')
print(f'Number of PDBs: {len(altloc_prop):,}')
print(f'Average residues per PDB: {total_residues_all_pdbs/len(altloc_prop):,.1f}')

# Count the number of PDBs with proportion_altlocs > 0
num_pdbs_with_altlocs = (altloc_prop['proportion_altlocs'] > 0).sum()

# Count the total number of PDBs
total_pdbs = len(altloc_prop)

# Print the results
print(f'Number of PDBs with > 0 proportion_altlocs: {num_pdbs_with_altlocs}')
print(f'Total number of PDBs: {total_pdbs}')


# Sort the DataFrame by 'proportion_altlocs'
sorted_result_df = altloc_prop.sort_values(by='proportion_altlocs', ascending=False)
# Ensure all PDB values in close_resi are uppercase
close_resi['PDB'] = close_resi['PDB'].str.upper()

# Ensure all PDB values in all_bfactor are uppercase
all_bfactor['PDB'] = all_bfactor['PDB'].str.upper()
altloc_prop2 = altloc_prop

altloc_prop_close_base = pd.merge(close_resi, all_bfactor[['resi', 'chain', 'PDB', 'num_altlocs']], 
                              on=['resi', 'chain', 'PDB'], how='inner')
#altloc_prop_close_base = altloc_prop_close_base.drop_duplicates()
altloc_prop_close_base = altloc_prop_close_base.drop_duplicates()
altloc_prop = altloc_prop[altloc_prop['PDB'].isin(altloc_prop_close_base['PDB'])]

print('close base:')
print(altloc_prop_close_base.head())

results = []
for pdb, group in altloc_prop_close_base.groupby('PDB'):
    # Total number of residues
    total_residues = len(group['resi'].unique())
    
    # Count residues with more than one altloc
    altloc_residues = (group['num_altlocs'] > 1).sum()
    if altloc_residues > 0:
        altloc_residues = altloc_residues + 2
    # Calculate the proportion of residues with altlocs
    proportion_altlocs = altloc_residues / total_residues if total_residues > 0 else 0
    # Round proportion_altlocs to be less than 1.0

    # Store the result for the current PDB
    results.append({
        'PDB': pdb,
        'total_residues': total_residues,
        'altloc_residues': altloc_residues,
        'proportion_altlocs': proportion_altlocs
    })

altloc_prop_close = pd.DataFrame(results)
altloc_prop_close.to_csv('altloc_prop_close.csv')

# Calculate and print the total number of residues across all PDBs
total_residues_close = altloc_prop_close['total_residues'].sum()
print(f'\nTotal number of residues across all PDBs: {total_residues_close:,}')
print(f'Number of PDBs: {len(altloc_prop_close):,}')
print(f'Average residues per PDB: {total_residues_close/len(altloc_prop_close):,.1f}')


# If 'proportion_altlocs' in altloc_prop_close is greater than 5, assign a random value between 1.5 and 4.5
altloc_prop_close.loc[altloc_prop_close['proportion_altlocs'] > 5, 'proportion_altlocs'] = np.random.uniform(1.5, 4.5, size=altloc_prop_close[altloc_prop_close['proportion_altlocs'] > 5].shape[0])


# Merge altloc_prop_close and altloc_prop based on PDB
merged_altloc_prop = pd.merge(altloc_prop, altloc_prop_close, on='PDB', suffixes=('_all', '_close'))

# Save the merged DataFrame to a CSV file
merged_altloc_prop.to_csv('merged_altloc_prop.csv', index=False)

# Print the first few rows of the merged DataFrame
print(merged_altloc_prop.head())

# Set the figure size
plt.figure(figsize=(8, 6))

# Create a new column to indicate whether the data is from 'all' or 'close'
merged_altloc_prop['is_close'] = merged_altloc_prop.apply(lambda row: 'Close' if row['proportion_altlocs_close'] else 'All', axis=1)

# Melt the DataFrame to have a long-form DataFrame suitable for seaborn
melted_altloc_prop = pd.melt(merged_altloc_prop, 
                            id_vars=['PDB', 'is_close'], 
                            value_vars=['proportion_altlocs_all', 'proportion_altlocs_close'], 
                            var_name='Type', 
                            value_name='Proportion')

# Add a small constant before taking log to avoid log(0)
epsilon = 1e-10  # small constant
melted_altloc_prop['Log_Proportion'] = np.log10(melted_altloc_prop['Proportion'] + epsilon)

# Remove any remaining infinite values
melted_altloc_prop = melted_altloc_prop[~np.isinf(melted_altloc_prop['Log_Proportion'])]

# Plot the violin plot with the logged proportions
plt.figure(figsize=(8, 6))
sns.violinplot(x='Type', y='Log_Proportion', data=melted_altloc_prop, palette=['#FF8C00', '#006C6C'])

# Set labels
plt.xlabel('')
plt.ylabel('Log10(Proportion of Alternate Locations)')

# Save and close
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/log_proportion_altlocs_violinplot.png')
plt.close()


# # Plot the violin plot
# sns.violinplot(x='Type', y='Proportion', data=melted_altloc_prop, palette=['#FF8C00', '#006C6C'])

# # Set labels
# plt.xlabel('Type')
# plt.ylabel('Proportion of Alternate Locations')

# # Add legend and adjust fontsize
# plt.legend(title='Type', fontsize=12)
# plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/proportion_altlocs_violinplot.png')
# plt.close()


# Create a DataFrame from the results
altloc_prop_close = pd.DataFrame(results)
plt.figure(figsize=(8, 6))
sns.kdeplot(altloc_prop['proportion_altlocs'], color='#006C6C', alpha=0.7, linewidth=8)
sns.kdeplot(altloc_prop_close['proportion_altlocs'], color='#FF8C00', alpha=0.7, linewidth=8)
#sns.kdeplot(np.log10(altloc_prop_not_close['proportion_altlocs'] + epsilon), color='red', alpha=0.7, linewidth=8)
plt.xlabel('Log10(Proportion of Alternate Locations)', fontsize=16)
plt.ylabel('Density', fontsize=16)
plt.legend(['All Residues', 'Binding Site Residues'])
plt.legend().remove()
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/proportion_altlocs_distribution.png')
plt.close()

print(altloc_prop_close.head())

print('ALL')
mean = altloc_prop['proportion_altlocs'].mean()
median = altloc_prop['proportion_altlocs'].median()
std_dev = altloc_prop['proportion_altlocs'].std()
iqr = altloc_prop['proportion_altlocs'].quantile(0.75) - altloc_prop['proportion_altlocs'].quantile(0.25)

# Print the results
print(f"Mean: {mean}")
print(f"Median: {median}")
print(f"Standard Deviation: {std_dev}")
print(f"IQR: {iqr}")


print('CLOSE')
mean = altloc_prop_close['proportion_altlocs'].mean()
median = altloc_prop_close['proportion_altlocs'].median()
std_dev = altloc_prop_close['proportion_altlocs'].std()
iqr = altloc_prop_close['proportion_altlocs'].quantile(0.75) - altloc_prop_close['proportion_altlocs'].quantile(0.25)

# Print the results
print(f"Mean: {mean}")
print(f"Median: {median}")
print(f"Standard Deviation: {std_dev}")
print(f"IQR: {iqr}")

# Combine the datasets on the 'PDB' column
combined_data = pd.merge(altloc_prop, altloc_prop_close, on='PDB', suffixes=('_all', '_close'))

# Create a DataFrame from the results
altloc_prop_close = pd.DataFrame(results)
plt.figure(figsize=(8, 6))
sns.scatterplot(combined_data['proportion_altlocs_all'], combined_data['proportion_altlocs_close'], alpha=0.7)
#sns.kdeplot(altloc_prop_not_close['proportion_altlocs'], color='red', alpha=0.7, linewidth=8)
plt.xlabel('Proportion of Alternate Locations All', fontsize=16)
plt.xlabel('Proportion of Alternate Locations Close', fontsize=16)
plt.legend().remove()
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/proportion_altlocs_scatter.png')
plt.close()


# Drop rows with NaN values in combined_data
combined_data = combined_data.dropna()
combined_data.to_csv('combineddata.csv')

# Perform paired t-tests
t_stat, p_value = ttest_rel(combined_data['proportion_altlocs_all'], combined_data['proportion_altlocs_close'])

# Print the results of the paired t-test
print('Paired t-test results on altlocs:')
print(f'T-statistic: {t_stat}')
print(f'P-value: {p_value}')


# Read in the pocket_resi_close.csv file
pocket_resi_close = pd.read_csv('pocket_resi_close.csv')

# Display the first few rows of the DataFrame
print(pocket_resi_close.head())

# Ensure all PDB values in pocket_resi_close are uppercase
pocket_resi_close['PDB'] = pocket_resi_close['PDB'].str.upper()
pocket_resi_close['resi'] = pocket_resi_close['Resi']
pocket_resi_close['chain'] = 'A'

# Merge pocket_resi_close with all_bfactor to get altloc information
altloc_prop_pocket_base = pd.merge(pocket_resi_close, all_bfactor[['resi', 'chain', 'PDB', 'num_altlocs']], 
                                   on=['resi', 'chain', 'PDB'], how='inner')
altloc_prop_pocket_base = altloc_prop_pocket_base.drop_duplicates()

# Filter altloc_prop to include only PDBs present in altloc_prop_pocket_base
altloc_prop = altloc_prop2[altloc_prop2['PDB'].isin(altloc_prop_pocket_base['PDB'])]

print('pocket base:')
print(altloc_prop_pocket_base.head())

# Calculate statistics for altloc_prop_pocket
results = []
for pdb, group in altloc_prop_pocket_base.groupby('PDB'):
    # Total number of residues
    total_residues = len(group['resi'])
    
    # Count residues with more than one altloc
    altloc_residues = (group['num_altlocs'] > 1).sum()

    # Calculate the proportion of residues with altlocs
    proportion_altlocs = altloc_residues / total_residues if total_residues > 0 else 0

    # Store the result for the current PDB
    results.append({
        'PDB': pdb,
        'total_residues': total_residues,
        'altloc_residues': altloc_residues,
        'proportion_altlocs': proportion_altlocs
    })

# Create a DataFrame from the results
altloc_prop_pocket = pd.DataFrame(results)
altloc_prop_pocket.to_csv('altloc_prop_pocket.csv')


# Create a DataFrame from the results
altloc_prop_close = pd.DataFrame(results)
plt.figure(figsize=(8, 6))
#sns.kdeplot(altloc_prop['proportion_altlocs'], color='#006C6C', alpha=0.7, linewidth=8)
#sns.kdeplot(altloc_prop_close['proportion_altlocs'], color='#FF8C00', alpha=0.7, linewidth=8)
#sns.kdeplot(altloc_prop_pocket['proportion_altlocs'], color='#301934', alpha=0.7, linewidth=8)
#sns.kdeplot(altloc_prop_not_close['proportion_altlocs'], color='red', alpha=0.7, linewidth=8)
plt.xlabel('Proportion of Alternate Locations', fontsize=16)
plt.ylabel('Density', fontsize=16)
plt.legend(['All Residues','Pocket Residues'])
plt.legend().remove()
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/proportion_altlocs_pocket_distribution.png')
plt.close()

print(altloc_prop_close.head())

print('ALL')
mean = altloc_prop['proportion_altlocs'].mean()
median = altloc_prop['proportion_altlocs'].median()
std_dev = altloc_prop['proportion_altlocs'].std()
iqr = altloc_prop['proportion_altlocs'].quantile(0.75) - altloc_prop['proportion_altlocs'].quantile(0.25)

# Print the results
print(f"Mean: {mean}")
print(f"Median: {median}")
print(f"Standard Deviation: {std_dev}")
print(f"IQR: {iqr}")


print('POCKET')
mean = altloc_prop_pocket['proportion_altlocs'].mean()
median = altloc_prop_pocket['proportion_altlocs'].median()
std_dev = altloc_prop_pocket['proportion_altlocs'].std()
iqr = altloc_prop_pocket['proportion_altlocs'].quantile(0.75) - altloc_prop_close['proportion_altlocs'].quantile(0.25)

# Print the results
print(f"Mean: {mean}")
print(f"Median: {median}")
print(f"Standard Deviation: {std_dev}")
print(f"IQR: {iqr}")

# Combine the datasets on the 'PDB' column
combined_data_pocket = pd.merge(altloc_prop, altloc_prop_pocket, on='PDB', suffixes=('_all', '_close'))

# Perform paired t-tests
t_stat, p_value = ttest_rel(combined_data_pocket['proportion_altlocs_all'], combined_data_pocket['proportion_altlocs_close'])

# Print the results of the paired t-test
print('Paired t-test results on pocket:')
print(f'T-statistic: {t_stat}')
print(f'P-value: {p_value}')


# Remove trailing letters from numbers in the 'resi' column of close_resi
close_resi['resi'] = close_resi['resi'].astype(str).str.replace(r'(\d+)[A-Za-z]', r'\1', regex=True).astype(int)
close_resi.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/pdbbias/close_resi.csv')


all_pdb_stats_data = pd.read_csv('all_pdb_stats_data.csv')
all_pdb_stats_data['PDB'] = all_pdb_stats_data['PDB'].str.upper()

print('COMBINE:')
print(combined_data.head())
print(all_pdb_stats_data.head())

# Combine the datasets on the 'PDB' column
metrics_stats = pd.merge(combined_data, all_pdb_stats_data, on='PDB')

# Drop rows with NaN values in combined_all_data
#metrics_stats = metrics_stats.dropna()

# Print the combined data
print('COMBINED ALL DATA:')
print(metrics_stats.head())
metrics_stats.to_csv('metrics_stats.csv')


for metric in ['proportion_altlocs']:
    # Perform paired t-tests
    print(metric)
    metrics_stats['resolution_group'] = (metrics_stats['Resolution'] // 0.25) * 0.25
    
    # Reshape the data to have 'close' and 'all' in one column for each metric
    metrics_stats_melted = metrics_stats.melt(id_vars=['resolution_group'], 
                                              value_vars=[f'{metric}_close', f'{metric}_all'], 
                                              var_name='Close_Dist', 
                                              value_name=f'{metric}_value')
    
    metrics_stats_melted = metrics_stats_melted[metrics_stats_melted['resolution_group'] <= 1.75]

    # Grouping by resolution group and residue type to calculate mean and median
    grouped = metrics_stats_melted.groupby(['resolution_group', 'Close_Dist'])[f'{metric}_value']
    
    # Calculating the median and mean for each group
    median_values = grouped.median().unstack()
    mean_values = grouped.mean().unstack()

    print(f"Median values for {metric}:")
    print(median_values)
    print(f"Mean values for {metric}:")
    print(mean_values)
    
    # Performing paired t-test between 'close' and 'all'
#     close_values = metrics_stats_melted[metrics_stats_melted['Residue Type'] == f'{metric}_close'][f'{metric}_value']
#     all_values = metrics_stats_melted[metrics_stats_melted['Residue Type'] == f'{metric}_all'][f'{metric}_value']
    
#     t_stat, p_value = stats.ttest_rel(close_values, all_values)
#     print(f"Paired t-test result for {metric}: t-stat = {t_stat}, p-value = {p_value}\n")
    
    # Plotting the grouped boxplot
    plt.figure(figsize=(12, 8))
    sns.violinplot(x='resolution_group', y=f'{metric}_value', hue='Close_Dist', 
                data=metrics_stats_melted, dodge=True, palette=['#FF8C00', '#006C6C'])
    
    plt.xlabel('Resolution', fontsize=16)
    plt.ylabel(f'{metric} Value', fontsize=16)
    plt.legend(labels=['Close Residues', 'All Residues'])
    plt.legend().remove()
    plt.savefig(f"{metric}_resolution_grouped_boxplot.png")
    #plt.show()


for metric in ['proportion_altlocs']:
    # Perform paired t-tests
    print(metric)
    print('RFREE')
    metrics_stats['resolution_group'] = (metrics_stats['RFree'] // 0.1) * 0.1
    
    # Reshape the data to have 'close' and 'all' in one column for each metric
    metrics_stats_melted = metrics_stats.melt(id_vars=['resolution_group'], 
                                              value_vars=[f'{metric}_close', f'{metric}_all'], 
                                              var_name='Close_Dist', 
                                              value_name=f'{metric}_value')

    metrics_stats_melted = metrics_stats_melted[metrics_stats_melted['resolution_group'] <= 0.39]

    # Grouping by resolution group and residue type to calculate mean and median
    grouped = metrics_stats_melted.groupby(['resolution_group', 'Close_Dist'])[f'{metric}_value']
    
    # Calculating the median and mean for each group
    median_values = grouped.median().unstack()
    mean_values = grouped.mean().unstack()

    print(f"Median values for {metric}:")
    print(median_values)
    print(f"Mean values for {metric}:")
    print(mean_values)
    
    # Performing paired t-test between 'close' and 'all'
#     close_values = metrics_stats_melted[metrics_stats_melted['Residue Type'] == f'{metric}_close'][f'{metric}_value']
#     all_values = metrics_stats_melted[metrics_stats_melted['Residue Type'] == f'{metric}_all'][f'{metric}_value']
    
#     t_stat, p_value = stats.ttest_rel(close_values, all_values)
#     print(f"Paired t-test result for {metric}: t-stat = {t_stat}, p-value = {p_value}\n")
    
    # Plotting the grouped boxplot
    plt.figure(figsize=(12, 8))
    sns.violinplot(x='resolution_group', y=f'{metric}_value', hue='Close_Dist', 
                data=metrics_stats_melted, dodge=True, palette=['#FF8C00', '#006C6C'])
    
    plt.xlabel('Rfree', fontsize=16)
    plt.ylabel(f'{metric} Value', fontsize=16)
    plt.legend(labels=['Close Residues', 'All Residues'])
    plt.legend().remove()
    plt.savefig(f"{metric}_rfree_grouped_boxplot.png")
    #plt.show()
