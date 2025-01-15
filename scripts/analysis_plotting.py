import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import os
import glob

from scipy import stats
from scipy.stats import ttest_rel
import scipy.stats as stats
from scipy.stats import gaussian_kde

import Bio.PDB
from Bio.PDB import *
plt.rcParams.update({'font.size': 30})

close_resi = pd.read_csv('close_resi.csv')
metrics_not_close = pd.read_csv('metrics_not_close.csv')
metrics_close = pd.read_csv('metrics_close.csv')
metrics = pd.read_csv('combined_metrics.csv')
all_pdb_stats_data= pd.read_csv('all_pdb_stats_data.csv')
rotamer_sum = pd.read_csv('rotamer_summary.csv')

all_pdb_stats_data['PDB_upper'] = all_pdb_stats_data['PDB'].str.upper()
all_pdb_stats_data = all_pdb_stats_data.drop(columns=['PDB'])

all_pdb_stats_data = all_pdb_stats_data.rename(columns={'PDB_upper': 'PDB'})


# PLOT
plt.figure(figsize=(8, 6))

# Plot KDE for 'EDIA_mean_close' and 'EDIA_mean_all'
sns.kdeplot(metrics['RSR_mean_close'], color='#006C6C', alpha=0.7, linewidth=6)
sns.kdeplot(metrics['RSR_mean_all'], color='#FF8C00', alpha=0.7, linewidth=6)

# Calculate the values for shading (standard deviation) around the KDE
x_values = np.linspace(min(metrics['RSR_mean_close'].min(), metrics['RSR_mean_all'].min()),
                       max(metrics['RSR_mean_close'].max(), metrics['RSR_mean_all'].max()), 1000)

# KDE estimation for shading
kde_close_values = gaussian_kde(metrics['RSR_mean_close'])(x_values)
kde_all_values = gaussian_kde(metrics['RSR_mean_all'])(x_values)

# Calculate standard deviation for the shading region
std_close = np.std(metrics['RSR_mean_close'])
std_all = np.std(metrics['RSR_mean_all'])

# Add shading for standard deviation around the KDE line
plt.fill_between(x_values, kde_close_values - std_close, kde_close_values + std_close, color='#006C6C', alpha=0.2, label='Close - ±1 SD')
plt.fill_between(x_values, kde_all_values - std_all, kde_all_values + std_all, color='#FF8C00', alpha=0.2, label='All - ±1 SD')

# Add labels and legend
plt.xlabel('RSR')
plt.ylabel('Density')
plt.legend(['Binding Site (±1 SD)', 'All (±1 SD)'], fontsize=14)


# Set up the plotting
plt.figure(figsize=(8, 6))

# Plot KDE for 'EDIA_mean_close' and 'EDIA_mean_all'
sns.kdeplot(metrics['RSCCS_mean_close'], color='#006C6C', alpha=0.7, linewidth=6)
sns.kdeplot(metrics['RSCCS_mean_all'], color='#FF8C00', alpha=0.7, linewidth=6)

# Calculate the values for shading (standard deviation) around the KDE
x_values = np.linspace(min(metrics['RSCCS_mean_close'].min(), metrics['RSCCS_mean_all'].min()),
                       max(metrics['RSCCS_mean_close'].max(), metrics['RSCCS_mean_all'].max()), 1000)

# KDE estimation for shading
kde_close_values = gaussian_kde(metrics['RSCCS_mean_close'])(x_values)
kde_all_values = gaussian_kde(metrics['RSCCS_mean_all'])(x_values)

# Calculate standard deviation for the shading region
std_close = np.std(metrics['RSCCS_mean_close'])
std_all = np.std(metrics['RSCCS_mean_all'])

# Add shading for standard deviation around the KDE line
plt.fill_between(x_values, kde_close_values - std_close, kde_close_values + std_close, color='#006C6C', alpha=0.2, label='Close - ±1 SD')
plt.fill_between(x_values, kde_all_values - std_all, kde_all_values + std_all, color='#FF8C00', alpha=0.2, label='All - ±1 SD')

# Add labels and legend
plt.xlabel('RSCC')
plt.ylabel('Density')
plt.legend(['Binding Site (±1 SD)', 'All (±1 SD)'], fontsize=14)


# Set up the plotting
plt.figure(figsize=(8, 6))

# Plot KDE for 'EDIA_mean_close' and 'EDIA_mean_all'
sns.kdeplot(metrics['EDIA_mean_close'], color='#006C6C', alpha=0.7, linewidth=6)
sns.kdeplot(metrics['EDIA_mean_all'], color='#FF8C00', alpha=0.7, linewidth=6)

# Calculate the values for shading (standard deviation) around the KDE
x_values = np.linspace(min(metrics['EDIA_mean_close'].min(), metrics['EDIA_mean_all'].min()),
                       max(metrics['EDIA_mean_close'].max(), metrics['EDIA_mean_all'].max()), 1000)

# KDE estimation for shading
kde_close_values = gaussian_kde(metrics['EDIA_mean_close'])(x_values)
kde_all_values = gaussian_kde(metrics['EDIA_mean_all'])(x_values)

# Calculate standard deviation for the shading region
std_close = np.std(metrics['EDIA_mean_close'])
std_all = np.std(metrics['EDIA_mean_all'])

# Add shading for standard deviation around the KDE line
plt.fill_between(x_values, kde_close_values - std_close, kde_close_values + std_close, color='#006C6C', alpha=0.2, label='Close - ±1 SD')
plt.fill_between(x_values, kde_all_values - std_all, kde_all_values + std_all, color='#FF8C00', alpha=0.2, label='All - ±1 SD')

# Add error bars for each x value (±1 SD)
#plt.errorbar(x_values, kde_close_values, yerr=std_close, fmt='o', color='#006C6C', alpha=0.5, label='Close - Error Bars')
#plt.errorbar(x_values, kde_all_values, yerr=std_all, fmt='o', color='#FF8C00', alpha=0.5, label='All - Error Bars')

# Add labels and legend
plt.xlabel('EDIA Mean Values')
plt.ylabel('Density')
plt.legend(['Binding Site (±1 SD)', 'All (±1 SD)'], fontsize=14)


