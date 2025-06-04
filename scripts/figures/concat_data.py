import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys


all_files = glob.glob("qfit/*_rfree.csv") #read in full protein files
li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',')
    li.append(df)
qFit_rfree = pd.concat(li, axis=0, ignore_index=True)
qFit_rfree = redo_rfree.rename(columns={'Rfree': 'Rfree_qFit'})

all_files = glob.glob("redo/*_rfree.csv") #read in full protein files
li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',')
    li.append(df)
redo_rfree = pd.concat(li, axis=0, ignore_index=True)

#rename columns of redo_rfree
redo_rfree = redo_rfree.rename(columns={'Rfree': 'Rfree_redo'})
print(redo_rfree.head())

#MATCH PDBS
merge_redo_qfit = redo_rfree.merge(qFit_rfree, on='PDB')
print(merge_redo_qfit.head())



