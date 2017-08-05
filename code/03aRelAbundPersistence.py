###############################################################################
# relAbund.py
# Copyright (c) 2017, Joshua J Hamilton
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Given a BIOM table of unique sequences abundances:
#   Rarefy all samples to a given depth
#   Compute relatibe abundances
#   Discard all sequences which aren't present at a specific abundance
# These relative abundance tables serve as the basis for modeling
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import re

#%%#############################################################################
### Define folder structure
################################################################################

deblurDir = '../results/deblur' # BIOM file from deblurring

#%%#############################################################################
#### Import otu table and visualize distribution of relative abundances
#################################################################################

# Import table
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)

# Compute counts to relative abundances
relAbundTable = otuTable / otuTable.sum(axis=0)

# Convert to a single series and remove 0 values
relAbundSeries = relAbundTable.unstack()
relAbundSeries = relAbundSeries[relAbundSeries!=0]

# Convert to log10 and plot a histogram with integer bin size
logRelAbundSeries = np.log10(relAbundSeries)
plt.figure(0)
plt.hist(logRelAbundSeries, bins = range(math.floor(logRelAbundSeries.min()), math.ceil(logRelAbundSeries.max())+1), normed=True)
plt.xlabel('Log10 relative abundance of OTU')
plt.ylabel('% Total OTUs')
plt.savefig(deblurDir+'/Total OTUs vs. Abundance.png')

#%%#############################################################################
### Now let's do the same but with the average of sample replicates
################################################################################

# Import table and sort
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)
otuTable = otuTable.reindex_axis(sorted(otuTable.columns), axis=1)

# Compute counts to relative abundances
avgRelAbundTable = otuTable / otuTable.sum(axis=0)

# Compute average for each sample (every 2 columns)
avgRelAbundTable = avgRelAbundTable.groupby(np.arange(len(avgRelAbundTable.columns))//2, axis=1).mean()

# Retrive the list of sample names and remove the '.R2'
sampleNames = list(otuTable.columns[::2])
sampleNames = [re.sub('.R1', '', sample) for sample in sampleNames]
avgRelAbundTable.columns = sampleNames

# Convert to a single series and remove 0 values
avgRelAbundSeries = avgRelAbundTable.unstack()
avgRelAbundSeries = avgRelAbundSeries[avgRelAbundSeries!=0]

# Convert to log10 and plot a histogram with integer bin size
logAvgRelAbundSeries = np.log10(avgRelAbundSeries)
plt.figure(1)
plt.hist(logAvgRelAbundSeries, bins = range(math.floor(logAvgRelAbundSeries.min()), math.ceil(logAvgRelAbundSeries.max())+1), normed=True)
plt.xlabel('Log10 average relative abundance of OTU')
plt.ylabel('% Total OTUs')
plt.savefig(deblurDir+'/Total OTUs vs. Average Abundance.png')

#%%#############################################################################
### Now let's plot the fraction of total abundance in each bin
################################################################################

# Create a dataframe showing the % total abundance accounted for by each OTU in each sample

logRelAbundSeriesTable  = pd.DataFrame(data=0, index=logRelAbundSeries.index, columns=['Log Rel Abund', 'Total Seqs'])
logRelAbundSeriesTable['Log Rel Abund'] = logRelAbundSeries
logRelAbundSeriesTable['Total Seqs'] = otuTable.unstack()
logRelAbundSeriesTable['% Total Seqs'] = logRelAbundSeriesTable['Total Seqs'] / logRelAbundSeriesTable['Total Seqs'].sum(axis=0)

# Use these data do calculate the fraction of sequences within each relative abundance window
histTotalSeqDF = pd.DataFrame(data=0, index = range(math.floor(logRelAbundSeries.min()), math.ceil(logRelAbundSeries.max())), columns=['% Total Seqs'])

for floor in histTotalSeqDF.index:
    histTotalSeqDF.loc[floor] = logRelAbundSeriesTable[(logRelAbundSeriesTable['Log Rel Abund'] >= floor) & (logRelAbundSeriesTable['Log Rel Abund'] < floor+1)]['% Total Seqs'].sum()

# Plot as a bar chart
plt.figure(2)
plt.bar(histTotalSeqDF.index, histTotalSeqDF['% Total Seqs'], align='edge')
plt.xlabel('Log10 relative abundance of OTU')
plt.ylabel('% Total Sequences')
plt.savefig(deblurDir+'/Total Sequences vs. Abundance.png')

#%%#############################################################################
### And the same but again with sample replicates
################################################################################

# First create a dataframe giving the total number of sequences in each set of replicates
summedOtuTable = otuTable.groupby(np.arange(len(otuTable.columns))//2, axis=1).sum()
sampleNames = list(otuTable.columns[::2])
sampleNames = [re.sub('.R1', '', sample) for sample in sampleNames]
summedOtuTable.columns = sampleNames

# Create a dataframe showing the % total abundance accounted for by each OTU in each sample
logAvgRelAbundSeriesTable  = pd.DataFrame(data=0, index=logAvgRelAbundSeries.index, columns=['Log Rel Abund', 'Total Seqs'])
logAvgRelAbundSeriesTable['Log Rel Abund'] = logAvgRelAbundSeries
logAvgRelAbundSeriesTable['Total Seqs'] = summedOtuTable.unstack()
logAvgRelAbundSeriesTable['% Total Seqs'] = logAvgRelAbundSeriesTable['Total Seqs'] / logAvgRelAbundSeriesTable['Total Seqs'].sum(axis=0)

# Use these data do calculate the fraction of sequences within each relative abundance window
histTotalSeqDF = pd.DataFrame(data=0, index = range(math.floor(logAvgRelAbundSeries.min()), math.ceil(logAvgRelAbundSeries.max())), columns=['% Total Seqs'])

for floor in histTotalSeqDF.index:
    histTotalSeqDF.loc[floor] = logAvgRelAbundSeriesTable[(logAvgRelAbundSeriesTable['Log Rel Abund'] >= floor) & (logAvgRelAbundSeriesTable['Log Rel Abund'] < floor+1)]['% Total Seqs'].sum()

# Plot as a bar chart
plt.figure(3)
plt.bar(histTotalSeqDF.index, histTotalSeqDF['% Total Seqs'], align='edge')
plt.xlabel('Log10 average relative abundance of OTU')
plt.ylabel('% Total Sequences')
plt.savefig(deblurDir+'/Total Sequences vs. Average Abundance.png')

#%%#############################################################################
#### Import otu table and visualize distribution of persistence
################################################################################

# Define a function to obtain a list of fractions
def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step

# Import table
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)

# Compute persistence
otuTable['Persistence'] = otuTable.astype(bool).sum(axis=1) / len(otuTable.columns)

## Plot a histogram with integer bin size
plt.figure(4)
# Assign weights to the sum of all bins equals one
# Necessary b/c bins are of width < 1
# See https://stackoverflow.com/questions/3866520/plotting-histograms-whose-bar-heights-sum-to-1-in-matplotlib/16399202#16399202)
weights = np.ones_like(otuTable['Persistence'])/float(len(otuTable['Persistence']))

plt.hist(otuTable['Persistence'], bins=list(frange(0, 1, 0.1)), weights=weights)
plt.xlabel('Persistence')
plt.ylabel('Fraction of OTUs')
plt.savefig(deblurDir+'/Total OTUs vs. Persistence.png')

#%%#############################################################################
#### Repeat for average of replicates
################################################################################
# Import table and sort
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)
otuTable = otuTable.reindex_axis(sorted(otuTable.columns), axis=1)

# Compute counts to relative abundances
avgRelAbundTable = otuTable / otuTable.sum(axis=0)

# Compute average for each sample (every 2 columns)
avgRelAbundTable = avgRelAbundTable.groupby(np.arange(len(avgRelAbundTable.columns))//2, axis=1).mean()

# Retrive the list of sample names and remove the '.R2'
sampleNames = list(otuTable.columns[::2])
sampleNames = [re.sub('.R1', '', sample) for sample in sampleNames]
avgRelAbundTable.columns = sampleNames

# Compute persistence
avgRelAbundTable['Persistence'] = avgRelAbundTable.astype(bool).sum(axis=1) / len(avgRelAbundTable.columns)

## Plot a histogram with integer bin size
plt.figure(5)
# Assign weights to the sum of all bins equals one
# Necessary b/c bins are of width < 1
# See https://stackoverflow.com/questions/3866520/plotting-histograms-whose-bar-heights-sum-to-1-in-matplotlib/16399202#16399202)
weights = np.ones_like(avgRelAbundTable['Persistence'])/float(len(avgRelAbundTable['Persistence']))

plt.hist(avgRelAbundTable['Persistence'], bins=list(frange(0, 1, 0.1)), weights=weights)
plt.xlabel('Average Persistence of Replicates')
plt.ylabel('Fraction of OTUs')
plt.savefig(deblurDir+'/Total OTUs vs. Average Persistence.png')

#%%#############################################################################
### Import otu table and visualize fraction ot total sequences in each
### persistence bin
################################################################################

# Define a function to obtain a list of fractions
def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step

# Import table
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)

# Compute total seqs and persistence
otuTable['Total Seqs'] = otuTable.sum(axis=1)
otuTable['Rel Abund'] = otuTable['Total Seqs'] / otuTable['Total Seqs'].sum()
otuTable['Persistence'] = otuTable.astype(bool).sum(axis=1) / len(otuTable.columns)

# Create a dataframe showing the % total abundance accounted for by each OTU in each sample
# Use these data do calculate the fraction of sequences within each relative abundance window
histPersistDF = pd.DataFrame(data=0, index = list(frange(0, 1, 0.05)), columns=['Rel Abund'])

for floor in histPersistDF.index:
    histPersistDF.loc[floor] = otuTable[(otuTable['Persistence'] >= floor) & (otuTable['Persistence'] < floor+0.05)]['Rel Abund'].sum()

# Plot as a bar chart
plt.figure(6)
plt.bar(histPersistDF.index, histPersistDF['Rel Abund'], align='edge')
plt.xlim([0, 1])
plt.xlabel('Persistence')
plt.ylabel('% Total Sequences')
plt.savefig(deblurDir+'/Total Sequences vs. Persistence.png')

#%%#############################################################################
### Repeat for average of each replicate
################################################################################

# Define a function to obtain a list of fractions
def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step

# Import table and sort
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)
otuTable = otuTable.reindex_axis(sorted(otuTable.columns), axis=1)

# Compute counts to relative abundances
avgRelAbundTable = otuTable / otuTable.sum(axis=0)

# Compute average for each sample (every 2 columns)
avgRelAbundTable = avgRelAbundTable.groupby(np.arange(len(avgRelAbundTable.columns))//2, axis=1).mean()

# Retrive the list of sample names and remove the '.R2'
sampleNames = list(otuTable.columns[::2])
sampleNames = [re.sub('.R1', '', sample) for sample in sampleNames]
avgRelAbundTable.columns = sampleNames
    
# Compute relative abundance and persistence
relAbundSeries = avgRelAbundTable.sum(axis=1) / avgRelAbundTable.sum(axis=1).sum()
persistSeries = avgRelAbundTable.astype(bool).sum(axis=1) / len(avgRelAbundTable.columns)

# Add to DF
avgRelAbundTable['Rel Abund'] = relAbundSeries
avgRelAbundTable['Persistence'] = persistSeries
    
# Create a dataframe showing the % total abundance accounted for by each OTU in each sample
# Use these data do calculate the fraction of sequences within each relative abundance window
histAvgPersistDF = pd.DataFrame(data=0, index = list(frange(0, 1, 0.05)), columns=['Rel Abund'])

for floor in histAvgPersistDF.index:
    histAvgPersistDF.loc[floor] = avgRelAbundTable[(avgRelAbundTable['Persistence'] >= floor) & (avgRelAbundTable['Persistence'] < floor+0.05)]['Rel Abund'].sum()

# Plot as a bar chart
plt.figure(7)
plt.bar(histPersistDF.index, histPersistDF['Rel Abund'], align='edge')
plt.xlim([0, 1])
plt.xlabel('Persistence of Averaged Replicates')
plt.ylabel('% Total Sequences')
plt.savefig(deblurDir+'/Total Sequences vs. Average Persistence.png')