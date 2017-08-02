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

#%%#############################################################################
### Define folder structure
################################################################################

deblurDir = '../data/deblur' # BIOM file from deblurring

#%%#############################################################################
### Import otu table and visualize frequency distribution
################################################################################

# Import table
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)

# Compute counts to frequencies
freqTable = otuTable / otuTable.sum(axis=0)

# Convert to a single series and remove 0 values
freqSeries = freqTable.unstack()
freqSeries = freqSeries[freqSeries!=0]

# Convert to log10 and plot a histogram with integer bin size
logSeries = np.log10(freqSeries)
fig1, ax1 = plt.subplots()
fig1 = plt.hist(logSeries, bins = range(math.floor(logSeries.min()), math.ceil(logSeries.max())+1), normed=True)
ax1.set_xlabel('Log10 relative abundance of OTU')
ax1.set_ylabel('% Total OTUs')
plt.savefig(deblurDir+'/Total OTUs vs. Abundance.png')

#%%#############################################################################
### This plot suggests there are no unusually abundant or unusally infrequent taxa
### Now let's plot the fraction of total abundance in each bin
################################################################################

# Create a dataframe showing the % total abundance accounted for by each OTU in each sample

freqAbundTable  = pd.DataFrame(data=0, index=logSeries.index, columns=['Log Rel Abund', 'Total Seqs'])
freqAbundTable['Log Rel Abund'] = logSeries
freqAbundTable['Total Seqs'] = otuTable.unstack()
freqAbundTable['% Total Seqs'] = freqAbundTable['Total Seqs'] / freqAbundTable['Total Seqs'].sum(axis=0)

# Use these data do calculate the fraction of sequences within each relative abundance window
histTotalSeqDF = pd.DataFrame(data=0, index = range(math.floor(logSeries.min()), math.ceil(logSeries.max())), columns=['% Total Seqs'])

for floor in histTotalSeqDF.index:
    histTotalSeqDF.loc[floor] = freqAbundTable[(freqAbundTable['Log Rel Abund'] >= floor) & (freqAbundTable['Log Rel Abund'] < floor+1)]['% Total Seqs'].sum()

# Plot as a bar chart
fig2, ax2 = plt.subplots()
fig2 = plt.bar(histTotalSeqDF.index, histTotalSeqDF['% Total Seqs'], align='edge')
ax2.set_xlabel('Log10 relative abundance of OTU')
ax2.set_ylabel('% Total Sequences')
plt.savefig(deblurDir+'/Total Sequences vs. Abundance.png')

#%%#############################################################################
### Taken together, these data suggest that a relative abundance threshold of 
### 0.001 will capture > 90% of the total sequences. Around 10% of all sequence
### variants account for this 90%.
###
### Let's simplify the OTU table and write to file
################################################################################

## Absolute abundances
# Replace all relative abundances below the threshhold and drop all-0 rows
simpleOtuTable = otuTable.copy()
simpleOtuTable[freqTable < 0.001] = 0
simpleOtuTable = simpleOtuTable.loc[~(simpleOtuTable==0).all(axis=1)]

## Relative abundances
# Replace all relative abundances below the threshhold and drop all-0 rows
simpleFreqTable = freqTable.copy()
simpleFreqTable[simpleFreqTable < 0.001] = 0
simpleFreqTable = simpleFreqTable.loc[~(simpleFreqTable==0).all(axis=1)]

# Write to file
simpleOtuTable.to_csv(deblurDir+'/simpleOtuTable.csv')
simpleFreqTable.to_csv(deblurDir+'/simpleFreqTable.csv')