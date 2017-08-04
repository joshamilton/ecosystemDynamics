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

deblurDir = '../results/deblur' # BIOM file from deblurring

#%%#############################################################################
### Import otu table and visualize distribution of relative abundances
################################################################################

# Import table
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)

# Compute counts to frequencies
relAbundTable = otuTable / otuTable.sum(axis=0)

# Convert to a single series and remove 0 values
relAbundSeries = relAbundTable.unstack()
relAbundSeries = relAbundSeries[relAbundSeries!=0]

# Convert to log10 and plot a histogram with integer bin size
logRelAbundSeries = np.log10(relAbundSeries)
fig1, ax1 = plt.subplots()
fig1 = plt.hist(logRelAbundSeries, bins = range(math.floor(logRelAbundSeries.min()), math.ceil(logRelAbundSeries.max())+1), normed=True)
ax1.set_xlabel('Log10 relative abundance of OTU')
ax1.set_ylabel('% Total OTUs')
plt.savefig(deblurDir+'/Total OTUs vs. Abundance.png')

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
fig2, ax2 = plt.subplots()
fig2 = plt.bar(histTotalSeqDF.index, histTotalSeqDF['% Total Seqs'], align='edge')
ax2.set_xlabel('Log10 relative abundance of OTU')
ax2.set_ylabel('% Total Sequences')
plt.savefig(deblurDir+'/Total Sequences vs. Abundance.png')