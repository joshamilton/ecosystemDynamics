###############################################################################
# compareReplicates.py
# Copyright (c) 2017, Joshua J Hamilton
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Given a set of replicates, find the minimum relative abundance which is 
# detected in both replicates.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import numpy as np
import pandas as pd
import re

#%%#############################################################################
### Define folder structure
################################################################################

deblurDir = '../results/deblur'

#%%#############################################################################
### Import otu table and calculate minimum relative abundance in each pair
### of replicates
################################################################################

# Import table
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)

# Sort on the basis of column names
otuTable = otuTable.reindex_axis(sorted(otuTable.columns), axis=1)

# Compute counts to relative abundances
relAbundTable = otuTable / otuTable.sum(axis=0)

# Compute minimum abundance of each set of replicates
minRelAbundTable = relAbundTable.groupby(np.arange(len(relAbundTable.columns))//2, axis=1).min()

# Retrieve the list of sample names and remove the '.R2'
sampleNames = list(otuTable.columns[::2])
sampleNames = [re.sub('.R1', '', sample) for sample in sampleNames]
minRelAbundTable.columns = sampleNames

# For each sample, find the minimal non-zero value
minInBoth = minRelAbundTable[minRelAbundTable>0].min(axis=0).min()

#%%#############################################################################
### Import otu table and calculate maximum relative abundance in each pair
### of replicates where one sample failed to detect the OTU
################################################################################

# However, this doesn't guarantee that an OTU with abundance above
# that value will get detected. Instead, we need to know the "worst case"
# relative abundance that is detected in one sample and not the other

# List to store the results
maxZeroList = []

# For each pair of replicates:
for curIndex in range(0, len(relAbundTable.columns)//2):
    tempAbundTable = relAbundTable[relAbundTable.columns[2*curIndex:2*(curIndex+1)]]

    # First drop all OTUs which are zero in both replicates
    tempAbundTable = tempAbundTable.loc[~(tempAbundTable==0).all(axis=1)]
    
    # Drop all OTUs which are non-zero in both replicates
    tempAbundTable = tempAbundTable.loc[~ (tempAbundTable != 0).all(axis=1)]
    
    # This will leave just those samples where an OTU was detected in on3
    # replicate but not the other
    
    # Of the remaining OTUs, return the max relative abundance
    maxZeroList.append(tempAbundTable.max().max())
    
# Convert into a pandas series and make a plot
maxZeroSeries = pd.Series(maxZeroList)

maxZeroSeries.plot(kind='box', title='Maximum OTU abundance', label='TBE')