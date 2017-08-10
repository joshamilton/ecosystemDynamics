###############################################################################
# relAbund.py
# Copyright (c) 2017, Joshua J Hamilton
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Given a set of filtering criterion, create a simplified OTU table
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

deblurDir = '../results/deblur' # BIOM file from deblurring


### Cutoffs. he previous analyses suggest that a relative abundance threshold 
### of  0.003 and a persistence of 0.1 will capture > 80% of all sequences, 
### while retaining over 50% of sequences in each sample. There will be ~ 100
### OTUs.
persist = 0.1
abund = 0.003

#%%#############################################################################
### Import otu table and calculate relative abundances
################################################################################

# Import table
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)

# Compute counts to frequencies
relAbundTable = otuTable / otuTable.sum(axis=0)

#%%#############################################################################
### Compute a filtered table of average relative abundances
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

# Simplify average relative abundances
simpleAvgRelAbundTable = avgRelAbundTable.copy()
simpleAvgRelAbundTable[simpleAvgRelAbundTable < abund] = 0
simpleAvgRelAbundTable = simpleAvgRelAbundTable.loc[~(simpleAvgRelAbundTable==0).all(axis=1)]

## Further filter based on persistence
simpleAvgRelAbundTable['Persistence'] = simpleAvgRelAbundTable.astype(bool).sum(axis=1) / len(simpleAvgRelAbundTable.columns)
simpleAvgRelAbundTable = simpleAvgRelAbundTable[simpleAvgRelAbundTable.Persistence > persist]
simpleAvgRelAbundTable = simpleAvgRelAbundTable.drop('Persistence', axis=1)

# Write to file
simpleAvgRelAbundTable.to_csv(deblurDir+'/simpleAvgRelAbundTable.csv')