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

#%%#############################################################################
### Import otu table and calculate relative abundances
################################################################################

# Import table
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)

# Compute counts to frequencies
relAbundTable = otuTable / otuTable.sum(axis=0)

#%%#############################################################################
### Taken together, the previous analyses suggest that a relative abundance
### threshold of  0.001 will capture > 90% of the total sequences and a
### persistence threshold of 0.7 will capture > 80% of the sequences.
###
### Let's simplify the OTU table and write to file
################################################################################

## Absolute abundances
# Replace all relative abundances below the threshhold and drop all-0 rows
simpleOtuTable = otuTable.copy()
simpleOtuTable[relAbundTable < 0.001] = 0
simpleOtuTable = simpleOtuTable.loc[~(simpleOtuTable==0).all(axis=1)]

## Relative abundances
# Replace all relative abundances below the threshhold and drop all-0 rows
simpleRelAbundTable = relAbundTable.copy()
simpleRelAbundTable[simpleRelAbundTable < 0.001] = 0
simpleRelAbundTable = simpleRelAbundTable.loc[~(simpleRelAbundTable==0).all(axis=1)]

## Further filter based on persistence
simpleOtuTable['Persistence'] = simpleOtuTable.astype(bool).sum(axis=1) / len(simpleOtuTable.columns)
simpleOtuTable = simpleOtuTable[simpleOtuTable.Persistence < 0.7]
simpleOtuTable = simpleOtuTable.drop('Persistence', axis=1)

simpleRelAbundTable['Persistence'] = simpleRelAbundTable.astype(bool).sum(axis=1) / len(simpleRelAbundTable.columns)
simpleRelAbundTable = simpleRelAbundTable[simpleRelAbundTable.Persistence > 0.7]
simpleRelAbundTable = simpleRelAbundTable.drop('Persistence', axis=1)

# Write to file
simpleOtuTable.to_csv(deblurDir+'/simpleOtuTable.csv')
simpleRelAbundTable.to_csv(deblurDir+'/simpleRelAbundTable.csv')

#%%#############################################################################
### Taken together, the previous analyses suggest that a relative abundance
### threshold of  0.001 will capture > 90% of the total sequences. Around 10% 
### of all sequence variants account for this 90%.
###
### Also compute a table of average relative abundances
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
simpleAvgRelAbundTable[simpleAvgRelAbundTable < 0.001] = 0
simpleAvgRelAbundTable = simpleAvgRelAbundTable.loc[~(simpleAvgRelAbundTable==0).all(axis=1)]

## Further filter based on persistence
simpleAvgRelAbundTable['Persistence'] = simpleAvgRelAbundTable.astype(bool).sum(axis=1) / len(simpleAvgRelAbundTable.columns)
simpleAvgRelAbundTable = simpleAvgRelAbundTable[simpleAvgRelAbundTable.Persistence > 0.7]
simpleAvgRelAbundTable = simpleAvgRelAbundTable.drop('Persistence', axis=1)

# Write to file
simpleAvgRelAbundTable.to_csv(deblurDir+'/simpleAvgRelAbundTable.csv')