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
### threshold of  0.001 will capture > 90% of the total sequences. Around 10% 
### of all sequence variants account for this 90%.
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

# Write to file
simpleAvgRelAbundTable.to_csv(deblurDir+'/simpleAvgRelAbundTable.csv')
