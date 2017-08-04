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

import pandas as pd

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
freqTable = otuTable / otuTable.sum(axis=0)

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
