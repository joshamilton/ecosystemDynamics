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

import pandas as pd

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

# Identify the sample with the fewest replicates
minSample = otuTable.sum().idxmin()

# Compute relative abundance of the otuTable
relAbundTable = otuTable / otuTable.sum(axis=0)

# Return the smallest (non-zero) relative abundance in the shallowest sample
cutoff =  relAbundTable[relAbundTable > 0][minSample].min()