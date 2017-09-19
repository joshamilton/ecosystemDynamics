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

import pandas as pd

#%%#############################################################################
### Define folder structure
################################################################################

deblurDir = '../results/deblur' # BIOM file from deblurring
cutoff = 0.001

#%%#############################################################################
### Import otu table and calculate relative abundances
################################################################################

# Import table
avgRelAbundTable = pd.read_csv(deblurDir+'/avgRelAbundTable.csv', sep=',', index_col=0)

#%%#############################################################################
### Filter and renormalize
################################################################################

# Simplify average relative abundances
simpleAvgRelAbundTable = avgRelAbundTable.copy()
simpleAvgRelAbundTable[simpleAvgRelAbundTable < cutoff] = 0
simpleAvgRelAbundTable = simpleAvgRelAbundTable.loc[~(simpleAvgRelAbundTable==0).all(axis=1)]

# Renormalize
renormAvgRelAbundTable = simpleAvgRelAbundTable / simpleAvgRelAbundTable.sum(axis=0)

# Write to file
renormAvgRelAbundTable.to_csv(deblurDir+'/simpleAvgRelAbundTable-renorm.csv')

#%%#############################################################################
### Calculate recovery
################################################################################

recoveryDF = pd.DataFrame(0.0, index=avgRelAbundTable.columns, columns=['Reads', 'OTUs'])
recoveryDF['Reads'] = simpleAvgRelAbundTable.sum(axis=0)
recoveryDF['OTUs'] = renormAvgRelAbundTable.astype(bool).sum(axis=0) / avgRelAbundTable.astype(bool).sum(axis=0)