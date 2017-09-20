###############################################################################
# scanPersistence.py
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

#%%#############################################################################
### Define folder structure
################################################################################

deblurDir = '../results/deblur'

#%%#############################################################################
### Import renormalized otu table
################################################################################

# Import table
otuTable = pd.read_csv(deblurDir+'/simpleAvgRelAbundTable-renorm.csv', sep=',', index_col=0)

#%%#############################################################################
### Consider a range of persistences
### Calculate the % of total sequences retained, as well as
### the minimum % retained across all samples
################################################################################

persistRange = list(range(0, 105, 5))
persistRange = [persist / 100 for persist in persistRange]

resultsDF = pd.DataFrame(0.0, index=persistRange, columns=['Frac Reads', 'Worst Reads', 'OTUs'])

# Loop over all combinations
# For each, filter the dataset based first on abundance, then persistence
for persist in persistRange:

    # Make a copy of master OTU table
    tempOtuTable = otuTable.copy()
    
    # Filter on persistence
    tempOtuTable['Persistence'] = tempOtuTable.astype(bool).sum(axis=1) / len(tempOtuTable.columns)
    tempOtuTable = tempOtuTable[tempOtuTable.Persistence > persist]
    tempOtuTable = tempOtuTable.drop('Persistence', axis=1)
    
    resultsDF.set_value(persist, 'Frac Reads', tempOtuTable.sum().sum() / len(tempOtuTable.columns))
    resultsDF.set_value(persist, 'Worst Reads', tempOtuTable.sum().min())
    resultsDF.set_value(persist, 'OTUs', len(tempOtuTable))
    
resultsDF.to_csv(deblurDir+'/persistenceRecovery.csv')

#%%#############################################################################
### Consider a range of relative abundances
### Retain all OTUs which have that abundance in at least one sample
### Calculate the % of total sequences retained, as well as
### the minimum % retained across all samples
################################################################################

abundRange = list(range(0, 105, 5))
abundRange = [abund / 100 for abund in abundRange]

resultsDF = pd.DataFrame(0.0, index=abundRange, columns=['Frac Reads', 'Worst Reads', 'OTUs'])
#abundRange = [abundRange[2]]
# Loop over all combinations
# For each, filter the dataset based first on abundance, then persistence
for abund in abundRange:

    # Make a copy of master OTU table
    tempOtuTable = otuTable.copy()
    
    # Filter on maximum abundance
    tempOtuTable['Abundance'] = tempOtuTable.max(axis=1)
    tempOtuTable = tempOtuTable[tempOtuTable.Abundance > abund]
    tempOtuTable = tempOtuTable.drop('Abundance', axis=1)
    
    resultsDF.set_value(abund, 'Frac Reads', tempOtuTable.sum().sum() / len(tempOtuTable.columns))
    resultsDF.set_value(abund, 'Worst Reads', tempOtuTable.sum().min())
    resultsDF.set_value(abund, 'OTUs', len(tempOtuTable))
    
resultsDF.to_csv(deblurDir+'/abundanceRecovery.csv')

#%%#############################################################################
### Consider a range of persistences and abundances
### Implement a double filter, in which taxa satisfying EITHER of the following
### criteria are retained:
###     The OTU persists across N samples
###     The OTU has abundance > 1 sample
### Calculate the % of total sequences retained, as well as
### the minimum % retained across all samples
################################################################################

persistRange = list(range(0, 105, 5))
persistRange = [persist / 100 for persist in persistRange]

abundRange = persistRange

recoveryTable = pd.DataFrame(0, index=abundRange, columns=persistRange)
worstTable = pd.DataFrame(0, index=abundRange, columns=persistRange)
totalOtuTable = pd.DataFrame(0, index=abundRange, columns=persistRange)
log10TotalOtuTable = pd.DataFrame(0, index=abundRange, columns=persistRange)

# Loop over all combinations
# For each, filter the dataset based first on abundance, then persistence

persistRange = [0.2]
abundRange = [0.1]
for abund in abundRange:
    for persist in persistRange:

        # Make a copy of master OTU table
        tempOtuTable = otuTable.copy()
        
        # Filter on persistence
        tempOtuTable['Persistence'] = tempOtuTable.astype(bool).sum(axis=1) / len(tempOtuTable.columns)
        tempOtuTable = tempOtuTable[(tempOtuTable.Persistence > persist) | (tempOtuTable.drop('Persistence', axis=1).max(axis=1) > abund)]
        tempOtuTable = tempOtuTable.drop('Persistence', axis=1)
        
        # Update survival data frame
        recoveryTable.loc[abund, persist] = tempOtuTable.sum().sum() / len(tempOtuTable.columns)
        worstTable.loc[abund, persist] = tempOtuTable.sum(axis=0).min()
        totalOtuTable.loc[abund, persist] = len(tempOtuTable.index)
        log10TotalOtuTable.loc[abund, persist] = np.log10(len(tempOtuTable.index))
        
# Write to file
recoveryTable.to_csv(deblurDir+'/persistAbund-TotalRecovery.csv')
worstTable.to_csv(deblurDir+'/persistAbund-WorstRecovery.csv')
totalOtuTable.to_csv(deblurDir+'/persistAbund-TotalOTUs.csv')