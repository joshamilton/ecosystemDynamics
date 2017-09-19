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