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

deblurDir = '../results/deblur'

#%%#############################################################################
### Import otu table and calculate average relative abundances
################################################################################

# Import table
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', sep=',', index_col=0)

# Convert counts to relative abundances
relAbundTable = otuTable / otuTable.sum(axis=0)

# Convert to a table of average relative abundances
otuTable = otuTable.reindex_axis(sorted(otuTable.columns), axis=1)

# Compute counts to relative abundances
avgRelAbundTable = otuTable / otuTable.sum(axis=0)

# Compute average for each sample (every 2 columns)
avgRelAbundTable = avgRelAbundTable.groupby(np.arange(len(avgRelAbundTable.columns))//2, axis=1).mean()

# Retrive the list of sample names and remove the '.R2'
sampleNames = list(otuTable.columns[::2])
sampleNames = [re.sub('.R1', '', sample) for sample in sampleNames]
avgRelAbundTable.columns = sampleNames

#%%#############################################################################
### Consider a range of relative abundances and persistences
### For each pairing, calculate the % of total sequences retained, as well as
### the minimum % retained across all samples
### Construct a 3-D plot of these data to help choose cutoffs
################################################################################

zero = 0.001
abundRange = np.linspace(0.001, 0.01, 10, endpoint=True)
persistRange = np.linspace(0, 0.5, 11, endpoint=True)
bloomRange = np.linspace(0, 0.15, 16, endpoint=True)

#abundPersistIndex = pd.MultiIndex.from_product([abundRange, persistRange], names=['Abundance', 'Persistence'])
recoveryTable = pd.DataFrame(0, index=abundRange, columns=persistRange)
worstTable = pd.DataFrame(0, index=abundRange, columns=persistRange)
totalOtuTable = pd.DataFrame(0, index=abundRange, columns=persistRange)
log10TotalOtuTable = pd.DataFrame(0, index=abundRange, columns=persistRange)

# Loop over all combinations
# For each, filter the dataset based first on abundance, then persistence
for bloom in bloomRange:
    for abund in abundRange:
        for persist in persistRange:
    
            # Make a copy of master OTU table
            simpleAvgRelAbundTable = avgRelAbundTable.copy()
            
            # Filter on zero threshold
            simpleAvgRelAbundTable[simpleAvgRelAbundTable < zero] = 0
            
            # Calculate max abundance
            simpleAvgRelAbundTable['Bloom'] = simpleAvgRelAbundTable.max(axis=1)
            
            # Calculate persistence (based on abundance)
            simpleAvgRelAbundTable['Persistence'] = (simpleAvgRelAbundTable.drop(['Bloom'], axis=1) > abund).sum(axis=1) / len(simpleAvgRelAbundTable.drop(['Bloom'], axis=1).columns)

            # Filter on either of these criteria
            simpleAvgRelAbundTable = simpleAvgRelAbundTable[(simpleAvgRelAbundTable.Persistence > persist) | (simpleAvgRelAbundTable.Bloom > bloom)]
            simpleAvgRelAbundTable = simpleAvgRelAbundTable.drop('Persistence', axis=1)
            simpleAvgRelAbundTable = simpleAvgRelAbundTable.drop('Bloom', axis=1)
            
            # Update survival data frame
            recoveryTable.loc[abund, persist] = simpleAvgRelAbundTable.sum().sum() / len(simpleAvgRelAbundTable.columns)
            worstTable.loc[abund, persist] = simpleAvgRelAbundTable.sum(axis=0).min()
            totalOtuTable.loc[abund, persist] = len(simpleAvgRelAbundTable.index)
            log10TotalOtuTable.loc[abund, persist] = np.log10(len(simpleAvgRelAbundTable.index))
        
    # Write to file
    recoveryTable.to_csv(deblurDir+'/persistAbund-TotalRecovery-'+str(bloom)+'.csv')
    worstTable.to_csv(deblurDir+'/persistAbund-WorstRecovery-'+str(bloom)+'.csv')
    totalOtuTable.to_csv(deblurDir+'/persistAbund-TotalOTUs-'+str(bloom)+'.csv')