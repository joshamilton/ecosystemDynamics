###############################################################################
# createOtuTables.py
# Copyright (c) 2017, Joshua J Hamilton
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Create OTU tables which serve as the basis for analysis.
# First, a simplified OTU table is created in which OTUs that don't satisfy
# our filtering criteria are removed.
# Then, the OTU table is split apart by year, and OTU tables for each year
# are created.
# Finally, we create log-transformed OTU tables. For this, two tables are
# created. One in which zeros have been replaced with small values, and a 
# second in which the data have been log-transformed.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from datetime import datetime
import numpy as np
import os
import pandas as pd
import re
import skbio.stats.composition as ssm

#%%#############################################################################
### Define folder structure
################################################################################

deblurDir = '../results/deblur'
forecastDir = '../results/forecasting'

if not os.path.exists(forecastDir):
    os.makedirs(forecastDir)
      
sampleList = ['TBE07', 'TBE08']

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

# Retrieve the list of sample names and remove the '.R2'
sampleNames = list(otuTable.columns[::2])
sampleNames = [re.sub('.R1', '', sample) for sample in sampleNames]
avgRelAbundTable.columns = sampleNames

#%%#############################################################################
### Filter based on persistence, abundance, and max. abundance
### Create dummy OTU
################################################################################

# Filter table based on persistence, and create a dummy OTU to represent the 
# remaining sequences
zero = 0.001
persist = 0.15
abund = 0.003
bloom = 0.1

# Calculate max abundance
avgRelAbundTable['Bloom'] = avgRelAbundTable.max(axis=1)

# Calculate persistence (based on abundance)
avgRelAbundTable['Persistence'] = (avgRelAbundTable.drop(['Bloom'], axis=1) > abund).sum(axis=1) / len(avgRelAbundTable.drop(['Bloom'], axis=1).columns)

# Filter on either of these criteria
avgRelAbundTable = avgRelAbundTable[(avgRelAbundTable.Persistence > persist) | (avgRelAbundTable.Bloom > bloom)]
avgRelAbundTable = avgRelAbundTable.drop('Persistence', axis=1)
avgRelAbundTable = avgRelAbundTable.drop('Bloom', axis=1)

# Add a dummy OTU containing the remaining sequences
# This operation reindexes the dataframe, so retain original mapping
mapDict = {}
curIter = 0
for otu in avgRelAbundTable.index:
    mapDict[curIter] = otu
    curIter += 1
avgRelAbundTable = avgRelAbundTable.append(1 - avgRelAbundTable.sum(), ignore_index=True)
mapDict[len(mapDict)] = 'Other'

# Then restore the original index
indexList = avgRelAbundTable.index
indexList = [mapDict[index] for index in indexList]
avgRelAbundTable.index = indexList

# Write to file
avgRelAbundTable.to_csv(deblurDir+'/simpleOtuTable.csv')

#%%#############################################################################
### Link samples with dates and update the otuTable to include times instead 
### of dates
### Create a dataframe linking sample names to dates and days since the
### earliest time
################################################################################

sortedOtuTable = avgRelAbundTable.copy()
timeDF = pd.DataFrame(0, index=sortedOtuTable.columns, columns=['Date', 'Time'])

# Convert sample names to dates
for curDate in timeDF.index:
    timeDF.set_value(curDate, 'Date', datetime.strptime(curDate[3:], '%d%b%y').date())

startTime = timeDF['Date'].min()

for curDate in timeDF.index:
    timeDF.set_value(curDate, 'Time', (timeDF.loc[curDate, 'Date'] - startTime).days)

# Rearrange the columns in the OTU table to reflect temporal ordering
columnList = sortedOtuTable.columns
columnList = [timeDF.loc[column, 'Time'] for column in columnList]
sortedOtuTable.columns = columnList
sortedOtuTable = sortedOtuTable.sort_index(axis=1)

# Write to file
sortedOtuTable.to_csv(forecastDir+'/simpleOtuTable-sorted.csv')

## Finally, replace zeros with small values

# Filter on zero threshold - when performing log transformations, these values
# will be replaced with non-zeros from a proper distribution
logOtuTable = sortedOtuTable.copy()
logOtuTable[logOtuTable < zero] = 0

# Multiplicative replacement destroys the dataframe, so before performing the
# operation, capture the index and columns
indexList = logOtuTable.index
columnList = logOtuTable.columns
# Replacment
logOtuTable = (ssm.multiplicative_replacement(logOtuTable.T, delta=0.00065)).T
# Restore the dataframe
logOtuTable = pd.DataFrame(logOtuTable, index=indexList, columns=columnList)

# Write to file
logOtuTable.to_csv(forecastDir+'/simpleOtuTable-sorted-NoZeros.csv')

# Compute log relative abundance and write to file
logOtuTable = np.log10(logOtuTable)
logOtuTable.to_csv(forecastDir+'/simpleOtuTable-sorted-NoZeros-Log.csv')

#%%#############################################################################
### Split otuTable into separate tables for each year (using a dict). Then we 
### need to link samples with dates and update the otuTable to include times 
### instead of dates.
### Create a dataframe linking sample names to dates and days since the
### earliest time
################################################################################

otuTableDict = {}

for sample in sampleList:
    lake = sample[0:3]
    year = sample[3:5]
    tempOtuTable = avgRelAbundTable.filter(regex=(lake+'\d{2}[A-Z]{3}'+year))

    timeDF = pd.DataFrame(0, index=tempOtuTable.columns, columns=['Date', 'Time'])
    
    # Write to file
    tempOtuTable.to_csv(deblurDir+'/simpleOtuTable-'+sample+'.csv')
    
    # Convert sample names to dates
    for curDate in timeDF.index:
        timeDF.set_value(curDate, 'Date', datetime.strptime(curDate[3:], '%d%b%y').date())
    
    startTime = timeDF['Date'].min()
    
    for curDate in timeDF.index:
        timeDF.set_value(curDate, 'Time', (timeDF.loc[curDate, 'Date'] - startTime).days)
    
    # Rearrange the columns in the OTU table to reflect temporal ordering
    columnList = tempOtuTable.columns
    columnList = [timeDF.loc[column, 'Time'] for column in columnList]
    tempOtuTable.columns = columnList
    tempOtuTable = tempOtuTable.sort_index(axis=1)
    
    otuTableDict[sample] = tempOtuTable

    # Write to file
    tempOtuTable.to_csv(forecastDir+'/simpleOtuTable-'+sample+'-sorted.csv')

    ## Finally, replace zeros with small values
    # Filter on zero threshold - when performing log transformations, these values
    # will be replaced with non-zeros from a proper distribution
    avgRelAbundTable[avgRelAbundTable < zero] = 0
    # Multiplicative replacement destroys the dataframe, so before performing the
    # operation, capture the index and columns
    indexList = tempOtuTable.index
    columnList = tempOtuTable.columns
    # Replacment
    tempOtuTable = (ssm.multiplicative_replacement(tempOtuTable.T, delta=0.00065)).T
    # Restore the dataframe
    tempOtuTable = pd.DataFrame(tempOtuTable, index=indexList, columns=columnList)
    
    # Write to file
    tempOtuTable.to_csv(forecastDir+'/simpleOtuTable-'+sample+'-sorted-NoZeros.csv')
    
    # Compute log relative abundance and write to file
    tempOtuTable = np.log10(tempOtuTable)
    tempOtuTable.to_csv(forecastDir+'/simpleOtuTable-'+sample+'-sorted-NoZeros-Log.csv')