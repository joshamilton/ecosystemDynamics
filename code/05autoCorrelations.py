###############################################################################
# transformationStationarity.py
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

from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import skbio.stats.composition as ssm
import statsmodels.graphics.tsaplots as sgt

#%%#############################################################################
### Define folder structure
################################################################################

deblurDir = '../results/deblur'
corrDir = '../results/auto-correlations'

if not os.path.exists(corrDir):
    os.makedirs(corrDir)
      
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

# Retrive the list of sample names and remove the '.R2'
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

## Finally, replace zeros with small values
# Filter on zero threshold - when performing log transformations, these values
# will be replaced with non-zeros from a proper distribution
avgRelAbundTable[avgRelAbundTable < zero] = 0
# Multiplicative replacement destroys the dataframe, so before performing the
# operation, capture the index and columns
indexList = avgRelAbundTable.index
columnList = avgRelAbundTable.columns
# Replacment
avgRelAbundTable = (ssm.multiplicative_replacement(avgRelAbundTable.T, delta=0.00065)).T
# Restore the dataframe
avgRelAbundTable = pd.DataFrame(avgRelAbundTable, index=indexList, columns=columnList)

# Write to file
avgRelAbundTable.to_csv(deblurDir+'/simpleOtuTable.csv')

#%%#############################################################################
### We're doing a time-series analysis, need to split otuTable into separate
### tables for each year (using a dict). Then we need to link samples with 
### dates and update the otuTable to include times instead of dates
### Create a dataframe linking sample names to dates and days since the
### earliest time
################################################################################

otuTableDict = {}

for sample in sampleList:
    lake = sample[0:3]
    year = sample[3:5]
    tempOtuTable = avgRelAbundTable.filter(regex=(lake+'\d{2}[A-Z]{3}'+year))

    timeDF = pd.DataFrame(0, index=tempOtuTable.columns, columns=['Date', 'Time'])
    
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
    tempOtuTable.to_csv(deblurDir+'/otuTable-'+sample+'.csv')

#%%#############################################################################
### Plot auto-correlation and partial auto-correlation
### This requires interpolating between time-points
### A variety of interpolation methods exist, and there are two methods for
### testing stationarity.
### Interpolation: Linear, Quadratic, Cubic, etc
### Stationarity: ADF, KPSS
################################################################################

# Skip 'quadratic', 'cubic', and 'krogh' interp b/c their non-monotonic and result in negative OTU values
interpList = ['linear', 'slinear', 'pchip', 'pad', 'ignore']

for interp in interpList:    
    for sample in sampleList:

        # Create the otu table with missing values
        otuTable = otuTableDict[sample]
        otuTable = otuTable.reindex(columns=list(range(otuTable.columns[0], otuTable.columns[-1]+1, 1)))
        otuList = otuTable.index

        # Interpolate within the otu table
        if interp == 'pad': # fill with most recent value
            otuTable = otuTable.fillna(method='pad', axis=1)
        elif interp == 'ignore': # dorp all 'nan' columns
            otuTable = otuTable.dropna(axis=1)
            otuTable.columns = list(range(0, len(otuTable.columns)))
        else: # otherwise interpolate usign the specified method
            otuTable = otuTable.interpolate(method=interp, axis=1)

        for otu in otuList:
            # Extract the time-series of interest and log-transform
            timeSeries = np.log(otuTable.loc[otu])
            
            # Perform first-differencing
            firstSeries = np.diff(timeSeries, n=1)
        
            try: 
                # Plot correlation coefficients
                # Of the undifferenced time-series
                plt.figure(otu)
                sgt.plot_acf(timeSeries, zero=False, lags=10)
                plt.xlabel('Lag')
                plt.ylabel('Auto-correlation')
                plt.savefig(corrDir+'/'+str(otu)+'-'+sample+'-'+interp+'-acf.png')
                plt.close()
                
                # Of the differenced time-series
                plt.figure(otu)
                sgt.plot_acf(firstSeries, zero=False, lags=10)
                plt.xlabel('Lag')
                plt.ylabel('Auto-correlation')
                plt.savefig(corrDir+'/'+str(otu)+'-'+sample+'-'+interp+'-diffenced-acf.png')
                plt.close()
            
                # Plot partialcorrelation coefficients
                # Of the undifferenced time-series
                plt.figure(otu)
                sgt.plot_pacf(timeSeries, zero=False, lags=10)
                plt.xlabel('Lag')
                plt.ylabel('Auto-correlation')
                plt.savefig(corrDir+'/'+str(otu)+'-'+sample+'-'+interp+'-pacf.png')
                plt.close()
                
                # Of the differenced time-series
                plt.figure(otu)
                sgt.plot_pacf(firstSeries, zero=False, lags=10)
                plt.xlabel('Lag')
                plt.ylabel('Partial auto-correlation')
                plt.savefig(corrDir+'/'+str(otu)+'-'+sample+'-'+interp+'-diffenced-pacf.png')
                plt.close()

            except:
                break
            