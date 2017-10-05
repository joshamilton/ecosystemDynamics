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
import math
import mannKendall as mk
import numpy as np
import os
import pandas as pd
import re
import scipy.stats as ss
import skbio.stats.composition as ssm
import statsmodels.tsa.stattools as sts

#%%#############################################################################
### Define folder structure
################################################################################

deblurDir = '../results/deblur'
resultsDir = '../results/stationarity'

if not os.path.exists(resultsDir):
    os.makedirs(resultsDir)
      
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
### Link samples with  dates and update the otuTable to include times instead 
### of dates
### Create a dataframe linking sample names to dates and days since the
### earliest time
################################################################################

tempOtuTable = avgRelAbundTable.copy()
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

# Write to file
tempOtuTable.to_csv(deblurDir+'/simpleOtuTable-sorted.csv')

#%%#############################################################################
### Perform tests of stationarity on each time-series
### This requires interpolating between time-points
### A variety of interpolation methods exist, and there are two methods for
### testing stationarity.
### Interpolation: Linear, Quadratic, Cubic, etc
### Stationarity: ADF, KPSS
################################################################################

# Skip 'quadratic', 'cubic', and 'krogh' interp b/c their non-monotonic and result in negative OTU values
interpList = ['linear', 'slinear', 'pchip']
testList = ['KM-Trend', 'ADF-Mean', 'ADF-Trend', 'KPSS-Mean', 'KPSS-Trend']

# Create the data frame to store results (p-values from stat tests)
# Index is all (sample, otu) pairs. Columns are all statistical tests pairs
resultsDF = pd.DataFrame(0.0, index=tempOtuTable.index, columns=testList)

for interp in interpList:    

    # Create the otu table with missing values
    tempOtuTable = tempOtuTable.reindex(columns=list(range(tempOtuTable.columns[0], tempOtuTable.columns[-1]+1, 1)))

    # Interpolate within the otu table
    tempOtuTable = tempOtuTable.interpolate(method=interp, axis=1)
            
    for curOtu in resultsDF.index:
        # Extract the time-series of interest and log-transform
        timeSeries = np.log(tempOtuTable.loc[curOtu])
        
        # Compute the first differences
        timeSeries = np.diff(timeSeries, n=1)
        
        # Plot the "time-series" of first differences
        plt.figure(otu)
        plt.plot(timeSeries)
        plt.xlabel('Time')
        plt.ylabel('First Difference')
        plt.savefig(resultsDir+'/'+str(curOtu)+'.png')
        plt.close()
    
        # Construct and plot a histogram of the first differences
        ## Specify number of bins using Freedman-Diaconis rule
        ## https://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule
        binSize = 2 * ss.iqr(timeSeries) / len(timeSeries)**(1/3)
        numBins = np.ceil((timeSeries.max() - timeSeries.min()) / binSize)
        if numBins == np.inf or math.isnan(numBins) or numBins > 100:
            numBins = 100
        numBins = int(numBins)

        plt.figure()
        n, bins, patches = plt.hist(timeSeries, numBins)
        plt.xlabel('First Difference')
        plt.ylabel('Count')
        plt.savefig(resultsDir+'/'+str(curOtu)+'-Hist.png')        
        plt.close()
    
        # Perform statistical tests
        # Mann-Kendall test: presence of a trend
        # Unavailable in Python, use function defined here
        # http://michaelpaulschramm.com/simple-time-series-trend-analysis/
        [trend, h, p, z] = mk.mk_test(timeSeries, alpha=0.05)
        resultsDF.set_value(curOtu, 'KM-Trend', p)
        
        # Stationary around a mean: ADF
        [adf, pval, lags, nobs, crit] = sts.adfuller(timeSeries, regression='c', autolag=None, store=False)
        resultsDF.set_value(curOtu, 'ADF-Mean', adf)
        
        # Stationary around a mean: KPSS
        [kpss, pval, lags, crit] = sts.kpss(timeSeries, regression='c', lags=None, store=False)
        resultsDF.set_value(curOtu, 'KPSS-Mean', kpss)
        
        # Stationary around a trend: ADF
        [adf, pval, lags, nobs, crit] = sts.adfuller(timeSeries, regression='ct', autolag=None, store=False)
        resultsDF.set_value(curOtu, 'ADF-Trend', adf)
        
        # Stationary around a trend: KPSS
        [kpss, pval, lags, crit] = sts.kpss(timeSeries, regression='ct', lags=None, store=False)
        resultsDF.set_value(curOtu, 'KPSS-Trend', kpss)
    
    resultsDF.to_csv(resultsDir+'/trend-stationarity-'+interp+'Interp.csv')