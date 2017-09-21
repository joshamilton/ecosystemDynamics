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
import os
import mannKendall as mk
import pandas as pd
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
### Import renormalized otu table. Filter based on persistence and abundance, 
### and create dummy OTU.
################################################################################

# Import table
otuTable = pd.read_csv(deblurDir+'/simpleAvgRelAbundTable-renorm.csv', sep=',', index_col=0)

# Filter table based on persistence, and create a dummy OTU to represent the 
# remaining sequences
persist = 0.15
abund = 0.1
otuTable['Persistence'] = otuTable.astype(bool).sum(axis=1) / len(otuTable.columns)
otuTable = otuTable[(otuTable.Persistence > persist) | (otuTable.drop('Persistence', axis=1).max(axis=1) > abund)]
otuTable = otuTable.drop('Persistence', axis=1)

# Add a dummy OTU containing the remaining sequences
# This operation reindexes the dataframe, so retain original mapping
mapDict = {}
curIter = 0
for otu in otuTable.index:
    mapDict[curIter] = otu
    curIter += 1
otuTable = otuTable.append(1 - otuTable.sum(), ignore_index=True)
mapDict[len(mapDict)] = 'Other'

# Then restore the original index
indexList = otuTable.index
indexList = [mapDict[index] for index in indexList]
otuTable.index = indexList

# Write to file
otuTable.to_csv(deblurDir+'/simpleAvgRelAbundTable-renorm-filterPersistAbund.csv')

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
    tempOtuTable = otuTable.filter(regex=(lake+'\d{2}[A-Z]{3}'+year))

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
### Construct the time-series for each OTU and plot it
### Perform statistical tests of stationarity for each time-series
################################################################################

for sample in sampleList:
    tempOtuTable = otuTableDict[sample]
    for curOtu in tempOtuTable.index:
        
        plt.figure(curOtu)
        plt.plot(tempOtuTable.loc[curOtu])
        plt.xlabel('Time')
        plt.ylabel('Relative Abundance')
        plt.savefig(resultsDir+'/'+str(curOtu)+'-'+sample+'.png')
        plt.close()

#%%#############################################################################
### Perform tests of stationarity on each time-series
### This requires interpolating between time-points
### A variety of interpolation methods exist, and there are two methods for
### testing stationarity.
### Interpolation: Linear, Quadratic, Cubic, etc
### Stationarity: ADF, KPSS
################################################################################
# Test samples TBE07 OTUs 1 and 4

interpList = ['linear', 'slinear', 'quadratic', 'cubic', 'krogh', 'pchip']
testList = ['KM-Trend', 'ADF-Mean', 'ADF-Trend', 'KPSS-Mean', 'KPSS-Trend']
otuList = otuTable.index[0:10]

for sample in sampleList:

    # Create the otu table with missing values
    otuTable = otuTableDict[sample]
    otuTable = otuTable.reindex(columns=list(range(otuTable.columns[0], otuTable.columns[-1]+1, 1)))

    for interp in interpList:
        # Interpolate within the otu table
        otuTable = otuTable.interpolate(method=interp, axis=1)
        
        # Create the data frame to store results (p-values from stat tests)
        # Index is all (sample, otu) pairs. Columns are all statistical tests pairs
        multiIndex = pd.MultiIndex.from_product([sampleList, otuList])
        resultsDF = pd.DataFrame(0.0, index=multiIndex, columns=testList)
        
        for otu in otuList:
            # Extract the time-series of interest
            timeSeries = otuTable.loc[otu]
        
            # Perform statistical tests

            # Mann-Kendall test: presence of a trend
            # Unavailable in Python, use function defined here
            # http://michaelpaulschramm.com/simple-time-series-trend-analysis/
            
            [trend, h, p, z] = mk.mk_test(timeSeries, alpha=0.05)
            resultsDF.set_value((sample, otu), 'KM-Trend', p)
            
            # Stationary around a mean: ADF
            [adf, pval, lags, nobs, crit] = sts.adfuller(timeSeries, regression='c', autolag=None, store=False)
            resultsDF.set_value((sample, otu), 'ADF-Mean', adf)
            
            # Stationary around a mean: KPSS
            [kpss, pval, lags, crit] = sts.kpss(timeSeries, regression='c', lags=None, store=False)
            resultsDF.set_value((sample, otu), 'KPSS-Mean', kpss)
            
            # Stationary around a trend: ADF
            [adf, pval, lags, nobs, crit] = sts.adfuller(timeSeries, regression='ct', autolag=None, store=False)
            resultsDF.set_value((sample, otu), 'ADF-Trend', adf)
            
            # Stationary around a trend: KPSS
            [kpss, pval, lags, crit] = sts.kpss(timeSeries, regression='ct', lags=None, store=False)
            resultsDF.set_value((sample, otu), 'KPSS-Trend', kpss)
        
        resultsDF.to_csv(resultsDir+'/trend-stationarity-'+interp+'Interp.csv')