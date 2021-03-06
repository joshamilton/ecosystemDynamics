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

import matplotlib.pyplot as plt
import mannKendall as mk
import numpy as np
import os
import pandas as pd
import statsmodels.tsa.stattools as sts

#%%#############################################################################
### Define folder structure
################################################################################

forecastDir = '../results/forecasting'
resultsDir = '../results/stationarity'

if not os.path.exists(resultsDir):
    os.makedirs(resultsDir)
      
sampleList = ['TBE07', 'TBE08']

#%%#############################################################################
### Construct the time-series for each OTU and plot it
################################################################################

for sample in sampleList:
    
    otuTable = pd.read_csv(forecastDir+'/simpleOtuTable-'+sample+'-sorted-NoZeros-Log.csv', index_col=0)
    
    for curOtu in otuTable.index:
        
        plt.figure(curOtu)
        plt.plot(otuTable.loc[curOtu])
        plt.xlabel('Time')
        plt.ylabel('Log(Relative Abundance)')
        plt.savefig(resultsDir+'/'+str(curOtu)+'-'+sample+'.png')
        plt.close()
    
#%%#############################################################################
### Construct a histogram of the relative abundance data
################################################################################

for sample in sampleList:

    otuTable = pd.read_csv(forecastDir+'/simpleOtuTable-'+sample+'-sorted-NoZeros-Log.csv', index_col=0)
    
    for curOtu in otuTable.index:

        # Compute the time-series
        timeSeries = otuTable.loc[curOtu]
        
        # Construct a histogram using bins of size one
        binMin = int(np.floor(timeSeries.min()))
        binMax = int(np.ceil(timeSeries.max()))
        numBins = abs(binMin - binMax)

        plt.figure()
        plt.show()
        n, bins, patches = plt.hist(timeSeries, bins=numBins, range=(binMin, binMax))
        plt.xlabel('Log(Relative Abundance)')
        plt.ylabel('Count')
        plt.savefig(resultsDir+'/'+str(curOtu)+'-'+sample+'-Hist.png')        
        plt.close()
        
#%%#############################################################################
### Perform tests of stationarity on each time-series
### This requires interpolating between time-points
### A variety of interpolation methods exist, and there are two methods for
### testing stationarity.
### Interpolation: Linear, Quadratic, Cubic, etc
### Stationarity: ADF, KPSS
################################################################################

# Skip 'quadratic', 'cubic', and 'krogh' interp b/c their non-monotonic and result in negative OTU values
interpList = ['linear', 'slinear', 'pchip', 'pad', 'ignore']
testList = ['KM-Trend', 'ADF-Mean', 'ADF-Trend', 'KPSS-Mean', 'KPSS-Trend']

# Create the data frame to store results (p-values from stat tests)
# Index is all (sample, otu) pairs. Columns are all statistical tests pairs
otuList = otuTable.index
multiIndex = pd.MultiIndex.from_product([sampleList, otuList])
resultsDF = pd.DataFrame(0.0, index=multiIndex, columns=testList)

for interp in interpList:    
    for sample in sampleList:

        # Create the otu table with missing values
        otuTable = pd.read_csv(forecastDir+'/simpleOtuTable-'+sample+'-sorted-NoZeros-Log.csv', index_col=0)
        otuTable.columns = otuTable.columns.astype(int)
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