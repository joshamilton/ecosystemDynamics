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
import math
import mannKendall as mk
import numpy as np
import os
import pandas as pd
import re
import scipy.stats as ss
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


for interp in interpList:    
    
    # Create the otu table with missing values
    otuTable = pd.read_csv(forecastDir+'/simpleOtuTable-sorted-NoZeros-Log.csv', index_col=0)
    otuTable.columns = otuTable.columns.astype(int)
    otuTable = otuTable.reindex(columns=list(range(otuTable.columns[0], otuTable.columns[-1]+1, 1)))
    otuList = otuTable.index

    # Create the data frame to store results (p-values from stat tests)
    # Index is all otus. Columns are all statistical tests pairs    
    resultsDF = pd.DataFrame(0.0, index=otuList, columns=testList)

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
        timeSeries = otuTable.loc[otu]
        
        # Compute the first differences
        timeSeries = np.diff(timeSeries, n=1)
        
        # Plot the "time-series" of first differences
        plt.figure(otu)
        plt.plot(timeSeries)
        plt.xlabel('Time')
        plt.ylabel('First Difference')
        plt.savefig(resultsDir+'/'+str(otu)+'-'+interp+'.png')
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
        plt.savefig(resultsDir+'/'+str(otu)+'-'+interp+'-Hist.png')        
        plt.close()
    
        # Perform statistical tests
        # Mann-Kendall test: presence of a trend
        # Unavailable in Python, use function defined here
        # http://michaelpaulschramm.com/simple-time-series-trend-analysis/
        [trend, h, p, z] = mk.mk_test(timeSeries, alpha=0.05)
        resultsDF.set_value(otu, 'KM-Trend', p)
        
        # Stationary around a mean: ADF
        [adf, pval, lags, nobs, crit] = sts.adfuller(timeSeries, regression='c', autolag=None, store=False)
        resultsDF.set_value(otu, 'ADF-Mean', adf)
        
        # Stationary around a mean: KPSS
        [kpss, pval, lags, crit] = sts.kpss(timeSeries, regression='c', lags=None, store=False)
        resultsDF.set_value(otu, 'KPSS-Mean', kpss)
        
        # Stationary around a trend: ADF
        [adf, pval, lags, nobs, crit] = sts.adfuller(timeSeries, regression='ct', autolag=None, store=False)
        resultsDF.set_value(otu, 'ADF-Trend', adf)
        
        # Stationary around a trend: KPSS
        [kpss, pval, lags, crit] = sts.kpss(timeSeries, regression='ct', lags=None, store=False)
        resultsDF.set_value(otu, 'KPSS-Trend', kpss)
    
    resultsDF.to_csv(resultsDir+'/trend-stationarity-'+interp+'Interp.csv')