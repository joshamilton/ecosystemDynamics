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
import numpy as np
import os
import pandas as pd
import statsmodels.graphics.tsaplots as sgt

#%%#############################################################################
### Define folder structure
################################################################################

forecastDir = '../results/forecasting'
corrDir = '../results/auto-correlations'

if not os.path.exists(corrDir):
    os.makedirs(corrDir)
      
sampleList = ['TBE07', 'TBE08']

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
            # Extract the time-series of interest and log-transform
            timeSeries = otuTable.loc[otu]
            
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
            