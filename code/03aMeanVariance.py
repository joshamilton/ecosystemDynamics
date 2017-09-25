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

import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import statsmodels.stats.diagnostic as smd

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

# Retrieve the list of sample names and remove the '.R2'
sampleNames = list(otuTable.columns[::2])
sampleNames = [re.sub('.R1', '', sample) for sample in sampleNames]
avgRelAbundTable.columns = sampleNames

# Write to file
avgRelAbundTable.to_csv(deblurDir+'/avgRelAbundTable.csv')

#%%#############################################################################
### Check mean-variance scaling: visualization
################################################################################

meanSeries = np.log10(avgRelAbundTable.mean(axis=1))
varianceSeries = np.log10(avgRelAbundTable.var(axis=1))

# Drop all OTUs w/ zero abundance
meanSeries = meanSeries.loc[meanSeries != -np.inf]
varianceSeries = varianceSeries.loc[varianceSeries != -np.inf]

plt.figure(0) # create the figure
plt.scatter(meanSeries, varianceSeries)
plt.xlabel('Log10 Mean(OTU)')
plt.ylabel('Log10 Variance(OTU)')

plt.savefig(deblurDir+'/Mean-Variance Scaling.png')

#%%#############################################################################
### Check mean-variance scaling: statistical testing
################################################################################

# Null hypothesis is that the data are homoscedastic. That is, the variance is
# constant. This can be tested by looking at the residuals: their standard
# deviation should be constant and not depend on the x-value.

# The  Breusch–Pagan test can be used to test for heteroskedasticity given
# residuals from a fit. It is implemented in the statsmodels package.

# Consider a range of cutoff values
cutoffList = list(np.logspace(-3, -2, 21, base=10))
pDF = pd.DataFrame(0.0, index=cutoffList, columns=['P-value'])

for cutoff in cutoffList:

    # Make a copy of the table
    tempAvgRelAbundTable = copy.deepcopy(avgRelAbundTable)
    
    # Filter table, removing all OTUs with abundance < cutoff
    tempAvgRelAbundTable[tempAvgRelAbundTable < cutoff] = 0
    
    # Compute mean and variance of each OTU
    meanSeries = np.log10(tempAvgRelAbundTable.mean(axis=1))
    varianceSeries = np.log10(tempAvgRelAbundTable.var(axis=1))
    
    # Drop all OTUs w/ zero abundance
    meanSeries = meanSeries.loc[meanSeries != -np.inf]
    varianceSeries = varianceSeries.loc[varianceSeries != -np.inf]
    
    # Regress variance on the mean
    slope, intercept = np.polyfit(meanSeries, varianceSeries, 1)
    
    # Compute residuals
    fitDF = pd.DataFrame(0, index=meanSeries.index, columns=['Mean', 'Variance', 'Estimate', 'Residual'])
    fitDF['Mean'] = meanSeries
    fitDF['Variance'] = varianceSeries
    fitDF['Estimate'] = slope*fitDF['Mean'] + intercept
    fitDF['Residual'] = fitDF['Estimate'] - fitDF['Variance']
    
    plt.figure(cutoffList.index(cutoff)) # create the figure
    plt.scatter(fitDF['Mean'], fitDF['Residual'])
    plt.xlabel('Log10 Mean(OTU)')
    plt.ylabel('Residual')
    fileName = 'Residuals of Mean-Variance Scaling-'+str(cutoff)+'.png'
    plt.savefig(deblurDir+'/'+fileName)
    plt.close()
    
    lm, lm_pvalue, fvalue, f_pvalue = smd.het_breuschpagan(fitDF['Residual'], fitDF['Mean'].to_frame())
    
    pDF.set_value(cutoff, 'P-value', f_pvalue)
    
pDF.to_csv(deblurDir+'/Heteroskedasticity.csv')