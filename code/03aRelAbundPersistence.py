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
import math
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import statsmodels.stats.diagnostic as smd

#%%#############################################################################
### Define folder structure
################################################################################

deblurDir = '../results/deblur' # BIOM file from deblurring

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
### Consider a range of relative abundances and persistences
### For each pairing, calculate the % of total sequences retained, as well as
### the minimum % retained across all samples
### Construct a 3-D plot of these data to help choose cutoffs
################################################################################

abundRange = list(np.logspace(-8, 0, 17, base=10))
persistRange = list(range(0, 105, 5))
persistRange = [persist / 100 for persist in persistRange]

#abundPersistIndex = pd.MultiIndex.from_product([abundRange, persistRange], names=['Abundance', 'Persistence'])
recoveryTable = pd.DataFrame(0, index=abundRange, columns=persistRange)
worstTable = pd.DataFrame(0, index=abundRange, columns=persistRange)
totalOtuTable = pd.DataFrame(0, index=abundRange, columns=persistRange)
log10TotalOtuTable = pd.DataFrame(0, index=abundRange, columns=persistRange)

# Loop over all combinations
# For each, filter the dataset based first on abundance, then persistence

for abund in abundRange:
    for persist in persistRange:

        # Make a copy of master OTU table
        simpleAvgRelAbundTable = avgRelAbundTable.copy()
        
        # Filter on abundance
        simpleAvgRelAbundTable[simpleAvgRelAbundTable < abund] = 0
        simpleAvgRelAbundTable = simpleAvgRelAbundTable.loc[~(simpleAvgRelAbundTable==0).all(axis=1)]
        
        # Filter on persistence
        simpleAvgRelAbundTable['Persistence'] = simpleAvgRelAbundTable.astype(bool).sum(axis=1) / len(simpleAvgRelAbundTable.columns)
        simpleAvgRelAbundTable = simpleAvgRelAbundTable[simpleAvgRelAbundTable.Persistence > persist]
        simpleAvgRelAbundTable = simpleAvgRelAbundTable.drop('Persistence', axis=1)
        
        # Update survival data frame
        recoveryTable.loc[abund, persist] = simpleAvgRelAbundTable.sum().sum() / len(simpleAvgRelAbundTable.columns)
        worstTable.loc[abund, persist] = simpleAvgRelAbundTable.sum(axis=0).min()
        totalOtuTable.loc[abund, persist] = len(simpleAvgRelAbundTable.index)
        log10TotalOtuTable.loc[abund, persist] = np.log10(len(simpleAvgRelAbundTable.index))
        
# Write to file
recoveryTable.to_csv(deblurDir+'/persistAbund-TotalRecovery.csv')
worstTable.to_csv(deblurDir+'/persistAbund-WorstRecovery.csv')
totalOtuTable.to_csv(deblurDir+'/persistAbund-TotalOTUs.csv')

#%%#############################################################################
### Visualize with contour plots
################################################################################

xRange = list(np.logspace(-8, 0, 9, base=10)) # steps for x axis. Use persist range for y axis.

plt.figure(0) # create the figure
plt.figure(figsize=(11,8.5))
plt.subplots_adjust(hspace=0.25)

### Fraction of sequence retained
plt.subplot(221) # first subplot
X, Y = np.meshgrid(recoveryTable.index, recoveryTable.columns)
CS = plt.contourf(X, Y, recoveryTable.values.transpose(), levels=persistRange, cmap=cm.gray)
plt.xscale('log') # log scale
plt.xticks(xRange) # ticks every order of magnitude
plt.yticks(persistRange) # ticks every 0.05
plt.title('Fraction Sequences Retained')
plt.xlabel('Relative Abundance')
plt.ylabel('Persistence')
# make a colorbar for the contour lines
cbar = plt.colorbar(CS, ticks=persistRange)


### Worst case retained
plt.subplot(222) # second subplot
X, Y = np.meshgrid(worstTable.index, worstTable.columns)
CS = plt.contourf(X, Y, worstTable.values.transpose(), levels=persistRange, cmap=cm.gray)
plt.xscale('log') # log scale
plt.xticks(xRange) # ticks every order of magnitude
plt.yticks(persistRange) # ticks every 0.05
plt.title('Frac. Seq. Retained (Worst)')
plt.xlabel('Relative Abundance')
plt.ylabel('Persistence')
# make a colorbar for the contour lines
cbar = plt.colorbar(CS, ticks=persistRange)


### Fraction of OTUs retained
fracOtuTable = np.log10(totalOtuTable / len(otuTable.index))
cmapRange = list(range(math.floor(fracOtuTable[fracOtuTable != fracOtuTable.loc[1, 0]].min().min()), 1, 1))
plt.subplot(223) # third subplot
X, Y = np.meshgrid(fracOtuTable.index, fracOtuTable.columns)
CS = plt.contourf(X, Y, fracOtuTable.values.transpose(), levels=cmapRange, cmap=cm.gray)
plt.xscale('log') # log scale
plt.xticks(xRange) # ticks every order of magnitude
plt.yticks(persistRange) # ticks every 0.05
plt.title('Fraction Remaining OTUs')
plt.xlabel('Relative Abundance')
plt.ylabel('Persistence')
# make a colorbar for the contour lines
cbar = plt.colorbar(CS, ticks=cmapRange)


### Total OTUs remaining
yRange = list(range(0, math.ceil(np.log10(len(otuTable.index)))+1)) # log scale
plt.subplot(224) # fourth subplot
X, Y = np.meshgrid(log10TotalOtuTable.index, log10TotalOtuTable.columns)
CS = plt.contourf(X, Y, log10TotalOtuTable.values.transpose(), levels=yRange, cmap=cm.gray)
plt.xscale('log') # log scale
plt.xticks(xRange) # ticks every order of magnitude
plt.yticks(persistRange) # ticks every 0.05
plt.title('Log10 Remaining OTUs')
plt.xlabel('Relative Abundance')
plt.ylabel('Persistence')
# make a colorbar for the contour lines
cbar = plt.colorbar(CS, ticks=yRange)

plt.savefig(deblurDir+'/Filter on Persist and Rel Abund.png')

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
    
    lm, lm_pvalue, fvalue, f_pvalue = smd.het_breuschpagan(fitDF['Residual'], fitDF['Mean'].to_frame())
    
    pDF.set_value(cutoff, 'P-value', f_pvalue)
    
pDF.to_csv(deblurDir+'/Heteroskedasticity.csv')