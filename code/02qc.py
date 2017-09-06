###############################################################################
# Quality Control.py
# Copyright (c) 2017, Joshua J Hamilton
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# For deblurred sequences, count the number of sequences in each file
# To look for outliers
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import os
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#%%#############################################################################
### Define folder structure
################################################################################
    
deblurDir = '../results/deblur'
outputDir = '../results/qc'

if not os.path.exists(outputDir):
    os.makedirs(outputDir)

#%%#############################################################################
### Construct a frequency distribution of sequence counts to determine if
### any samples are outliers.
################################################################################

# Read in the otuTable and count total number of sequences per sample
otuTable = pd.read_csv(deblurDir+'/otuTable.csv', index_col=0)

# Create a new series with total sequence counts and write to file
seqCountSeries = otuTable.sum()
seqCountSeries.to_csv(outputDir+'/numSeqPerSample.csv')

# Convert to log10 and plot a histogram with integer bin size
logSeries = np.log10(list(seqCountSeries))
plt.hist(logSeries, bins = range(math.floor(logSeries.min()), math.ceil(logSeries.max())+1), normed=True)
plt.xlabel('Log10 Sequence Count')
plt.ylabel('Frequency of Samples')
plt.savefig(outputDir+'/numSeqPerSample.png')

#%%#############################################################################
### These data suggest that most samples have > 1000 sequences
### Exclude any sample and its replicates with < 1000 sequences
### Create and store a list of these filenames for exclusion from further analysis
### Also exlude from otuTable
################################################################################

excludeList = ['FBE02AUG07.R1', 'HKE28JUN07.R1', 'HKE28JUN07.R2', 
               'MAH06JUL07.R1', 'MAH06JUL07.R2', 'NSH03JUL09.R1', 
               'NSH03JUL09.R2', 'SSH27JUL07.R1', 'SSH27JUL07.R2', 
               'TBE13JUN07.R1', 'TBE13JUN07.R2', 'WSE02OCT07.R1', 
               'WSE02OCT07.R2', 'WSH03JUL07.R1', 'WSH03JUL07.R2', 
               'WSH25JUN07.R1', 'WSH25JUN07.R2']

with open(outputDir+'/excludeMe.csv', 'w') as outFile:
    for sample in excludeList:
        outFile.write(sample+'\n')

otuTable = otuTable.drop(excludeList, axis=1)
otuTable = otuTable.reindex_axis(sorted(otuTable.columns), axis=1)
otuTable.to_csv(deblurDir+'/otuTable.csv')

#%%#############################################################################
### Split otutable by layer - epi and hypo
### For each layer, drop all 0 entries and save to file
################################################################################

# Extract epi and hypo columns
epiSamples = ['CBE','FBE','HKE','MAE','NSE','NSE','SSE','TBE','TBE','WSE']
epiColumns = [col for col in otuTable.columns if col[:3] in epiSamples]

hypoSamples = ['CBH','FBH','HKH','MAH','NSH','NSH','SSH','TBH','TBH','WSH']
hypoColumns = [col for col in otuTable.columns if col[:3] in hypoSamples]

# Subset dataframe
epiOtuTable = otuTable[epiColumns]
hypoOtuTable = otuTable[hypoColumns]

# Drop OTUs missing from that layer
epiOtuTable = epiOtuTable.loc[~(epiOtuTable==0).all(axis=1)]
hypoOtuTable = hypoOtuTable.loc[~(hypoOtuTable==0).all(axis=1)]

# Write to file
epiOtuTable.to_csv(deblurDir+'/epiOtuTable.csv')
hypoOtuTable.to_csv(deblurDir+'/hypoOtuTable.csv')
