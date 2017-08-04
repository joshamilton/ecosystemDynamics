###############################################################################
# demultiplexAndQC.py
# Copyright (c) 2017, Joshua J Hamilton
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Demultiplex and QC raw FASTA files
# Count number of reads in each FASTA file
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
import os
import math
import matplotlib.pyplot as plt
import numpy as np

#%%#############################################################################
### Define folder structure
################################################################################
    
fastqDir = '../data/rawSequences' # raw fastq files
qcDir = '../data/demultiplexedAndQCd' # demultiplexed and QCd
outputDir = '../results/qc'

if not os.path.exists(outputDir):
    os.makedirs(outputDir)

#%%#############################################################################
### Demultiplex and QC samples (new samples only)
################################################################################

#%%#############################################################################
### Now all samples are in qcDir
### Construct a frequency distribution of sequence counts to determine if
### any samples are outliers.
################################################################################

sampleList = []
for sample in os.listdir(qcDir):
    if sample.endswith('.fastq'):
       sampleList.append(sample)
sampleList = [sample.replace('.fastq', '') for sample in sampleList]

# Construct a dictionary of sample sizes
sampleSizeDict = {}
for sample in sampleList:
    print('Processing sample '+str(sampleList.index(sample) + 1)+' of '+str(len(sampleList)))
    recordList = list(SeqIO.parse(qcDir+'/'+sample+'.fastq', 'fastq'))
    sampleSizeDict[sample] = len(recordList)

# Write the dictionary to file
with open(outputDir+'/numSeqPerSample.csv', 'w') as outFile:
    outFile.write('Sample,NumSeqs\n')
    for sample in sampleSizeDict.keys():
        outFile.write(sample+','+str(sampleSizeDict[sample])+'\n')

# Convert to log10 and plot a histogram with integer bin size
logSeries = np.log10(list(sampleSizeDict.values()))
plt.hist(logSeries, bins = range(math.floor(logSeries.min()), math.ceil(logSeries.max())+1), normed=True)
plt.xlabel('Log10 Sequence Count')
plt.ylabel('Frequency of Samples')
plt.savefig(outputDir+'/numSeqPerSample.png')

#%%#############################################################################
### These data suggest that two samples ought to be excluded:
### TBE13JUN07.R1 and TBE13JUN07.R2
### Create and store a list of these filenames for exclusion from further analysis
################################################################################

excludeList = ['TBE13JUN07.R1', 'TBE13JUN07.R2']

with open(outputDir+'/excludeMe.csv', 'w') as outFile:
    for sample in excludeList:
        outFile.write(sample+'\n')