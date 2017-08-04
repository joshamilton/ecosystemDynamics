###############################################################################
# demultiplexAndQC.py
# Copyright (c) 2017, Joshua J Hamilton
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Demultiplex and QC raw FASTA files
# Compute read distribution
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
import os
import matplotlib.pyplot as plt
import pandas as pd

#%%#############################################################################
### Define folder structure
################################################################################

qcDir = '../data/demultiplexedAndQCd' # demultiplexed and QCd
outputDir = '../results/qc'
seqLenDir = outputDir+'/seqLength'
excludeList = outputDir+'/excludeMe.csv' # sampeles to exclude

# Read in all samples
sampleList = []
for sample in os.listdir(qcDir):
    if sample.endswith('.fastq'):
       sampleList.append(sample)
sampleList = [sample.replace('.fastq', '') for sample in sampleList]

# Read in samples to exclude
excludeList = []
with open(outputDir+'/excludeMe.csv', 'r') as inFile:
    excludeList = inFile.readlines()
excludeList = [sample.strip() for sample in excludeList]

# Define a new sample list
sampleList = [sample for sample in sampleList if sample not in excludeList]

#%%#############################################################################
### Construct a frequency distribution of sequence length in each file
### Because deblurring requires all sequences to be the same length, use the
### length distribution to find the shortest sequence length
################################################################################

lengthCountDF = pd.DataFrame(data=0, index=range(1, 151), columns=sampleList)

for sample in sampleList:
    print('Processing sample '+str(sampleList.index(sample) + 1)+' of '+str(len(sampleList)))
    for record in list(SeqIO.parse(qcDir+'/'+sample+'.fastq', 'fastq')):
        lengthCountDF.loc[len(record.seq), sample] = lengthCountDF.loc[len(record.seq), sample] + 1
    
lengthCountDF = lengthCountDF.loc[~(lengthCountDF==0).all(axis=1)]

lengthCountDF['Total'] = lengthCountDF.sum(axis=1)
lengthCountDF['Fraction'] = lengthCountDF['Total'] / lengthCountDF['Total'].sum()

# Plot for each sample to check for outliers
for sample in sampleList:
    lengthCountDF['Fraction'] = lengthCountDF[sample] / lengthCountDF[sample].sum()
    plt.bar(lengthCountDF.index, lengthCountDF['Fraction'])
    plt.xlabel('Sequence Length')
    plt.ylabel('Fraction of Sequences')
    plt.savefig(seqLenDir+'/'+sample+'.png')
    plt.clf()

# Plot frequency distribution for all samples
lengthCountDF['Total'] = lengthCountDF[sampleList].sum(axis=1)
lengthCountDF['Fraction'] = lengthCountDF['Total'] / lengthCountDF['Total'].sum()

plt.bar(lengthCountDF.index, lengthCountDF['Fraction'])
plt.xlabel('Sequence Length')
plt.ylabel('Fraction of Sequences')
plt.savefig(outputDir+'/seqLength.png')

lengthCountDF.to_csv(outputDir+'/seqLengthPerSample.csv')

