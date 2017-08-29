###############################################################################
# deblur.py
# Copyright (c) 2017, Joshua J Hamilton
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# After fastq files have been demultiplexed and QC'd, the deblur algorithm
# is used to identify and count unique sequences in each sample.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
import biom
import os
import pandas as pd
import subprocess

#%%#############################################################################
### Define folder structure
################################################################################

qcDir = '../data/demultiplexedAndQCd' # demultiplexed and QCd fasta files
outputDir = '../results/deblur'

minLength = 150 # trim sequences to this length

if not os.path.exists(outputDir):
    os.makedirs(outputDir)

#%%#############################################################################
### Call deblur
################################################################################

# Call deblur with the following options
#  deblur workflow --seqs-fp ../data/demultiplexedAndQCd
#  --output-dir ../data/deblur/
#  --trim-length 114
#  --overwrite
#  --min-reads 1
#  --min-size 1
#  --threads-per-sample 1

deblurCommand = 'deblur workflow --seqs-fp '+qcDir+' --output-dir '+outputDir+' --trim-length '+str(minLength)+' --overwrite --min-reads 1 --min-size 1 --threads-per-sample 1'
subprocess.call(deblurCommand, shell=True)

#%%#############################################################################
### Convert BIOM file to a human-readable OTU table
################################################################################

# Read in the BIOM table
otuTable = biom.load_table(outputDir+'/all.biom')

# Convert to a TSV
otuTableString = otuTable.to_tsv() # returns a string
with open(outputDir+'/otuTable.tsv', 'w') as outFile:
    outFile.write(otuTableString)

# Read back into a pandas dataframe and reindex
otuDF = pd.read_csv(outputDir+'/otuTable.tsv', sep='\t', skiprows=1)
index = list(range(1, len(otuDF)+1)) # Create an index
index = [str(seqInt).zfill(6) for seqInt in index] # Convert to strings of equal length
otuDF['IndexCol'] = index
otuDF.index = otuDF['IndexCol'] 
otuDF = otuDF.drop(labels=['#OTU ID', 'IndexCol'], axis =1)
otuDF.to_csv(outputDir+'/otuTable.csv')
os.remove(outputDir+'/otuTable.tsv')

# Also update the FASTA file, replacing sequence strings with the appropriate numeric ID
seqInt = 1
with open(outputDir+'/all.seqs.numIDs.fasta', 'w') as outFile:
    for seqRecord in SeqIO.parse(outputDir+'/all.seqs.fa', 'fasta'):
        seqRecord.id = str(seqInt).zfill(6)
        outFile.write('>'+seqRecord.id+'\n')
        outFile.write(str(seqRecord.seq)+'\n')
        seqInt += 1
