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

import os
import subprocess

#%%#############################################################################
### Define folder structure
################################################################################

qcDir = '../data/demultiplexedAndQCd' # demultiplexed and QCd fasta files
outputDir = '../data/deblur'

minLength = 114 # trim sequences to this length

if not os.path.exists(outputDir):
    os.makedirs(outputDir)
    
#%%#############################################################################
### Call deblur
################################################################################

# Call deblur with the following options
#  deblur workflow --seqs-fp ../data/demultiplexedAndQCd/2008
#  --output-dir ../data/deblur/
#  --trim-length 114
#  --overwrite
#  --min-reads 1
#  --min-size 1
#  --threads-per-sample 1

deblurCommand = 'deblur workflow --seqs-fp '+qcDir+'/2008 --output-dir '+outputDir+' --trim-length '+str(minLength)+' --overwrite --min-reads 1 --min-size 1 --threads-per-sample 1'
subprocess.call(deblurCommand, shell=True)