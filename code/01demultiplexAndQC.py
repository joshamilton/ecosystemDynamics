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

import os
import subprocess

#%%#############################################################################
### Define folder structure
################################################################################
    
fastqDir = '../data/rawSequences' # raw fastq files
qcDir = '../data/demultiplexedAndQCd' # demultiplexed and QCd

if not os.path.exists(qcDir):
    os.makedirs(qcDir)

#%%#############################################################################
### Demultiplex and QC samples (new samples only)
################################################################################

#%%#############################################################################
### Move samples to appropriate folder (new samples only)
################################################################################

#%%#############################################################################
### Now all samples are in qcDir, with one folder per year
### Use awk to construct a frequency distribution of sequence length
################################################################################

awkCommand  = "awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' "+qcDir+"/*/*.fastq >"+qcDir+"/seqDist.txt"
subprocess.call(awkCommand, shell=True)
    

