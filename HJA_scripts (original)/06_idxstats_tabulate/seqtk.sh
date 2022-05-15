#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to calculate total number of reads per sample
#######################################################################################
#######################################################################################

# Usage: run the commmands interactively

# run this on hpc (i ran these four commands in parallel in multiple terminal windows)
# opt-cmd-return to send to Terminal console
ssh hpc
interactive
PATH=$PATH:~/scripts/seqkitdir
seqkit # chk that command is working

cd ~/_Oregon/HJAdryad/2.trimmeddata/ # go to root folder and run these commands sequentially
# run in separate windows because each takes around 2 hours to run

# HJAdryad/2.trimmeddata
nohup seqkit stats BWA*/*/*fq.gz > fastq_read_counts.txt &
# download to local folder
