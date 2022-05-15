#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# create a directory of the sequence files from which working files have been removed
# this was run on macOS
#######################################################################################
#######################################################################################

# log into ada.uea.ac.uk
interactive
cd ~/_Oregon/
# mkdir HJAdryad/2.trimmeddata/ # receiving directory
# cp sequence files to receiving directory
sbatch _2_copy_bwa_dirs.bsub

# remove non-sequence directories and files (e.g. filterreadsoutput/, *.bsub, *.err, *.out, *.sh, *.sub)
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.bsub #  | wc -l # first check the command
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.err #  | wc -l # first check the command
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.out #  | wc -l # first check the command
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.sh #  | wc -l # first check the command
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.sub #  | wc -l # first check the command
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.md5 #  | wc -l # first check the command
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.txt #  | wc -l # first check the command
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/filterreadsoutput/ #  | wc -l # first check the command

rm -Iv ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.bsub #  | wc -l # -Iv is interactive and verbose, which is safer
rm -Iv ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.err #  | wc -l # -Iv is interactive and verbose, which is safer
rm -Iv ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.out #  | wc -l # -Iv is interactive and verbose, which is safer
rm -Iv ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.sh #  | wc -l # -Iv is interactive and verbose, which is safer
rm -Iv ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/*.sub #  | wc -l # -Iv is interactive and verbose, which is safer
rm -dIv ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/filterreadsoutput/ #  | wc -l # -Iv is interactive and verbose, which is safer

# re-run the above ls commands
# confirm the number of copied directories
ls -d ~/_Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/* | wc -l # 242 directories, which is correct
