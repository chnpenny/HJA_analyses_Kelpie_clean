#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to loop through a set of fastq files, run FilterReads, and move output files to a single folder
##################################################################################################
##################################################################################################

# Usage: bash _parallel_FilterReads_YYYYMMDD.sh

# use the find/array/parallel method to FilterReads all fastq files and save outputs in a single directory

# upload _parallel_FilterReads_YYYYMMDD.sub and _parallel_FilterReads_YYYYMMDD.sh *into* ~/_Oregon/2019Sep_shotgun/2.trimmeddata/
# use launch_FilterReads.sh to copy into each BWA folder

# use the find/array/parallel method to filter all the fastq files and save outputs in a single directory
     # read in folder list and make a bash array
     # find ./ -maxdepth 3 -mindepth 1 -type d -iname "*.fq.gz" -exec basename {} \; | sort > fastqlist.txt
find ./ -maxdepth 3 -mindepth 1 -type d -exec basename {} \; | sort > folderlist.txt
sed -i '/minimap2_outputs/d' ./folderlist.txt # remove minimap2 folder from folderlist.txt
sed -i '/filterreadsoutput/d' ./folderlist.txt # remove minimap2 folder from folderlist.txt
cat folderlist.txt

# make array of folder names
sample_info=folderlist.txt # put folderlist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array

# FilterReads
parallel -k -j 4 "gzip -d < {1}/{1}_1_val_1.fq.gz > {1}/{1}_1_val_1.fq; gzip -d < {1}/{1}_2_val_2.fq.gz > {1}/{1}_2_val_2.fq; FilterReads -r COI -qt 30 -fasta -t 6 +f /gpfs/home/b042/GenbankCOI_24-9-2019/GenBank_24-9-19_COI_C99_20.mer 25pct {1}/{1}_?_val_?.fq; rm {1}/{1}_1_val_1.fq; rm {1}/{1}_2_val_2.fq" ::: "${sample_names[@]}"   # output files are: {1}/{1}_?_val_?_COI.fa

# mv all *_COI.fa files to a single folder where i can run kelpie on them
# this runs within each BWA folder
if [ ! -d filterreadsoutput ] # if directory filterreadsoutput does not exist.
then
	mkdir filterreadsoutput
fi
mv */*_COI.fa filterreadsoutput/
