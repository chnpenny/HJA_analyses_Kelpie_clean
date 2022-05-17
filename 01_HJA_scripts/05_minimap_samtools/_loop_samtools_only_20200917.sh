#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to loop through a set of bam files and run samtools, bedtools
#######################################################################################
#######################################################################################

# Usage: bash _loop_samtools_only_20200917.sh

# SCRIPT OUTLINE:
	# make list of folders
  # loop through folder name to get $FOLDER
  # get name for $SAMPLEname
  # move fastq files to a working folder
  # run the samtools and bedtools commands
  # move the result to minimap2_outputs
	# rm working folder
  # next folder

PIPESTART=$(date)

# upload _loop_samtools_only_20200917.sh to BWA** folder that contains the bam files that i want to filter
# I don't place the shell script in the minimap2_outputs because it's harder to manage the files
HOMEFOLDER=$(pwd) # this sets working directory to that folder.

echo "Home folder is ${HOMEFOLDER}"
cd "${HOMEFOLDER}"/minimap2_outputs/ || exit  # cd into minimap2_outputs folder

# read in folder list and make a bash array
# find * -maxdepth 0 -type d -name "Sample_*" > folderlist.txt  # find all folders (-type d) starting with "Sample*"
sample_info=folderlist.txt # put folderlist.txt into variable. folderlist.txt should already be in the minimap2_outputs folder from minimap2 script
sample_names=($(cut -f 1 "$sample_info" | uniq)) # convert variable to array this way
# echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "folders that will be processed." # echo number of elements in the array

INDEX=1 # set variable

for sample in "${sample_names[@]}"  # ${sample_names[@]} is the full bash array
do
     echo "Start of Sample" ${INDEX} of ${#sample_names[@]}

     FOLDER="${sample}" # this sets the folder name to the value in the bash array (which is a list of the folders)

     echo "**** Working folder is" $FOLDER

     #### samtools filtering pipeline using -F 2308 -q 1 ####
          # samtools view -F 2308 # exclude UNMAP,SECONDARY,SUPPLEMENTARY # exclude unmapped, secondary, and supplementary alignments
          # see https://broadinstitute.github.io/picard/explain-flags.html;  see also https://www.biostars.org/p/256448/
          # samtools view -q 1 # to remove multiply mapped reads = q 0
     # samtools view -b -F 2308 -q 1 ${sample}_sorted.bam -o ${sample}_F2308_q1_sorted.bam
     # samtools index ${sample}_F2308_q1_sorted.bam # creates ${sample}_F2308_q1_sorted.bam.bai file
     # samtools idxstats ${sample}_F2308_q1_sorted.bam > ${sample}_F2308_q1_sorted.bam_idxstats.txt # calculate idxstats
     # bedtools genomecov -d -ibam ${sample}_F2308_q1_sorted.bam | gzip > ${sample}_F2308_q1_sorted_genomecov_d.txt.gz # calc genome cov stats

     #### samtools filtering pipeline using -F 2308 -q 48 ####
          # samtools view -F 2308 # exclude UNMAP,SECONDARY,SUPPLEMENTARY # exclude unmapped, secondary, and supplementary alignments
          # see https://broadinstitute.github.io/picard/explain-flags.html;  see also https://www.biostars.org/p/256448/
          # samtools view -q 48 # based on viewing the distribution of q scores in a couple of samples
     samtools view -b -F 2308 -q 48 ${sample}_F0x4_sorted.bam -o ${sample}_F2308_q48_sorted.bam
     samtools index ${sample}_F2308_q48_sorted.bam # creates ${sample}_F2308_q48_sorted.bam.bai file
     samtools idxstats ${sample}_F2308_q48_sorted.bam > ${sample}_F2308_q48_sorted.bam_idxstats.txt # calculate idxstats
     bedtools genomecov -d -ibam ${sample}_F2308_q48_sorted.bam | gzip > ${sample}_F2308_q48_sorted_genomecov_d.txt.gz # calc genome cov stats

     #### samtools filtering pipeline using -F 2308 -f 0x2 -q 48 ####
          # samtools view -F 2308 # exclude UNMAP,SECONDARY,SUPPLEMENTARY # exclude unmapped, secondary, and supplementary alignments
          # samtools view -f 0x2 # include PROPER_PAIR
          # see https://broadinstitute.github.io/picard/explain-flags.html;  see also https://www.biostars.org/p/256448/
          # samtools view -q 48 # based on viewing the distribution of q scores in a couple of samples
     samtools view -b -F 2308 -f 0x2 -q 48 ${sample}_F0x4_sorted.bam -o ${sample}_F2308_f0x2_q48_sorted.bam
     samtools index ${sample}_F2308_f0x2_q48_sorted.bam # creates ${sample}_F2308_q48_sorted.bam.bai file
     samtools idxstats ${sample}_F2308_f0x2_q48_sorted.bam > ${sample}_F2308_f0x2_q48_sorted.bam_idxstats.txt # calculate idxstats
     bedtools genomecov -d -ibam ${sample}_F2308_f0x2_q48_sorted.bam | gzip > ${sample}_F2308_f0x2_q48_sorted_genomecov_d.txt.gz # calc genome cov stats

     # #### samtools filtering pipeline using -F 3332 -f 0x2 -q 1 ####
     #      # samtools view -f 0x2 # keep PROPER_PAIR   .. each segment properly aligned according to the aligner
     #      # samtools view -F 3332 # exclude DUP,UNMAP,SECONDARY,SUPPLEMENTARY # exclude PCR duplicates, unmapped, secondary, and supplementary alignments
     #      # see https://broadinstitute.github.io/picard/explain-flags.html;  see also https://www.biostars.org/p/256448/
     #      # samtools view -q 1 # to remove multiply mapped reads = q 0
     # samtools view -b -F 3332 -f 0x2 -q 1 ${sample}_sorted.bam -o ${sample}_F3332_f0x2_q1_sorted.bam
     # samtools index ${sample}_F3332_f0x2_q1_sorted.bam # creates ${sample}_F3332_f0x2_q1_sorted.bam.bai file
     # samtools idxstats ${sample}_F3332_f0x2_q1_sorted.bam > ${sample}_F3332_f0x2_q1_sorted.bam_idxstats.txt # calculate idxstats
     # bedtools genomecov -d -ibam ${sample}_F3332_f0x2_q1_sorted.bam | gzip > ${sample}_F3332_f0x2_q1_sorted_genomecov_d.txt.gz # calc genome cov stats
     #
     # #### samtools filtering pipeline using -F 3332 -f 0x2 -q 60, more stringent ####
     #      # samtools view -f 0x2 # keep PROPER_PAIR   .. each segment properly aligned according to the aligner
     #      # samtools view -F 3332 # exclude DUP,UNMAP,SECONDARY,SUPPLEMENTARY # exclude PCR duplicates, unmapped, secondary, and supplementary alignments
     #      # see https://broadinstitute.github.io/picard/explain-flags.html;  see also https://www.biostars.org/p/256448/
     #      # samtools view -q 60 # to remove multiply mapped reads = q 0, and keep only high quality mappings
     # samtools view -b -F 3332 -f 0x2 -q 60 ${sample}_sorted.bam -o ${sample}_F3332_f0x2_q60_sorted.bam
     # samtools index ${sample}_F3332_f0x2_q60_sorted.bam # creates ${sample}_F3332_f0x2_q60_sorted.bam.bai file
     # samtools idxstats ${sample}_F3332_f0x2_q60_sorted.bam > ${sample}_F3332_f0x2_q60_sorted.bam_idxstats.txt # calculate idxstats
     # bedtools genomecov -d -ibam ${sample}_F3332_f0x2_q60_sorted.bam | gzip > ${sample}_F3332_f0x2_q60_sorted_genomecov_d.txt.gz # calc genome cov stats

     echo "End of Sample" ${INDEX} of ${#sample_names[@]}
     INDEX=$((INDEX+1))

done

echo "Pipeline started at $PIPESTART"
NOW=$(date)
echo "Pipeline ended at $NOW"
