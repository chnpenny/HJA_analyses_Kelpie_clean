#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to copy idx and genomcov files into a directory to be downloaded to my computer
# run after running the samtools script, does not have to be uploaded to hpc
#######################################################################################
#######################################################################################

ada
interactive
# to use parallel without a pathname in bsub scripts
PATH=$PATH:~/scripts/parallel-20200922/bin/ # GNU Parallel


####### set up output folder to hold everything before running
OUTPUTFOLDER="minimap2_outputs"
cd ~/_Oregon/HJAdryad/2.trimmeddata/; ls # production
if [ ! -d ${OUTPUTFOLDER} ] # if directory minimap2_outputs does not exist.
then
     mkdir ${OUTPUTFOLDER}
fi
ls
"ls" -l BWA*/samtl*${PLATE}*.out # check if all BWA folders have a current samtl*.out file, showing that samtools ran correctly.  Check the sizes of the samtools.out file.  they should all be about the same.

####### code to copy samtools output files to: outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}.tar.gz #######
####### This folder is then downloaded to my laptop to process with R:  idxstats_tabulate_macOS_Plates*.Rmd

# set filters and minimum mapping quality scores
FILTER1="F2308_f0x2" # filter 1
FILTER2="F2308" # filter 2
echo $FILTER1
echo $FILTER2
QUAL1=""; echo $QUAL1 # set to "" if i don't want to use this variable for, say, q1
QUAL2=48; echo $QUAL2

echo $OUTPUTFOLDER

# copy output files into minimap2_outputs_Plates${PLATE}/
# e.g. Sample_IPO3916_C5_F2308_f0x2_q60_sorted.bam_idxstats.txt
# check before copying
parallel ls -l BWA*/${OUTPUTFOLDER}/*_{1}_q{2}_sorted.bam_idxstats.txt ::: ${FILTER1} ${FILTER2} ::: ${QUAL1} ${QUAL2}
parallel ls BWA*/${OUTPUTFOLDER}/*_{1}_q{2}_sorted.bam_idxstats.txt ::: ${FILTER1} ${FILTER2} ::: ${QUAL1} ${QUAL2} | wc -l # 484 (2 X 242)
parallel cp BWA*/${OUTPUTFOLDER}/*_{1}_q{2}_sorted.bam_idxstats.txt ${OUTPUTFOLDER}/ ::: ${FILTER1} ${FILTER2} ::: ${QUAL1} ${QUAL2}
# check after copying
ls ${OUTPUTFOLDER}/
ls ${OUTPUTFOLDER}/*_q*_sorted.bam_idxstats.txt | wc -l # 484
# copy flagstat files
cp BWA*/${OUTPUTFOLDER}/*_sorted.bam.flagstat.txt ${OUTPUTFOLDER}/
ls ${OUTPUTFOLDER}/
ls ${OUTPUTFOLDER}/*_sorted.bam.flagstat.txt | wc -l # half the idxstats values

# e.g. Sample_IPO3916_C5_F2308_f0x2_q1_sorted_genomecov_d.txt.gz
# check
parallel ls -l BWA*/${OUTPUTFOLDER}/*_{1}_q{2}_sorted_genomecov_d.txt.gz ::: ${FILTER1} ${FILTER2} ::: ${QUAL1} ${QUAL2}
parallel ls -l BWA*/${OUTPUTFOLDER}/*_{1}_q{2}_sorted_genomecov_d.txt.gz ::: ${FILTER1} ${FILTER2} ::: ${QUAL1} ${QUAL2} | wc -l # 484
# copy
parallel cp BWA*/${OUTPUTFOLDER}/*_{1}_q{2}_sorted_genomecov_d.txt.gz ${OUTPUTFOLDER}/ ::: ${FILTER1} ${FILTER2} ::: ${QUAL1} ${QUAL2}  # takes longer than above
# check again
parallel ls ${OUTPUTFOLDER}/*_{1}_q{2}_sorted_genomecov_d.txt.gz ::: ${FILTER1} ${FILTER2} ::: ${QUAL1} ${QUAL2}
parallel ls ${OUTPUTFOLDER}/*_{1}_q{2}_sorted_genomecov_d.txt.gz ::: ${FILTER1} ${FILTER2} ::: ${QUAL1} ${QUAL2} | wc -l # 484

# rename, tar, and gzip for download
MAPDATE="20200929"
# TARGET="kelpie_20200916_LERAY_derep_filter3_vsearch97_rmdup_spikes"
TARGET="kelpie_20200927_BF3BR2_derep_filter3_vsearch97_rmdup_spikes" #
du -sh ${OUTPUTFOLDER}/ # 787M
# set filename to something that i can understand after download
mv ${OUTPUTFOLDER} outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}
ls # outputs_F2308_f0x2_q48_minimap2_outputs_20191219_kelpie_20200916_BF3BR2_derep_filter3_vsearch97_rmdup_spikes
# tar gzip for download, takes about 25 mins
nohup tar -czvf outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}.tar.gz outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}/
# filename format:  outputs_F2308_f0x2_q48_minimap2_outputs_20191219_kelpie_20200916_BF3BR2_derep_filter3_vsearch97_rmdup_spikes.tar.gz
ls

# uncomment and run when ready
# rm -rf outputs_${FILTER1}_q${QUAL2}_${OUTPUTFOLDER}_${MAPDATE}_${TARGET}/
ls

# download *.gz file to local computer where the file can be processed by 1.1_idxstats_tabulate.Rmd

######## remove the minimap2_output/ folders after i've finished the mapping jobs. these are the bam, bam.bai, idxstats, flagstats, and genomecov files left behind in the BWA folders
PATH=$PATH:~/scripts/parallel-20200922/bin/ # GNU Parallel
# FILTER="F2308"
# OUTPUTFOLDER="minimap2_outputs"

cd ~/_Oregon/HJAdryad/2.trimmeddata/; ls # production
ls
parallel "ls BWA{}/minimap2_outputs/" ::: 01 02 03 04 05 06 07 08 09 10
# parallel "rm -rf BWA{}/minimap2_outputs/" ::: 01 02 03 04 05 06 07 08 09 10
parallel "ls BWA{}/minimap2_outputs/" ::: 01 02 03 04 05 06 07 08 09 10
# ls: should say:  cannot access BWA*/minimap2_outputs/: No such file or directory
du -sh ~/_Oregon # should be around 1.5 T
