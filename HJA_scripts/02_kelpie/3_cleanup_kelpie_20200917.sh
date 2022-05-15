#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to interactively clean up after running 2_parallel_kelpie_20200917.sh
##################################################################################################
##################################################################################################

# on ada
ssh ada
interactive # 24 cores

# set paths and add modules
PATH=$PATH:~/scripts/vsearch-2.15.0-linux-x86_64/bin/ # downloaded 12 Jul 2020 from github
PATH=$PATH:~/scripts/parallel-20200922/bin/ # GNU Parallel
PATH=$PATH:~/scripts/seqkitdir/ # downloaded 12 Dec 2019 from github
# clustering methods not used
    # PATH=$PATH:~/scripts/swarm-3.0.0-linux-x86_64/bin # downloaded 16 Dec 2019 from github
    # module add Sumaclust/1.0.34 # not avail on ada

# go to correct folder
HOMEFOLDER="/gpfs/home/userid/_Oregon/HJAdryad/" # should be ~/_Oregon/HJAdryad/
TARGETFOLDER="2.trimmeddata/" # testkelpie
echo "Home folder is ${HOMEFOLDER}${TARGETFOLDER}" # ~/_Oregon/HJAdryad/2.trimmeddata
cd ${HOMEFOLDER}${TARGETFOLDER} # || exit

#### After running 2_parallel_kelpie_20200917.sh, there will be fasta files in two output folders
# kelpieoutputneighbors/
# kelpieoutputindiv/
# Now, i concatenate, deduplicate names, and clean up

# concatenate the fasta files in both output folders into a single large fasta file
primer="BF3BR2" # BF3BR2, LERAY
minlen=400 # 400 for BF3BR2, 300 for LERAY
timestamp="20200927_${primer}" # e.g. 20200715_LERAYFOL, 20200916_BF3BR2
echo $timestamp

if [ ! -d kelpie_output_${primer} ] # if directory kelpie_output_${primer} does not exist.
then
     mkdir kelpie_output_${primer}
fi
ls

# concatenate the fasta files from the kelpieoutputneighbors/ directory
cat kelpieoutputneighbors/*${primer}.fas > kelpie_output_${primer}/kelpie_${timestamp}.fas
ls kelpie_output_${primer}/
seqkit stats kelpie_output_${primer}/kelpie_${timestamp}.fas # 626,764 BF3BR2 seqs 2.0.10, 347,349 LERAY seqs 2.0.10

# concatenate the fasta files from the kelpieoutputindiv/ directory
cat kelpieoutputindiv/*${primer}.fas >> kelpie_output_${primer}/kelpie_${timestamp}.fas
ls kelpie_output_${primer}/

cd kelpie_output_${primer}/
seqkit stats kelpie_${timestamp}.fas # 1,233,694 BF3BR2 2.0.10, 708,938 LERAY 2.0.10 # 1,244,402 BF3BR2 seqs 2.0.8, 714,508 LERAY seqs 2.0.8
seqkit seq -m "${minlen}" kelpie_${timestamp}.fas -o kelpie_${timestamp}_min${minlen}.fas
seqkit stats kelpie_${timestamp}_min${minlen}.fas kelpie_${timestamp}.fas # 1,233,585 seqs BF3BR2 2.0.10, 708,786 LERAY 2.0.10
mv kelpie_${timestamp}_min${minlen}.fas kelpie_${timestamp}.fas
seqkit stats kelpie_${timestamp}.fas

# a problem is that the amplicon names are reused across multiple samples (e.g. >R1 is used in each sample fasta)
# this causes a problem after clustering if the names collide (although usually avoided because sizes are different)
# so i use seqkit rename to rename duplicated names
# deduplicate header names
grep "R1$" kelpie_${timestamp}.fas # there are duplicate headers, because i concatenated multiple fasta files
grep "R1$" kelpie_${timestamp}.fas | wc -l # should equal to number of samples + number of nearest neighbor groups: 335
seqkit rename kelpie_${timestamp}.fas > kelpie_${timestamp}_rename.fas # renames duplicate headers
mv kelpie_${timestamp}_rename.fas kelpie_${timestamp}.fas # replace old with new, deduplicated version
grep "R1$" kelpie_${timestamp}.fas # duplicates have new names, plus the original name info
seqkit stats kelpie_${timestamp}.fas # 1,233,585 seqs BF3BR2 2.0.10

# dereplicate the fasta file
    # orig command
    # vsearch --derep_fulllength kelpie_${timestamp}.fas --sizeout --sortbysize --threads 0 --output kelpie_${timestamp}_derep.fas
      # --threads 0 # zero to use all avail cores
      # --relabel_sha1 # hash each sequence as the name (to compare sequences)
      # --fastq_width 0 # no line breaks in sequence
      # --sizeout # include size= information
# command formatted for input to swarm 3.0.0
vsearch --derep_fulllength kelpie_${timestamp}.fas --sizeout --fasta_width 0 --threads 0 --output kelpie_${timestamp}_derep.fas # --relabel_sha1
head kelpie_${timestamp}_derep.fas
tail kelpie_${timestamp}_derep.fas
seqkit stats kelpie_${timestamp}_derep.fas # 5,560 uniq seqs BF3BR2 2.0.10, 4,152 uniq seqs LERAY 2.0.10 # 5,351 uniq BF3BR2 seqs 2.0.8; 4,227 uniq Leray seqs 2.0.8

# uncomment when ready to run
ls
# rm -f kelpie_${timestamp}.fas
ls

# kelpie_${timestamp}_derep.fas is what i use as input to GBIF, after which i cluster into 97% OTUs and do a translation align to remove false OTUs

# remove intermediate subdirectories
HOMEFOLDER="/gpfs/home/userid/_Oregon/HJAdryad/" # should be ~/_Oregon/HJAdryad/
TARGETFOLDER="2.trimmeddata/" # testkelpie
echo "Home folder is ${HOMEFOLDER}${TARGETFOLDER}" # ~/_Oregon/HJAdryad/2.trimmeddata
cd ${HOMEFOLDER}${TARGETFOLDER} # || exit

rm -rf kelpieoutputneighbors/
rm -rf kelpieoutputindiv/
rm -rf allfilterreadsoutput/


# concatenate the discard fasta files into a single large fasta file
echo $timestamp # e.g. 20200715_LERAYFOL, 20200916_BF3BR2, 20200916_LERAY
cd ~/_Oregon/HJAdryad/2.trimmeddata/; ls
cat kelpieoutputindiv/*_discards.fas > kelpie_output_${primer}/kelpie_${timestamp}_discards.fas
cat kelpieoutputneighbors/*_discards.fas >> kelpie_output_${primer}/kelpie_${timestamp}_discards.fas
cd kelpie_output_${primer}/
seqkit stats kelpie_${timestamp}_discards.fas # 16,450 seqs BF3BR2 2.0.10, 9,554 seqs LERAY 2.0.10
head kelpie_${timestamp}_discards.fas
grep "D1" kelpie_${timestamp}_discards.fas # duplicates have new names, plus the original name info
seqkit rename kelpie_${timestamp}_discards.fas > kelpie_${timestamp}_discards_rename.fas # renames duplicate headers
grep "D1" kelpie_${timestamp}_discards_rename.fas # duplicates have new names, plus the original name info
mv kelpie_${timestamp}_discards_rename.fas kelpie_${timestamp}_discards.fas
seqkit stats kelpie_${timestamp}_discards.fas # 16,450 seqs BF3BR2 2.0.10 # 12,883 BF3BR2 seqs kelpie 2.0.8; 10,962 LERAY seqs 2.0.8

# END




#### vsearch clustering code
# vsearch --cluster_fast kelpie_${timestamp}_derep.fas --sizein --sizeout --id 0.97 --centroids kelpie_${timestamp}_vsearch97.fas # --uc kelpie_${timestamp}_clusters.uc
# vsearch --sortbysize kelpie_${timestamp}_vsearch97.fas --output kelpie_${timestamp}_vsearch97_allseqs.fas
# seqkit stats kelpie_${timestamp}_vsearch97_allseqs.fas # 1,155 OTUs 2.0.4, 1,285 OTUs 2.0.6
# swarm results in multiple OTUs that get assigned to the same species in BOLD, which means that they will not receive mappings due to competing, similar sequences
# sumaclust does not have an option to output only the representative sequences
# five species only appeared as OTUs with swarm and pure dereplication (e.g. Helina troene, Hybomitra affinis, Hybomitra rhombica, Laphria columbica, Phaonia errans). I will see if they reappear when i have a larger dataset


#### swarm clustering code
# # swarm clustering
# swarm -v # Swarm 3.0.0
# # -d 1 # default and recommended
# # -f # fastidious
# # -z # use usearch size= for abundance
# # -s filename # statistics file
# # -w filename # OTU representative fasta file
# # > /dev/null # discard std output
# swarm -t 16 -d 1 -f -z -s OTU_stats_${timestamp}.txt -w kelpie_${timestamp}_swarmOTUs.fas kelpie_${timestamp}_derep.fas > /dev/null
# #
# seqkit stats kelpie_${timestamp}_swarmOTUs.fas # 439 OTUs
# # for download to use on bold
# # split representative OTU set into subsets of 100 seqs for upload to BOLDSystems (provides the BOLD consensus)
# seqkit split -s 100 kelpie_${timestamp}_swarmOTUs.fas


# split representative OTU set into subsets of 100 seqs for upload to BOLDSystems (provides the BOLD consensus)
# not needed, use https://www.gbif.org/tools/sequence-id
# seqkit split -s 100 kelpie_${timestamp}_centroids_sort.fas

# sumaclust clustering
# sumaclust -t 0.97 -e kelpie_${timestamp}.fas > kelpie_${timestamp}_sumaclust97.fas
# tail kelpie_${timestamp}_sumaclust97.fas
