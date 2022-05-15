#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# This script navigates through all the bioinformatic scripts from sequences to input tables to sjSDM
#######################################################################################
#######################################################################################

# Download the sequence files from DataDryad and save on a Linux OS server
# On my server, the pathnames are ~/Oregon/HJAdryad/2.trimmeddata/BWA{01,02,03,04,05,06,07,08,09,10}/

# 1. run code in 02_kelpie/1_launch_FilterReads.sh
  # FilterReads reduces the sequence files to those that match COI sequences
  # on ada, i run filterreads on all 10 BWA folders in parallel, requiring ~3 hrs wallclock time
  # (as opposed to a sequential runtime of 10*3 = 30 hrs wallclock)

# 2. run code in 02_kelpie/2_parallel_kelpie_20200917.sh
  # Kelpie carries out in-silico PCR on the outputs of 1_launch_FilterReads.sh

# 3. run code in 02_kelpie/3_cleanup_kelpie_20200917.sh
  # concatenate all kelpie outputs, deduplicate names, and clean up

# 4. run code in 05_minimap_samtools/1_launch_minimap2_samtools_20191217.sh
  # run minimap2 to map reads against HJAdryad/reference_seqs/kelpie_20200927_BF3BR2_derep_filter3_vsearch97_rmdup_spikes.fas
  # run samtools and bedtools to count the number and percent coverage of reads mapped to each OTU representative sequence

# 5. run code in 2_launch_postsamtools_copy_idx_genomecov_files.sh
  # copy idx and genomcov files into a directory to be downloaded to my computer

# 6. run code in 06_idxstats_tabulate/1_idxstats_tabulate.Rmd
  # This file takes the idxstats.txt outputs of samtools and bedtools (idxstats, genomecov) and combines them with the sample metadata and environmental covariates, fixes the metadata, does sanity checks, removes failed samples, and creates sampleXspecies tables

# 7. run code in 06_idxstats_tabulate/2_FSL_and_lysis_correction.Rmd
  # Two additional corrections: one for the proportion of the total lysis buffer (lysis_ratio) that was used from each sample and one for the weighted mean COI DNA spike-in reads (COISpike_wt_mean).

# 8. run code in 06_idxstats_tabulate/3_sjsdm_dataprep.Rmd
  # Prepares datasets for sjSDM analysis



########## remove code below before publishing

# upload HJAdryad/2.trimmeddata/ to S3
ada
interactive
cd _Oregon/
ls

# tar HJAdryad/
sbatch 3_tar_HJAdryad.bsub

# upload HJAdryad.tar to s3://amazon-dryad-douglasyu
module add s3cmd/2.2.0 # makes s3cmd available
s3cmd ls s3://amazon-dryad-douglasyu
# sbatch 4_s3cmd_upload_sequence_data.bsub
squeue -u b042
ls
# manually check md5s of HJAdryad.tar on S3 and ada

# delete folder and file
ls ~/_Oregon/HJAdryad/
# rm -rf ~/_Oregon/HJAdryad/
# rm  ~/_Oregon/HJAdryad.tar

# set file public on Amazon S3
s3cmd setacl --acl-public s3://amazon-dryad-douglasyu/HJAdryad.tar.partaa
s3cmd setacl --acl-public s3://amazon-dryad-douglasyu/HJAdryad.tar.partab
s3cmd setacl --acl-public s3://amazon-dryad-douglasyu/HJAdryad.tar.partac
s3cmd setacl --acl-public s3://amazon-dryad-douglasyu/HJAdryad.tar.partad
# should get this message back from Amazon:
	# "s3://amazon-dryad-douglasyu/HJAdryad.tar: ACL set to Public  [1 of 1]"

# the public URL then follows this format:
	# https://<bucket-name>.s3.amazonaws.com/<object or key name>
https://amazon-dryad-douglasyu.s3.amazonaws.com/HJAdryad.tar
# or from Transmit
https://amazon-dryad-douglasyu.s3-eu-west-2.amazonaws.com/HJAdryad.tar.partaa
https://amazon-dryad-douglasyu.s3-eu-west-2.amazonaws.com/HJAdryad.tar.partab
https://amazon-dryad-douglasyu.s3-eu-west-2.amazonaws.com/HJAdryad.tar.partac
https://amazon-dryad-douglasyu.s3-eu-west-2.amazonaws.com/HJAdryad.tar.partad

# 1 TB tar file too big for datadryad
# download tar file and split
sbatch 5_download_HJAdryad.bsub
squeue -u b042
ls
# compare MD5s
module add s3cmd/2.2.0 # makes s3cmd available
s3cmd ls --list-md5 s3://amazon-dryad-douglasyu > amazon-dryad-douglasyu_20220206.md5
sbatch 6_md5_HJAdryad.bsub

# split tar file
split -b 290G HJAdryad.tar "HJAdryad.tar.part"

cat HJAdryad.tar.partaa HJAdryad.tar.partab HJAdryad.tar.partac HJAdryad.tar.partad > HJAdryad.tar



# upload HJAdryad.tar to s3://amazon-dryad-douglasyu
module add s3cmd/2.2.0 # makes s3cmd available
s3cmd ls s3://amazon-dryad-douglasyu
sbatch 8_s3cmd_upload_sequence_data.bsub
squeue -u b042
ls
# manually check md5s of HJAdryad.tar on S3 and ada

# delete folder and file
ls ~/_Oregon/HJAdryad/
# rm ~/_Oregon/HJAdryad.tar.parta{a,b,c,d}
# rm  ~/_Oregon/HJAdryad.tar

# set file public on Amazon S3
s3cmd setacl --acl-public s3://amazon-dryad-douglasyu/HJAdryad.tar
# should get this message back from Amazon:
	# "s3://amazon-dryad-douglasyu/HJAdryad.tar: ACL set to Public  [1 of 1]"

# the public URL then follows this format:
	# https://<bucket-name>.s3.amazonaws.com/<object or key name>
https://amazon-dryad-douglasyu.s3.amazonaws.com/HJAdryad.tar
# or from Transmit
https://amazon-dryad-douglasyu.s3-eu-west-2.amazonaws.com/HJAdryad.tar


# initiate download by datadryad
# upload to datadryad the binaries for filterreads and kelpie
# README
# download all parts, cat to rejoin, and extract
cat HJAdryad.tar.part* > HJAdryad.tar
tar -xvf HJAdryad.tar
