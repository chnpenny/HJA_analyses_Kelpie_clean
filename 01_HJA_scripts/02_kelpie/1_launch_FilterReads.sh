#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to launch FilterReads in 10 different BWA directories
#######################################################################################
#######################################################################################

# All the scripts are run on the UEA ADA supercluster, which uses the Slurm Workload Manager
# We have 242 directories containing the shotgun sequence files
# so to speed things up, i have manually parallelised the script by
# manually distributing the 242 directories into 10 directories, named BWA01 to BWA10
# the reason for using "BWA" is that i originally was going to use bwa-mem, but i later changed to using minimap2

# this script copies _parallel_FilterReads_20200915.sub and _parallel_FilterReads_20200915.sh into the 10
# BWA directories and then launches them on the server
# pathname to these directories on ADA is ~/Oregon/HJAdryad/2.trimmeddata/.  Adjust your script accordingly

############# edit kelpie script #############
# log into ADA
# path to GNU parallel
PATH=$PATH:~/scripts/parallel-20200922/bin/ # GNU Parallel

############# copy the FilterReads shell and sub scripts into each BWA folder and edit the jobIDs
cd ~/_Oregon/HJAdryad/2.trimmeddata/; ls
# upload _parallel_FilterReads_20200915.sub and _parallel_FilterReads_20200915.sh into HJAdryad/2.trimmeddata/
FILTREAD2_SUB="_parallel_FilterReads_20200915.sub"; head -n30 ${FILTREAD2_SUB}
FILTREAD2_SH="_parallel_FilterReads_20200915.sh"; head -n60 ${FILTREAD2_SH}

parallel cp ${FILTREAD2_SUB} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${FILTREAD2_SH} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel ls -lrt BWA{} ::: 01 02 03 04 05 06 07 08 09 10
# ls BWA{01,02,03,04,05,06,07,08,09,10}

# edit the bsub files so that the correct jobID will show up (i suppose i could have instead run a job array...)
cd ~/_Oregon/HJAdryad/2.trimmeddata/; ls

parallel "sed 's/filtrd/filtrd{}/g' BWA{}/${FILTREAD2_SUB} > BWA{}/${FILTREAD2_SUB}_tmp" ::: 01 02 03 04 05 06 07 08 09 10
parallel "mv BWA{}/${FILTREAD2_SUB}_tmp BWA{}/${FILTREAD2_SUB}" ::: 01 02 03 04 05 06 07 08 09 10
head -n14 BWA{01,02,03,04,05,06,07,08,09,10}/${FILTREAD2_SUB} # check.  should be filtrd{01,02,03,...}
     # check that i'm using #SBATCH -p compute-24-96
tail -n2 BWA{01,02,03,04,05,06,07,08,09,10}/${FILTREAD2_SUB} # check.  should be the correct shell file

ls # BWA* folders should now sort to bottom

####### launch the 10 FilterReads scripts #######
# cd into each BWA folder and submit bsub job
# each job takes about an hour to run, 4 jobs are allowed to run simultaneously on ada
cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA01; ls
echo ${FILTREAD2_SUB}
sbatch ${FILTREAD2_SUB}
squeue -u userid
ls

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA02; ls
sbatch ${FILTREAD2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA03; ls
sbatch ${FILTREAD2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA04; ls
sbatch ${FILTREAD2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA05; ls
sbatch ${FILTREAD2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA06; ls
sbatch ${FILTREAD2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA07; ls
sbatch ${FILTREAD2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA08; ls
sbatch ${FILTREAD2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA09; ls
sbatch ${FILTREAD2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA10; ls
sbatch ${FILTREAD2_SUB}
squeue -u userid

# check for filterreads output files in one of the earlier files to be filtered
# look for: HOBO-357-M1-S2_BDSW190603169-1a_1_val_1_COI.fa
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA01/HOBO-357-M1-S2_BDSW190603169-1a/
squeue -u userid
ls
