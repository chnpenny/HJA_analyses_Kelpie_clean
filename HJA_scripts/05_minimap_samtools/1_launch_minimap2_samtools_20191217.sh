#!/bin/bash
set -e
set -u
set -o pipefail
#######################################################################################
#######################################################################################
# a shell script to launch bsub files.
#######################################################################################
#######################################################################################

# upload the _loop_minimap2_only_20200917.{sh,sub} and _loop_samtools_only_20200917.{sh,sub} to _Oregon/2.trimmeddata

############# edit minimap2 and samtools scripts #############
ada
interactive
# to use parallel without a pathname in bsub scripts
PATH=$PATH:~/scripts/parallel-20200922/bin/ # GNU Parallel

cd ~/_Oregon/HJAdryad/2.trimmeddata/; ls

############# copy the minimap and samtools shell and bsub scripts into each BWA folder and edit the jobIDs
MINIMAP2_SUB="_loop_minimap2_only_20200917.sub"; echo ${MINIMAP2_SUB}
MINIMAP2_SH="_loop_minimap2_only_20200917.sh"; echo ${MINIMAP2_SH}
SAMTOOLS_SUB="_loop_samtools_only_20200917.sub"; echo ${SAMTOOLS_SUB}
SAMTOOLS_SH="_loop_samtools_only_20200917.sh"; echo ${SAMTOOLS_SH}

cd ~/_Oregon/HJAdryad/2.trimmeddata/; ls # 2.trimmeddata
parallel cp ${MINIMAP2_SUB} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${MINIMAP2_SH} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${SAMTOOLS_SUB} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
parallel cp ${SAMTOOLS_SH} BWA{} ::: 01 02 03 04 05 06 07 08 09 10
ls BWA{01,02,03,04,05,06,07,08,09,10}

# edit the bsub files so that the correct jobID will show up (i suppose i could have instead run a job array...)
cd ~/_Oregon/HJAdryad/2.trimmeddata/; ls # 2.trimmeddata

parallel "sed 's/mnmp01/mnmp{}/g' BWA{}/${MINIMAP2_SUB} > BWA{}/${MINIMAP2_SUB}_tmp" ::: 01 02 03 04 05 06 07 08 09 10
parallel "mv BWA{}/${MINIMAP2_SUB}_tmp BWA{}/${MINIMAP2_SUB}" ::: 01 02 03 04 05 06 07 08 09 10
head -n7 BWA{01,02,03,04,05,06,07,08,09,10}/${MINIMAP2_SUB} # check.  should be mnmp{01,02,03,...}
tail -n2 BWA{01,02,03,04,05,06,07,08,09,10}/${MINIMAP2_SUB} # check.  should be the correct shell file

parallel "sed 's/samtl01/samtl{}/g' BWA{}/${SAMTOOLS_SUB} > BWA{}/${SAMTOOLS_SUB}_tmp" ::: 01 02 03 04 05 06 07 08 09 10
parallel "mv BWA{}/${SAMTOOLS_SUB}_tmp BWA{}/${SAMTOOLS_SUB}" ::: 01 02 03 04 05 06 07 08 09 10
head -n 7 BWA{01,02,03,04,05,06,07,08,09,10}/${SAMTOOLS_SUB} # check.  should smtl{01,02,03,...}
tail -n 1 BWA{01,02,03,04,05,06,07,08,09,10}/${SAMTOOLS_SUB} # check.  should have the correct samtools shell filename

ls # BWA* folders should now sort to bottom

####### launch the minimap2 scripts #######
# cd into each BWA folder and submit bsub job
cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA01; ls # 2.trimmeddata
echo "${MINIMAP2_SUB}" # if incorrect or NULL, go up to top of script
sbatch ${MINIMAP2_SUB}
squeue -u userid
ls # should show working folder

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA02; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u userid
ls

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA03; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA04; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA05; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA06; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA07; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA08; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA09; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA10; ls # 2.trimmeddata
sbatch ${MINIMAP2_SUB}
squeue -u userid

squeue -u userid
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA02/ # check if there is a working folder
tail ~/_Oregon/HJAdryad/2.trimmeddata/BWA02/mnmp02.out
tail ~/_Oregon/HJAdryad/2.trimmeddata/BWA02/mnmp02.err



######  WAIT FOR THE MINIMAP2 JOBS TO FINISH BEFORE LAUNCHING THE SAMTOOLS SCRIPTS
ada
interactive

####### launch samtools scripts #######
cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA01; ls # 2.trimmeddata
echo "${SAMTOOLS_SUB}" # if incorrect or NULL, go up to top of script
sbatch ${SAMTOOLS_SUB}
squeue -u userid
ls minimap2_outputs/

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA02; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA03; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA04; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA05; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA06; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA07; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA08; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA09; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u userid

cd ~/_Oregon/HJAdryad/2.trimmeddata/BWA10; ls # 2.trimmeddata
sbatch ${SAMTOOLS_SUB}
squeue -u userid

# dilution series samples
# cd ~/_Oregon/HJAdryad/2.trimmeddata/_BWA11_dilution_series_only; ls # 2.trimmeddata
# sbatch ${SAMTOOLS_SUB}
# squeue -u userid

ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA01/minimap2_outputs/ # should show new genomecov files
ls ~/_Oregon/HJAdryad/2.trimmeddata/BWA02/minimap2_outputs/ # should show new genomecov files
squeue -u userid
