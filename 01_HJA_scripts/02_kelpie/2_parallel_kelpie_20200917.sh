#!/bin/bash
set -e
set -u
set -o pipefail
##################################################################################################
##################################################################################################
# a shell script to loop through a set of fasta files, run kelpie, and move output files to a single folder
##################################################################################################
##################################################################################################

# I usually run kelpie in an interactive session because it is fast when run on filtered datasets, but it could be run
# as a batch job without needing editing, using sbatch
ada
interactive
PATH=$PATH:~/scripts/vsearch-2.15.0-linux-x86_64/bin/ # downloaded 12 Jul 2020 from github
PATH=$PATH:~/scripts/parallel-20200922/bin/ # GNU Parallel
PATH=$PATH:~/scripts/WorkingDogs/Kelpie_v2/ubuntu-16.04/ # v 2.0.10

cd ~/_Oregon/HJAdryad/2.trimmeddata/ || exit
ls

#### copy all outputs from FilterReads run into allfilterreadsoutput/
# create folder
if [ ! -d allfilterreadsoutput ] # if directory allfilterreadsoutput/ does not exist.
then
     mkdir allfilterreadsoutput
fi
ls
# mv files
mv BWA*/filterreadsoutput/*_COI.fa ./allfilterreadsoutput
ls allfilterreadsoutput
find allfilterreadsoutput -type f -iname "*_COI.fa" | wc -l # should be 484 COI.fa files


#### run kelpie on indiv fasta files and save to kelpieoutputindiv/
# create folder to hold
if [ ! -d kelpieoutputindiv ] # if directory kelpieoutputindiv does not exist.
then
     mkdir kelpieoutputindiv
fi
ls
# assumption: all FilterReads outputs will be in one folder filterreadsoutput, one up from the BWA folders
find ./allfilterreadsoutput -type f -iname "*_COI.fa" -exec basename {} \; > fastalist.txt
	# remove _{1,2}_val_{1,2}_COI.fa from filenames 076361-M1-S1_BDSW190602952-1a_1_val_1_COI.fa
sed -i 's/_[1,2]_val_[1,2]_COI.fa//g' ./fastalist.txt
cat fastalist.txt
cat fastalist.txt | wc -l # 484 files
cat fastalist.txt | sort | uniq | wc -l # 242 unique names

# make array of fasta files
sample_info=fastalist.txt # put fastalist.txt into variable
sample_names=($(cut -f 1 "$sample_info" | sort | uniq)) # convert variable to array this way
echo "${sample_names[@]}" # echo all array elements
echo "There are" ${#sample_names[@]} "files that will be processed." # 242, echo number of elements in the array

# run Kelpie in parallel on each pair of COI.fa files. -j n means n samples at a time, -k keep same order as in array, --dryrun see the generated commands
# BF3BR2, -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA, -max 500
	nohup parallel -k -j 3 "Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -primers -filtered -min 400 -max 500 allfilterreadsoutput/{1}_?_val_?_COI.fa kelpieoutputindiv/{1}_BF3BR2.fas" ::: "${sample_names[@]}" &
# # takes about 70 mins
ls kelpieoutputindiv/*BF3BR2.fas
ls kelpieoutputindiv/*BF3BR2.fas | wc -l # 242

# alternative code for running kelpie with Leray Fol-degen-rev primers
# Leray Fol-degen-rev, -f GGWACWGGWTGAACWGTWTAYCCYCC -r TANACYTCNGGRTGNCCRAARAAYCA, -min 300 -max 400
	# nohup parallel -k -j 3 "Kelpie_v2 -f GGWACWGGWTGAACWGTWTAYCCYCC -r TANACYTCNGGRTGNCCRAARAAYCA -primers -filtered -min 300 -max 400 allfilterreadsoutput/{1}_?_val_?_COI.fa kelpieoutputindiv/{1}_LERAY.fas" ::: "${sample_names[@]}" &
# takes about 24 mins as a batch job
# ls kelpieoutputindiv/*LERAY.fas
# ls kelpieoutputindiv/*LERAY.fas | wc -l # 242



#### run kelpie on nearest-neighbor sets of files (each sample + five nearest neighbors)
     # read in each line of neighbors_20191204_wide.csv, index i
     # 96 rows
     # save the fields into separate sitename variables
     # run find command on all the sitename and cat the output files into a large filtered.fas
     # run kelpie on it, name the output by the index i, save to kelpieoutputneighbors folder
# http://www.compciv.org/topics/bash/loops/
# https://codefather.tech/blog/bash-loop-through-lines-file/
# https://www.cyberciti.biz/faq/unix-howto-read-line-by-line-from-file/
# https://bash.cyberciti.biz/guide/While_loop#Reading_A_Text_File_With_Separate_Fields

# allfilterreadsoutput/ should already exist from above
cd ~/_Oregon/HJAdryad/2.trimmeddata/ || exit

# upload 03_reference_sequences_datasets/neighbors_20191204_wide.csv into ~/_Oregon/HJAdryad/2.trimmeddata/

if [ ! -d kelpieoutputneighbors ] # if directory kelpieoutputneighbors does not exist.
then
     mkdir kelpieoutputneighbors
fi
ls

# expected to take ~1.2 mins per file (~ 2 hrs total)
i=0 # set index to 0
echo $i
# input each line separately
while IFS=, read -r f1 f2 f3 f4 f5 f6
  do
      ((i=i+1)) # increment i by 1
      echo "starting i should be 1: $i"
      # find any files with these sitenames and concatenate into a single fasta file
      echo "creating kelpieinput_${i}.fa from these files:"
      find ./allfilterreadsoutput -type f -iname "$f1*" -o -iname "$f2*" -o -iname "$f3*" -o -iname "$f4*" -o -iname "$f5*" -o -iname "$f6*"
      echo "$f1, $f2, $f3, $f4, $f5, $f6"

	 # cat the *_1_COI.fa files together and cat the *_2_COI.fa files together
      find ./allfilterreadsoutput -type f -iname "$f1*_1_COI.fa" -o -iname "$f2*_1_COI.fa" -o -iname "$f3*_1_COI.fa" -o -iname "$f4*_1_COI.fa" -o -iname "$f5*_1_COI.fa" -o -iname "$f6*_1_COI.fa" -exec cat {} + > kelpieinput_${i}_1.fa
      find ./allfilterreadsoutput -type f -iname "$f1*_2_COI.fa" -o -iname "$f2*_2_COI.fa" -o -iname "$f3*_2_COI.fa" -o -iname "$f4*_2_COI.fa" -o -iname "$f5*_2_COI.fa" -o -iname "$f6*_2_COI.fa" -exec cat {} + > kelpieinput_${i}_2.fa

      if [ ! -s kelpieinput_${i}_1.fa ] # if kelpieinput_${i}_1.fa has filesize==0, delete it and exit loop
      then
           echo "deleting kelpieinput_${i}_1.fa"
           rm -f kelpieinput_${i}_1.fa || exit
      fi

	  if [ ! -s kelpieinput_${i}_2.fa ] # if kelpieinput_${i}_2.fa has filesize==0, delete it and exit loop
       then
            echo "deleting kelpieinput_${i}_2.fa"
            rm -f kelpieinput_${i}_2.fa || exit
       fi

      if [ -s kelpieinput_${i}_1.fa ] && [ -s kelpieinput_${i}_2.fa ] # if kelpieinput_$i_{1,2}.fa exist and have filesizes > 0
      then # check that the commands inside the then fi statement are preceded only by spaces, no tabs!
           echo "running kelpie on line ${i}"
           # BF3BR2
               Kelpie_v2 -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA -primers -filtered -min 400 -max 500 kelpieinput_${i}_?.fa kelpieoutputneighbors/${i}_BF3BR2.fas
           # Leray-FolDegenRev
               # Kelpie_v2 -f GGWACWGGWTGAACWGTWTAYCCYCC -r TANACYTCNGGRTGNCCRAARAAYCA -primers -filtered -min 300 -max 400 kelpieinput_${i}_?.fa kelpieoutputneighbors/${i}_LERAY.fas
           echo "deleting kelpieinput_${i}_1.fa and kelpieinput_${i}_2.fa"
           rm -f kelpieinput_${i}_{1,2}.fa || exit
      fi
done < <(tail --lines=+2 neighbors_20191204_wide.csv)


#### other primers, with kelpie syntax

# fwhF2 Fol-degen-rev, -f CCHGAYATRGCHTTYCCHCG -r TCDGGRTGNCCRAARAAYCA, -max 500 (because some OTUs seem to get insertions, which i can remove later)

# 12S rRNA (MT-RNR1) (leech project, 82-150 bp), -f ACTGGGATTAGATACCCC -r YRGAACAGGCTCCTCTAG

# 16Smam rRNA (MT-RNR2, 81-117 bp), -f CGGTTGGGGTGACCTCGGA -r GCTGTTATCCCTAGGGTAACT

# trnL_P6 -f GGGCAATCCTGAGCCAA -r CCATYGAGTCTCTGCACCTATC

# 16S rRNA V4 (Earth Microbiome Project) 515F (Parada)–806R (Apprill) -f GTGYCAGCMGCCGCGGTAA -r GGACTACNVGGGTWTCTAAT
