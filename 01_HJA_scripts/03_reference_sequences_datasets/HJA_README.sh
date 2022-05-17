# FilterReads, which reduces the sequence files to those that match COI sequences
# run code in 02_kelpie/1_launch_FilterReads.sh

# Kelpie
# run code in 2_parallel_kelpie_20200917.sh



# test run scripts in 02_kelpie/ on HJAdryad/2.trimmeddata/BWA*
# if successful, upload HJAdryad/2.trimmeddata/ to S3
  # tar HJAdryad/2.trimmeddata/
  # create new bucket on Amazon S3
  # upload HJAdryad/2.trimmeddata.tar to Amazon S3
  # make public on S3
  # initiate download by datadryad
  # upload filterreads, kelpie, and other used binaries to datadryad




#######################################################################################
#######################################################################################
# set a file public on Amazon S3
# this was run on macOS
#######################################################################################
#######################################################################################

s3cmd ls s3://amazon-oregon-douglasyu
# s3://amazon-oregon-douglasyu/2019Sep_shotgun_2.trimmeddata.tar
s3cmd setacl --acl-public s3://amazon-oregon-douglasyu/2019Sep_shotgun_2.trimmeddata.tar

# message back from Amazon:
	# "s3://amazon-oregon-douglasyu/2019Sep_shotgun_2.trimmeddata.tar: ACL set to Public  [1 of 1]"

# the public URL then follows this format:
	# https://<bucket-name>.s3.amazonaws.com/<object or key name>
https://amazon-oregon-douglasyu.s3.amazonaws.com/2019Sep_shotgun_2.trimmeddata.tar
# or from Transmit
https://amazon-oregon-douglasyu.s3.eu-west-1.amazonaws.com/2019Sep_shotgun_2.trimmeddata.tar


s3cmd ls --list-md5 s3://amazon-oregon-douglasyu/
