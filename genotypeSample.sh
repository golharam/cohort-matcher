#!/bin/bash

#$ -cwd
#$ -j y

# These environment variables will be passed in to script:
# SAMPLE
# BAMFILE
# SCRATCH_DIR

echo Genotyping $SAMPLE
echo BAMFILE $BAMFILE
echo SCRATCH_DIR $SCRATCH_DIR
echo

if [[ $BAMFILE == s3* ]]; then
	echo "Download $BAMFILE -> $SCRATCH_DIR"
	aws s3 cp $BAMFILE $SCRATCH_DIR
fi


