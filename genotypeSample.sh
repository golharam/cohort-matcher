#!/bin/bash
#$ -S /bin/bash
#$ -cwd y
#$ -j y
#$ -pe orte 1

echo Genotyping $BAMFILE
export PATH=/ngs/apps/Python-2.7.8/bin:$PATH
/ngs/apps/bam-matcher/cohort-matcher.py
