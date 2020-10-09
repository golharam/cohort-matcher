#!/bin/bash

PROJECTID=<insert project id here>

# Pre-reqs:
# pre-1. Establish your local working directory
# WORKDIR=~/workspace/NGS
# mkdir -p $WORKDIR
# cd $WORKDIR

# pre-2. Clone the BMS-NGS-Python package
# git clone https://biogit.pri.bms.com/NGS/BMS.git
# cd BMS

# pre-3. Create a python 3 virtual environment 
# python3 -m venv env

# pre-4. Active virtual environment
# source env/bin/activate

# pre-5. Install BMS-NGS-Python package into virtual environment
# python setup.py install

# pre-6. Clone cohort-matcher repo
# cd ..
# git clone https://biogit.pri.bms.com/NGS/cohort-matcher.git

# pre-7. Make sure you have Rscript installed and available in your PATH (beyond the scope of this doc).

# Now all pre-reqs are satisified.   Run the steps below on any projects you want to run cohort-matcher on.

# 1. Parse the manifest -> bamsheet.txt, sampleToSubject.txt

source ~/workspace/NGS/BMS/env/bin/activate
manifest download -p $PROJECTID -o manifest.csv
deactivate

Rscript --vanilla ~/workspace/NGS/cohort-matcher/analysisScripts/parseManifest.R --manifest manifest.csv

# 2. Genotype samples
source ~/workspace/NGS/cohort-matcher/env/bin/activate
~/workspace/NGS/cohort-matcher/src/genotypeSamples.py -b bamsheet.txt -o s3://bmsrd-ngs-results/$PROJECTID/cohort-matcher

# 3. Compare sample genotypes (generate per-sample meltedResults)
~/workspace/NGS/cohort-matcher/src/compareSamples.py -b bamsheet.txt -CD s3://bmsrd-ngs-results/$PROJECTID/cohort-matcher

# 4. Merge results
~/workspace/NGS/cohort-matcher/src/mergeResults.py -b bamsheet.txt -CD s3://bmsrd-ngs-results/$PROJECTID/cohort-matcher

# 5. Find swaps
Rscript --vanilla ~/workspace/NGS/cohort-matcher/analysisScripts/findSwaps.R

