#!/bin/bash

PROJECTID=<insert project id here>

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

