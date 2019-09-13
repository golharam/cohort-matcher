#!/bin/bash

# 3.
~/NGS/cohort-matcher/src/compareSamples.py -b bamsheet.txt -CD s3://bmsrd-ngs-results/$PROJECTID/cohort-matcher

# 4.
~/NGS/cohort-matcher/src/mergeResults.py -b bamsheet.txt -CD s3://bmsrd-ngs-results/$PROJECTID/cohort-matcher

# 5.
Rscript --vanilla ~/NGS/cohort-matcher/analysisScripts/findSwaps.R
