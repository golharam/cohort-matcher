#!/bin/bash

# 3.
~/NGS/cohort-matcher/compareSamples.py -b bamsheet.txt -CD s3://bmsrd-ngs-results/P-20170601-0007/cohort-matcher

# 4.
~/NGS/cohort-matcher/mergeResults.py -b bamsheet.txt -CD s3://bmsrd-ngs-results/P-20170601-0007/cohort-matcher

# 5.
Rscript --vanilla ~/NGS/cohort-matcher/analysisScripts/findSwaps.R
