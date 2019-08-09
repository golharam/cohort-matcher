#!/bin/bash

# TODO: Provide the manifest file for the project, since this is currently the
# only place to get the sample to subject mapping

MANIFEST=BMS-Manifest-P-20170601-0007-CA209-358-EA-20181204-1.csv

Rscript --vanilla  ~/NGS/cohort-matcher/analysisScripts/parseManifest.R \
  --manifest $MANIFEST

