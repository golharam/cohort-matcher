#!/usr/bin/env Rscript
###
### This script simply makes a sample to subject mapping file from the manifest
### Inputs:
### 1.  Manifest file
### Outputs:
### 1.  Sample to Subject file

### Read arguments
library(argparser, quietly=TRUE)
p <- arg_parser("cohort-matcher prereq")
p <- add_argument(p, "--manifest", help="Manifest File", default="BMS-Manifest.csv")
p <- add_argument(p, "--sampleToSubject", help="Sample to patient mapping", default="sampleToSubject.txt")
argv <- parse_args(p)

manifestFile <- argv$manifest
sampleToSubjectFile <- argv$sampleToSubject
###

### Read manifest
manifest <- read.csv(manifestFile)
manifest <- manifest[manifest$VRUNIDSUFFIX==1,]
sample_to_subject <- data.frame(manifest$USUBJID, paste(manifest$VENDORNAME, manifest$VRUNID, sep=""))

### Write sample to subject file
write.table(sample_to_subject, file=sampleToSubjectFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

