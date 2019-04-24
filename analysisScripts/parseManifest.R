#!/usr/bin/env Rscript
###
### This script simply makes a sample to subject mapping file from the manifest
### Inputs:
### 1.  Manifest file
### Outputs:
### 1.  Sample to Subject file

### Read arguments
#install.packages(c('argparser'))
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

### Create the sample to subject file
paste("Writing", sampleToSubjectFile)
sample_to_subject <- data.frame(paste(manifest$VENDORNAME, manifest$VRUNID, sep=""), manifest$USUBJID)
write.table(sample_to_subject, file=sampleToSubjectFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

### Create the RNA-Seq bamsheet
paste("Writing", "rna_bamsheet.txt")
rna_samples <- manifest[manifest$ASSAYMETHOD=="RNA-Seq", ]
rna_samples$sample_name <- paste(rna_samples$VENDORNAME, rna_samples$VRUNID, sep="")
rna_samples$bamfile <- paste("s3://bmsrd-ngs-results/", rna_samples$BMSPROJECTID, "/bam/",
                             rna_samples$sample_name, ".GRCh37ERCC-ensembl75.decontaminated.genome.bam", sep="")
rna_bamsheet <- data.frame(sample=rna_samples$sample_name, bamfile=rna_samples$bamfile, reference="GRCh37ERCC")
write.table(rna_bamsheet, file="rna_bamsheet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

### Create the WES bamsheet
# TODO: Sujaya, please fill in code here
#wes_bamsheet <- data.frame(sample=..., bamfile=..., reference="hg19")
wes_bamsheet <- read.table("wes_bamsheet.txt", header=FALSE, sep="\t")
colnames(wes_bamsheet) <- colnames(rna_bamsheet)

# Merge rna_bamsheet + wes_bamsheet and write out master bamsheet
bamsheet <- rbind(wes_bamsheet, rna_bamsheet)
write.table(bamsheet, file="bamsheet.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

