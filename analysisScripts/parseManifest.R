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
p <- add_argument(p, "--manifest", help="[Input] Manifest File", default="BMS-Manifest.csv")
p <- add_argument(p, "--bamsheet", help="[Output] BAM sheet", default="bamsheet.txt")
p <- add_argument(p, "--sampleToSubject", help="[Output] Sample to patient mapping", default="sampleToSubject.txt")
argv <- parse_args(p)

manifestFile <- argv$manifest
sampleToSubjectFile <- argv$sampleToSubject
bamsheetFile <- argv$bamsheet
###

### Read manifest file
manifest <- read.csv(manifestFile, header = TRUE)

# Select 1 row for each sample 
manifest <- manifest[manifest$VRUNIDSUFFIX==1,]

# Get WGS Samples
wgsSamples <- subset(manifest, ASSAYMETHOD=="WGS", 
                     select=c("USUBJID", "SPECTYPE", "ASSAYMETHOD", "VRUNID", "VENDORNAME", "BMSPROJECTID"))
wgsNormal <- subset(wgsSamples, grepl("NORMAL", SPECTYPE), select=c("USUBJID", "VRUNID"))
wgsTumor <- subset(wgsSamples, grepl("TUMOR", SPECTYPE), select=c("USUBJID", "VRUNID"))

# Get WES Samples
wesSamples <- subset(manifest, ASSAYMETHOD=="WES", 
                     select=c("USUBJID", "SPECTYPE", "ASSAYMETHOD", "VRUNID", "VENDORNAME", "BMSPROJECTID"))
wesNormal <- subset(wesSamples, grepl("NORMAL", SPECTYPE), select=c("USUBJID", "VRUNID"))
wesTumor <- subset(wesSamples, grepl("TUMOR", SPECTYPE), select=c("USUBJID", "VRUNID"))

# Get RNA-Seq Samples
rnaSamples <- subset(manifest, ASSAYMETHOD=="RNA-Seq", 
                     select=c("USUBJID", "SPECTYPE", "ASSAYMETHOD", "VRUNID", "VENDORNAME", "BMSPROJECTID"))
rnaNormal <- subset(rnaSamples, grepl("NORMAL", SPECTYPE), select=c("USUBJID", "VRUNID"))
rnaTumor <- subset(rnaSamples, grepl("TUMOR", SPECTYPE), select=c("USUBJID", "VRUNID"))

# Match WES Tumor - Normal
#wesMatchedSamples <- merge(wesTumor, wesNormal, by="USUBJID", all = TRUE)
#colnames(wesMatchedSamples) <- c("USUBJID", "WES.Tumor.VRUNID", "WES.Normal.VRUNID")
  
# Match RNA-Seq Tumor - Normal
#rnaSeqMatchedSamples <- merge(rnaTumor, rnaNormal, by="USUBJID", all = TRUE)
#colnames(rnaSeqMatchedSamples) <- c("USUBJID", "RNASeq.Tumor.VRUNID", "RNASeq.Normal.VRUNID")

# Print out sample counts
message("WGS Tumor: ", dim(wgsTumor)[1])
message("WGS Normal: ", dim(wgsNormal)[1])
message("WGS Samples: ", dim(wgsSamples)[1])

message("WES Tumor: ", dim(wesTumor)[1])
message("WES Normal: ", dim(wesNormal)[1])
message("WES Samples: ", dim(wesSamples)[1])

message("RNA Tumor: ", dim(rnaTumor)[1])
message("RNA Normal: ", dim(rnaNormal)[1])
message("RNA Samples: ", dim(rnaSamples)[1])

message("Total Samples: ", dim(manifest)[1])

# Create the sample to subject file
sample_to_subject <- data.frame(sample=paste(manifest$VENDORNAME, manifest$VRUNID, sep=""), subject=manifest$USUBJID)
paste("Writing", sampleToSubjectFile)
write.table(sample_to_subject, file=sampleToSubjectFile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

bamsheet <- data.frame(sample=character(), bamfile=character(), reference=character())
if (dim(wgsSamples)[1] > 0) {
  wgsSamples$sample_name <- paste(wgsSamples$VENDORNAME, wgsSamples$VRUNID, sep="")
  wgsSamples$bamfile <- paste("s3://bmsrd-ngs-results/", wgsSamples$BMSPROJECTID, "/WGS/hg19/BAM/",
                               wgsSamples$sample_name, ".sorted.dedup.realigned.recal.hg19.bam", sep="")
  wgs_bamsheet <- data.frame(sample=wgsSamples$sample_name, bamfile=wgsSamples$bamfile, reference="hg19")
  bamsheet <- rbind(bamsheet, wgs_bamsheet)
}

# Create the WES bamsheet
if (dim(wesSamples)[1] > 0) {
  wesSamples$sample_name <- paste(wesSamples$VENDORNAME, wesSamples$VRUNID, sep="")
  wesSamples$bamfile <- paste("s3://bmsrd-ngs-results/", wesSamples$BMSPROJECTID, "/WES/hg19/BAM/",
                               wesSamples$sample_name, ".sorted.dedup.realigned.recal.hg19.bam", sep="")
  wes_bamsheet <- data.frame(sample=wesSamples$sample_name, bamfile=wesSamples$bamfile, reference="hg19")
  bamsheet <- rbind(bamsheet, wes_bamsheet)
}

# Create the RNA-Seq bamsheet and merge to wes bamsheet
if (dim(rnaSamples)[1] > 0) {
  rnaSamples$sample_name <- paste(rnaSamples$VENDORNAME, rnaSamples$VRUNID, sep="")
  rnaSamples$bamfile <- paste("s3://bmsrd-ngs-results/", rnaSamples$BMSPROJECTID, "/bam/",
                               rnaSamples$sample_name, ".GRCh37ERCC-ensembl75.decontaminated.genome.bam", sep="")
  rna_bamsheet <- data.frame(sample=rnaSamples$sample_name, bamfile=rnaSamples$bamfile, reference="GRCh37ERCC")
  bamsheet <- rbind(bamsheet, rna_bamsheet)
}

# Write out the master bamsheet
paste("Writing", bamsheetFile)
write.table(bamsheet, file=bamsheetFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
