#!/usr/bin/env Rscript
### Author: Ryan Golhar <ryan.golhar@bms.com>

### This script finds sample swaps based on genotyping from cohort-matcher. 
### Inputs:
### 1.  Cohort-matcher melted results
### 2.  Sample to Subject (This should be a master list of all samples mapping to subject)
###

### Install required packages to run
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("argparser", "logging"))
library(logging)
basicConfig()
###

### Get the command line arguments
library(argparser, quietly=TRUE)
p <- arg_parser("findSampleSwaps")

p <- add_argument(p, "--installPackages", help="Install packages (Default: False)", flag=TRUE)

p <- add_argument(p, "--meltedResults", help="Cohort-matcher melted results", default="meltedResults.txt")
p <- add_argument(p, "--sampleToSubject", help="Sample to subject mapping", default="sampleToSubject.txt")
p <- add_argument(p, "--outputFile", help="Output file to list swaps", default="swaps.txt")

p <- add_argument(p, "--circosFile", help="Output file for circos plot", default="circos.pdf")
p <- add_argument(p, "--showSampleLabels", help="Show sample labels in circos plot (Default: False)", flag=TRUE)
p <- add_argument(p, "--showSubjectLabels", help="Show subject labels in circos plot (Default: False)", flag=TRUE)
p <- add_argument(p, "--canvasXpress", help="Make canvasXpress plot (Default: False)", flag=TRUE)

argv <- parse_args(p)
rm(p)
###

if (argv$installPackages) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(c("circlize", "devtools"))
  if (argv$canvasXpress) {
    devtools::install_github('neuhausi/canvasXpress')
  }
}

# Read cohort matcher results
loginfo(paste("Reading", argv$meltedResults, sep=" "))
cm <- read.table(argv$meltedResults, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Read sample to Subject map
loginfo(paste("Reading", argv$sampleToSubject, sep=" "))
sample_to_subject <- read.table(argv$sampleToSubject, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                                na.strings = "")
colnames(sample_to_subject) <- c("sample", "subject")

loginfo( paste("Read", length(unique(sample_to_subject$subject)), "subjects and",
             length(unique(sample_to_subject$sample)), "samples", sep=" ") )

# Check that the sample set in cohort matcher results and sampleToSubject map match
samples <- unique(c(cm$Sample1, cm$Sample2))
if (setequal(samples, unique(sample_to_subject$sample)) == FALSE) {
  logerror("Samples in meltedResults and samplesToSubject map don't match")
}
rm(samples)

# Plot the comparison results
# TODO: Since we aren't loading a matrix anymore this code doesn't work.   I need to find a way to
# salvage plotNumSNPsCompared.  That was a useful plot.
#source("plots.R")
#pdf(pdfFile)
#plotSampleSimilarity(as.matrix(table))
#plotNumSNPsCompared(snpsCompared)
#dev.off()

print("Matching samples...")
sink(argv$outputFile)
# For each sample in the list of samples, get the best match
matches <- data.frame(from_sample=character(), from_subject=character(),
                      to_sample=character(), to_subject=character(), swap=character(),
                      stringsAsFactors = FALSE)
for (sample in sample_to_subject$sample) {
  cat(sample, sep="\n")
  sample_matches <- cm[ which((cm$Sample1==sample | cm$Sample2==sample) & cm$Judgement=="SAME"), ]
  # if sample_matches is Empty, then sample doesn't match to another subject
  # for all sample_matches, make sure the sample matches to the same subject
  writeLines(paste("  - Found", nrow(sample_matches), "match(es):", sep=" "))
  if (nrow(sample_matches) > 0) {
    cat(paste("     Found", nrow(sample_matches), "match(s)", sep=" "), sep="\n")
    for (row_index in 1:nrow(sample_matches)) {
      if (sample_matches[row_index, "Sample1"] == sample) {
        matched_sample <- sample_matches[row_index, "Sample2"]
      } else {
        matched_sample <- sample_matches[row_index, "Sample1"]
      }
      cat(paste("     Matched sample is", matched_sample, sep=" "), sep="\n")
      # sample is the sample of interest
      # matched_sample is the matched sample
      # Now, make sure they are from the same subject
      subject1 <- sample_to_subject[sample_to_subject$sample==sample, 'subject']
      subject2 <- sample_to_subject[sample_to_subject$sample==matched_sample, 'subject']
      if (subject1 == subject2) {
        swap <- ""
      } else {
        writeLines(paste("    - Sample ", sample, " (USUBJID: ", subject1, ") matches to ", 
                         matched_sample, " (USUBJID: ", subject2, ")", sep=""))
        swap <- "x"
      }
      matches <- rbind(matches, data.frame(from_sample=sample, from_subject=subject1,
                                           to_sample=matched_sample, to_subject=subject2,
                                           swap=swap))
    }
  }
}
sink()

# Clean up
rm(matched_sample, row_index, sample, subject1, subject2, swap)
loginfo(paste("Writing results to", argv$outputFile, sep=" "))
write.table(matches, argv$outputFile, quote=FALSE, sep="\t", row.names=FALSE)

# Fix the data.frame to remove levels (for circlize)
matches <- data.frame(lapply(matches, as.character), stringsAsFactors = FALSE)

# matches contain duplicate rows of (from_sample, to_sample) and (to_sample, from_sample).   
# We need to remove these duplicates before plotting,
# https://stackoverflow.com/questions/25297812/pair-wise-duplicate-removal-from-dataframe
cols = c(1, 2, 3, 4)
newdf = matches[,cols]
for (i in 1:nrow(matches)){
  newdf[i, ] = sort(matches[i,cols])
}
matches <- matches[!duplicated(newdf),]
rm(cols, newdf)

## Construct the circos plot
loginfo(paste("Saving circos plot to", argv$circosFile, sep=" "))
library(circlize)

subjects = c(structure(matches$from_subject, names=matches$from_sample),
             structure(matches$to_subject, names=matches$to_sample))
subjects = subjects[!duplicated(names(subjects))]
subjects = subjects[order(subjects, names(subjects))]

subject_color = structure(2:(length(unique(subjects))+1), names=unique(subjects))
sample_color = structure(2:(length(names(subjects))+1), names=names(subjects))

pdf(argv$circosFile)

# Set up the circos parameters:
# 1. Put a small gap between each subject
if (length(subjects) < 100) {
  # If we don't have a massive # of samples (not sure what massive is yet), this code works
  gap.after = do.call("c", lapply(table(subjects), function(i) c(rep(2, i-1),
                      (length(names(subjects))+1))))
  circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))
} else {
  # else use this code
  circos.par(gap.degree = 0.1, cell.padding = c(0, 0, 0, 0))
}

# Draw the samples and connect them
chordDiagram(matches[, c(1, 3)],
             order = names(subjects),
             grid.col = sample_color,
             directional = 1,
             annotationTrack = "grid",
             preAllocateTracks = list(
               list(track.height = 0.02)
             )
             #preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(matches)))))
)

# Add sample labels (only worth doing if we don't have lots of samples)
if (argv$showSampleLabels) {
  #circos.track(track.index = 2, panel.fun = function(x, y) {
  #  xlim = get.cell.meta.data("xlim")
  #  ylim = get.cell.meta.data("ylim")
  #  sector.index = get.cell.meta.data("sector.index")
  #  circos.text(mean(xlim), mean(ylim), sector.index, col = "white", cex = 0.6, 
  #              facing = "inside", niceFacing = TRUE)
  #}, bg.border = NA)

  # Track index needs to be added here to move sample labels to the outside
  circos.track(ylim = c(0, 1), 
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim, CELL_META$sector.index,
                             facing = "clockwise", niceFacing = TRUE, adj = c(1, 0.5), cex = 0.6)
               },
               track.height = 0.01, bg.border = NA)
}

# Add the subject track
for(b in unique(subjects)) {
  model = names(subjects[subjects == b])
  if (argv$showSubjectLabels) {
    highlight.sector(sector.index = model, track.index = 1, col = subject_color[b],
                     niceFacing = TRUE,
                     facing = "bending.inside", text = b, text.vjust = -1)
  } else {
    highlight.sector(sector.index = model, track.index = 1, col = subject_color[b],
                     niceFacing = TRUE)
  }
}
circos.clear()
dev.off()

# CanvasXpress
if (argv$canvasXpress) {
  library(canvasXpress)
  smpAnnot = as.data.frame(subjects)
  smpAnnot$sample = names(subjects)
  connections = apply(as.matrix(cbind(matches[,1], matches[,3])), 1, as.list)
  cxData = data.frame(var1=rep.int(1, length(smpAnnot$sample)))
  rownames(cxData) = smpAnnot$sample
  result <- canvasXpress(data = t(cxData), 
             smpAnnot = smpAnnot, 
             graphType = "Circular",
             connections = connections,
             smpOverlays = c("subject", "sample"),
             ringsOrder = c("subject", "sample"),
             segregateSamplesBy = list("subject"),
             showLegend = FALSE,
             arcSegmentsSeparation = 1,
             xAxisShow = FALSE)
  htmlwidgets::saveWidget(result, file = "canvasXpress.html")
}
