#!/usr/bin/env Rscript
### Author: Ryan Golhar <ryan.golhar@bms.com>

### This script finds sample swaps based on genotyping from cohort-matcher. 
### Inputs:
### 1.  Cohort-matcher melted results
### 2.  Sample to Subject (This should be a master list of all samples mapping to subject)
###

### Install required packages to run
source("http://bioconductor.org/biocLite.R")
biocLite(c("argparser", "circlize", "logging"))
#biocLite("canvasXpress")
library(logging)
basicConfig()
###

### Get the command line arguments
library(argparser, quietly=TRUE)
p <- arg_parser("findSampleSwaps")
p <- add_argument(p, "--meltedResults", help="Cohort-matcher melted results", default="meltedResults.txt")
p <- add_argument(p, "--sampleToSubject", help="Sample to subject mapping", default="sampleToSubject.txt")
p <- add_argument(p, "--showSampleLabels", help="Show sample labels in circos plot (Default: False)", flag=TRUE)
p <- add_argument(p, "--showSubjectLabels", help="Show subject labels in circos plot (Default: False)", flag=TRUE)
p <- add_argument(p, "--outputFile", help="Output file to list swaps", default="swaps.txt")
argv <- parse_args(p)
###

# Read cohort matcher results
print(paste("Reading", argv$meltedResults, sep=" "))
cm <- read.table(argv$meltedResults, header=TRUE, sep="\t", stringsAsFactors=FALSE)
cm.same <- cm[ which(cm$Judgement=="SAME"), ]

# Read sample to Subject map
print(paste("Reading", argv$sampleToSubject, sep=" "))
sample_to_subject <- read.table(argv$sampleToSubject, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                                na.strings = "")
colnames(sample_to_subject) <- c("sample", "subject")

print( paste("Read", length(unique(sample_to_subject$subject)), "subjects and",
             length(unique(sample_to_subject$sample)), "samples", sep=" ") )

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
fileConn <- file(argv$outputFile)
# For each sample in the list of samples, get the best match
matches <- data.frame(from_sample=character(), from_subject=character(),
                      to_sample=character(), to_subject=character(),
                      stringsAsFactors = FALSE)
for (sample in sample_to_subject$sample) {
  writeLines(sample)
  sample_matches <- cm[ which((cm$Sample1==sample | cm$Sample2==sample) & cm$Judgement=="SAME"), ]
  # if sample_matches is Empty, then sample doesn't match to another subject
  # for all sample_matches, make sure the sample matches to the same subject
  if (nrow(sample_matches) > 0) {
    writeLines(paste("     Found", nrow(sample_matches), "match(s)", sep=" "))
    for (row_index in 1:nrow(sample_matches)) {
      if (sample_matches[row_index, "Sample1"] == sample) {
        matched_sample <- sample_matches[row_index, "Sample2"]
      } else {
        matched_sample <- sample_matches[row_index, "Sample1"]
      }
      writeLines(paste("     Matched sample is", matched_sample, sep=" "))
      # sample is the sample of interest
      # matched_sample is the matched sample
      # Now, make sure they are from the same subject
      subject1 <- sample_to_subject[sample_to_subject$sample==sample, 'subject']
      subject2 <- sample_to_subject[sample_to_subject$sample==matched_sample, 'subject']
      if (subject1 != subject2) {
        writeLines(paste("     Sample ", sample, " (USUBJID: ", subject1, ") matches to ", matched_sample,
                    " (USUBJID: ", subject2, ")", sep=""))
      }
      matches <- rbind(matches, data.frame(from_sample=sample, from_subject=subject1, to_sample=matched_sample, to_subject=subject2))
    }
  }
}
close(fileConn)

# Fix the data.frame to remove levels
matches <- data.frame(lapply(matches, as.character), stringsAsFactors = FALSE)

# Construct the circos plot
subject = c(structure(matches$from_subject, names=matches$from_sample), structure(matches$to_subject, names=matches$to_sample))
subject = subject[!duplicated(names(subject))]
subject = subject[order(subject, names(subject))]

subject_color = structure(2:(length(unique(subject))+1), names=unique(subject))
sample_color = structure(2:(length(names(subject))+1), names=names(subject))

library(circlize)
# If we don't have a massive # of samples (not sure what massive is yet), this code works
#gap.after = do.call("c", lapply(table(subject), function(i) c(rep(2, i-1), (length(names(subject))+1))))
#circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))
# else use this code
circos.par(gap.degree = 0.1, cell.padding = c(0, 0, 0, 0))

# Draw the samples and connect them
chordDiagram(matches[, c(1, 3)], order = names(subject), grid.col = sample_color,
             directional = 1, annotationTrack = "grid", preAllocateTracks = list(
               list(track.height = 0.02))
)

# Add sample labels (only worth doing if we don't have lots of samples)
if (argv$showSampleLabels) {
  circos.track(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), mean(ylim), sector.index, col = "white", cex = 0.6, facing = "inside", niceFacing = TRUE)
  }, bg.border = NA)
}

# Add the subject track
for(b in unique(subject)) {
  model = names(subject[subject == b])
  if (argv$showSubjectLabels) {
    highlight.sector(sector.index = model, track.index = 1, col = subject_color[b],
                     text = b, text.vjust = -1, niceFacing = TRUE)
  } else {
    highlight.sector(sector.index = model, track.index = 1, col = subject_color[b], niceFacing = TRUE)
  }
}

circos.clear()
print(paste("Output written to", argv$outputFile))
print("Plot saved to Rplots.pdf")

# CanvasXpress
#library(canvasXpress)
#smpAnnot = as.data.frame(subject)
#smpAnnot$sample = names(subject)
#connections = apply(as.matrix(cbind(matches[,1], matches[,3])), 1, as.list)
#cxData = data.frame(var1=rep.int(1, length(smpAnnot$sample)))
#rownames(cxData) = smpAnnot$sample
#result <- canvasXpress(data = t(cxData), 
#             smpAnnot = smpAnnot, 
#             graphType = "Circular",
#             connections = connections,
#             smpOverlays = c("subject", "sample"),
#             ringsOrder = c("subject", "sample"),
#             segregateSamplesBy = list("subject"),
#             showLegend = FALSE,
#             arcSegmentsSeparation = 1,
#             xAxisShow = FALSE)
#htmlwidgets::saveWidget(result, file = "canvasXpress.html")
