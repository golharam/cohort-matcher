#!/usr/bin/env Rscript
### Author: Ryan Golhar <ryan.golhar@bms.com>

### This script finds sample swaps based on genotyping from cohort-matcher. 
### Inputs:
### 1.  Cohort-matcher melted results
### 2.  Sample to Patient (This should be a master list of all samples mapping to patient)
###

### Get the command line arguments
library(argparser, quietly=TRUE)
p <- arg_parser("findSampleSwaps")
p <- add_argument(p, "--meltedResults", help="Cohort-matcher melted results", default="meltedResults.txt")
p <- add_argument(p, "--sampleToPatient", help="Sample to patient mapping", default="sampleToPatient.txt")
p <- add_argument(p, "--showSampleLabels", help="Show sample labels in circos plot (Default: False)", flag=TRUE)
p <- add_argument(p, "--showPatientLabels", help="Show patient labels in circos plot (Default: False)", flag=TRUE)
argv <- parse_args(p)
###

# Read cohort matcher results
print(paste("Reading", argv$meltedResults, sep=" "))
cm <- read.table(argv$meltedResults, header=TRUE, sep="\t", stringsAsFactors=FALSE)
cm.same <- cm[ which(cm$Judgement=="SAME"), ]

# Read sample to Patient map
print(paste("Reading", argv$sampleToPatient, sep=" "))
sample_to_patient <- read.table(argv$sampleToPatient, header=FALSE, sep="\t", stringsAsFactors=FALSE, 
                                na.strings = "")
colnames(sample_to_patient) <- c("sample", "patient")

# TODO: Verify the list of samples in meltedResults matches the samples in sampleToPatient

print( paste("Read", length(unique(sample_to_patient$patient)), "patients and", 
             length(unique(sample_to_patient$sample)), "samples", sep=" ") )

# For each sample in the list of samples, get the best match
matches <- data.frame(from_sample=character(), from_patient=character(),
                      to_sample=character(), to_patient=character(), 
                      stringsAsFactors = FALSE)
for (sample in sample_to_patient$sample) {
  print(sample)
  sample_matches <- cm[ which((cm$Sample1==sample | cm$Sample2==sample) & cm$Judgement=="SAME"), ]
  # if sample_matches is Empty, then sample doesn't match to another patient
  # for all sample_matches, make sure the sample matches to the same patient
  if (nrow(sample_matches) > 0) {
    print(paste("     Found", nrow(sample_matches), "match(s)", sep=" "))
    for (row_index in 1:nrow(sample_matches)) {
      if (sample_matches[row_index, "Sample1"] == sample) {
        matched_sample <- sample_matches[row_index, "Sample2"]
      } else {
        matched_sample <- sample_matches[row_index, "Sample1"]
      }
      print(paste("     Matched sample is", matched_sample, sep=" "))
      # sample is the sample of interest
      # matched_sample is the matched sample
      # Now, make sure they are from the same patient
      patient1 <- sample_to_patient[sample_to_patient$sample==sample, 'patient']
      patient2 <- sample_to_patient[sample_to_patient$sample==matched_sample, 'patient']
      print(paste("     ", patient1, "-", patient2, sep=" "))
      if (patient1 != patient2) {
        print(paste("     Sample ", sample, " (USUBJID: ", patient1, ") matches to ", matched_sample, 
                    " (USUBJID: ", patient2, ")", sep=""))
      }
      matches <- rbind(matches, data.frame(from_sample=sample, from_patient=patient1, to_sample=matched_sample, to_patient=patient2))
    }
    
  }
}
# Fix the data.frame to remove levels
matches <- data.frame(lapply(matches, as.character), stringsAsFactors = FALSE)

# Construct the circos plot
patient = c(structure(matches$from_patient, names=matches$from_sample), structure(matches$to_patient, names=matches$to_sample))
patient = patient[!duplicated(names(patient))]
patient = patient[order(patient, names(patient))]

patient_color = structure(2:(length(unique(patient))+1), names=unique(patient))
sample_color = structure(2:(length(names(patient))+1), names=names(patient))

# If we don't have a massive # of samples (not sure what massive is yet), this code works
library(circlize)
gap.after = do.call("c", lapply(table(patient), function(i) c(rep(2, i-1), (length(names(patient))+1))))
circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))
# else use this code
#circos.par(gap.degree = 0.1, cell.padding = c(0, 0, 0, 0))

# Draw the samples and connect them
chordDiagram(matches[, c(1, 3)], order = names(patient), grid.col = sample_color,
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

# Add the patient track
for(b in unique(patient)) {
  model = names(patient[patient == b])
  if (argv$showPatientLabels) {
    highlight.sector(sector.index = model, track.index = 1, col = patient_color[b], 
                     text = b, text.vjust = -1, niceFacing = TRUE)
  } else {
    highlight.sector(sector.index = model, track.index = 1, col = patient_color[b], niceFacing = TRUE)
  }
}

circos.clear()
