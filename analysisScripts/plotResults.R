source("plots.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cohort_matcher_results <- "meltedResults.txt"
#  total_compared_file <- "cohort-matcher-results.total_compared.txt"
} else {
  cohort_matcher_results <- args[1]
#  total_compared_file <- args[2]
}
paste("Reading", cohort_matcher_results, sep=" ")
table <- read.table(cohort_matcher_results, header=TRUE, row.names=1)
#paste("Reading", total_compared_file, sep=" ")
#total_compared <- read.table(total_compared_file, header=TRUE, row.names=1)

reportTopMatches <- function(x, topMatchesFile, ...) {
  f <- file(topMatchesFile, "w")
  writeLines(paste("sample", "match1", "score1", "match2", "score2", "match3", "score3", "match4", "score4", "match5", "score5", sep="\t"), f)
  for (sample in rownames(x)) {
    sample_matches <- sort(x[sample,], decreasing=TRUE)
    writeLines(paste(sample,
                     colnames(sample_matches[1]), sample_matches[1,1],
                     colnames(sample_matches[2]), sample_matches[1,2],
                     colnames(sample_matches[3]), sample_matches[1,3],
                     colnames(sample_matches[4]), sample_matches[1,4],
                     colnames(sample_matches[5]), sample_matches[1,5],
                     sep="\t"),
               f)
  }
  close(f)
}

# ----- END plot function ----- #

topMatchesFile <- gsub(".txt", ".topmatches.txt", cohort_matcher_results)
paste("Writing top matches in", topMatchesFile, sep=" ")
reportTopMatches(table, topMatchesFile)

pdfFile <- gsub(".txt", ".pdf", cohort_matcher_results)
paste("Writing", pdfFile, sep=" ")
pdf(pdfFile)
plotSampleSimilarity(as.matrix(table))
plotNumSNPsCompared(total_compared)
dev.off()

plotFile <- gsub(".txt", ".plot1.tiff", cohort_matcher_results)
tiff(plotFile)
plotSampleSimilarity(as.matrix(table))
dev.off()
plotFile <- gsub(".txt", ".plot2.tiff", cohort_matcher_results)
tiff(plotFile)
plotNumSNPsCompared(total_compared)
dev.off()

