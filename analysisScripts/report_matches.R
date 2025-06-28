#!/usr/bin/env Rscript

### Get the command line arguments
library(argparser, quietly=TRUE)
p <- arg_parser("plot_snps.R")

p <- add_argument(p, "--meltedResults", help = "Cohort-matcher melted results", default = "meltedResults.txt")
p <- add_argument(p, "--bamsheet", help = "Cohort-matcher bamsheet", default = "bamsheet.txt")

argv <- parse_args(p)
rm(p)
###

bamsheet <- argv$bamsheet
samples <- read.table(bamsheet, header = TRUE, sep = "\t")

melted_results <- argv$meltedResults
cm <- read.table(melted_results, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#
get_sample_matches <- function(sample_id, cm) {
    # Get all the rows where the sample is Sample1
    matches1 <- cm[cm$Sample1 == sample_id, ]
    # Get all the rows where the sample is Sample2
    matches2 <- cm[cm$Sample2 == sample_id, ]
    # Swap Sample2 and Sample1 
    colnames(matches2) <- c("Subject2", "Sample2", "Subject1", "Sample1", "n_S2", "n_S1", "SNPs_Compared", "Fraction_Match", "Binomial_PV", "Fraction_Match_Plus", "PV", "Judgement", "Swap")

    # Fix the subset string
    # Step 1: Temporarily replace first string with a placeholder
    matches2$Judgement <- gsub(
    "LIKELY SAME\\. \\(BAM2 is subset of BAM1\\)",
    "TEMP_SWAP_TEXT",
    matches2$Judgement
    )

    # Step 2: Replace second string with the first
    matches2$Judgement <- gsub(
    "LIKELY SAME\\. \\(BAM1 is subset of BAM2\\)",
    "LIKELY SAME. (BAM2 is subset of BAM1)",
    matches2$Judgement
    )

    # Step 3: Replace placeholder with the second
    matches2$Judgement <- gsub(
    "TEMP_SWAP_TEXT",
    "LIKELY SAME. (BAM1 is subset of BAM2)",
    matches2$Judgement
    )

    # Merge the rows
    matches <- rbind(matches1, matches2)

    # Select rows where the Judgement column contains "SAME"
    get_sample_matches <- matches[grepl("SAME|INCONCLUSIVE", matches$Judgement, ignore.case = TRUE), ]
}


outfile <- "report.txt"
if (file.exists(outfile)) {
  file.remove(outfile)
}
for (sample_id in samples$Sample_ID) {
    matches <- get_sample_matches(sample_id, cm)
    if (length(matches) > 0) {
        write.table(matches, file = outfile, sep = "\t", append = TRUE, quote = FALSE, row.names = FALSE)
        cat("\n", file = outfile, append = TRUE)
    }
}
