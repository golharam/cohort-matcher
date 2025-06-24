#!/usr/bin/env Rscript

### Get the command line arguments
library(argparser, quietly=TRUE)
p <- arg_parser("plot_snps.R")

p <- add_argument(p, "--meltedResults", help = "Cohort-matcher melted results", default = "meltedResults.txt")
p <- add_argument(p, "--bamsheet", help = "Cohort-matcher bamsheet", default = "bamsheet.txt")

argv <- parse_args(p)
rm(p)
###

melted_results <- argv$meltedResults
cm <- read.table(melted_results, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Generate a histogram of fraction match
library(ggplot2)
p <- ggplot(cm, aes(x = Fraction_Match, fill = Judgement)) +
       geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6, color = "black") +
       labs(title = "Histogram of Fraction_Match by Judgement",
         x = "Fraction_Match", y = "Count") +
       theme_minimal()
ggsave("fraction_match_histogram.png", plot = p, width = 6, height = 4, dpi = 300)


# Generate a plot of the callable SNPs per sample
# Combine Sample1 and Sample2 into one vector with corresponding counts
samples_s1 <- data.frame(Sample = cm$Sample1, Count = cm$n_S1, stringsAsFactors = FALSE)
samples_s2 <- data.frame(Sample = cm$Sample2, Count = cm$n_S2, stringsAsFactors = FALSE)

# Combine both data frames
all_samples <- rbind(samples_s1, samples_s2)

# Remove duplicates, keeping the first occurrence (i.e., from Sample1 if duplicated)
unique_samples <- all_samples[!duplicated(all_samples$Sample), ]
unique_samples$Sample <- factor(unique_samples$Sample)

p <- ggplot(unique_samples, aes(x = Sample, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
  labs(title = "Counts per Sample", x = "Sample", y = "Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),  # smaller text size
    axis.ticks.length = unit(0.2, "cm")
  )
write.table(unique_samples, "callable_snps_per_sample.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Save plot with increased width
ggsave("callable_snps_per_sample.png", plot = p, width = 16, height = 6, dpi = 300)
