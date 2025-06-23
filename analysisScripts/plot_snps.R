#!/usr/bin/env Rscript

meltedResults <- "meltedResults.txt"
cm <- read.table(meltedResults, header=TRUE, sep="\t", stringsAsFactors=FALSE)

hist(cm$Fraction_Match,
     main = "Histogram of Fraction_Match",
     xlab = "Fraction_Match",
     col = "skyblue",
     border = "white")

library(ggplot2)
p <- ggplot(cm, aes(x = Fraction_Match, fill = Judgement)) +
       geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6, color = "black") +
       labs(title = "Histogram of Fraction_Match by Judgement",
         x = "Fraction_Match", y = "Count") +
       theme_minimal()
ggsave("fraction_match_histogram.png", plot = p, width = 6, height = 4, dpi = 300)

