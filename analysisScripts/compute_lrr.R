if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation")

library(VariantAnnotation)

vcf_file_A <- "sampleA.vcf.gz"
vcf_file_B <- "sampleB.vcf.gz"
allele_freq_file <- "allele_frequencies.tsv"  # CHROM, POS, REF, ALT, FREQ_REF

freqs <- read.table(allele_freq_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(freqs) <- paste(freqs$CHROM, freqs$POS, sep = ":")

vcf_A <- readVcf(vcf_file_A)
vcf_B <- readVcf(vcf_file_B)

geno_prob <- function(geno, p) {
  q <- 1 - p
  if (geno == "AA") return(p^2)
  if (geno == "AB") return(2 * p * q)
  if (geno == "BB") return(q^2)
  return(0)
}

recode_geno <- function(gt, ref, alt) {
  g <- unlist(strsplit(gt, "[|/]"))
  alleles <- c(ref, alt)
  sorted_alleles <- sort(alleles[as.integer(g) + 1])
  geno_str <- paste(sorted_alleles, collapse = "")
  if (geno_str == refref <- paste0(ref, ref)) return("AA")
  if (geno_str == paste0(ref, alt)) return("AB")
  if (geno_str == paste0(alt, alt)) return("BB")
  return(NA)
}

compute_llr <- function(vcfA, vcfB, freqs, error_rate = 0.01) {
  total_llr <- 0

  for (i in 1:nrow(vcfA)) {
    chrom <- as.character(seqnames(vcfA)[i])
    pos <- start(vcfA)[i]
    key <- paste(chrom, pos, sep = ":")

    if (!(key %in% rownames(freqs))) next
    p <- freqs[key, "FREQ_REF"]

    gt_A <- geno(vcfA)$GT[i,1]
    gt_B <- geno(vcfB)$GT[i,1]
    ref <- as.character(ref(vcfA)[i])
    alt <- as.character(unlist(alt(vcfA)[i]))

    geno_A <- recode_geno(gt_A, ref, alt)
    geno_B <- recode_geno(gt_B, ref, alt)
    if (is.na(geno_A) || is.na(geno_B)) next

    L_same <- if (geno_A == geno_B) 1 - error_rate else error_rate
    L_diff <- geno_prob(geno_A, p) * geno_prob(geno_B, p)
    if (L_diff == 0) next

    total_llr <- total_llr + log10(L_same / L_diff)
  }

  return(total_llr)
}

llr <- compute_llr(vcf_A, vcf_B, freqs)
cat(sprintf("Total LLR: %.4f\n", llr))
