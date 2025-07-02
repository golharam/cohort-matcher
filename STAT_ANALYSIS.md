# Goal:

Determine whether two samples are likely from the same individual or different individuals, 
based on the number of matching genotypes out of 91 SNPs.

## Method 1: Binomial Test for Identity

Assume:

1. If the two samples are from the same individual, the expected match rate = ~1 (barring technical errors).

2. If the samples are from different individuals, then the expected match rate depends on population allele frequencies.

Binomial Model:

Let:

n = number of SNPs
k = number of matching SNPs
p0 = expected match probability under null hypothesis

If you assume:

H0: samples are unrelated individuals -> match rate p0 ~ 0.5 (random genotypes)
H1: sample are identical -> match rate p1 ~ 1.0

Then you can compute:

```{r}
binom.test(k, n, p = 0.5, alternative = "greater")
```

This tests whether the number of matches k is significantly higher than expected under the unrelated hypothesis.

## Method 2: Likelihood Ratio Test (Preferred for Identity by State)

You compute:

Lsame: likelihood of the observed genotypes assuming same individual
Ldiff: likelihood assuming unrelated individuals

Then compute the log likelihood ratio (LLR):

LRR = log10 (Lsame/Ldiff)

If the LLR is high (e.g. >10), that strongly supports identity.

This requires:

1. Knowing allele frequencies for each SNP in the population
2. Genotype call confidence (optional)
3. Assuming independence between SNPs

This approach is used in tools like KING, PLINK IBS, and relatedness estimators.

### Concept

You are comparing the likelihood of observing a set of genotypes between two samples under two hypotheses:

1. Hsame: The two samples come from the same individual -> expect identical genotypes (except for errors).
2. Hdiff: The samples come from unrelated individuals -> genotype concordance depends on population allele frequencies.

For each SNP, you compute:

- Lsample = SUM(i=1, n) P(G_Ai, G_Bi | Hsame)
- Ldiff = SUM(i=1, n) P(G_Ai, G_Bi | Hdiff)

Then: LRR = log10(Lsame/Ldiff)

### What you need

For each SNP i:

- The observed genotype for both samples (e.g. AA, AG, GG)
- The allele frequencies in the population (e.g. A: 0.7, G: 0.3)

### How to compute each likelihood

1. P(G_Ai, G_Bi | Hsame)

If samples are from the same individual, assume:

- With no error, genotypes should match → probability = 1 if equal, else 0.
- With error model (e.g., 1% chance of mismatch), you can use:
  P(G_A, G_B | Hsame) =
    1 - e if G_A = G_B
    e     if G_A != G_B

For example, if error rate e = 0.01, then:

- Matching genotype → likelihood = 0.99
- Mismatch → likelihood = 0.01

Or more finely, incorporate genotype likelihoods if available.

2. P(G_Ai, G_Bi | Hdiff)

Assume independent draws from population:

P(G_A, G_B | Hdiff) = P(G_A) x P(G_B)

where: P(G) is computed using Hardy-Weinberg equilibrium:

If allele A has frequency p, G has frequency:

- AA: p^2
- AG: 2pg
- GG: q^2, where q=1-p

Example Calculation:

Let's say:

- SNP1 Allele A (p=0.7), G(q=0.3)
- Sample A: AA
- Sample B: AA

Compute:

Under same-individual (H_same):

- Genotypes match -> P = 0.99

Under unrelated (H_diff):

- P(AA) = 0.7^2 = 0.49
- So P(G_A, G_B | H_diff) = 0.49 x 0.49 = 0.2401

Log likelihood ratio (base 10): log10(0.99 / 0.2401) ~ log10(4.12) ~ 0.615

You do this for each SNP and sum the log10 LLRs:

Total LLR = SUM(i=1,91) log10( P(G_Ai, G_Bi | Hsame) / P(G_Ai, G_Bi | Hdiff) )


## Method 3: Match Score / Hamming Distance

Match score = k / n

Compare this to an empirical distribution of match scores from known related/unrelated individuals

## Example

Let's say:

1. You observe 88/91 SNPs match
2. Assume 0.5 match probability under random individuals

```{r}
binom.test(88, 91, p = 0.5, alternative = "greater")
```

This gives a very low p-value → strong evidence of identity.

## Summary

| Approach              | Assumptions                | Requires Allele Frequencies | Output                |
| -------------------   | -----------------------    | --------------------------- | ------                |
| Binomial test         | Equal random match prob    | No                          | p-value               |
| Likelihood ratio test | Genotype probability model | Yes                         | log-likelihood ratio  |
| Match score           | Empirical                  | No                          | % match or similarity |

# Code

## Assumptions

You have a table of SNPs with:

1. genotype_A: genotype of sample A (e.g. "AA", "AG", "GG")
2. genotype_B: genotype of sample B
3. allele1 and allele2: e.g. "A" and "G"
4. freq1: frequency of allele1 (e.g. 0.7)

You assume a small error rate for mismatches under the same individual (default: e= 0.01)

## Python Script

```{python}

import math
import pandas as pd

# Sample input DataFrame
# You should replace this with your own input file or structure
data = pd.DataFrame({
    'SNP': ['rs1', 'rs2', 'rs3'],
    'allele1': ['A', 'C', 'G'],
    'allele2': ['G', 'T', 'C'],
    'freq1': [0.7, 0.6, 0.4],
    'genotype_A': ['AA', 'CT', 'GG'],
    'genotype_B': ['AA', 'TT', 'GG']
})

def genotype_prob(geno, p):
    q = 1 - p
    if geno == 'AA':
        return p ** 2
    elif geno == 'AB' or geno == 'BA':
        return 2 * p * q
    elif geno == 'BB':
        return q ** 2
    else:
        return 0.0

def simplify(geno):
    # Sort alleles so AG == GA
    return ''.join(sorted(geno))

def compute_llr(data, error_rate=0.01):
    total_llr = 0.0

    for _, row in data.iterrows():
        gA = simplify(row['genotype_A'])
        gB = simplify(row['genotype_B'])
        p = row['freq1']

        # Likelihood under H_same
        if gA == gB:
            L_same = 1 - error_rate
        else:
            L_same = error_rate

        # Likelihood under H_diff
        P_A = genotype_prob(gA, p)
        P_B = genotype_prob(gB, p)
        L_diff = P_A * P_B

        # Avoid log(0)
        if L_diff == 0:
            continue

        llr = math.log10(L_same / L_diff)
        total_llr += llr

    return total_llr

llr = compute_llr(data)
print(f"Total Log Likelihood Ratio: {llr:.4f}")

```

## R script

```{r}
# Sample input data.frame
data <- data.frame(
  SNP = c("rs1", "rs2", "rs3"),
  allele1 = c("A", "C", "G"),
  allele2 = c("G", "T", "C"),
  freq1 = c(0.7, 0.6, 0.4),
  genotype_A = c("AA", "CT", "GG"),
  genotype_B = c("AA", "TT", "GG"),
  stringsAsFactors = FALSE
)

genotype_prob <- function(geno, p) {
  q <- 1 - p
  geno <- paste(sort(strsplit(geno, "")[[1]]), collapse = "")
  if (geno == "AA") {
    return(p^2)
  } else if (geno == "AB") {
    return(2 * p * q)
  } else if (geno == "BB") {
    return(q^2)
  } else {
    return(0)
  }
}

simplify_genotype <- function(g) {
  paste(sort(unlist(strsplit(g, ""))), collapse = "")
}

compute_llr <- function(data, error_rate = 0.01) {
  total_llr <- 0

  for (i in 1:nrow(data)) {
    gA <- simplify_genotype(data$genotype_A[i])
    gB <- simplify_genotype(data$genotype_B[i])
    p <- data$freq1[i]

    # Under same individual
    L_same <- if (gA == gB) 1 - error_rate else error_rate

    # Under unrelated
    P_A <- genotype_prob(gA, p)
    P_B <- genotype_prob(gB, p)
    L_diff <- P_A * P_B

    if (L_diff > 0) {
      total_llr <- total_llr + log10(L_same / L_diff)
    }
  }

  return(total_llr)
}

llr <- compute_llr(data)
cat(sprintf("Total Log Likelihood Ratio: %.4f\n", llr))
```
