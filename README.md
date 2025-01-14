<!-- confluence-page-id: 462094455 -->
<!-- confluence-space-key: KSR -->
<!-- generated-by: Do not edit. This page was auto-generated from the pipeline repo. -->

# cohort-matcher

A workflow for comparing a cohort of [BAM files](https://samtools.github.io/hts-specs/SAMv1.pdf) to determine if they
contain reads sequenced from the same samples or patients/subjects by counting genotype matches at common SNPs.
Cohort-matcher is an efficient, cloud-enabled variation of BAM-matcher.

## Steps to run pipeline

The basic workflow consists of:

1. Construct a bamsheet, using analysisScripts/parseManifest.R

   - This will create a bamsheet and sample-to-subject file, both required inputs to cohort-matcher.

2. Genotype all the samples to be compared.

   - script: genotypeSamples.py

3. Compare the genotypes of each sample against the genotypes of all the other samples.

   - script: compareSamples.py uses compareGenotypes.py to compare a sample to remaining cohort of samples

4. Merge the results of the sample comparisons

   - script: mergeResults.py

5. Evalulate matches to determine observed sample to subject association.

6. Generate plots based on results and known patient-to-sample association.

## Matching Algorithm

A minimum of 20 comparable SNPs are required to make a reliable comparison between samples (based on several
peer-reviewed papers, references to follow).  We normally observe 100-2000 comparable SNPs between samples, due to
sequencing coverage.

Based on empirical evidence, two samples that originate from the same subject have a % similarity > 70% are considered
a match.  We typically see >70% similarity for RNA-RNA or RNA-DNA, and > 90% (for DNA-DNA).

Two samples originating from different subjects have a % similarity < 50% (RNA-DNA or DNA-DNA). Anything in between is
representative a cross-contamination.

If a sample has 1 match, and the match comes from the same subject (as annotated in the sample to subject map), there is
no swap. If a sample as more than 1 match, then:

- If all the matches are from the same subject, there is no swap.
- If any of the matches are from a different subject, then a swap exists.

## How to run

Refer to the example/ directory

## Output

mergeResults.py creats meltedResults.txt, which contains the sample-to-sample comparisons.
findSwaps.R creates a circos plot of subjects, samples, and swaps, to give a quick visual representation of analysis.
It also creates matches.txt listing the best matched sample.

To quickly find the samples that are swapped, load matches.txt in R and filter on swap column:

```{r}
matches <- read.table("matches.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
swaps <- matches[matches$swap=="x",]
swaps
```

## Genome Reference

cohort-matcher currently works on cohorts of samples mapped to hg19 and/or GRCh37.

Other combinations of references will not work.  In version 2, the chromosome map has been eliminated, and the VCF to
TSV process removes the 'chr' chromosome prefix, if one exists, allowing all VCFs to be compared against each other.

Reference/Target Paths for GRCh37ERCC:

- s3://bmsrd-ngs-repo/cohort-matcher/GRCh37ERCC.tar.bz2
- s3://bmsrd-ngs-repo/cohort-matcher/GRCh37ERCC.cohort-matcher.bed

Reference/Target Paths for hg19:

- s3://bmsrd-ngs-repo/cohort-matcher/hg19.tar.bz2
- s3://bmsrd-ngs-repo/cohort-matcher/hg19.cohort-matcher.bed

## Optimization Notes

In order to efficiently, some steps are parallelized to reduce runtime.  Specifically:

1. Genotype each sample independently of each other
2. Compare a sample's genotype against all other samples (to create a sample's meltedResults file)

## Variant Callers

(Require at least one)

- GATK (requires Java)
- VarScan2 (requires Java and Samtools)
- Freebayes

Note: Cohort-matcher only supports Freebayes at this time.  I haven't tested with GATK or VarScan2.

## Installation

```{bash}
git clone https://github.com/golharam/cohort-matcher
pip install -r cohort-matcher/requirements.txt
```

The repository includes 3 VCF files which can be used for comparing human data (hg19/GRCh37).

These VCF files also contain variants extracted from 1000 Genomes project which are all exonic and have high likelihood
of switching between REF and ALT alleles (global allele frequency between 0.45 and 0.55). The only difference between
them is the number of variants contained within.

The repository also includes several BAM files which can be used for testing (under test_data directory), as well as the
expected results for various settings.

Cohort-matcher adds unit tests to test the python code.

## LICENSE

The code is released under the [Creative Commons by Attribution licence](http://creativecommons.org/licenses/by/4.0/).
You are free to use and modify it for any purpose (including commercial), so long as you include appropriate
attribution.

## Citation

*cohort-matcher* - in prep

## Contact

[Ryan Golhar](ryan.golhar@bms.com)
