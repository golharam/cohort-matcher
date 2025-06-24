# cohort-matcher #

A workflow for comparing a cohort of samples [BAM files](https://samtools.github.io/hts-specs/SAMv1.pdf) to determine if they contain reads sequenced from the same samples or patients/subjects by counting genotype matches at common SNPs. Cohort-matcher is an efficient, cloud-enabled variation of BAM-matcher.

# Steps to run pipeline #

The basic workflow consists of:
1. Construct a bamsheet (using analysisScripts/parseManifest.R)
   - This will create a "bamsheet", which consists of a subject id, sample id, bamfile path, genome reference.
   - This file contains a header line.

2. Genotype all the samples to be compared.
  - script: genotypeSamples.py
  - If genome reference is not provided, a s3 path must be provided on the command line, and all samples must be mapped to provided reference.
  - If genome reference is provided in the bamsheet, then samples mapped to different versions of the same genome can be used (e.g. hg38 or GRCh38)

3. Compare the genotypes of each sample against the genotypes of all the other samples.
  - script: compareSamples.py uses compareGenotypes.py to compare a sample to remaining cohort of samples

4. Merge the results of the sample comparisons
  - script: mergeResults.py

5. Generate circos plot of subjects & samples

a. Implementation based on circlize: https://github.com/jokergoo/circlize
b. Implementation based on pyCirclize: https://github.com/moshi4/pyCirclize
  - script: work/cohort-matcher/analysisScripts/pycircos.py


6. Evalulate matches to determine observed sample to subject association.

7. Generate plots based on results and known patient-to-sample association.

# Matching Algorithm #

A minimum of 20 comparable SNPs are required to make a reliable comparison between samples (based on several peer-reviewed papers, references to follow).  We normally observe 100-2000 comparable SNPs between samples, due to sequencing coverage.

Based on empirical evidence, two samples that originate from the same subject have a % similarity > 70% are considered a match.  We typically see >70% similarity for RNA-RNA or RNA-DNA, and > 90% (for DNA-DNA).

Two samples originating from different subjects have a % similarity < 50% (RNA-DNA or DNA-DNA). Anything in between is representative a cross-contamination.

If a sample has 1 match, and the match comes from the same subject (as annotated in the sample to subject map), there is no swap.
If a sample as more than 1 match, then:
  - If all the matches are from the same subject, there is no swap.
  - If any of the matches are from a different subject, then a swap exists.

## Installation ##

```
git clone https://github.com/golharam/cohort-matcher
pip install -r cohort-matcher/requirements.txt
```

The repository includes 3 VCF files which can be used for comparing human data (hg19/GRCh37). 

These VCF files also contain variants extracted from 1000 Genomes project which are all exonic and have high likelihood of switching between REF and ALT alleles (global allele frequency between 0.45 and 0.55). The only difference between them is the number of variants contained within.

The repository also includes several BAM files which can be used for testing (under test_data directory), as well as the expected results for various settings.

Cohort-matcher adds unit tests to test the python code.

# How to run #

Pre-req:  Make input bamsheet

Construct a single 3 column tab-delimited text file consisting of sampleName, S3 path to the sample bamfile, and reference sample is mapped to (hg19 or GRCh37ERCC) for all the samples. For example:

bamsheet.txt:

| sample  | s3 path to bamfile | reference |
| ------------- | ------------- | ----- |
| sample1 | s3://<s3bucket>/<path_to_sample_bam> | GRCh37ERCC |
| sample2 | s3://<s3bucket>/<path_to_sample_bam> | hg19 |


1.  Call genotypeSamples.py

```
genotypeSamples.py -b bamsheet.txt -o s3://<s3bucket/<path-for-cohort-matcher-cache-and-results>
```

2.  Call compareSamples.py

```
compareSamples.py -b bamsheet.txt -CD s3://<s3bucket/<path-for-cohort-matcher-cache-and-results>
```

3.  Call mergeResults.py

```
mergeResults.py -b bamsheet.txt -CD s3://<s3bucket/<path-for-cohort-matcher-cache-and-results>
```

4.  Call findSwaps.R
```
Rscript analysisScripts/findSwaps.R
```
or via Docker
```
docker run -ti --rm -v $PWD:/work -w /work -v /home/ec2-user/NGS/cohort-matcher:/cohort-matcher cohort-matcher-r Rscript /cohort-matcher/analysisScripts/findSwaps.R
```

## Output ##

mergeResults.py creats meltedResults.txt, which contains the sample-to-sample comparisons.
findSwaps.R creates a circos plot of subjects, samples, and swaps, to give a quick visual representation of analysis.  It also creates matches.txt listing the best matched sample.

To quickly find the samples that are swapped, load matches.txt in R and filter on swap column:
```
matches <- read.table("matches.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
swaps <- matches[matches$swap=="x",]
swaps
```

# Genome Reference #

cohort-matcher currently works on cohorts of samples mapped to hg19 and/or GRCh37.

Other combinations of references will not work.  In version 2, the chromosome map has been eliminated, and the VCF to TSV process removes the 'chr' chromosome prefix, if one exists, allowing all VCFs to be compared against each other.

Reference/Target Paths for GRCh37ERCC:
  - s3://bucket/cohort-matcher/GRCh37ERCC.tar.bz2
  - s3://bucket/cohort-matcher/GRCh37ERCC.cohort-matcher.bed
  
Reference/Target Paths for hg19:
  - s3://bucket/cohort-matcher/hg19.tar.bz2
  - s3://bucket/cohort-matcher/hg19.cohort-matcher.bed

## Optimization Notes ##

In order to efficiently, some steps are parallelized to reduce runtime.  Specifically:
1.  Genotype each sample independently of each other
3.  Compare a sample's genotype against all other samples (to create a sample's meltedResults file)

## Variant Callers ##

(Require at least one)

* GATK (requires Java)
* VarScan2 (requires Java and Samtools)
* Freebayes

Note: Cohort-matcher only supports Freebayes at this time.  I haven't tested with GATK or VarScan2.

# LICENSE #

The code is released under the Creative Commons by Attribution licence (http://creativecommons.org/licenses/by/4.0/). You are free to use and modify it for any purpose (including commercial), so long as you include appropriate attribution. 

# Citation #

*cohort-matcher* - in prep

# Contact #

Ryan Golhar (ryan.golhar@bms.com)
