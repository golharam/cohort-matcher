[![Build Status](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher/statusbadges-build/icon)](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher)
[![Code Grade](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher/statusbadges-grade/icon)](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher)
[![Coverage](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher/statusbadges-coverage/icon)](https://jenkins-ci.pri.bms.com:8443/job/cohort-matcher)

# cohort-matcher #

A workflow for comparing multiple cohorts of [BAM files](https://samtools.github.io/hts-specs/SAMv1.pdf) to determine if they contain reads sequenced from the same samples or patients by counting genotype matches at common SNPs.  Cohort-matcher is an efficient, cloud-enabled variation of BAM-matcher.

# Algorithm #

The basic workflow consists of:
1. Genotype all the samples to be compared. (genotypeSamples.py)
2. Compare the genotypes of each sample against the genotypes of all the other samples. (compareSamples.py which in turn uses compareGenotypes.py to compare a sample to reamining cohort of samples)
3. Merge the results of the sample comparisons
4. Generate plots based on results and known patient-to-sample assocation.

In order to efficiently, some steps are parallelized to reduce runtime.  Specifically:
1.  Genotype each sample independently of each other
2.  Compare a sample's genotype against all other samples (to create a sample's meltedResults file)

# Genome Reference #

The focus of cohort-matcher v2 is on human (hg19 / GRCh37, and hg38 / GRCh38). 
Samples must be mapped against either: 

1) hg19 or GRCh37

OR

2) hg38 or GRCh38.

Other combinations of references will not work.  In version 2, the chromosome map has been eliminated, and the VCF to TSV process removes the 'chr' chromosome prefix, if one exists, allowing all VCFs to be compared against each other.

# How to run #

1.  Make input bamsheets

Construct a single 2 column tab-delimited file consisting of sampleName and S3 path to the sample's bamfile, for each set of samples mapped to a specific reference. For instance, if comparing a set of samples all mapped to GRCh38, create 1 TSV file.  If compare two sets of samples mapped to two different references (hg39 and GRCh38), create two TSV files, one for each set of samples.

set1.txt:
sample1 s3://bmsrd-ngs-results/P-12345678-1234/RNA-Seq/bam/sample1.GRCh38ERCC-ensembl91.bam

set2.txt:
sample2 s3://bmsrd-ngs-results/P-12345678-1234/WES/bam/sample2.hg38.bam

## Variant Callers ##

(Require at least one)

* GATK (requires Java)
* VarScan2 (requires Java and Samtools)
* Freebayes

Note: Cohort-matcher only supports Freebayes at this time.

## Installation ##

```
git clone https://github.com/golharam/cohort-matcher
pip install -r cohort-matcher/requirements.txt
```

The repository includes 3 VCF files which can be used for comparing human data (hg19/GRCh37). 

These VCF files also contain variants extracted from 1000 Genomes project which are all exonic and have high likelihood of switching between REF and ALT alleles (global allele frequency between 0.45 and 0.55). The only difference between them is the number of variants contained within.

The repository also includes several BAM files which can be used for testing (under test_data directory), as well as the expected results for various settings.

Cohort-matcher adds unit tests to test the python code.

# LICENSE #

The code is released under the Creative Commons by Attribution licence (http://creativecommons.org/licenses/by/4.0/). You are free to use and modify it for any purpose (including commercial), so long as you include appropriate attribution. 

# Citation #

*cohort-matcher* - in prep

# Contact #

Ryan Golhar (ryan.golhar@bms.com)
