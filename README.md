# cohort-matcher #

A workflow for comparing multiple cohorts of [BAM files](https://samtools.github.io/hts-specs/SAMv1.pdf) to determine if they contain reads sequenced from the same samples or patients by counting genotype matches at common SNPs.  Cohort-matcher is an efficient, cloud-enabled variation of BAM-matcher.

# Algorithm #

The basic workflow consists of:
1. Genotype all the samples to be compared. (genotypeSamples.py)
2. Compare the genotypes of each sample against the genotypes of all the other samples. (compareSamples.py which in turn uses compareGenotypes.py to compare a sample to reamining cohort of samples)
3. Merge the results of the sample comparisons (mergeResults.py)
4. Generate plots based on results and known patient-to-sample assocation.

In order to efficiently, some steps are parallelized to reduce runtime.  Specifically:
1.  Genotype each sample independently of each other
2.  Compare a sample's genotype against all other samples (to create a sample's meltedResults file)

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

mergeResults.py created meltedResults.txt, which contains the sample-to-sample comparisons.

# Genome Reference #

The focus of cohort-matcher v2 is on human (hg19 / GRCh37, and hg38 / GRCh38). 
Samples must be mapped against either: 

1) hg19 or GRCh37

OR

2) hg38 or GRCh38

Other combinations of references will not work.  In version 2, the chromosome map has been eliminated, and the VCF to TSV process removes the 'chr' chromosome prefix, if one exists, allowing all VCFs to be compared against each other.

Reference/Target Paths for GRCh37ERCC:
  - s3://bucket/cohort-matcher/GRCh37ERCC.tar.bz2
  - s3://bucket/cohort-matcher/GRCh37ERCC.cohort-matcher.bed
  
Reference/Target Paths for hg19:
  - s3://bucket/cohort-matcher/hg19.tar.bz2
  - s3://bucket/cohort-matcher/hg19.cohort-matcher.bed

## Variant Callers ##

(Require at least one)

* GATK (requires Java)
* VarScan2 (requires Java and Samtools)
* Freebayes

Note: Cohort-matcher only supports Freebayes at this time.

# LICENSE #

The code is released under the Creative Commons by Attribution licence (http://creativecommons.org/licenses/by/4.0/). You are free to use and modify it for any purpose (including commercial), so long as you include appropriate attribution. 

# Citation #

*cohort-matcher* - in prep

# Contact #

Ryan Golhar (ryan.golhar@bms.com)
