#!/bin/bash

aws s3 cp s3://bmsrd-ngs-repo/cohort-matcher/GRCh37ERCC.cohort-matcher.bed .

mkdir vcfs
for vcf in `aws s3 ls s3://bmsrd-ngs-results/P-20170601-0007/cohort-matcher/ | grep vcf | awk '{print "s3://bmsrd-ngs-results/P-20170601-0007/cohort-matcher/"$4}'`; do
    aws s3 cp $vcf vcfs/
done
ls vcfs/* > vcffiles.txt

source ~/NGS/cohort-matcher/env/bin/activate
~/NGS/cohort-matcher/constructGenotypeFrequencyTable -B GRCh37ERCC.cohort-matcher.bed -L vcffiles.txt > genotypeFrequencyTable.txt
aws s3 cp genotypeFrequencyTable.txt s3://bmsrd-ngs-results/P-20170601-0007/cohort-matcher/ --sse
