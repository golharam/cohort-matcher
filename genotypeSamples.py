#!/usr/bin/env python
'''
This script calls freebayes for each sample.
'''
import argparse
import boto3
import logging
import os
import sys

from common import find_bucket_key, listFiles, readSamples

__appname__ = 'genotypeSamples'
__version__ = "0.1"

logger = logging.getLogger(__appname__)

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logger.info("%s v%s" % (__appname__, __version__))
    logger.info(args)

    batch = boto3.client('batch')

    samples = readSamples(args.bamsheet)
    if samples is False:
        return -1

    vcfFiles = listFiles(args.outputfolder_s3_path, suffix='.vcf')
    for sample in samples:
        sampleName = sample['name']
        vcfFile = "%s/%s.vcf" % (args.outputfolder_s3_path, sampleName)
        if vcfFile not in vcfFiles:
            if args.dryRun:
                logger.info("Would genotype %s", sampleName)
            else:
                logger.info("Genotype %s", sampleName)
                response = batch.submit_job(jobName='freebayes-%s' % sampleName,
                                            jobQueue='ngs-spot-job-queue',
                                            jobDefinition='freebayes:1',
                                            containerOverrides={
                                                'memory': 4096,
                                                'command': ["--sample_name", sampleName,
                                                            "--bamfile_s3_path", sample['bam'],
                                                            "--targets_s3_path", args.targets_s3_path,
                                                            "--reference_s3_path", args.reference_s3_path,
                                                            "--s3_output_folder_path", args.outputfolder_s3_path]
                                            })
                logger.debug(response)
        else:
            logger.info("%s already genotyped.  Skipping...", sampleName)

def parseArguments(argv):
    ''' Parse Arguments '''
    parser = argparse.ArgumentParser(description='Genotype a set of samples')
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    parser.add_argument('-d', '--dryRun', default=False, action="store_true",
                        help="Simulate job submission")
    parser.add_argument('-b', '--bamsheet', required=True, help="Bamsheet")

    parser.add_argument('-r', '--reference_s3_path', required=True,
                        help='S3 Path to Reference Fastq File')
    parser.add_argument('-t', '--targets_s3_path', required=True,
                        help="S3 Path to Targets BED File")

    parser.add_argument('-o', '--outputfolder_s3_path', required=True, 
                        help="Specify S3 path for cached VCF/TSV files")

    args = parser.parse_args(argv)
    return args

if __name__ == '__main__':
    main(sys.argv[1:])
