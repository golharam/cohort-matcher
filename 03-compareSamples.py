#!/usr/bin/env python
'''
Script to submit compareGenotypes.py, if <sample>.meltedResults.txt is not present.
'''
import argparse
import boto3
import logging
import os
import sys

from common import find_bucket_key, listFiles, readSamples

__appname__ = 'compareSamples'
__version__ = "0.1"

logger = logging.getLogger(__appname__)

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logger.info("%s v%s" % (__appname__, __version__))
    logger.info(args)

    batch = boto3.client('batch')

    # Get a list of VCF and meltedResults files
    files = listFiles(args.s3_cache_folder)
    samples = []
    meltedResultsFiles = []
    for f in files:
        if f.endswith('.vcf'):
            sampleid = os.path.basename(f).split('.')[0]
            samples.append(sampleid)
        elif f.endswith('.meltedResults.txt'):
            meltedResultsFiles.append(f)

    jobCount = 0
    for sample in samples:
        meltedResults = "%s/%s.meltedResults.txt" % (args.s3_cache_folder, sample)
        if meltedResults not in meltedResultsFiles:
            logger.info("Comparing genotype of %s to other samples", sample)
            if not args.dry_run:
                response = batch.submit_job(jobName='compareGenotypes-%s' % sample,
                                        jobQueue=args.aws_batch_job_queue,
                                        jobDefinition='cohort-matcher:2',
                                        containerOverrides={
                                            'vcpus': 1,
                                            'command': ['/compareGenotypes.py', '-s', sample,
                                                        '--s3_cache_folder', args.s3_cache_folder]
                                        })
            jobCount += 1

    if args.dry_run:
        logger.info("Would have submitted %s jobs", jobCount)
    else:
        logger.info("Submitted %s jobs", jobCount)

def parseArguments(argv):
    ''' Parse arguments '''
    parser = argparse.ArgumentParser(description='Compare a sample to a set of samples')
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument('--dry-run', default=False, action="store_true",
                        help="Simulates everything, except for actually submitting a job")
    parser.add_argument("--s3_cache_folder", "-CD", required=True, 
                        help="Specify S3 path for cached VCF/TSV files")

    parser.add_argument('-q', "--aws-batch-job-queue", required=True,
                        help="Specific the job queue to submit jobs to")

    args = parser.parse_args(argv)
    return args

if __name__ == '__main__':
    main(sys.argv[1:])

