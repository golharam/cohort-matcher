#!/usr/bin/env python
'''
Script to compare genotypes 

Not every sample is going to have a meltedResults.txt.   For a list of samples (sample1, ..., sampleN),
an all v all nieve approach would require O(n^2) comparisons, but really, we can do this in O(n log(n)). 

For example, 5 samples, sample1, ..., sample 5:

	s1	s2	s3	s4	s5
s1
s2	x
s3	x	x
s4	x	x	x
s5	x	x	x	x

Instead of visiting every cell, we only need to visit the ones with a X because the matrix is symmetrical
about the axis
'''
import argparse
import boto3
import logging
import os
import time
import subprocess
import sys

from common import find_bucket_key, listFiles, readSamples, uploadFile

__appname__ = 'compareSamples'
__version__ = "0.2"

logger = logging.getLogger(__appname__)

def _run(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    p.wait()
    if p.returncode != 0:
        return err
    return out

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
    # We don't need the last sample in the list so let's remove it
    samples.pop()
 
    # Get a list of meltedResults files
    meltedResultsFiles = listFiles(args.s3_cache_folder, suffix='.meltedResults.txt')

    # Upload the bamsheet to the cache directory (because each job will need to determine what samples 
    #to compare against based on its order in the bamsheet)
    if not args.dry_run:
        uploadFile(args.bamsheet, "%s/bamsheet.txt" % args.s3_cache_folder)

    jobs = []
    for sample in samples:
        meltedResults = "%s/%s.meltedResults.txt" % (args.s3_cache_folder, sample['name'])
        if meltedResults not in meltedResultsFiles:
            logger.info("Comparing genotype of %s to other samples", sample['name'])
            if args.local:
                cmd = ["%s/compareGenotypes.py" % os.path.dirname(__file__),
                       "-s", sample['name'],
                       "--s3_cache_folder", args.s3_cache_folder]
                if args.dry_run:
                    logger.info("Would call %s", cmd)
                else:
                    response = _run(cmd)
            else:
                if args.dry_run:
                    logger.info("Would call batch.submit_job: compareGenotypes.py -s %s --s3_cache_folder %s",
                                sample['name'], args.s3_cache_folder)
                else:
                    response = batch.submit_job(jobName='compareGenotypes-%s' % sample['name'],
                                                jobQueue=args.job_queue,
                                                jobDefinition=args.job_definition,
                                                containerOverrides={'vcpus': 1,
                                                                    'command': ['/compareGenotypes.py',
                                                                                '-s', sample['name'],
                                                                                '--s3_cache_folder',
                                                                                args.s3_cache_folder]})
                    jobId = response['jobId']
                    jobs.append(jobId)
                logger.debug(response)

    logger.info("Submitted %s jobs", len(jobs))
    completed_jobs = []
    failed_jobs = []
    while jobs:
        logger.info("Sleeping 60 secs")
        time.sleep(60)
        logger.info("Checking job %s", jobs[0])
        response = batch.describe_jobs(jobs=[jobs[0]])
        logger.info("Job %s state is %s", jobs[0], response['jobs'][0]['status'])
        if response['jobs'][0]['status'] == 'SUCCEEDED':
            completed_jobs.append(jobs.pop())
        elif response['jobs'][0]['status'] == 'FAILED':
            failed_jobs.append(jobs.pop())
    logger.info("Successed: %s", len(completed_jobs))
    logger.info("Failed: %s", len(failed_jobs))

def parseArguments(argv):
    ''' Parse arguments '''
    parser = argparse.ArgumentParser(description='Compare a set of samples')
    parser.add_argument('-l', '--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument('-d', '--dry-run', default=False, action="store_true",
                        help="Simulate job submission")

    required_args = parser.add_argument_group("Required")
    required_args.add_argument('-b', '--bamsheet', required=True, help="Bamsheet")
    required_args.add_argument("-CD", "--s3_cache_folder", required=True, 
                               help="Specify S3 path for cached VCF/TSV files")

    job_args = parser.add_argument_group("AWS Batch Job Settings")
    job_args.add_argument('--local', action="store_true", default=False,
                          help="Run locally instead of in AWS Batch")
    job_args.add_argument('-q', "--job-queue", action="store", default="ngs-job-queue",
                          help="AWS Batch Job Queue")
    job_args.add_argument('-j', '--job-definition', action="store", default="cohort-matcher:2",
                          help="AWS Batch Job Definition")
    args = parser.parse_args(argv)
    return args

if __name__ == '__main__':
    main(sys.argv[1:])

