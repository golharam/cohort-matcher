#!/usr/bin/env python
'''
Script to compare genotypes

Not every sample is going to have a meltedResults.txt. For a list of samples
(sample1, ..., sampleN), an all v all nieve approach would require O(n^2) comparisons,
but really, we can do this in O(n log(n)).

For example, 5 samples, sample1, ..., sample 5:

	s1	s2	s3	s4	s5
s1
s2	x
s3	x	x
s4	x	x	x
s5	x	x	x	x

Instead of visiting every cell, we only need to visit the ones with a X because the matrix is
symmetrical about the axis
'''
import argparse
import logging
import os
import time
import subprocess
import sys

import boto3
import botocore

from genotypeSamples import read_samples
from common import listFiles, uploadFile

__appname__ = 'compareSamples'
__version__ = "0.3"

def _run(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    p.wait()
    if p.returncode != 0:
        return p.returncode, err
    return p.returncode, out

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logging.info("%s v%s", __appname__, __version__)
    logging.info(args)

    batch = boto3.client('batch')

    samples = read_samples(args.bamsheet)
    if samples is False:
        return -1

    # Make sure VCF files are present before submitting jobs
    vcf_files = listFiles(args.s3_cache_folder, suffix='.vcf')
    ok = True
    for sample in samples:
        vcf_file = "%s/%s.vcf" % (args.s3_cache_folder, sample['sample_id'])
        if vcf_file not in vcf_files:
            logging.error("%s not found.", vcf_file)
            ok = False
    if not ok:
        return -1

    # We don't need the last sample in the list so let's remove it
    samples.pop()

    # Get a list of meltedResults files
    if args.force is True:
        melted_results_files = []
    else:
        melted_results_files = listFiles(args.s3_cache_folder, suffix='.meltedResults.txt')

    # Upload the bamsheet to the cache directory (because each job will need to determine
    # what samples to compare against based on its order in the bamsheet)
    if not args.dry_run:
        uploadFile(args.bamsheet, "%s/bamsheet.txt" % args.s3_cache_folder)

    response = None
    jobs = []
    for idx, sample in enumerate(samples):
        melted_results = "%s/%s.meltedResults.txt" % (args.s3_cache_folder, sample['sample_id'])
        if melted_results not in melted_results_files:
            logging.info("[%d/%d] Comparing genotype of %s to other samples", idx+1, len(samples), sample['sample_id'])
            if args.local:
                cmd = ["%s/compareGenotypes.py" % os.path.dirname(__file__),
                       "-s", sample['sample_id'],
                       "--s3_cache_folder", args.s3_cache_folder,
                       '--allele-freq', args.allele_freqs,
                       "--working_dir", args.working_dir]
                if args.consider_alleles:
                    cmd.append("--consider_alleles")
                if args.dry_run:
                    logging.info(cmd)
                else:
                    ret_code, response = _run(cmd)
                    if ret_code != 0:
                        logging.error(response)
                        return -1
            else:
                if args.dry_run:
                    logging.info("Would call batch.submit_job: compareGenotypes.py -s %s "
                                 "--s3_cache_folder %s", sample['sample_id'], args.s3_cache_folder)
                else:
                    cmd = ['/compareGenotypes.py',
                            '-s', sample['sample_id'],
                            '--s3_cache_folder',
                            args.s3_cache_folder,
                            '--allele-freq', args.allele_freqs]
                    if args.consider_alleles:
                        cmd.append("--consider_alleles")
                    response = batch.submit_job(jobName='compareGenotypes-%s' % sample['sample_id'],
                                                jobQueue=args.job_queue,
                                                jobDefinition=args.job_definition,
                                                containerOverrides={'vcpus': 1,
                                                                    'command': cmd})
                    jobId = response['jobId']
                    jobs.append(jobId)
                logging.debug(response)

    logging.info("Submitted %s jobs", len(jobs))
    completed_jobs = []
    failed_jobs = []
    for counter, jobid in enumerate(jobs):
        status = ''
        while status != 'SUCCEEDED' and status != 'FAILED':
            logging.info("[%d/%d] Checking job %s", counter+1, len(jobs), jobid)
            try:
                response = batch.describe_jobs(jobs=[jobid])
                status = response['jobs'][0]['status']
            except botocore.exceptions.ClientError as err:
                logging.error(err.response)
                pass

            logging.info("Job %s state is %s", jobid, status)
            if status == 'SUCCEEDED':
                completed_jobs.append(jobid)
            elif status == 'FAILED':
                failed_jobs.append(jobid)
            else:
                logging.info("Sleeping 60 secs")
                time.sleep(60)

    logging.info("Successed: %s", len(completed_jobs))
    logging.info("Failed: %s", len(failed_jobs))

def parseArguments(argv):
    ''' Parse arguments '''
    parser = argparse.ArgumentParser(description='Compare a set of samples')
    parser.add_argument('-l', '--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument('-d', '--dry-run', default=False, action="store_true",
                        help="Simulate job submission")
    parser.add_argument('-f', '--force', action="store_true", default=False,
                          help="Force re-run")

    required_args = parser.add_argument_group("Required")
    required_args.add_argument('-b', '--bamsheet', required=True, help="Bamsheet")
    required_args.add_argument("-CD", "--s3_cache_folder", required=True,
                               help="Specify S3 path for cached VCF/TSV files")

    optional_args = parser.add_argument_group("Allele Specific Expression")
    optional_args.add_argument('--consider-alleles', action="store_true", default=False, 
                               help="Consider allele-specific expression")

    optional_args = parser.add_argument_group("Log-Odds Ratio")
    optional_args.add_argument('-af', '--allele-freqs', required=True, help="Path to allele frequencies TSV file")
    optional_args.add_argument('--error-rate', type=float, default=0.01, help="Error rate for genotype comparison")

    job_args = parser.add_argument_group("AWS Batch Job Settings")
    job_args.add_argument('-q', "--job-queue", action="store", default="ngs-job-queue",
                          help="AWS Batch Job Queue")
    job_args.add_argument('-j', '--job-definition', action="store", default="cohort-matcher:2",
                          help="AWS Batch Job Definition")

    job_args = parser.add_argument_group("Local Execution Settings")
    job_args.add_argument('--local', action="store_true", default=False,
                          help="Run locally instead of in AWS Batch")
    job_args.add_argument('--working_dir', type=str, default='/scratch')

    args = parser.parse_args(argv)
    return args

if __name__ == '__main__':
    main(sys.argv[1:])
