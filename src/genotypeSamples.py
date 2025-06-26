#!/usr/bin/env python
'''
This script calls freebayes for each sample.
'''
import argparse
import logging
import sys
import time

import boto3

from common import listFiles, exists

__appname__ = 'genotypeSamples'
__version__ = "0.3"

def read_samples(bamsheet):
    # bamsheet is a tab-separated file with columns: Subject_ID, Sample_ID, bam_path, [reference]
    ''' Read Samples from Bamsheet '''
    samples = []
    with open(bamsheet, 'r') as f:
        # First line is a header line,
        header = f.readline().strip().split('\t')

        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                logging.error("Invalid line in bamsheet: %s", line.strip())
                continue
            sample = {
                'subject_id': parts[0],
                'sample_id': parts[1],
                's3_bam_path': parts[2],
                'reference': parts[3] if len(parts) > 3 else None
            }

            samples.append(sample)
    return samples

def main(argv):
    ''' Main Entry Point '''
    args = parse_arguments(argv)
    logging.basicConfig(level=args.log_level)
    logging.info("%s v%s", __appname__, __version__)
    logging.info(args)

    batch = boto3.client('batch')

    samples = read_samples(args.bamsheet)
    logging.info("Read %d samples from %s", len(samples), args.bamsheet)

    vcf_files = listFiles(args.s3_cache_folder, suffix='.vcf')

    genotyping_jobs = []
    for idx, sample in enumerate(samples):
        sample_id = sample['sample_id']
        vcf_file = "%s/%s.vcf" % (args.s3_cache_folder, sample_id)
        if vcf_file in vcf_files:
            logging.info("[%d/%d] Skipping %s, VCF file already exists: %s", idx+1, len(samples), sample_id, vcf_file)
            continue
        if vcf_file not in vcf_files:
            logging.info("[%d/%d] Genotyping %s", idx+1, len(samples), sample_id)
            if sample['reference'] == 'hg19':
                s3_reference_path = "s3://bmsrd-ngs-repo/cohort-matcher/hg19.tar.bz2"
                s3_targets_path = "s3://bmsrd-ngs-repo/cohort-matcher/hg19.cohort-matcher.bed"
            elif sample['reference'] == 'GRCh37ERCC':
                s3_reference_path = "s3://bmsrd-ngs-repo/cohort-matcher/GRCh37ERCC.tar.bz2"
                s3_targets_path = "s3://bmsrd-ngs-repo/cohort-matcher/GRCh37ERCC.cohort-matcher.bed"
            elif sample['reference'] == 'hg38':
                s3_reference_path = "s3://bmsrd-ngs-repo/cohort-matcher/hg38.tar.bz2"
                s3_targets_path = "s3://bmsrd-ngs-repo/NGSCheckMate/SNP_GRCh38_hg38_wChr.sorted.bed"
            elif sample['reference'] == 'GRCh38ERCC':
                s3_reference_path = "s3://bmsrd-ngs-repo/cohort-matcher/GRCh38ERCC.tar.bz2"
                s3_targets_path = "s3://bmsrd-ngs-repo/NGSCheckMate/SNP_GRCh38_hg38_woChr.sorted.bed"
            elif args.s3_reference_path and args.s3_targets_path:
                s3_reference_path = args.s3_reference_path
                s3_targets_path = args.s3_targets_path
            else:
                logging.warning("Skipping %s due to unknown reference", sample_id)
                continue

            if args.dryRun:
                logging.info("Would call batch.submit_job")
            else:
                response = batch.submit_job(jobName='freebayes-%s' % sample_id,
                                            jobQueue=args.job_queue,
                                            jobDefinition=args.job_definition,
                                            retryStrategy={
                                                'attempts': 2
                                            },
                                            containerOverrides={
                                                'memory': args.memory,
                                                'command': ["--sample_name", sample_id,
                                                            "--bamfile_s3_path", sample['s3_bam_path'],
                                                            "--targets_s3_path", s3_targets_path,
                                                            "--reference_s3_path", s3_reference_path,
                                                            "--s3_output_folder_path", args.s3_cache_folder]
                                            })
                logging.debug(response)
                jobId = response['jobId']
                genotyping_jobs.append(jobId)

    logging.info("Submitted %s jobs", len(genotyping_jobs))
    completed_jobs = []
    failed_jobs = []
    for counter, jobid in enumerate(genotyping_jobs):
        status = ''
        while status != 'SUCCEEDED' and status != 'FAILED':
            logging.info("[%d/%d] Checking job %s", counter+1, len(genotyping_jobs), jobid)
            response = batch.describe_jobs(jobs=[jobid])
            status = response['jobs'][0]['status']
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

def parse_arguments(argv):
    ''' Parse Arguments '''
    parser = argparse.ArgumentParser(description='Genotype a set of samples')
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument('-d', '--dryRun', default=False, action="store_true",
                        help="Simulate job submission")

    required_args = parser.add_argument_group("Required")
    required_args.add_argument('-b', '--bamsheet', required=True, help="Bamsheet")
    required_args.add_argument("-CD", "--s3_cache_folder", required=True,
                               help="Specify S3 path for VCF cache")

    parser.add_argument('-r', '--s3_reference_path', required=False,
                        help='S3 Path to Reference Fastq File (bzip2)')
    parser.add_argument('-t', '--s3_targets_path', required=False,
                        help="S3 Path to Targets BED File")

    job_args = parser.add_argument_group("AWS Batch Job Settings")
    job_args.add_argument('-q', '--job-queue', action="store",
                          default="freeBayesJobQueue-e36f01a082d3ba2",
                          help="AWS Batch Job Queue")
    job_args.add_argument('-j', '--job-definition', action="store", default="freebayes:1",
                          help="AWS Batch Job Definition")
    job_args.add_argument('-m', '--memory', action="store", default=4096, type=int,
                          help="Memory Required (MB)")

    args = parser.parse_args(argv)
    return args

if __name__ == '__main__':
    main(sys.argv[1:])
