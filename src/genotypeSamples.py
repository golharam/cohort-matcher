#!/usr/bin/env python
'''
This script calls freebayes for each sample.
'''
import argparse
import logging
import sys
import time

import boto3

from common import listFiles, readSamples, exists

__appname__ = 'genotypeSamples'
__version__ = "0.2"

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logging.info("%s v%s", __appname__, __version__)
    logging.info(args)

    batch = boto3.client('batch')

    samples = readSamples(args.bamsheet)
    if samples is False:
        return -1

    vcfFiles = listFiles(args.outputfolder_s3_path, suffix='.vcf')

    genotypingJobs = []
    for sample in samples:
        sampleName = sample['name']
        vcfFile = "%s/%s.vcf" % (args.outputfolder_s3_path, sampleName)
        if vcfFile not in vcfFiles:
            if not exists(sample['bam']):
                logging.error("%s does not exist.", sample['bam'])
                continue
            logging.info("Genotyping %s", sampleName)
            if sample['reference'] == 'hg19':
                reference_s3_path = "s3://bmsrd-ngs-repo/cohort-matcher/hg19.tar.bz2"
                targets_s3_path = "s3://bmsrd-ngs-repo/cohort-matcher/hg19.cohort-matcher.bed"
            elif sample['reference'] == 'GRCh37ERCC':
                reference_s3_path = "s3://bmsrd-ngs-repo/cohort-matcher/GRCh37ERCC.tar.bz2"
                targets_s3_path = "s3://bmsrd-ngs-repo/cohort-matcher/GRCh37ERCC.cohort-matcher.bed"
            elif sample['reference'] == 'hg38':
                reference_s3_path = "s3://bmsrd-ngs-repo/cohort-matcher/hg38.tar.bz2"
                targets_s3_path = "s3://bmsrd-ngs-repo/NGSCheckMate/SNP_GRCh38_hg38_wChr.sorted.bed"
            elif sample['reference'] == 'GRCh38ERCC':
                reference_s3_path = "s3://bmsrd-ngs-repo/cohort-matcher/GRCh38ERCC.tar.bz2"
                targets_s3_path = "s3://bmsrd-ngs-repo/NGSCheckMate/SNP_GRCh38_hg38_woChr.sorted.bed"
            else:
                logging.warn("Skipping %s due to unknown reference", sampleName)
                continue

            if args.dryRun:
                logging.info("Would call batch.submit_job")
            else:
                response = batch.submit_job(jobName='freebayes-%s' % sampleName,
                                            jobQueue=args.job_queue,
                                            jobDefinition=args.job_definition,
                                            retryStrategy={
                                                'attempts': 2
                                            },
                                            containerOverrides={
                                                'memory': args.memory,
                                                'command': ["--sample_name", sampleName,
                                                            "--bamfile_s3_path", sample['bam'],
                                                            "--targets_s3_path", targets_s3_path,
                                                            "--reference_s3_path", reference_s3_path,
                                                            "--s3_output_folder_path", args.outputfolder_s3_path]
                                            })
                logging.debug(response)
                jobId = response['jobId']
                genotypingJobs.append(jobId)

    logging.info("Submitted %s jobs", len(genotypingJobs))
    completed_jobs = []
    failed_jobs = []
    for counter, jobid in enumerate(genotypingJobs):
        status = ''
        while status != 'SUCCEEDED' and status != 'FAILED':
            logging.info("[%d/%d] Checking job %s", counter+1, len(genotypingJobs), jobid)
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

def parseArguments(argv):
    ''' Parse Arguments '''
    parser = argparse.ArgumentParser(description='Genotype a set of samples')
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument('-d', '--dryRun', default=False, action="store_true",
                        help="Simulate job submission")

    required_args = parser.add_argument_group("Required")
    required_args.add_argument('-b', '--bamsheet', required=True, help="Bamsheet")
    required_args.add_argument('-o', '--outputfolder_s3_path', required=True,
                               help="Specify S3 path for cached VCF/TSV files")

    #parser.add_argument('-r', '--reference_s3_path', required=True,
    #                    help='S3 Path to Reference Fastq File (bzip2)')
    #parser.add_argument('-t', '--targets_s3_path', required=True,
    #                    help="S3 Path to Targets BED File")

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
