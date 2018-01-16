#!/usr/bin/env python
from __future__ import print_function
''' This is a wrapper script for Docker '''

import logging
import multiprocessing
import os
import shlex
import subprocess
from argparse import ArgumentParser

from common_utils.s3_utils import download_file, upload_file, download_folder, upload_folder
from common_utils.job_utils import generate_working_dir, delete_working_dir, uncompress

__version__ = "1.0.3-alpha"
logger = logging.getLogger(__name__)

def download_reference(s3_path, working_dir):
    """
    Downloads reference folder that has been configured to run with Isaac
    :param s3_path: S3 path that the folder resides in
    :param working_dir: working directory
    :return: local path to the folder containing the reference
    """

    reference_folder = os.path.join(working_dir, 'reference')

    try:
        os.mkdir(reference_folder)
    except Exception as e:
        pass

    download_folder(s3_path, reference_folder)

    # Update sorted reference
    update_sorted_reference(reference_folder)

    return reference_folder


def download_fastq_files(fastq1_s3_path, fastq2_s3_path, working_dir):
    """
    Downlodas the fastq files
    :param fastq1_s3_path: S3 path containing FASTQ with read1
    :param fastq2_s3_path: S3 path containing FASTQ with read2
    :param working_dir: working directory
    :return: local path to the folder containing the fastq
    """
    fastq_folder = os.path.join(working_dir, 'fastq')

    try:
        os.mkdir(fastq_folder)
    except Exception as e:
        pass

    local_fastq1_path = download_file(fastq1_s3_path, fastq_folder)
    local_fastq2_path = download_file(fastq2_s3_path, fastq_folder)

    # Isaac requires the fastqs to be symlinked as lane1_read1.fastq.gz and lane1_read2.fastq.gz
    os.symlink(local_fastq1_path, os.path.join(fastq_folder, 'lane1_read1.fastq.gz'))
    os.symlink(local_fastq2_path, os.path.join(fastq_folder, 'lane1_read2.fastq.gz'))

    return fastq_folder


def upload_bam(bam_s3_path, local_folder_path):
    """
    Uploads results folder containing the bam file (and associated output)
    :param bam_s3_path: S3 path to upload the alignment results to
    :param local_folder_path: local path containing the alignment results
    """

    upload_folder(bam_s3_path, local_folder_path)

def run_cohort_matcher(log_level, bam_sheet1, bam_sheet2, reference1, reference2, working_dir,
    output_prefix, max_jobs):
    """
    Runs Cohort-matcher
    :param bam_sheet1: local path to bam sheet1
    :param bam_sheet2: local path to bam sheet2
    :param reference1: hg19 or GRCh37
    :param reference2: hg19 or GRCh37
    :param working_dir: working directory
    :param output_prefix: output prefix
    :return: path to results
    """

    os.chdir(working_dir)
    cache_dir = os.path.join(working_dir, 'cache')
    try:
        os.mkdir(cache_dir)
    except Exception as e:
        pass
   
    ref_dir = working_dir + '/reference'
    if reference1 == 'hg19':
        ref = os.path.join(ref_dir, 'hg19.fa')
        vcf = os.path.join(ref_dir, 'hg19.exome.highAF.7550.vcf')
        if reference2 == 'hg19':
            ref2 = ''
        else:
            ref2 = '-R2 %s -V2 %s -CM %s' % (os.path.join(ref_dir, 'GRCh37.fa'),
                                             os.path.join(ref_dir, '1kg.exome.highAF.7550.vcf'),
                                             os.path.join(ref_dir, 'hg19.chromosome_map'))
    else:
        ref = os.path.join(ref_dir, 'GRCh37.fa')
        vcf = os.path.join(ref_dir, '1kg.exome.highAF.7550.vcf')
        if reference2 == 'GRCh37':
            ref2 = ''
        else:
            ref2 = '-R2 %s -V2 %s -CM %s' % (os.path.join(ref_dir, 'hg19.fa'),
                                             os.path.join(ref_dir, 'hg19.exome.highAF.7550.vcf'),
                                             os.path.join(ref_dir, 'hg19.chromosome_map'))

    cmd = '/cohort-matcher/cohort_matcher.py --log-level %s --set1 %s --set2 %s --cache-dir %s ' \
          '--scratch-dir %s --caller freebayes --max-jobs %d -R %s -V %s %s ' \
          '--freebayes-path /usr/local/bin/freebayes --aws /usr/bin/aws ' \
          '--Rscript /usr/bin/Rscript --samtools /usr/local/bin/samtools -O %s --output-prefix %s' % \
          (log_level, bam_sheet1, bam_sheet2, cache_dir, working_dir, max_jobs, ref, vcf, ref2,
           working_dir, output_prefix)
    logger.info("Running: %s", cmd)
    subprocess.check_call(shlex.split(cmd))
    return working_dir

def parseArguments():
    argparser = ArgumentParser()
    argparser.add_argument('-l', '--log-level', help="Prints warnings to console by default",
                           default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    file_path_group = argparser.add_argument_group(title='File paths')
    file_path_group.add_argument('--set1_s3_path', type=str, 
                                 help='S3 path to first set of samples', required=True)
    file_path_group.add_argument('--set2_s3_path', type=str, 
                                 help='S3 path to second set of samples', required=True)
    file_path_group.add_argument('--set1_reference', type=str, choices=['hg19', 'GRCh37'],
                                 help='Reference genome for set1', required=True)
    file_path_group.add_argument('--set2_reference', type=str, choices=['hg19', 'GRCh37'],
                                 help='Reference genome for set2', required=True)
    file_path_group.add_argument('--s3_output_folder_path', type=str, 
                                 help='S3 path for output files', required=True)

    run_group = argparser.add_argument_group(title='Run command args')
    run_group.add_argument('--output_prefix', type=str, default='cohort-matcher-results',
                           help='Output prefix')

    argparser.add_argument('--working_dir', type=str, default='/scratch', 
                           help="Working directory (default: /scratch)")
    argparser.add_argument('--max_jobs', type=int, default=None,
                           help="Maximum # of parallel genotyping jobs (default: all cores")

    return argparser.parse_args()

def main():
    args = parseArguments()
    logging.basicConfig(level=args.log_level)
    logger.info("Run cohort-matcher Docker CLI v%s", __version__)
    logger.info(args)

    working_dir = generate_working_dir(args.working_dir)
    logger.info("Working in %s", working_dir)

    # Download fastq files and reference files
    logger.info('Downloading bam sheets')
    set1_bamsheet = download_file(args.set1_s3_path, working_dir)
    set2_bamsheet = download_file(args.set2_s3_path, working_dir)

    # Download reference bundles
    refdir = working_dir + '/reference'
    os.mkdir(refdir)
    if args.set1_reference == 'hg19' or args.set2_reference == 'hg19':
        logger.info("Downloading hg19 reference bundle")
        download_file('s3://bmsrd-ngs-repo/reference/hg19-cohort-matcher.tar.bz2', refdir)
        logger.info("Uncompressing hg19 reference bundle")
        uncompress(os.path.join(refdir, 'hg19-cohort-matcher.tar.bz2'), refdir)
        os.remove(os.path.join(refdir, 'hg19-cohort-matcher.tar.bz2'))

    if args.set2_reference == 'GRCh37' or args.set2_reference == 'GRCh37':
        logger.info("Downloading GRCh37 reference bundle")
        download_file('s3://bmsrd-ngs-repo/reference/GRCh37-cohort-matcher.tar.bz2', refdir)
        logger.info("Uncompressing GRCh37 reference bundle")
        uncompress(os.path.join(refdir, 'GRCh37-cohort-matcher.tar.bz2'), refdir)
        os.remove(os.path.join(refdir, 'GRCh37-cohort-matcher.tar.bz2'))

    # Run cohort-matcher
    if args.max_jobs is None:
        max_jobs = multiprocessing.cpu_count()
    else:
        max_jobs = args.max_jobs
    logger.info('Running cohort-matcher using %d threads.', max_jobs)
    output_folder_path = run_cohort_matcher(args.log_level, set1_bamsheet, set2_bamsheet,
                                            args.set1_reference, args.set2_reference,
                                            working_dir, args.output_prefix, max_jobs)
    delete_working_dir(refdir)
    logger.info('Uploading results to %s', args.s3_output_folder_path)
    upload_folder(args.s3_output_folder_path, working_dir)
    logger.info('Cleaning up working dir')
    delete_working_dir(working_dir)
    logger.info('Completed')

if __name__ == '__main__':
    main()
