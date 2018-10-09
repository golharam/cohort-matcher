#!/usr/bin/env python
'''
Script to merge results of compareSamples.py/compareGenotypes.py.
'''
import argparse
import boto3
import logging
import os
import sys

from common import find_bucket_key, listFiles, downloadFile

__appname__ = 'mergeResults'
__version__ = "0.1"

logger = logging.getLogger(__appname__)

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logger.info("%s v%s" % (__appname__, __version__))
    logger.info(args)

    # Get a list of meltedResults files
    meltedResultsFiles = listFiles(args.s3_cache_folder, suffix=".meltedResults.txt")
    logger.info("Found %s melted result files", len(meltedResultsFiles))

    localFiles = []
    for i, meltedResultFile in enumerate(meltedResultsFiles):
        f = "%s/%s" % (args.working_dir, os.path.basename(meltedResultFile))
        if os.path.exists(f) is False:
            logger.info("[%d/%d] Downloading %s -> %s", i+1, len(meltedResultsFiles), meltedResultFile, f)
            downloadFile(meltedResultFile, f)
        localFiles.append(f)

    meltedResultsFile = "%s/meltedResults.txt" % args.working_dir
    with open(meltedResultsFile, 'w') as outfile:
        header_written = False
        for i, fname in enumerate(localFiles):
            logger.info("[%d/%d] Merging %s -> %s", i+1, len(localFiles), fname, meltedResultsFile)
            with open(fname) as infile:
                if header_written is False:
                    # Write complete file out
                    outfile.write(infile.read())
                    header_written = True
                else:
                    # Skip first line and write out rest
                    next(infile)
                    for line in infile:
                        outfile.write(line)

def parseArguments(argv):
    ''' Parse arguments '''
    parser = argparse.ArgumentParser(description='Compare a sample to a set of samples')
    parser.add_argument('-l', '--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument('-d', '--dry-run', default=False, action="store_true",
                        help="Simulates everything, except for actually submitting a job")
    parser.add_argument("-CD", "--s3_cache_folder", required=True, 
                        help="Specify S3 path for cached VCF/TSV files")
    parser.add_argument("-w", "--working_dir", required=False, default="/scratch",
                        help="Local working directory")

    args = parser.parse_args(argv)
    return args

if __name__ == '__main__':
    main(sys.argv[1:])

