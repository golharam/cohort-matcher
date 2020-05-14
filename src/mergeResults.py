#!/usr/bin/env python
'''
Script to merge results of compareSamples.py/compareGenotypes.py.
'''
import argparse
import logging
import sys
from common import downloadFile, readSamples, uploadFile, generate_working_dir

__appname__ = 'mergeResults'
__version__ = "0.2"

logging.getLogger(__appname__)

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logging.info("%s v%s", __appname__, __version__)
    logging.info(args)

    samples = readSamples(args.bamsheet)
    # We don't need the last sample in the list so let's remove it
    samples.pop()

    working_dir = generate_working_dir(args.working_dir)
    logging.info("Working in %s", working_dir)

    localMeltedResultsFiles = []
    for i, sample in enumerate(samples):
        meltedResultFile = "%s/%s.meltedResults.txt" % (args.s3_cache_folder, sample['name'])
        f = "%s/%s.meltedResults.txt" % (working_dir, sample['name'])
        logging.info("[%d/%d] Downloading %s -> %s", i+1, len(samples), meltedResultFile, f)
        downloadFile(meltedResultFile, f)
        localMeltedResultsFiles.append(f)

    meltedResultsFile = "meltedResults.txt"
    with open(meltedResultsFile, 'w') as outfile:
        header_written = False
        for i, fname in enumerate(localMeltedResultsFiles):
            logging.info("[%d/%d] Merging %s -> %s", i+1, len(localMeltedResultsFiles), fname,
                         meltedResultsFile)
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

    s3_path = "%s/meltedResults.txt" % args.s3_cache_folder
    logging.info("Uploading %s -> %s", meltedResultsFile, s3_path)
    uploadFile(meltedResultsFile, s3_path)
    logging.info("Done.")

def parseArguments(argv):
    ''' Parse arguments '''
    parser = argparse.ArgumentParser(description='Merge set of meltedResults')
    parser.add_argument('-l', '--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument('-d', '--dry-run', default=False, action="store_true",
                        help="Simulates everything, except for actually submitting a job")

    required_args = parser.add_argument_group("Required")
    required_args.add_argument('-b', '--bamsheet', required=True, help="Bamsheet")
    required_args.add_argument("-CD", "--s3_cache_folder", required=True,
                               help="Specify S3 path for cached meltedResults files")

    parser.add_argument("-w", "--working_dir", required=False, default="/scratch",
                        help="Local working directory")

    args = parser.parse_args(argv)
    return args

if __name__ == '__main__':
    main(sys.argv[1:])
