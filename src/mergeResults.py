#!/usr/bin/env python
'''
Script to merge results of compareSamples.py/compareGenotypes.py.
'''
import argparse
import logging
import sys
from genotypeSamples import read_samples
from common import downloadFile, uploadFile, generate_working_dir

__appname__ = 'mergeResults'
__version__ = "0.2"

logging.getLogger(__appname__)

def download_melted_results(samples, working_dir, s3_cache_folder):
    localMeltedResultsFiles = []
    for i, sample in enumerate(samples):
        meltedResultFile = "%s/%s.meltedResults.txt" % (s3_cache_folder, sample['sample_id'])
        f = "%s/%s.meltedResults.txt" % (working_dir, sample['sample_id'])
        logging.info("[%d/%d] Downloading %s -> %s", i+1, len(samples), meltedResultFile, f)
        downloadFile(meltedResultFile, f)
        localMeltedResultsFiles.append(f)
    return localMeltedResultsFiles

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logging.info("%s v%s", __appname__, __version__)
    logging.info(args)

    samples = read_samples(args.bamsheet)
    
    # Create a sample-to-subject map
    sample_to_subject = {}
    for sample in samples:
        sample_to_subject[sample['sample_id']] = sample['subject_id']

    # We don't need the last sample in the list so let's remove it
    samples.pop()

    working_dir = generate_working_dir(args.working_dir)
    logging.info("Working in %s", working_dir)

    localMeltedResultsFiles = download_melted_results(samples, working_dir, args.s3_cache_folder)

    meltedResultsFile = "meltedResults.txt"
    with open(meltedResultsFile, 'w') as outfile:
        # Write the header
        outfile.write('\t'.join(['Subject1', 'Sample1', 'Subject2', 'Sample2', 'n_S1', 'n_S2', 'SNPs_Compared', 'Fraction_Match', 'Binomial_PV', 'LRR', 'Fraction_Match_Plus', 'PV', 'Judgement', 'Swap']) + "\n")
        
        for i, fname in enumerate(localMeltedResultsFiles):
            logging.info("[%d/%d] Merging %s -> %s", i+1, len(localMeltedResultsFiles), fname, meltedResultsFile)
            with open(fname) as infile:
                # Skip the header line
                next(infile)
                for line in infile:
                    sample1, sample2, n_s1, n_s2, snps_compared, fraction_match, binomial_pv, lrr, fraction_match_plus, pv, judgement = line.strip().split('\t')
                    subject1 = sample_to_subject[sample1]
                    subject2 = sample_to_subject[sample2]
                    swap = 'x' if subject1 != subject2 and 'SAME' in judgement else ''
                    outfile.write('\t'.join([subject1, sample1, subject2, sample2, n_s1, n_s2, snps_compared, fraction_match, binomial_pv, lrr, fraction_match_plus, pv, judgement, swap]) + "\n")

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
