#!/usr/bin/env python
'''
This script parses Cohort-Matcher meltedResults output and compares
to the subject to sample map.
'''
import argparse
import logging
import os
import sys

__appname__ = os.path.basename(__file__).split('.py')[0]
__version__ = '0.01'

SAMPLE_TO_SUBJECT = dict()
SUBJECT_TO_SAMPLE = dict()
CM_SUBJECT_TO_SAMPLE = dict()
NGSCHECKMATE_SUBJECT_TO_SAMPLE = dict()

def read_cohortmatcher_results(melted_results_file):
    '''
    Read the melted results file and generate a dict
    of subject to matched samples
    :param meltedResultsFile: Cohort-matcher meltedResults.txt
    :param subjectToSample: subject to sample map
    :return: subject to sample map based on meltedResults.txt
    '''
    logging.info("Reading %s", melted_results_file)
    with open(melted_results_file, 'r') as fin:
        fin.readline()    # skip header line
        for line in fin:
            line = line.rstrip('\r\n')
            sample1, sample2, ns1, ns2, snps, f_match, judgement = line.split('\t')
            if 'DIFFERENT' in judgement:
                continue
            if 'INCONCLUSIVE' in judgement:
                continue

            subject1 = SAMPLE_TO_SUBJECT[sample1]
            subject2 = SAMPLE_TO_SUBJECT[sample2]

            if subject1 not in CM_SUBJECT_TO_SAMPLE:
                CM_SUBJECT_TO_SAMPLE[subject1] = []
            if subject2 not in CM_SUBJECT_TO_SAMPLE:
                CM_SUBJECT_TO_SAMPLE[subject2] = []

            if sample1 not in CM_SUBJECT_TO_SAMPLE[subject1]:
                CM_SUBJECT_TO_SAMPLE[subject1].append(sample1)
            if sample2 not in CM_SUBJECT_TO_SAMPLE[subject2]:
                CM_SUBJECT_TO_SAMPLE[subject2].append(sample2)

def read_ngscheckmate_output(ngscheckmate_matchedoutput):
    ''' Read NGS CheckMate output file '''
    if os.path.exists(ngscheckmate_matchedoutput) is False:
        logging.warn("%s not found.  Skipping", ngscheckmate_matchedoutput)
        return
    logging.info("Reading %s", ngscheckmate_matchedoutput)
    with open(ngscheckmate_matchedoutput, 'r') as fin:
        fin.readline()
        for line in fin:
            line = line.rstrip('\r\n')
            sample1, _, sample2, _, _ = line.split('\t')
            # s1 and s2 are the filename, so strip everything past the period
            sample1 = sample1[0:sample1.find('.')]
            sample2 = sample2[0:sample2.find('.')]

            subject1 = SAMPLE_TO_SUBJECT[sample1]
            subject2 = SAMPLE_TO_SUBJECT[sample2]

            if subject1 not in NGSCHECKMATE_SUBJECT_TO_SAMPLE:
                NGSCHECKMATE_SUBJECT_TO_SAMPLE[subject1] = []
            if subject2 not in NGSCHECKMATE_SUBJECT_TO_SAMPLE:
                NGSCHECKMATE_SUBJECT_TO_SAMPLE[subject2] = []

            if sample1 not in NGSCHECKMATE_SUBJECT_TO_SAMPLE[subject1]:
                NGSCHECKMATE_SUBJECT_TO_SAMPLE[subject1].append(sample1)
            if sample2 not in NGSCHECKMATE_SUBJECT_TO_SAMPLE[subject2]:
                NGSCHECKMATE_SUBJECT_TO_SAMPLE[subject2].append(sample2)

def read_sample_to_subject(sample_to_subject_file):
    '''
    Read the sample to subject file and generate a dict
    of sample->subject and subject->sample
    :param sampleToSubjectFile: sample to subject file
    '''
    logging.info("Reading %s", sample_to_subject_file)
    with open(sample_to_subject_file, 'r') as fin:
        fin.readline()    # skip header line
        for line in fin:
            line = line.rstrip('\r\n')
            sample, subject = line.split('\t')
            SAMPLE_TO_SUBJECT[sample] = subject
            if subject not in SUBJECT_TO_SAMPLE:
                SUBJECT_TO_SAMPLE[subject] = []
            if sample in SUBJECT_TO_SAMPLE[subject]:
                logging.error("%s already in map", sample)
            SUBJECT_TO_SAMPLE[subject].append(sample)
    logging.debug("Subject to Sample Map: %s", SUBJECT_TO_SAMPLE)

def write_subject_to_sample(map_description, subject_to_sample_map):
    ''' Write subject to sample map '''
    if not subject_to_sample_map:
        return
    txt = "%s.txt" % map_description
    logging.info("Writing %s", txt)
    with open(txt, 'w') as fout:
        for subject in sorted(subject_to_sample_map.keys()):
            fout.write("%s: %s\n" % (subject, sorted(subject_to_sample_map[subject])))

def do_work(args):
    ''' Main work function '''
    # 1. Read the sample to subject mapping and build a graph of
    # subject -> [samples]
    read_sample_to_subject(args.sample_to_subject)

    # 2. Read meltedResults and build/write similar graph
    read_cohortmatcher_results(args.melted_results)
    # Samples that do not match to anything are not reported, hence not included.

    # 2b. Read NGSCheckMate results
    read_ngscheckmate_output(args.ngscheckmate_matchedoutput)

    # 3. Write maps
    write_subject_to_sample('expected', SUBJECT_TO_SAMPLE)
    write_subject_to_sample('cohort-matcher-actual', CM_SUBJECT_TO_SAMPLE)
    write_subject_to_sample('NGSCheckMate-actual', NGSCHECKMATE_SUBJECT_TO_SAMPLE)

def main(argv):
    ''' Main Entry Point '''
    args = parse_arguments(argv)
    logging.basicConfig(level=args.log_level)
    logging.info("%s v%s", __appname__, __version__)
    logging.info(args)
    do_work(args)

def parse_arguments(argv):
    ''' Parse command-line arguments '''
    parser = argparse.ArgumentParser()

    parser.add_argument('-l', '--log-level',
                        help="Prints warnings to console by default",
                        default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument('-v', '--version', action="version",
                        version=__version__)

    parser.add_argument('-m', '--melted-results',
                        default='meltedResults.txt',
                        help="Melted Results")
    parser.add_argument('-n', '--ngscheckmate-matchedoutput',
                        default="output_matched.txt",
                        help="NGS CheckMate Matched Output")
    parser.add_argument('-s', '--sample-to-subject',
                        default='sampleToSubject.txt',
                        help="Sample to Subject Map")
    return parser.parse_args(argv)

if __name__ == "__main__":
    main(sys.argv[1:])
