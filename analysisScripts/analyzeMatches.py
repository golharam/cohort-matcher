#!/usr/bin/env python
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

def readCohortMatcherMeltedResults(meltedResultsFile):
    '''
    Read the melted results file and generate a dict
    of subject to matched samples
    :param meltedResultsFile: Cohort-matcher meltedResults.txt
    :param subjectToSample: subject to sample map
    :return: subject to sample map based on meltedResults.txt
    '''
    logging.info("Reading %s", meltedResultsFile)
    with open(meltedResultsFile, 'r') as f:
        f.readline()    # skip header line
        for line in f:
            line = line.rstrip('\r\n')
            s1, s2, ns1, ns2, snps_cmp, fmatch, judgement = line.split('\t')
            if 'DIFFERENT' in judgement:
                continue
            if 'INCONCLUSIVE' in judgement:
                continue
            subject1 = SAMPLE_TO_SUBJECT[s1]
            subject2 = SAMPLE_TO_SUBJECT[s2]
            if subject1 not in CM_SUBJECT_TO_SAMPLE:
                CM_SUBJECT_TO_SAMPLE[subject1] = []
            CM_SUBJECT_TO_SAMPLE[subject1].append(s1)
            if subject2 not in CM_SUBJECT_TO_SAMPLE:
                CM_SUBJECT_TO_SAMPLE[subject2] = []
            CM_SUBJECT_TO_SAMPLE[subject2].append(s2)

def readNGSCheckMateOutput(ngscheckmate_matchedoutput):
    logging.info("Reading %s", ngscheckmate_matchedoutput)
    if os.path.exists(ngscheckmate_matchedoutput) is False:
        logging.warn("%s not found.  Skipping", ngscheckmate_matchedoutput)
        return
    with open(ngscheckmate_matchedoutput, 'r') as f:
        f.readline()
        for line in f:
            line = line.rstrip('\r\n')
            s1, status, s2, corr, depth = line.split('\t')
            # s1 and s2 are the filename, so strip everything past the period
            s1 = s1[0:s1.find('.')]
            s2 = s2[0:s2.find('.')]
            subject1 = SAMPLE_TO_SUBJECT[s1]
            subject2 = SAMPLE_TO_SUBJECT[s2]
            if subject1 not in NGSCHECKMATE_SUBJECT_TO_SAMPLE:
                NGSCHECKMATE_SUBJECT_TO_SAMPLE[subject1] = []
            NGSCHECKMATE_SUBJECT_TO_SAMPLE[subject1].append(s1)
            if subject2 not in NGSCHECKMATE_SUBJECT_TO_SAMPLE:
                NGSCHECKMATE_SUBJECT_TO_SAMPLE[subject2] = []
            NGSCHECKMATE_SUBJECT_TO_SAMPLE[subject2].append(s2)

def readSampleToSubject(sampleToSubjectFile):
    '''
    Read the sample to subject file and generate a dict
    of sample->subject and subject->sample
    :param sampleToSubjectFile: 
    '''
    logging.info("Reading %s", sampleToSubjectFile)
    with open(sampleToSubjectFile, 'r') as f:
        f.readline()    # skip header line
        for line in f:
            line = line.rstrip('\r\n')
            sample, subject = line.split('\t')
            SAMPLE_TO_SUBJECT[sample] = subject
            if subject not in SUBJECT_TO_SAMPLE:
                SUBJECT_TO_SAMPLE[subject] = []
            SUBJECT_TO_SAMPLE[subject].append(sample)
    logging.debug("Subject to Sample Map: %s", SUBJECT_TO_SAMPLE)

def writeSubjectToSample(mapDescription, subjectToSampleMap):
    txt = "%s.txt" % mapDescription
    logging.info("Writing %s", txt)
    with open(txt, 'w') as f:
        for subject in sorted(subjectToSampleMap.keys()):
            f.write("%s: %s\n" % (subject, subjectToSampleMap[subject]))

def doWork(args):
    # 1. Read the sample to subject mapping and build a graph of
    # subject -> [samples]
    readSampleToSubject(args.sample_to_subject)

    # 2. Read meltedResults and build/write similar graph
    readCohortMatcherMeltedResults(args.melted_results)

    # 2b. Read NGSCheckMate results
    readNGSCheckMateOutput(args.ngscheckmate_matchedoutput)

    # 3. Write maps
    writeSubjectToSample('expected', SUBJECT_TO_SAMPLE)
    writeSubjectToSample('cohort-matcher-actual', CM_SUBJECT_TO_SAMPLE)
    writeSubjectToSample('NGSCheckMate-actual', NGSCHECKMATE_SUBJECT_TO_SAMPLE)

def main(argv):
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logging.info("%s v%s", __appname__, __version__)
    logging.info(args)
    doWork(args)

def parseArguments(argv):
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
