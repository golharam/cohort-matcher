#!/usr/bin/env python
'''
findSwaps - Check samples pairs

1. Genotype all WES Tumor, Normal, RNA-Seq samples
2. Perform WES Tumor x WES Normal, RNA x WES Tumor or Normal.
3. Merge meltedResults
	a. head -1 normal_rnaseq.meltedResults.txt > meltedResults.txt
	b. tail -n +2 *.meltedResults.txt >> meltedResults.txt
4. Run findSwaps.py
'''
import argparse
import logging
import sys
from pandas import read_csv

logger = logging.getLogger("checkSamples")
__version__ = "1.0"

def readPatientToSample(patientToSampleFile):
    with open(patientToSampleFile, 'r') as f:
        patients = dict()
        for line in f:
            line = line.strip("\n")
            fields = line.split("\t")
            patientid = fields.pop(0)
            patients[patientid] = fields
    logger.info("Read %d patients", len(patients))
    return patients

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logger.info("findSwaps v%s" % __version__)
    logger.info(args)
    # Read in files
    cm = read_csv(args.meltedResults, sep="\t")
    logger.info("Read %d lines", len(cm))
    cm = cm[cm.SNPs_Compared >= args.threshold]
    logger.info("%d lines after filtering", len(cm))
    patients = readPatientToSample(args.patientToSample)
    # Scan each patient
    # For each sample, make sure the top match is another sample from the same patient
    badSamples = []
    for patient in patients:
        logger.info("%s", patient)
        samples = patients[patient]
        if len(samples) > 1:
            for sample in samples:
                logger.info("\t\t%s", sample)
                sample_matches = cm[cm.Sample1==sample]
                if len(sample_matches) > 0:
                    otherSample = 'Sample2'
                else:
                    sample_matches = cm[cm.Sample2==sample]
                    otherSample = 'Sample1'
                    if len(sample_matches) == 0:
                        logger.warn("Unable to find sample %s in cohort-matcher meltedResults. " \
                                     "It could have been filtered out due to low coverage.", sample)
                        badSamples.append(sample)
                        
                if len(sample_matches) > 0:
                    sample_matches = sample_matches.sort_values(by="Fraction_Match", ascending=False)
                    topMatch = sample_matches.iloc[0][otherSample]
                    if topMatch not in samples:
                        logger.warn("Sample %s does not match to samples from same patient (%s).  It matches to %s", sample, patient, topMatch)
    logger.warn("List of non-matched samples: %s", badSamples)
    
def parseArguments(argv):
    ''' Parse Arguments '''
    parser = argparse.ArgumentParser(description="Find swaps in cohort-matcher results")
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    parser_grp1 = parser.add_argument_group("Required")
    parser_grp1.add_argument("--meltedResults", required=True, help="Melted results file from cohort-matcher")
    parser_grp1.add_argument("--patientToSample", required=True, help="Patient to Sample Mapping file")
    parser_grp1.add_argument("--threshold", required=False, default=15, help="Min # of SNPS compared between samples")
    args = parser.parse_args(argv)
    return args

if __name__ == '__main__':
    main(sys.argv[1:])
