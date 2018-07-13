#!/usr/bin/env python
'''
matchSamples - Check samples pairs

1. Genotype all WES Tumor, Normal, RNA-Seq samples
2. Perform WES Tumor x WES Normal, RNA x WES Tumor or Normal.
3. Merge meltedResults
	a. head -1 normal_rnaseq.meltedResults.txt > meltedResults.txt
	b. tail -n +2 *.meltedResults.txt >> meltedResults.txt
4. Run matchSamples.py
'''
import argparse
import logging
import os
import sys
from pandas import read_csv
from common import downloadFile

__appname__ = "matchSamples"
logger = logging.getLogger(__appname__)
__version__ = "1.0"

def readPatientToSample(patientToSampleFile):
    '''
    Read a patient to sample mapping file.
    Each line has at least two columns, where the first two columns are patient_id and sample_id

    :param patientToSampleFile: tsv file
    :return dict of patient ids to sample id lists
    '''
    patients = dict()
    with open(patientToSampleFile, 'r') as f:
        for line in f:
            line = line.strip("\n")
            fields = line.split("\t")
            patientid = fields[0]
            sampleid = fields[1]
            if patientid in patients:
                patients[patientid].append(sampleid)
            else:
                patients[patientid] = [sampleid]
    logger.debug("Read %d patients", len(patients))
    return patients

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logging.getLogger("botocore").setLevel(logging.WARNING)
    logger.info("%s v%s" % (__appname__, __version__))
    logger.info(args)

    # Read in meltedResults
    patients = readPatientToSample(args.patientToSample)

    # Match samples to patients
    logger.info("Matching samples to patients")
    matchSamplesToPatients(patients, args.meltedResults_s3path, args.meltedResults_localpath,
                           args.threshold, args.resultsFile)
    logger.info("Done.")

def matchSamplesToPatients(patients, meltedResults_s3path, meltedResults_localpath, threshold, resultsFile):
    '''
    Match samples to patient, find swaps, etc.
    For each patient that has more than 1 sample
      i. For each sample
        1) Read <sample>.meltedResults.txt
        2) Check best matching sample (BMS)
        3) If BMS is from the same patient, all is okay, else report BMS.

    :param patients: list of patient -> samples
    :param meltedResults_s3path: S3 Path to meltedResults files (mutually exclusive with meltedResults_localpath)
    :param meltedResults_localpath: Local path to meltedResults files (m.e with meltedResults_s3path)
    :param threshold: Minimum # of SNPs compared against sample
    :param resultsFile: Output results files
    '''
    totalSamples = 0
    missingSamples = []
    badSamples = []
    mismatchedSamples = []
    matchedSamples = 0
    with open(resultsFile, 'w') as output:
      output.write("Sample1\tSample2\tn_S1\tn_S2\tSNPs_Compared\tFraction_Match\n")
      for patient in patients:
        samples = patients[patient]
        logger.debug("Patient: %s, Samples: %s", patient, samples)
        if len(samples) > 1:
            for sample in samples:
                # Get the sample meltedResults
                if meltedResults_s3path:
                    meltedResultsFile = "%s.meltedResults.txt" % sample
                    downloadFile("%s/%s.meltedResults.txt" % (meltedResults_s3path, sample), meltedResultsFile)
                    if not os.path.exists(meltedResultsFile):
                        logger.error("Unable to download %s", meltedResultsFile)
                        missingSamples.append(sample)
                        continue
                else:
                    meltedResultsFile = "%s/%s.meltedResults.txt" % (meltedResults_localpath, sample)
                    if not os.path.exists(meltedResultsFile):
                        logger.error("%s not found", meltedResultsFile)
                        missingSamples.append(sample)
                        continue
                cm = read_csv(meltedResultsFile, sep="\t")
                if meltedResults_s3path:
                    os.remove(meltedResultsFile)
                # Get the samples above threshold
                cm = cm[cm.SNPs_Compared >= threshold]
                snpsCompared = len(cm)
                if snpsCompared == 0:
                    logger.warn("%s -> Not enough SNPs", sample)
                    missingSamples.append(sample)
                else:
                    totalSamples += 1
                    samples_to_print = 1
                    cm = cm.sort_values(by="Fraction_Match", ascending=False)
                    topMatch = cm.iloc[0]['Sample2']
                    fractionMatch = cm.iloc[0]['Fraction_Match']
                    if fractionMatch < 0.7:
                        badSamples.append(sample)
                    else:
                        if topMatch in samples:
                            matchedSamples += 1
                            logger.debug("%s -> %s", sample, topMatch)
                        else:
                            logger.warn("%s -> %s:%s", sample, patient, topMatch)
                            mismatchedSamples.append(sample)
                            samples_to_print = 5
                    for i in range(0, samples_to_print):
                         output.write("%s\t%s\t%d\t%d\t%d\t%0.4f\n" % (cm.iloc[i]['Sample1'],
                                                                       cm.iloc[i]['Sample2'],
                                                                       cm.iloc[i]['n_S1'],
                                                                       cm.iloc[i]['n_S2'],
                                                                       cm.iloc[i]['SNPs_Compared'],
                                                                       cm.iloc[i]['Fraction_Match']))

    logger.info("Samples compared: %d", totalSamples)
    logger.info("Matched samples: %d", matchedSamples)
    logger.info("Mismatched samples: %d", len(mismatchedSamples))
    logger.info("Missing samples: %d (no meltedResults files or low coverage)", len(missingSamples))
    logger.info("No matching samples: %d (Fraction_Match < 0.7)", len(badSamples))

def parseArguments(argv):
    ''' Parse Arguments '''
    parser = argparse.ArgumentParser(description="Find swaps in cohort-matcher results")
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    parser_grp1 = parser.add_argument_group("Required")

    parser_mr = parser.add_mutually_exclusive_group(required=True)
    parser_mr.add_argument("--meltedResults_s3path", help="Melted results files from cohort-matcher")
    parser_mr.add_argument("--meltedResults_localpath", help="Melted results files from cohort-matcher")
    
    parser_grp1.add_argument("--patientToSample", required=True, help="Patient to Sample Mapping file")
    parser_grp1.add_argument("--threshold", default=15, help="Min # of SNPS compared between samples")
    parser_grp1.add_argument("--resultsFile", default="cohort-matcher_results.txt", help="Results of mismatched samples")
    args = parser.parse_args(argv)
    return args

if __name__ == '__main__':
    main(sys.argv[1:])
