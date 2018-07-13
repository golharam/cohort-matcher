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
    logger.info("Read %d patients", len(patients))
    return patients

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logging.getLogger("botocore").setLevel(logging.WARNING)
    logger.info("%s v%s" % (__appname__, __version__))
    logger.info(args)

    # Read in meltedResults
    #cm = read_csv(args.meltedResults, sep="\t")
    #logger.info("Read %d lines", len(cm))
    #cm = cm[cm.SNPs_Compared >= args.threshold]
    # Read in Patient to Sample mapping
    #logger.info("%d lines after filtering", len(cm))
    patients = readPatientToSample(args.patientToSample)

    # Match samples to patients
    matchSamplesToPatients(patients, args.meltedResults_s3path, args.threshold, args.resultsFile)

def matchSamplesToPatients(patients, meltedResults_s3path, threshold, resultsFile):
    # Scan each patient
    # For each sample, make sure the top match is another sample from the same patient
    missingSamples = []
    badSamples = []
    mismatchedSamples = 0
    with open(resultsFile, 'w') as output:
      output.write("Sample1\tSample2\tn_S1\tn_S2\tSNPs_Compared\tFraction_Match\n")
      for patient in patients:
        samples = patients[patient]
        logger.info("Patient: %s, Samples: %s", patient, samples)
        if len(samples) > 1:
            for sample in samples:
                meltedResultsFile = "%s.meltedResults.txt" % sample
                downloadFile("%s/%s.meltedResults.txt" % (meltedResults_s3path, sample), meltedResultsFile)
                if not os.path.exists(meltedResultsFile):
                    logger.error("Unable to download %s", meltedResultsFile)
                    continue
                cm = read_csv(meltedResultsFile, sep="\t")
                os.remove(meltedResultsFile)
                cm = cm[cm.SNPs_Compared >= threshold]
                sample_matches = cm[cm.Sample1==sample]
                if len(sample_matches) > 0:
                    otherSample = 'Sample2'
                else:
                    sample_matches = cm[cm.Sample2==sample]
                    otherSample = 'Sample1'
                    if len(sample_matches) == 0:
                        missingSamples.append(sample)

                if len(sample_matches) > 0:
                    sample_matches = sample_matches.sort_values(by="Fraction_Match", ascending=False)
                    topMatch = sample_matches.iloc[0][otherSample]
                    fractionMatch = sample_matches.iloc[0]['Fraction_Match']
                    if fractionMatch < 0.7:
                        badSamples.append(sample)
                    else:
                        if topMatch in samples:
                            logger.info("%s -> %s", sample, topMatch)
                        else:
                            logger.warn("%s -> %s:%s", sample, patient, topMatch)
                            mismatchedSamples += 1
                            for i in range(0, 5):
                                output.write("%s\t%s\t%d\t%d\t%d\t%0.4f\n" % (sample_matches.iloc[i]['Sample1'],
                                                                     sample_matches.iloc[i]['Sample2'],
                                                                     sample_matches.iloc[i]['n_S1'],
                                                                     sample_matches.iloc[i]['n_S2'],
                                                                     sample_matches.iloc[i]['SNPs_Compared'],
                                                                     sample_matches.iloc[i]['Fraction_Match']))
                else:
                    logger.warn("%s -> No matches", sample)

    logger.info("Mismatched Samples: %d" % mismatchedSamples)
    logger.info("There were %d samples not found in cohort-matcher meltedResults. " \
          "They could have been filtered out due to low coverage: %s" % (len(missingSamples), missingSamples))
    logger.info("There were %d samples that did not have a good match (Fraction_Match >= 0.7) to any samples: %s" %
          (len(badSamples), badSamples))

def parseArguments(argv):
    ''' Parse Arguments '''
    parser = argparse.ArgumentParser(description="Find swaps in cohort-matcher results")
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    parser_grp1 = parser.add_argument_group("Required")
    parser_grp1.add_argument("--meltedResults_s3path", required=True, help="Melted results files from cohort-matcher")
    parser_grp1.add_argument("--patientToSample", required=True, help="Patient to Sample Mapping file")
    parser_grp1.add_argument("--threshold", default=15, help="Min # of SNPS compared between samples")
    parser_grp1.add_argument("--resultsFile", default="cohort-matcher_results.txt", help="Results of mismatched samples")
    args = parser.parse_args(argv)
    return args

if __name__ == '__main__':
    main(sys.argv[1:])
