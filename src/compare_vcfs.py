'''
Script to test compareGenotypes.py by providing two local VCF files
'''
import argparse
import sys

from compareGenotypes import VCFtoTSV, get_tsv_variants, getIntersectingVariants, compareGenotypes

def parse_arguments(argv):
    parser = argparse.ArgumentParser(description='Compare a sample to a set of samples')
    parser.add_argument('--vcf1', required=True)
    parser.add_argument('--vcf2', required=True)
    return parser.parse_args(argv)

def main(argv):
    ''' Main Entry Point '''
    args = parse_arguments(argv)
    dp_threshold = 15

    VCFtoTSV(args.vcf1, 'vcf1.tsv')
    var_list1 = get_tsv_variants('vcf1.tsv', dp_threshold)
    VCFtoTSV(args.vcf2, 'vcf2.tsv')
    var_list2 = get_tsv_variants('vcf2.tsv', dp_threshold)

    intersection = getIntersectingVariants(var_list1, var_list2)
    results = compareGenotypes(var_list1, var_list2, intersection)

if __name__ == '__main__':
    main(sys.argv[1:])