'''
Script to test compareGenotypes.py by providing two local VCF files
'''
import argparse
import sys

from compareGenotypes import VCFtoTSV, get_tsv_variants, getIntersectingVariants, compareGenotypes

def read_allele_freqs(allele_freqs_file, sep="\t"):
    var_list = {}
    with open(allele_freqs_file, "r") as fin:
        for line in fin:
            if line.startswith('CHROM'):
                continue
            chrom, pos, ref, alt, af = line.strip("\n").split("\t")
            chrom_ = chrom.replace("chr", "")
            var_list[f"{chrom}\t{pos}"] = {
                'REF': ref, 'ALT': alt, 'AF': float(af)}
    return var_list

def parse_arguments(argv):
    parser = argparse.ArgumentParser(description='Compare a sample to a set of samples')
    parser.add_argument('--vcf1', required=True)
    parser.add_argument('--vcf2', required=True)
    parser.add_argument("-t", "--dp_threshold", default=15, help="Depth of Coverage Threshold")

    # Allele-Specific Expression
    optional_args = parser.add_argument_group("Allele-Specific Expression")
    optional_args.add_argument('--consider-alleles', action="store_true", default=False, 
                               help="Consider allele-specific expression")

    # LRR
    optional_args = parser.add_argument_group("Log-Odds Ratio")
    optional_args.add_argument('-af', '--allele-freqs', required=True, help="Path to allele frequencies TSV file")
    optional_args.add_argument('--error-rate', type=float, default=0.01, help="Error rate for genotype comparison")

    return parser.parse_args(argv)

def main(argv):
    ''' Main Entry Point '''
    args = parse_arguments(argv)
    dp_threshold = args.dp_threshold

    # Load allele frequencies
    allele_freqs = read_allele_freqs(args.allele_freqs)

    VCFtoTSV(args.vcf1, 'vcf1.tsv')
    var_list1 = get_tsv_variants('vcf1.tsv', dp_threshold)
    VCFtoTSV(args.vcf2, 'vcf2.tsv')
    var_list2 = get_tsv_variants('vcf2.tsv', dp_threshold)

    intersection = getIntersectingVariants(var_list1, var_list2)
    results = compareGenotypes(var_list1, var_list2, intersection, allele_freqs, args.error_rate, args.consider_alleles)

if __name__ == '__main__':
    main(sys.argv[1:])