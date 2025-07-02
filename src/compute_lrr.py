import argparse
import logging
import sys

import pysam
import pandas as pd
import math

def get_genotype(record):
    gt = record.samples[0]['GT']
    if None in gt:
        return None
    alleles = [record.alleles[i] for i in gt]
    return ''.join(sorted(alleles))

def genotype_prob(geno, p):
    q = 1 - p
    if geno == "AA":
        return p ** 2
    elif geno == "AB":
        return 2 * p * q
    elif geno == "BB":
        return q ** 2
    return 0.0

def simplify(geno):
    return ''.join(sorted(geno))

def compute_llr(vcf_A, vcf_B, error_rate, freqs):
    total_llr = 0.0

    for rec_A in vcf_A.fetch():
        chrom, pos = rec_A.contig, rec_A.pos
        try:
            rec_B = next(vcf_B.fetch(chrom, pos - 1, pos))
        except StopIteration:
            continue

        key = (int(chrom.replace('chr', '')), pos)
        if key not in freqs.index:
            continue

        #p = freqs.loc[key, "FREQ_REF"]
        #q = 1 - p
        q = freqs.loc[key, "AF"]
        p = 1 - q

        gA = get_genotype(rec_A)
        gB = get_genotype(rec_B)

        if gA is None or gB is None:
            continue

        gA = simplify(gA)
        gB = simplify(gB)

        # Recode to AA, AB, BB (assuming REF is allele1)
        def recode(geno, ref, alt):
            if geno == ref + ref:
                return "AA"
            elif geno == ref + alt or geno == alt + ref:
                return "AB"
            elif geno == alt + alt:
                return "BB"
            return None

        alt_A = rec_A.alts[0] if rec_A.alts and rec_A.alts[0] != '.' else None
        alt_B = rec_B.alts[0] if rec_B.alts and rec_B.alts[0] != '.' else None
        geno_A = recode(gA, rec_A.ref, alt_A)
        geno_B = recode(gB, rec_B.ref, alt_B)

        if geno_A is None or geno_B is None:
            continue

        L_same = 1 - error_rate if geno_A == geno_B else error_rate
        L_diff = genotype_prob(geno_A, p) * genotype_prob(geno_B, p)
        if L_diff == 0:
            continue
        total_llr += math.log10(L_same / L_diff)

    return total_llr

def parse_arguments(argv):
    ''' Parse arguments '''
    parser = argparse.ArgumentParser(description='Compute Log Likelihood Ratio from VCF files')
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument('--vcf1', required=True, help="Path to first VCF.GZ file")
    parser.add_argument('--vcf2', required=True, help="Path to second VCF.GZ file")
    parser.add_argument('--freqs', required=True, help="Path to allele frequencies TSV file")
    parser.add_argument('--error-rate', type=float, default=0.01, help="Error rate for genotype comparison")
    args = parser.parse_args(argv)
    return args

def main(argv):
    ''' Main Entry Point '''
    args = parse_arguments(argv)
    logging.basicConfig(level=args.log_level)
    logger = logging.getLogger(__name__)
    logger.info(args)

    # Configuration
    vcf_A_path = args.vcf1
    vcf_B_path = args.vcf2
    freq_file = args.freqs # "allele_frequencies.tsv"  # TSV with columns: CHROM, POS, REF, ALT, FREQ_REF

    error_rate = args.error_rate

    # Load allele frequencies
    freqs = pd.read_csv(freq_file, sep="\t")

    # Index by (CHROM, POS)
    freqs.set_index(["CHROM", "POS"], inplace=True)

    # Open VCFs
    vcf_A = pysam.VariantFile(vcf_A_path)
    vcf_B = pysam.VariantFile(vcf_B_path)

    llr = compute_llr(vcf_A, vcf_B, error_rate, freqs)
    print(f"Total Log Likelihood Ratio: {llr:.4f}")

if __name__ == '__main__':
    main(sys.argv[1:])
