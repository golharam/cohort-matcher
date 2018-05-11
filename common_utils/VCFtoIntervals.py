import argparse
import logging
import sys

import vcf

logger = logging.getLogger("VCFtoIntervals")
__version__ = "0.1"

'''
This script really does nothing more than convert a VCF to a BED file appropriate for the
variant caller used.  It can also create a BED interval window around a SNP.   The real
benefit of this is that it can convert a VCF from one reference to a BED file for another
reference using the chromosome map, eg hg19 - GRCh37 (chr to no-chr prefix).
'''

def get_chrom_data_from_map(chrom_map_file):
    ''' Get mapping of chr names '''
    chrom_ct = 0
    default_chroms = []
    alternate_chroms = []
    def_to_alt = {}
    alt_to_def = {}

    with open(chrom_map_file, 'r') as fin:
        fin.readline()
        for line in fin:
            if line.strip() == "":
                continue
            chrom_ct += 1
            bits = line.strip().split()
            default_chroms.append(bits[0])
            alternate_chroms.append(bits[1])
            def_to_alt[bits[0]] = bits[1]
            alt_to_def[bits[1]] = bits[0]
    return default_chroms, alternate_chroms, def_to_alt, alt_to_def

def main(argv):
    ''' Main Entry Point '''
    config = parseArguments(argv)
    logging.basicConfig(level=config.log_level)
    logger.info("cohort-matcher v%s" % __version__)
    logger.info(config)

    vcfToIntervals(config.vcf, config.bed, cmap=config.chromosome_map)

def parseArguments(argv):
    ''' Parse Arguments '''
    parser = argparse.ArgumentParser(description="Convert VCF to BED file")
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    parser.add_argument('-V', '--vcf', required=True, help="VCF File")
    parser.add_argument('-B', '--bed', required=True, help="Output BED file")
    parser.add_argument("--chromosome-map", "-CM", default=None,
                        help="Chromosome mapping, if two reference is different than VCF")

    args = parser.parse_args(argv)
    return args

def vcfToIntervals(vcfFile, bedFile, window=0, caller="freebayes", cmap=None):
    '''
    Convert a vcf file to a 3-column interval/bed file
    '''
    # Since bedFile is being written to cache, don't use the old one.
    #if os.path.exists(bedFile) is True:
    #    logger.info("%s already exists.", bedFile)
    #    return

    if cmap:
        default_chroms, alternate_chroms, \
        def_to_alt, alt_to_def = get_chrom_data_from_map(cmap)
    else:
        default_chroms, alternate_chroms, def_to_alt, alt_to_def = None, None, None, None

    vcf_read = vcf.Reader(open(vcfFile, "r"))
    fout = open(bedFile, "w")
    for var in vcf_read:
        # intervals format
        chrom_ = var.CHROM
        if cmap != None:
            if chrom_ in default_chroms:
                chrom_ = def_to_alt[chrom_]
            elif chrom_ in alternate_chroms:
                chrom_ = alt_to_def[chrom_]

        if caller == "gatk" or caller == "varscan":
            start_pos = var.POS - window
            end_pos = start_pos + window
            fout.write("%s:%d-%d\n" % (chrom_, start_pos, end_pos))
        # BED format
        elif caller == "bed" or caller == "freebayes":
            start_pos = var.POS - window - 1
            end_pos = start_pos + window + 1
            fout.write("%s\t%d\t%d\n" % (chrom_, start_pos, end_pos))
    fout.close()
    return

if __name__ == "__main__":
    main(sys.argv[1:])

