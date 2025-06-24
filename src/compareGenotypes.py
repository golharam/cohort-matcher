#!/usr/bin/env python
'''
Script to compare genotypes of samples
'''
from __future__ import division
import argparse
from math import log
import logging
import os
import sys
import csv
from fisher import pvalue
import vcf

from genotypeSamples import read_samples

from common import downloadFile, uploadFile, generate_working_dir, delete_working_dir

__appname__ = 'compareGenotypes'
__version__ = "0.3"

logger = logging.getLogger(__appname__)

def is_hom(gt):
    ''' Test for homozygosity '''
    gt_ = gt.split("/")
    return gt_[0] == gt_[1]

def is_same_gt(gt1, gt2):
    ''' Check is same genotype '''
    if gt1 == gt2:
        return True
    gt1_ = sorted(gt1.split("/"))
    gt2_ = sorted(gt2.split("/"))
    return gt1_ == gt2_

def is_subset(hom_gt, het_gt):
    ''' Test if hom_gt is subset of het_gt '''
    gt_hom = hom_gt.split("/")[0]
    return gt_hom in het_gt

def compareGenotypes(var_list, var_list2, intersection):
    ''' Compare genotypes of two samples '''
    ct_common = 0
    comm_hom_ct = 0
    comm_het_ct = 0
    ct_diff = 0
    diff_hom_ct = 0
    diff_het_ct = 0
    diff_1sub2_ct = 0
    diff_hom_het_ct = 0
    diff_2sub1_ct = 0
    diff_het_hom_ct = 0

    for pos_ in intersection:
        gt1 = var_list[pos_]['GT']
        #if alternate_chroms is not None:
        #    bits = pos_.split('\t')
        #    gt2 = var_list2[def_to_alt[bits[0]] + "\t" + bits[1]]['GT']
            #gt2 = def_to_alt[pos_]
        #else:
        gt2 = var_list2[pos_]['GT']

        # if genotypes are the same
        if is_same_gt(gt1, gt2):
            ct_common += 1
            if is_hom(gt1):
                comm_hom_ct += 1
            else:
                comm_het_ct += 1

            # P-Value calc
            bits = pos_.split('\t')
        else:
            ct_diff += 1
            # both are hom and different
            if is_hom(gt1) and is_hom(gt2):
                diff_hom_ct += 1
            # both are het and different
            elif is_hom(gt1) is False and is_hom(gt2) is False:
                diff_het_ct += 1
            # one is hom, one is het, test for subset
            elif is_hom(gt1):
                if is_subset(gt1, gt2):
                    diff_1sub2_ct += 1
                else:
                    diff_hom_het_ct += 1
            elif is_hom(gt2):
                if is_subset(gt2, gt1):
                    diff_2sub1_ct += 1
                else:
                    diff_het_hom_ct += 1
            else:
                print("WTF?")
                print(gt1, gt2)
                exit(1)

    total_compared = ct_common + ct_diff
    frac_common_plus = 0
    frac_common = 0
    if total_compared > 0:
        frac_common = float(ct_common)/total_compared
        frac_common_plus = float(ct_common + max(diff_2sub1_ct,
                                                 diff_1sub2_ct))/total_compared

    # test of allele-specific genotype subsets
    allele_subset = ""
    sub_sum = diff_1sub2_ct + diff_2sub1_ct
    # don't bother if fewer than 10
    if sub_sum > 10:
        pv_set = pvalue(diff_1sub2_ct, diff_2sub1_ct, ct_diff/2, ct_diff/2)
        pv_ = min(pv_set.left_tail, pv_set.right_tail)
        if pv_ < 0.05:
            if diff_1sub2_ct < diff_2sub1_ct:
                allele_subset = "2sub1"
                frac_common_plus = float(ct_common + diff_2sub1_ct) / total_compared
            else:
                allele_subset = "1sub2"
                frac_common_plus = float(ct_common + diff_1sub2_ct) / total_compared

    results = {}
    results['total_compared'] = total_compared
    results['ct_common'] = ct_common
    results['frac_common'] = frac_common
    results['frac_common_plus'] = frac_common_plus
    results['comm_hom_ct'] = comm_hom_ct
    results['comm_het_ct'] = comm_het_ct
    results['ct_diff'] = ct_diff
    results['diff_hom_ct'] = diff_hom_ct
    results['diff_het_ct'] = diff_het_ct
    results['diff_hom_het_ct'] = diff_hom_het_ct
    results['diff_het_hom_ct'] = diff_het_hom_ct
    results['diff_1sub2_ct'] = diff_1sub2_ct
    results['diff_2sub1_ct'] = diff_2sub1_ct
    results['allele_subset'] = allele_subset
    results['judgement'], results['short_judgement'] = makeJudgement(total_compared,
                                                                     frac_common,
                                                                     frac_common_plus,
                                                                     allele_subset)
    return results

def makeJudgement(total_compared, frac_common, frac_common_plus, allele_subset):
    ''' Make judgement of sample similarity based on genotype comparison '''
    A_BIT_LOW = """ the number of comparable genomic loci is a bit low. Try using a
                    different variants list (--VCF) file which have more appropriate
                    genomic positions for comparison."""

    if total_compared <= 20:
        judgement = "Inconclusive: Too few loci to compare"
        short_judgement = "INCONCLUSIVE"
    elif total_compared <= 100:
        # allow for 0.90 frac_common for low loci count
        if frac_common >= 0.9 or frac_common_plus >= 0.9:
            judgement = "LIKELY SAME SOURCE: %s" % A_BIT_LOW
            short_judgement = "LIKELY SAME"
            # if there is possible allele-specific genotype subset
            if allele_subset == "1sub2" or allele_subset == "2sub1":
                sub_ = allele_subset.split("sub")[0]
                over_ = allele_subset.split("sub")[1]
                judgement += """ BAM%s genotype appears to be a subset of BAM%s.
                                 Possibly BAM%s is RNA-seq data or BAM%s is 
                                 ontaminated.""" % (sub_, over_, sub_, over_)
                short_judgement += ". (BAM%s is subset of BAM%s)" % (sub_, over_)
        elif frac_common <= 0.6:
            judgement = "LIKELY FROM DIFFERENT SOURCES: %s" % A_BIT_LOW
            short_judgement = "LIKELY DIFFERENT"
        else:
            judgement = "INCONCLUSIVE: %s" % A_BIT_LOW
            short_judgement = "INCONCLUSIVE"
    # 3. >100 sites compared
    else:
        if frac_common >= 0.95:
            judgement = "BAM FILES ARE FROM THE SAME SOURCE"
            short_judgement = "SAME"
        elif frac_common_plus >= 0.95:
            judgement = "BAM FILES ARE VERY LIKELY FROM THE SAME SOURCE"
            short_judgement = "LIKELY SAME"
            if allele_subset == "1sub2" or allele_subset == "2sub1":
                sub_ = allele_subset.split("sub")[0]
                over_ = allele_subset.split("sub")[1]
                judgement += """, but with possible allele specific genotype. BAM%s
                genotype appears to be a subset of BAM%s. Possibly BAM%s is RNA-seq
                data or BAM%s is contaminated.""" % (sub_, over_, sub_, over_)
                short_judgement += ". (BAM%s is subset of BAM%s)" % (sub_, over_)
        elif frac_common <= 0.6:
            judgement = "BAM FILES ARE FROM DIFFERENT SOURCES"
            short_judgement = "DIFFERENT"
        elif frac_common >= 0.8:
            judgement = """LIKELY FROM THE SAME SOURCE. However, the fraction of sites
                           with common genotype is lower than expected. This can happen
                           with samples with low coverage."""
            short_judgement = "LIKELY SAME"
        else:
            judgement = "LIKELY FROM DIFFERENT SOURCES"
            short_judgement = "LIKELY DIFFERENT"

    return judgement, short_judgement

def getIntersectingVariants(var_list, var_list2):
    ''' Given two tsv list of variants, find the intersection '''
    intersection = []
    for var2 in var_list2:
        var_ = var2
        if var_ in var_list:
            intersection.append(var_)
    return intersection

def main(argv):
    ''' Main Entry Point '''
    args = parseArguments(argv)
    logging.basicConfig(level=args.log_level)
    logger.info("%s v%s", __appname__, __version__)
    logger.info(args)

    working_dir = generate_working_dir(args.working_dir)
    logger.info("Working in %s", working_dir)

    sample_name = args.sample

    # Download and read the bamsheet from the s3 cache directory
    downloadFile("%s/bamsheet.txt" % args.s3_cache_folder, "%s/bamsheet.txt" % working_dir)
    # Determine the current sample index in the list of samples
    samples = read_samples("%s/bamsheet.txt" % working_dir)
    sample_index = -1
    for i, sample in enumerate(samples):
        if sample['sample_id'] == sample_name:
            sample_index = i
            break
    if sample_index == -1:
        logger.error("Unable to locate sample in bamsheet")
        return -1

    # Get the list of variants in the reference sample
    s3_vcf_file = "%s/%s.vcf" % (args.s3_cache_folder, sample_name)
    vcf_file = "%s/%s.vcf" % (working_dir, sample_name)
    downloadFile(s3_vcf_file, vcf_file)
    if os.path.exists(vcf_file):
        tsvFile = "%s/%s.tsv" % (working_dir, sample_name)
        VCFtoTSV(vcf_file, tsvFile)
        if os.path.exists(tsvFile):
            var_list = get_tsv_variants(tsvFile, args.dp_threshold)
            os.remove(tsvFile)
        else:
            logger.error("Failed to convert VCF to TSV, %s -> %s", vcf_file, tsvFile)
            return -1
        os.remove(vcf_file)
    else:
        logger.error("Failed to download %s", s3_vcf_file)
        return -1

    meltedResultsFile = "%s/%s.meltedResults.txt" % (working_dir, sample_name)
    with open(meltedResultsFile, "w") as fout:
        fout.write("Sample1\tSample2\tn_S1\tn_S2\tSNPs_Compared\tFraction_Match\tJudgement\n")
        sample_index += 1
        while sample_index < len(samples):
            sample = samples[sample_index]
            logger.info("[%d/%d] Comparing %s - %s", sample_index+1, len(samples),
                        sample_name, sample['sample_id'])

            s3_vcf_file = "%s/%s.vcf" % (args.s3_cache_folder, sample['sample_id'])
            vcf_file = "%s/%s.vcf" % (working_dir, sample['sample_id'])
            downloadFile(s3_vcf_file, vcf_file)
            if not os.path.exists(vcf_file):
                logger.error("Error downloading %s.  Aborting.", s3_vcf_file)
                return -1
            tsvFile = vcf_file.replace('.vcf', '.tsv')
            VCFtoTSV(vcf_file, tsvFile)
            var_list2 = get_tsv_variants(tsvFile, args.dp_threshold)
            os.remove(vcf_file)
            os.remove(tsvFile)

            # compare the genotypes
            intersection = getIntersectingVariants(var_list, var_list2)
            results = compareGenotypes(var_list, var_list2, intersection)
            n1 = '%d' % len(var_list)
            n2 = '%d' % len(var_list2)
            fm = '%.4f' % results['frac_common']
            tc = '%d' % results['total_compared']
            j = results['short_judgement']
            fout.write('\t'.join([sample_name, sample['sample_id'], n1, n2, tc, fm, j]) + "\n")
            sample_index += 1
    logger.info("Uploading %s to %s", meltedResultsFile, args.s3_cache_folder)
    uploadFile(meltedResultsFile, "%s/%s" % (args.s3_cache_folder,
                                             os.path.basename(meltedResultsFile)))

    logger.info('Cleaning up working dir')
    os.remove(meltedResultsFile)
    delete_working_dir(working_dir)
    logger.info('Completed')

def parseArguments(argv):
    ''' Parse arguments '''
    parser = argparse.ArgumentParser(description='Compare a sample to a set of samples')
    parser.add_argument('--log-level', help="Prints warnings to console by default",
                        default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    required_args = parser.add_argument_group("Required")
    required_args.add_argument('-s', '--sample', required=True, help="Sample of interest")
    required_args.add_argument("--s3_cache_folder", "-CD", required=True,
                               help="Specify S3 path for cached VCF/TSV files")

    parser.add_argument('--working_dir', type=str, default='/scratch')
    parser.add_argument("-t", "--dp_threshold", default=15, help="Depth of Coverage Threshold")

    args = parser.parse_args(argv)
    return args

def get_tsv_variants(tsvFile, dp_threshold):
    ''' Return a list of variants that pass threshold '''
    var_list = {}
    with open(tsvFile, "r") as fin:
        for line in fin:
            if line.startswith("CHROM\t"):
                continue
            bits = line.strip("\n").split("\t")
            if bits[5] == "NA":
                continue
            elif int(bits[5]) < dp_threshold:
                continue
            else:
                var_list["\t".join(bits[:2])] = {'REF': bits[2], 'ALT': bits[3],
                                                 'DP': int(bits[5]), 'GT': bits[7]}
    return var_list

def VCFtoTSV(invcf, outtsv, caller="freebayes"):
    '''
    Convert a VCF to TSV
    '''
    logger.debug("%s -> %s", invcf, outtsv)
    fout = open(outtsv, "w")
    vcf_in = vcf.Reader(open(invcf, "r"))
    var_ct = 0
    if caller == "gatk" or caller == "varscan":
        fields_to_extract = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AD", "GT"]
    elif caller == "freebayes":
        fields_to_extract = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AO", "GT"]
    fout.write("%s\n" % "\t".join(fields_to_extract))
    try:
      for var in vcf_in:
        if var.is_indel:
            continue
        if var.is_snp is False and var.is_monomorphic is False:
            pass
        logger.debug("%s - %s:%s", invcf, var.CHROM, var.POS)
        chrom_ = var.CHROM.replace("chr", "")
        pos_ = str(var.POS)
        ref_ = var.REF
        qual_ = str(var.QUAL)
        if caller == "varscan":
            qual_ = str(var.samples[0]["GQ"])

        dp_ = "0"
        #ad_or_ao = "NA"
        ad_str = "NA"
        gt_ = "NA"
        alt_str = ""

        if var.samples[0].called is False:
            continue

        # usually need to bypass indels, however,
        # homozygous REF is considered indel by pyvcf... WTF?
        if var.is_monomorphic:
            alt_str = "."
            gt_ = "%s/%s" % (ref_, ref_)
            if caller == "freebayes":
                dp_ = var.samples[0].data.DP
                ro_ = var.samples[0].data.RO
                ad_str = str(dp_ - ro_)
            elif caller == "gatk":
                dp_ = var.samples[0].data.DP
                ad_str = "0"
            elif caller == "varscan":
                dp_ = var.INFO["ADP"]
            dp_ = str(dp_)
        else:
            alt_ = var.ALT[0]
            alt_str = "."
            if alt_ != None:
                alt_str = alt_.sequence
            if caller == "freebayes" or caller == "gatk":
                dp_ = str(var.INFO["DP"])
            else:
                dp_ = str(var.INFO["ADP"])

            gt_ = var.samples[0].gt_bases

            if caller == "gatk":
                ad_ = var.samples[0]["AD"]
                for a_ in ad_:
                    ad_str += ",%d" % a_
                ad_str = ad_str[1:]
            elif caller == "varscan":
                ad_str = str(var.samples[0]["AD"])
            else:
                ad_str = str(var.samples[0]["AO"])

        if ad_str == "NA":
            ad_str = "0"

        data_bits = [chrom_, pos_, ref_, alt_str, qual_, dp_, ad_str, gt_]
        fout.write("%s\n" % "\t".join(data_bits))
        var_ct += 1
      fout.close()
    except:
      os.remove(outtsv)
      logger.error("Could not convert %s -> %s", invcf, outtsv)
    return var_ct

if __name__ == '__main__':
    main(sys.argv[1:])
