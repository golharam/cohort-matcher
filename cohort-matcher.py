#!/bin/env python
import argparse
import ConfigParser
import json
import os
import pysam
import subprocess
import sys
import vcf
from fisher import pvalue

def getArgParser(genotyping = False):
	parser = argparse.ArgumentParser(description="Compare two vcf files to see if \
	they are from the same samples, using frequently occuring SNPs \
	reported in the 1000Genome database")

	parser_grp1 = parser.add_argument_group("Required")
	if genotyping == False:
		parser_grp1.add_argument("--set1", "-S1", required=True, help="First set of samples")
		parser_grp1.add_argument("--set2", "-S2", required=True, help="Second set of samples")

	parser_grp2 = parser.add_argument_group("Configuration")
	parser_grp2.add_argument("--cache-dir", "-CD", required=False,
							default="./cache",
							help="Specify directory for cached data. (Default: ./cache)")
	parser_grp2.add_argument("--scratch-dir", "-SD", required=False,
							default="/scratch",
							help="Specify scratch directory. (Default: /scratch)")
	
	parser_grp3 = parser.add_argument_group("Genotyper")
	parser_grp3.add_argument("--caller", "-CL", required=False,
							default="freebayes",
							choices=('gatk', 'freebayes', 'varscan'),
							help="Specify which caller to use (default = 'freebayes')")
	parser_grp3.add_argument("--dp-threshold",   "-DP", required=False,
							default=15, type=int, 
							help="Minimum required depth for comparing variants")
	parser_grp3.add_argument("--number_of_snps", "-N", required=False,
							type=int, help="Number of SNPs to compare.")
	
	parser_grp4 = parser.add_argument_group("Reference")
	parser_grp4.add_argument("--vcf", "-V",  required=True,
							 help="VCF file containing SNPs to check (default can be specified in config file instead)")
	parser_grp4.add_argument("--vcf2", "-V2",  default=None,
							 help="VCF file containing SNPs to check (default can be specified in config file instead) when using reference2")
	parser_grp4.add_argument("--reference", "-R", required=True,
							help="Reference FASTA File (indexed with samtools faidx) for set1 or both sets")
	parser_grp4.add_argument("--reference2", "-R2", default=None,
							help="Reference FASTA File (indexed with samtools faidx) for set2")
	parser_grp4.add_argument("--chromosome-map", "-CM", default=None,
							help="Chromosome mapping, if two reference are used")
	
	parser_grp5 = parser.add_argument_group("Freebayes Settings")
	parser_grp5.add_argument("--freebayes-path", required=False,
							default="freebayes",
							help="Path to freebayes binary (if not in PATH")
							
	parser_grp6 = parser.add_argument_group("Cluster")
	parser_grp6.add_argument("--cluster", required=False,
							default="local",
							choices=('local', 'sge'),
							help="Specify type of cluster architecture to use (Default: local)")
	
	parser_grp7 = parser.add_argument_group("Miscellaneous")
	parser_grp7.add_argument("--aws-path", required=False, default="/usr/bin/aws", help="Specify path to aws cli")

	parser_grp8 = parser.add_argument_group("Output")	
	parser_grp8.add_argument("--output-dir", "-O", required=False,
						default="./cohort-matcher-output",
						help="Specify output directory for sample comparisons (Default: ./cohort-matcher-output)")
	parser_grp8.add_argument("--short-output", "-so", required=False, default=False, action="store_true",
						help="Short output format (Default: False")
	'''
	# Freebayes options
	parser_grp3.add_argument("--fastfreebayes",  "-FF", required=False, default=False,
							 action="store_true", help="Use --targets option for Freebayes.")
	# GATK options
	parser_grp3.add_argument("--gatk-mem-gb" ,   "-GM", required=False,
							 type=int, help="Specify Java heap size for GATK (GB, int)")
	parser_grp3.add_argument("--gatk-nt" ,	   "-GT", required=False,
							 type=int, help="Specify number of threads for GATK UnifiedGenotyper (-nt option)")
	# VarScan Options
	parser_grp3.add_argument("--varscan-mem-gb", "-VM", required=False,
							 type=int, help="Specify Java heap size for VarScan2 (GB, int)")
	'''
	
	parser.add_argument("--verbose", "-v", required=False, action="store_true", help="Verbose (Default = False)")
	parser.add_argument("--debug", "-d", required=False, action="store_true", help="Verbose (Default = False)")
	return parser

def parseArgs(argv = None):
	# Do argv default this way, as doing it in the functional
	# declaration sets it at compile time.
	if argv is None:
		argv = sys.argv[1:]
	
	parser = getArgParser()
	args = parser.parse_args(argv)
	
	if args.verbose:
		print "Set1: {}".format(args.set1)
		print "Set2: {}".format(args.set2)
		print ""
		print "Cache Dir: {}".format(args.cache_dir)
		print "Scratch Dir: {}".format(args.scratch_dir)
		print "Output Dir: {}".format(args.output_dir)
		print "Short Report: {}".format(args.short_output)
		print ""
		print "Genotype Caller: {}".format(args.caller)
		if args.caller == "freebayes":
			print "Path to freebayes: {}".format(args.freebayes_path)			
		print "Depth Threshold: {}".format(args.dp_threshold)
		print "Number of SNPs: {}".format(args.number_of_snps)
		print ""
		print "VCF File: {}".format(args.vcf)
		if args.vcf2 is not None:
			print "VCF2 File: {}".format(args.vcf2)
		print "Reference: {}".format(args.reference)
		if args.reference2 is not None:
			print "Reference2: {}".format(args.reference2)
		if args.chromosome_map is not None:
			print "Chromosome Map: {}".format(args.chromosome_map)
		print ""
		print "Cluster Type: {}".format(args.cluster)
		print "AWS Path: {}".format(args.aws_path)
		print ""
		print "Working in " + os.getcwd()
	return args

def checkConfig(config):
	if os.path.isdir(config.cache_dir) == False:
		print "Cache dir (%s) is not a directory" % config.cache_dir
		exit(1)
		
	if os.path.isdir(config.scratch_dir) == False:
		print "Scratch dir (%s) is not a directory" % config.scratch_dir
		exit(1)

	if os.path.exists(config.vcf) == False:
		print "VCF file (%s) is not accessible" % config.vcf
		exit(1)

	if config.vcf2 is not None:
		if os.path.exists(config.vcf2) == False:
			print "VCF2 (%s) is not accessible" % config.vcf2
			exit(1)
			
	if os.path.exists(config.reference) == False:
		print "Reference FASTA file (%s) is not accessible" % config.reference
		exit(1)
	reference_index = config.reference + ".fai"
	if os.path.exists(reference_index) == False:
		print "Reference FASTA file (%s) is not indexed" % reference_index
		exit(1)

	if config.reference2 is not None:
		if os.path.exists(config.reference2) == False:
			print "Reference2 FASTA file (%s) is not accessible" % config.reference2
			exit(1)
		reference_index = config.reference2 + ".fai"
		if os.path.exists(reference_index) == False:
			print "Reference 2 FASTA file (%s) is not indexed" % reference_index
			exit(1)
	if config.chromosome_map is not None:
		if os.path.exists(config.chromosome_map) == False:
			print "Chromosome map (%s) could not be read" % config.chromosome_map
			exit(1)
	
	''' If reference2, vcf2, or chromosome-map is specified, make sure all are specified '''
	if config.reference2 is not None or config.vcf2 is not None or config.chromosome_map is not None:
		if config.reference2 is None:
			print "Reference2 must be specified if vcf2 or chromosome-map is specified"
			exit(1)
		if config.vcf2 is None:
			print "vcf2 must be speicfied if reference2 or chromosome-map is specified"
			exit(1)
		if config.chromosome_map is None:
			print "Chromosome map must be specified is referenc2 or vcf2 is specified"
			exit(1)
			
	if os.path.isdir(config.output_dir) == False:
		print "Output directory (%s) does not exist" %config.output_dir
		exit(1)
		
	if config.caller == "freebayes" and config.freebayes_path is None:
		print "Freebayes-path must be set when caller is Freebayes"
		exit(1)
		
def readSamples(sampleSheetFile, verbose):
    if os.path.isfile(sampleSheetFile) == False:
        print "%s does not exist" % sampleSheetFile
        exit(1)    
    
    if verbose:
    	print "Reading %s" % sampleSheetFile
    	
    samples = []        
    with open(sampleSheetFile, 'r') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if len(line) == 0 or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) != 2:
            	print "Expect 2 fields (sampleName, bamFile) but encountered %d" % len(fields)
            	exit(1)
            	 
            sampleBamfile = fields[1]
            sample = {"name": fields[0],
                      "bam": fields[1]}
            samples.append(sample)
    return samples

def vcfToIntervals(vcfFile, bedFile, window = 0, format="freebayes", cmap=None):
    vcf_read = vcf.Reader(open(vcfFile, "r"))
    fout = open(bedFile, "w")

    for var in vcf_read:
        # intervals format
        chrom_ = var.CHROM
        if cmap != None:
            if var.CHROM not in cmap:
                continue
            else:
                chrom_ = cmap[var.CHROM]

        if format == "gatk" or format == "varscan":
            start_pos = var.POS - window
            end_pos   = start_pos + window
            fout.write("%s:%d-%d\n" % (chrom_, start_pos, end_pos))
            
        # BED format
        elif format == "bed" or format == "freebayes":
            start_pos = var.POS - window - 1
            end_pos   = start_pos + window + 1
            fout.write("%s\t%d\t%d\n" % (chrom_, start_pos, end_pos))

    fout.close()
    return

def VCFtoTSV(invcf, outtsv, caller):
    fout = open(outtsv, "w")
    vcf_in = vcf.Reader(open(invcf, "r"))
    var_ct = 0
    if caller == "gatk" or caller == "varscan":
        fields_to_extract = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AD", "GT"]
    elif caller == "freebayes":
        fields_to_extract = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AO", "GT"]
    fout.write("%s\n" % "\t".join(fields_to_extract))
    for var in vcf_in:

        chrom_ = var.CHROM
        pos_   = str(var.POS)
        ref_   = var.REF
        qual_  = str(var.QUAL)
        if caller == "varscan":
            qual_ = str(var.samples[0]["GQ"])

        dp_    = "0"
        ad_or_ao = "NA"
        ad_str = "NA"
        gt_ = "NA"
        alt_str = ""

        if var.samples[0].called == False:
            continue

        # usually need to bypass indels, however, homozygous REF is considered indel by pyvcf... WTF?
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
            alt_   = var.ALT[0]
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
    return var_ct

def isFileInAmazon(srcFile, config):
	cmd = [config.aws_path, "s3", "ls", srcFile]
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	p.wait()
	lines = out.splitlines()
	for line in lines:
		fields = line.split()
		if fields[3] == os.path.basename(srcFile):
			return True
	return False
	
def downloadFileFromAmazon(srcFile, destDirectory, config):
	if len(config.aws_path) == 0:
		print "AWS Path not set"
		exit(1)

	if isFileInAmazon(srcFile, config) == False:
		print "File (%s) is not accessible" % srcFile
		return None
		
	cmd = [config.aws_path, "s3", "cp", srcFile, destDirectory]
	if config.verbose:
		print "Downloading file: {}".format(srcFile)
	p = subprocess.Popen(cmd)
	p.wait()
	if p.returncode != 0:
		print "Error downloading file {}".format(srcFile)
		return None
	
	localFile = os.path.join(destDirectory, os.path.basename(srcFile))
	if os.access(localFile, os.R_OK) == False:
		print "{} is not accessible.".format(localFile)
		return None
	 
	return localFile 
	
# Get list of chromosome names from the BAM file
def get_chrom_names_from_BAM(bam_file):
    chrom_list = []
    inbam = pysam.AlignmentFile(bam_file, "rb")
    header_sq = inbam.header["SQ"]
    for sq_ in header_sq:
        chrom_list.append(sq_["SN"])
    return chrom_list

def get_chrom_names_from_REF(ref_fasta):
	ref_idx = ref_fasta + ".fai"
	chrom_list = []
	fin = open(ref_idx, "r")
	for line in fin:
		if line.strip() == "":
			continue
		chrom_list.append(line.strip().split()[0])
	return chrom_list

def get_chrom_names_from_VCF(vcf_file):
	chrom_list = []
	with open(vcf_file, "r") as fin:
		vcf_reader = vcf.Reader(fin)
		for vcfRecord in vcf_reader:
			if chrom_list.count(vcfRecord.CHROM) == 0:
				chrom_list.append(vcfRecord.CHROM)
	return chrom_list

def genotypeSample(sample, bamFile, reference, vcf, config):
	if config.verbose:
		print "Genotyping {}".format(sample)

	deleteBam = False
	outputVcf = os.path.join(config.cache_dir, sample + ".vcf")
	if os.path.exists(outputVcf) == False:
		''' Download the BAM file and index if they are not local '''
		if bamFile.startswith("s3://"):
			localBamFile = os.path.join(config.scratch_dir, os.path.basename(bamFile))
			if os.path.exists(localBamFile):
				print "Using cached bam file: {}".format(localBamFile)
			else:
				if downloadFileFromAmazon(bamFile, config.scratch_dir, config) is None:
					print "File does not exist in Amazon."
					return None
			deleteBam = True
			
			''' If the index is already downloaded, use it '''
			bamIndex1 = bamFile.rstrip(".bam") + ".bai"
			bamIndex2 = bamFile + ".bai"
			localBamIndex1 = localBamFile.rstrip(".bam") + ".bai"
			localBamIndex2 = localBamFile + ".bai"
			localBamIndex = ''
		
			if os.path.exists(localBamIndex1):
				localBamIndex = localBamIndex1
			elif os.path.exists(localBamIndex2):
				localBamIndex = localBamIndex2
			else:
				''' Else, try to download index '''
				if downloadFileFromAmazon(bamIndex1, config.scratch_dir, config):
					localBamIndex = localBamIndex1
				else:
					''' Try to download index 2 '''
					if downloadFileFromAmazon(bamIndex2, config.scratch_dir, config):
						localBamIndex = localBamIndex2
					else:
						print "Could not find matching bam index.  Generating."
						if len(config.samtools_path) == 0:
							print "samtools path not specified"
							exit(1)
						cmd = [config.samtools_path, 'index', localBamFile]
						p = subprocess.Popen(cmd)
						p.wait()
						if p.returncode != 0:
							print "Unable to generated BAM index"
							exit(1)
		else:
			localBamFile = os.path.abspath(bamFile)

		''' Make sure BAM and reference have matching chromosomes '''
		bam_chroms = get_chrom_names_from_BAM(localBamFile)
		if config.debug:
			for chr in bam_chroms:
				print "\t" + chr
				
		if config.verbose:
			print "Checking reference " + reference
		REF_CHROMS = get_chrom_names_from_REF(reference)
		if config.debug:
			for chr in REF_CHROMS:
				print "\t" + chr
		
		if set(REF_CHROMS).issubset(set(bam_chroms)) == False:
			#bamREF_diff = set(bam_chroms).difference(set(REF_CHROMS))
			bamREF_dff = set(REF_CHROMS).difference(set(bam_chroms))
			#if len(bamREF_diff) >= (len(REF_CHROMS) / 2):
			print "Sample {} contains chromosomes not in reference {}:".format(sample, reference)
			for chr in bamREF_diff:
				print chr
			exit(1)
			
		''' Make sure the VCF file and the BAM file have matching chromosome names '''
		if config.verbose:
			print "Checking VCF " + vcf
		VCF_chroms = get_chrom_names_from_VCF(vcf)
		if config.debug:
			for chr in VCF_chroms:
				print "\t" + chr
		
		if set(VCF_chroms).issubset(set(bam_chroms)) == False:
			bamVCF_diff = set(VCF_chroms).difference(set(bam_chroms))
			#if len(bamVCF_diff) >= (len(VCF_chroms) / 2):
			print "Sample {} contains chromosomes not in VCF {}".format(sample, vcf)
			for chr in bamVCF_diff:
				print "\t" + chr
			exit(1)
				
		''' Make the intervals file for the BAM file and matching reference '''
		intervalsFile = os.path.join(config.scratch_dir, sample +".intervals")		
		if config.verbose:
			print "Generating intervals file: " + intervalsFile
		vcfToIntervals(vcf, intervalsFile)

		if config.caller == 'freebayes':
			if config.verbose:
				print "Calling freebayes"
			cmd = [config.freebayes_path, "--fasta-reference", reference, "--targets", intervalsFile, "--no-indels",
				"--min-coverage", str(config.dp_threshold), "--report-all-haplotype-alleles", "--report-monomorphic", 
				"--vcf", outputVcf, localBamFile]
		
			p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			p.wait()
			if p.returncode != 0:
				print "Error executing {}".format(' '.join(cmd))
				print err
				os.remove(outputVcf)
				exit(1)
			if os.path.exists(outputVcf) == False:
				print "Output VCF file, {}, could not be found in cache.".format(outputVcf)
				exit(1)

		os.remove(intervalsFile)

	''' Convert the vcf to tsv '''
	if config.verbose:
		print "Convert VCF to TSV"
	out_tsv = os.path.join(config.cache_dir, sample + ".tsv")
	if os.path.exists(out_tsv) == False:
		VCFtoTSV(outputVcf, out_tsv, config.caller)

	if deleteBam == True:
		os.remove(localBamFile)	
		os.remove(localBamIndex)
	return None

def submitSample(sample, reference, vcf, config):
	''' Only keep the variables we need for genotyping a sample '''
	settings = { "sample": { "name": sample["name"],
						     "bam": sample["bam"]
						   }, 
				"REFERENCE": reference,
				"VCF": vcf,

				"CACHE_DIR": config.cache_dir,
				"SCRATCH_DIR": config.scratch_dir,
				"CALLER": config.caller,
				"DP_THRESHOLD": config.dp_threshold,
				"NUM_SNPS": config.number_of_snps,
				"FREEBAYES_PATH": config.freebayes_path,
				"AWS_PATH": config.aws_path,
							
				"VERBOSE": config.verbose,
				"DEBUG": config.debug,
			}
	jsonFile = os.path.join(config.cache_dir, sample["name"] + ".json")
	with open(jsonFile, 'w') as f:
		json.dump(settings, f)
	
	cmd = ["qsub", "-N", sample["name"], "-cwd", "-S", sys.executable, "-v", "JSONFILE="+jsonFile, __file__]
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	p.wait()
	if p.returncode != 0:
		print "Error executing {}".format(' '.join(cmd))
		print err
		exit(1)
	print out
	exit(0)
	
def genotypeSamples(sampleSet, reference, vcf, config):
	jobs = []
	for sampleIndex in range(len(sampleSet)):
		sample = sampleSet[sampleIndex]
		tsvFile = os.path.join(config.cache_dir, sample["name"] + ".tsv")
		if os.path.exists(tsvFile) == False:
			if config.cluster == "local":
				genotypeSample(sample["name"], sample["bam"], reference, vcf, config)
			elif config.cluster == "sge":
				jobId = submitSample(sample, reference, vcf, config)
				jobs.append(jobId)
		
def get_chrom_data_from_map(chrom_map_file):
    chrom_ct = 0
    default_chroms = []
    alternate_chroms = []
    def_to_alt = {}
    alt_to_def = {}
	
    with open(chrom_map_file, 'r') as fin:
		header = fin.readline()
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
		
def same_gt(gt1, gt2):
	if gt1 == gt2:
		return True
	else:
		gt1_ = sorted(gt1.split("/"))
		gt2_ = sorted(gt2.split("/"))
		if gt1_ == gt2_:
			return True
		else:
			return False

def is_hom(gt):
    gt_ = gt.split("/")
    if gt_[0] == gt_[1]:
        return True
    else:
        return False

def is_subset(hom_gt, het_gt):
    gt_hom = hom_gt.split("/")[0]
    if gt_hom in het_gt:
        return True
    else:
        return False

def compareSamples(sampleSet1, sampleSet2, config):
	A_BIT_LOW = """the number of comparable genomic loci is a bit low.
Try using a different variants list (--VCF) file which have more appropriate genomic positions for comparison."""

	if config.chromosome_map:
		default_chroms, alternate_chroms, def_to_alt, alt_to_def = get_chrom_data_from_map(config.chromosome_map)
	else:
		default_chroms, alternate_chroms, def_to_alt, alt_to_def = None, None, None, None

	matrix = {}
	for sample1 in sampleSet1:
		for sample2 in sampleSet2:
			if config.verbose:
				print "Comparing {} - {}".format(sample1["name"], sample2["name"])
				
			var_list = {}
			
			''' Get a list of variants that pass in sample 1 '''
			tsv1 = os.path.join(config.cache_dir, sample1["name"] + ".tsv")
			with open(tsv1, "r") as fin:
				for line in fin:	
					if line.startswith("CHROM\t"):
						continue
					bits = line.strip("\n").split("\t")
					if bits[5] == "NA":
						continue
					elif int(bits[5]) < config.dp_threshold:
						continue
					else:
						var_list["\t".join(bits[:2])] = 1
						
			''' then parse second tsv file to get list of variants that passed in both samples '''
			tsv2 = os.path.join(config.cache_dir, sample2["name"] + ".tsv")
			with open(tsv2, "r") as fin:
				for line in fin:
					if line.startswith("CHROM\t"):
						continue
					bits = line.strip("\n").split("\t")
					if alternate_chroms is not None:
						var_ = alt_to_def[bits[0]] + "\t" + bits[1]
					else:
						var_ = "\t".join(bits[:2])

					if var_ in var_list:
						var_list[var_] = 2
			
			#-------------------------------------------------------------------------------
			# write out bam1 variants
			bam1_var = os.path.join(config.scratch_dir, "sample1.variants") 
			with open(bam1_var, "w") as fout:
				with open(tsv1, "r") as fin: 
					for line in fin:
					    if line.startswith("CHROM\t"):
					        continue
					    bits = line.strip("\n").split("\t")
					    out_line = "%s\t%s\t%s\t%s\t%s\n" % (bits[0], bits[1], bits[2], bits[3], bits[7])
					    var_ = "\t".join(bits[:2])
					    if var_ in var_list:
					        if var_list[var_] == 2:
					            fout.write(out_line)
			
			#-------------------------------------------------------------------------------
			# write out bam2 variants
			bam2_var = os.path.join(config.scratch_dir, "sample2.variants") 
			with open(bam2_var, "w") as fout:
				with open(tsv2, "r") as fin: 
					for line in fin:
						if line.startswith("CHROM\t"):
							continue
						bits = line.strip("\n").split("\t")
						out_line = "%s\t%s\t%s\t%s\t%s\n" % (bits[0], bits[1], bits[2], bits[3], bits[7])
						if alternate_chroms is not None:
							var_ = alt_to_def[bits[0]] + "\t" + bits[1]
						else:
							var_ = "\t".join(bits[:2])
						if var_list.has_key(var_) and var_list[var_] == 2:
							fout.write(out_line)

			''' at this point sample1.variants and sample2.variants should have the same number of lines ie variants '''
			''' get the genotypes from each sample '''
			bam1_gt = {}
			bam2_gt = {}
			pos_list = []
			with open(bam1_var, "r") as fin:
				for line in fin:
					bits = line.strip("\n").split("\t")
					pos_ = "_".join(bits[:2])
					geno = bits[4]
					pos_list.append(pos_)
					bam1_gt[pos_] = geno
			with open(bam2_var, "r") as fin:
				for line in fin:
				    bits = line.strip("\n").split("\t")
				    if alternate_chroms is not None:
				    	pos_ = alt_to_def[bits[0]] + "_" + bits[1]
				    else:
				    	pos_ = "_".join(bits[:2])
				    geno = bits[4]
				    bam2_gt[pos_] = geno
			
			''' compare the genotypes '''
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
			for pos_ in pos_list:
				gt1 = bam1_gt[pos_]
				gt2 = bam2_gt[pos_]
				# if genotypes are the same
				if same_gt(gt1, gt2):
					ct_common += 1
					if is_hom(gt1):
						comm_hom_ct += 1
					else:
						comm_het_ct += 1
				else:
					ct_diff += 1
					# both are hom and different
					if is_hom(gt1) and is_hom(gt2):
						diff_hom_ct += 1
					# both are het and different
					elif is_hom(gt1) == False and is_hom(gt2) == False:
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
						print "WTF?"
						print gt1, gt2
			total_compared = ct_common + ct_diff
			frac_common_plus = 0
			frac_common = 0
			if total_compared > 0:
				frac_common = float(ct_common)/total_compared
				frac_common_plus = float(ct_common + max(diff_2sub1_ct, diff_1sub2_ct))/total_compared
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
			if total_compared <= 20:
				judgement = "Inconclusive: Too few loci to compare"
			elif total_compared <= 100:
				# allow for 0.90 frac_common for low loci count
				if frac_common >= 0.9 or frac_common_plus >= 0.9:
					judgement = "LIKELY SAME SOURCE: %s" % A_BIT_LOW
					short_judgement = "LIKELY SAME"
					# if there is possible allele-specific genotype subset
					if allele_subset == "1sub2" or allele_subset == "2sub1":
						sub_ = allele_subset.split("sub")[0]
						over_ = allele_subset.split("sub")[1]
						judgement += """BAM%s genotype appears to be a subset of BAM%s. Possibly BAM%s is RNA-seq data or BAM%s is contaminated.""" % (sub_, over_, sub_, over_)
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
					if allele_subset == "1sub2" or allele_subset == "2sub1":
						sub_ = allele_subset.split("sub")[0]
						over_ = allele_subset.split("sub")[1]
						judgement += """, but with possible allele specific genotype.\nBAM%s genotype appears to be a subset of BAM%s. Possibly BAM%s is RNA-seq data or BAM%s is contaminated. """ % (sub_, over_, sub_, over_)
						short_judgement += ". (BAM%s is subset of BAM%s)" % (sub_, over_)
				elif frac_common <= 0.6:
					judgement = "BAM FILES ARE FROM DIFFERENT SOURCES"
					short_judgement = "DIFFERENT"
				elif frac_common >= 0.8:
					judgement = """LIKELY FROM THE SAME SOURCE. However, the fraction of sites with common genotype is lower than expected. This can happen with samples with low coverage."""
					short_judgement = "LIKELY SAME"
				else:
					judgement = "LIKELY FROM DIFFERENT SOURCES"
					short_judgement = "LIKELY DIFFERENT"

			# so pad numeric string to 6 spaces
			diff_hom = ("%d" % diff_hom_ct).rjust(5)
			diff_het = ("%d" % diff_het_ct).rjust(5)
			diff_het_hom = ("%d" % diff_het_hom_ct).rjust(5)
			diff_hom_het = ("%d" % diff_hom_het_ct).rjust(5)
			diff_1sub2 = ("%d" % diff_1sub2_ct).rjust(5)
			diff_2sub1 = ("%d" % diff_2sub1_ct).rjust(5)
			std_report_str = """sample1:\t%s
sample2:\t%s
variants:\t%s
depth threshold: %d
________________________________________

Positions with same genotype:   %d
    breakdown:    hom: %d
                  het: %d
________________________________________

Positions with diff genotype:   %d
	 breakdown:
                         SAMPLE 1
                  | het  | hom  | subset
            -------+------+------+------
            het   |%s |%s |%s |
            -------+------+------+------
SAMPLE 2    hom   |%s |%s |   -  |
            -------+------+------+------
            subset|%s |   -  |   -  |
________________________________________

Total sites compared: %d
Fraction of common: %f (%d/%d)
________________________________________
CONCLUSION:
%s"""  % (sample1["name"], sample2["name"], config.vcf, config.dp_threshold,
					ct_common, comm_hom_ct, comm_het_ct,
					ct_diff, diff_het, diff_hom_het, diff_1sub2, diff_het_hom, diff_hom, diff_2sub1,
					total_compared, frac_common, ct_common, total_compared, judgement)
			short_report_str = """# sample1\t sample2\t DP_thresh\t FracCommon\t Same\t Same_hom\t Same_het\t Different\t 1het-2het\t 1het-2hom\t 1het-2sub\t 1hom-2het\t 1hom-2hom\t 1sub-2het\t Conclusion
			%s\t%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s""" % (sample1["name"],
				   sample2["name"], config.dp_threshold, frac_common, ct_common, comm_hom_ct, comm_het_ct,
				   ct_diff, diff_het_ct, diff_het_hom_ct, diff_2sub1_ct, diff_hom_het_ct,
				   diff_hom_ct, diff_1sub2_ct, short_judgement)
			REPORT_PATH = "%s/%s-%s" % (config.output_dir, sample1["name"], sample2["name"])
			with open(REPORT_PATH, "w") as fout:
				if config.short_output:
					fout.write(short_report_str)
					if config.verbose:
						print short_report_str
				else:
					fout.write(std_report_str)
					if config.verbose:
						print std_report_str
			''' save to grand matrix '''
			if sample1["name"] not in matrix:
				matrix[sample1["name"]] = {}
			matrix[sample1["name"]][sample2["name"]] = frac_common
	''' print out grand matrix '''
	REPORT_PATH = "./cohort-matcher-results.txt"
	with open(REPORT_PATH, "w") as fout:
		for sample1 in sampleSet1:
			fout.write("\t" + sample1["name"])
		fout.write("\n")
		for sample2 in sampleSet2:
			fout.write(sample2["name"])
			for sample1 in sampleSet1:
				s = '%.4f' % matrix[sample1["name"]][sample2["name"]]
				fout.write("\t" + s)
			fout.write("\n")
			
def main(argv = None):
	config = parseArgs(argv)
	checkConfig(config)
	
	sampleSet1 = readSamples(config.set1, config.verbose)
	sampleSet2 = readSamples(config.set2, config.verbose)

	genotypeSamples(sampleSet1, config.reference, config.vcf, config)
	if config.reference2 is None:
		genotypeSamples(sampleSet2, config.reference, config.vcf, config)
	else:
		genotypeSamples(sampleSet2, config.reference2, config.vcf2, config)	

	compareSamples(sampleSet1, sampleSet2, config)
	
if __name__ == "__main__":
	if os.environ.has_key("JOB_ID"):
		jsonConfigFile = os.environ.get("JSONFILE")
		with open(jsonConfigFile, 'r') as f:
			jsonConfig = json.load(f)
			# TODO: Convert jsonConfig to Namespace
			config = argparse.Namespace()
			config.cache_dir = jsonConfig["CACHE_DIR"]
			config.scratch_dir = jsonConfig["SCRATCH_DIR"]
			config.caller = jsonConfig["CALLER"]
			config.dp_threshold = jsonConfig["DP_THRESHOLD"] 
			config.number_of_snps = jsonConfig["NUM_SNPS"] 
			config.freebayes_path = jsonConfig["FREEBAYES_PATH"]
			config.aws_path = jsonConfig["AWS_PATH"] 
						
			config.verbose = jsonConfig["VERBOSE"]
			config.debug = jsonConfig["DEBUG"]
			
			genotypeSample(jsonConfig['sample']['name'], jsonConfig['sample']['bam'], jsonConfig["REFERENCE"], jsonConfig["VCF"], config)
			os.remove(jsonConfigFile)
	else:
		sys.exit(main())
