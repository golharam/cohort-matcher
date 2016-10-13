# Name: cohort-matcher.py
import argparse
import ConfigParser
import os
import pysam
import subprocess
import sys
import vcf

def parseArgs(argv = None):
	# Do argv default this way, as doing it in the functional
	# declaration sets it at compile time.
	if argv is None:
		argv = sys.argv[1:]
	
	parser = argparse.ArgumentParser(description="Compare two vcf files to see if \
	they are from the same samples, using frequently occuring SNPs \
	reported in the 1000Genome database")

	parser_grp1 = parser.add_argument_group("Required")
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
	parser_grp7.add_argument("--aws-path", required=False, help="Specify path to aws cli")
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

	args = parser.parse_args(argv)
	
	if args.verbose:
		print "Set1: {}".format(args.set1)
		print "Set2: {}".format(args.set2)
		print ""
		print "Cache Dir: {}".format(args.cache_dir)
		print "Scratch Dir: {}".format(args.scratch_dir)
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

def downloadFileFromAmazon(srcFile, destDirectory, config):
	if len(config.aws_path) == 0:
		print "AWS Path not set"
		exit(1)

	cmd = [config.aws_path, "s3", "cp", srcFile, destDirectory]
	if config.verbose:
		print "Downloading file: {}".format(srcFile)
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
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

	''' Download the BAM file and index if they are not local '''
	deleteBam = False
	if bamFile.startswith("s3://"):
		localBamFile = os.path.join(config.scratch_dir, os.path.basename(bamFile))
		if os.path.exists(localBamFile):
			print "Using cached bam file: {}".format(localBamFile)
		else:
			downloadFileFromAmazon(bamFile, config.scratch_dir, config)
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
		localBamFile = bamFile

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
	
	outputVcf = os.path.join(config.cache_dir, sample + ".vcf")
	if os.path.exists(outputVcf) == False:
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
				print "Error executing {}".format(cmd)
				print err
				os.remove(outputVcf)
				exit(1)
			if os.path.exists(outputVcf) == False:
				print "Output VCF file, {}, could not be found in cache.".format(outputVcf)
				exit(1)

	''' Convert the vcf to tsv '''
	if config.verbose:
		print "Convert VCF to TSV"
	out_tsv = os.path.join(config.cache_dir, sample + ".tsv")
	if os.path.exists(out_tsv) == False:
		VCFtoTSV(outputVcf, out_tsv, config.caller)

	if deleteBam == True:
		os.remove(localBamFile)	
		os.remove(localBamIndex)
	os.remove(intervalsFile)
	return None

def genotypeSamples(sampleSet, reference, vcf, config):
	for sampleIndex in range(len(sampleSet)):
		sample = sampleSet[sampleIndex]
		tsvFile = os.path.join(config.cache_dir, sample["name"] + ".tsv")
		if os.path.exists(tsvFile) == False:
			if config.cluster == "local":
				if sample["bam"].startswith("s3://"):
					genotypeSample(sample["name"], sample["bam"], reference, vcf, config)
				else:
					genotypeSample(sample["name"], os.path.abspath(sample["bam"]), reference, vcf, config)

def compareSamples(sampleSet1, sampleSet2, config):
	for sample1 in sampleSet1:
		for sample2 in sampleSet2:
			var_list = {}
			tsv1 = os.path.join(config.cache_dir, sample1["name"] + ".tsv")
			with open(tsv1, "r") as fin:
				for line in fin:	
				    if line.startswith("CHROM\t"):
				        continue
				    bits = line.strip("\n").split("\t")
				    if bits[5] == "NA":
				        continue
				    elif int(bits[5]) < DP_THRESH:
				        continue
				    else:
				        # add variants to list
				        var_list["\t".join(bits[:2])] = 1
			
			# then parse second tsv file to get list of variants that passed in both bams
			tsv2 = os.path.join(config.cache_dir, sample2["name"] + ".tsv")
			with open(tsv2, "r") as fin:
				for line in fin:
				    if line.startswith("CHROM\t"):
				        continue
				    bits = line.strip("\n").split("\t")
				    var_ = "\t".join(bits[:2])
				
				    if var_ in var_list:
				        var_list[var_] = 2
			
			#-------------------------------------------------------------------------------
			# write out bam1 variants
			bam1_var = os.path.join(config.scratch_dir, sample1["name"] + ".variants") 
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
			bam2_var = os.path.join(config.scratch_dir, sample2["name"] + ".variants") 
			with open(bam2_var, "w") as fout:
				with open(tsv2, "r") as fin: 
					for line in fin:
					    if line.startswith("CHROM\t"):
					        continue
					    bits = line.strip("\n").split("\t")
					    out_line = "%s\t%s\t%s\t%s\t%s\n" % (bits[0], bits[1], bits[2], bits[3], bits[7])
					    var_ = "\t".join(bits[:2])
					    if var_ in var_list:
					        if var_list[var_] == 2:
					            fout.write(out_line)
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
	sys.exit(main())
