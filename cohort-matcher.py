# Name: cohort-matcher.py
import argparse
import ConfigParser
import os
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
	parser_grp4.add_argument("--reference", "-R", required=True,
							help="Reference FASTA File (indexed with samtools faidx)")
	
	parser_grp5 = parser.add_argument_group("Freebayes Settings")
	parser_grp5.add_argument("--freebayes-path", required=False,
							default="freebayes",
							help="Path to freebayes binary (if not in PATH")
							
	parser_grp6 = parser.add_argument_group("Cluster")
	parser_grp6.add_argument("--cluster", required=False,
							default="local",
							choices=('local', 'sge'),
							help="Specify type of cluster architecture to use (Default: local)")
	

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
		print "Reference: {}".format(args.reference)
		print ""
		print "Cluster Type: {}".format(args.cluster)
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

	if os.path.exists(config.reference) == False:
		print "Reference FASTQ file (%s) is not accessible" % config.reference
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

def genotypeSample(sample, bamFile, intervalsFile, config):
	if config.verbose:
		print "Genotyping {}: {}".format(sample, bamFile)
		
	''' Make sure BAM file is accessible '''
	if os.access(bamFile, os.R_OK) == False:
		print "{} is not accessible.".format(bamFile)
		exit(1)
		
	''' Make sure BAM file is indexed '''
	bam_index2 = bamFile.rstrip(".bam") + ".bai"
	bam_index1 = bamFile + ".bai"
	if (os.access(bam_index1, os.R_OK) == False and
		os.access(bam_index2, os.R_OK) == False):
		print "BAM file (%s) is either missing index or the index file is not readable." % bamFile
		exit(1)
	
	''' TODO: Make sure BAM header and reference have matching chromosome names '''
	
	''' TODO: Make sure the VCF file and the BAM file have matching chromosome names '''
	
	outputVcf = os.path.join(config.cache_dir, sample + ".vcf")
	if os.path.exists(outputVcf) == False:
		if config.caller == 'freebayes':
			cmd = [config.freebayes_path, "--fasta-reference", config.reference, "--targets", intervalsFile, "--no-indels",
				"--min-coverage", str(config.dp_threshold), "--report-all-haplotype-alleles", "--report-monomorphic", 
				"--vcf", outputVcf, bamFile]
		
			p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			p.wait()
				
			if p.returncode != 0:
				print "Error executing {}".format(cmd)
				print err
				return p.returncode
			
			if os.path.exists(outputVcf) == False:
				print "Output VCF file, {}, could not be found in cache.".format(outputVcf)
				return p.returncode
		
	''' Convert the vcf to tsv '''
	out_tsv = os.path.join(config.cache_dir, sample + ".tsv")
	if os.path.exists(out_tsv) == False:
		VCFtoTSV(outputVcf, out_tsv, config.caller)
	return None
	
def genotypeSamples(sampleSet, intervalsFile, config):
	for sampleIndex in range(len(sampleSet)):
		sample = sampleSet[sampleIndex]
		tsvFile = os.path.join(config.cache_dir, sample["name"] + ".tsv")
		if os.path.exists(tsvFile) == False:
			if config.cluster == "local":
				genotypeSample(sample["name"], os.path.abspath(sample["bam"]), intervalsFile, config)

def compareSamples(sampleSet1, sampleSet2, cache_dir, results_dir):
	compareJobs = []
	for sampleIndex in range(len(sampleSet1)):
		sample1 = sampleSet1[sampleIndex]
		sample2 = sampleSet2[sampleIndex]
		result = os.path.join(results_dir, sample1) + "-" + os.path.join(results_dir, sample2)
		if os.path.exist(result) == False:
			compareJobs.append(compareGenotypes(sample1, sample2))

	waitForCompareJobs(compareJobs)
	
	for sampleIndex in range(len(sampleSet1)):
		sample1 = sampleSet1[sampleIndex]
		sample2 = sampleSet2[sampleIndex]
		result = os.path.join(results_dir, sample1) + "-" + os.path.join(results_dir, sample2)
		pctIdentity = getCompareResult(result)
		print sample1 + "\t" + sample2 + "\t" + pctIdentity
		
def main(argv = None):
	config = parseArgs(argv)
	checkConfig(config)
	
	sampleSet1 = readSamples(config.set1, config.verbose)
	sampleSet2 = readSamples(config.set2, config.verbose)

	intervalsFile = os.path.join(config.scratch_dir, "genotype.intervals")
	if config.verbose:
		print "Generating intervals file: " + intervalsFile
	vcfToIntervals(config.vcf, intervalsFile)

	genotypeSamples(sampleSet1, intervalsFile, config)
	genotypeSamples(sampleSet2, intervalsFile, config)

	compareSamples(sampleSet1, sampleSet2, args.cache_dir, args.results_dir)
	
if __name__ == "__main__":
	sys.exit(main())