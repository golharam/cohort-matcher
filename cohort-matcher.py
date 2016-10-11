# Name: cohort-matcher.py
import argparse
import ConfigParser
import os
import subprocess
import sys

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
	parser_grp3.add_argument("--vcf", "-V",  required=True,
							 help="VCF file containing SNPs to check (default can be specified in config file instead)")
	
	parser_grp4 = parser.add_argument_group("Cluster")
	parser_grp4.add_argument("--cluster", required=False,
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
		print "Depth Threshold: {}".format(args.dp_threshold)
		print "Number of SNPs: {}".format(args.number_of_snps)
		print "VCF File: {}".format(args.vcf)
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
	
def genotypeSample(sample, bamFile, config):
	if config.cluster == "local":
		env = dict(os.environ)
		env['SAMPLE'] = sample
		env['BAMFILE'] = bamFile
		env['SCRATCH_DIR'] = config.scratch_dir
		
		cmd = "./genotypeSample.sh"
		if config.verbose:
			print "Genotyping {}".format(sample)

		p = subprocess.Popen(cmd, shell=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate()
		p.wait()
		p.returncode
		
		if p.returncode != 0:
			print "Error executing {}".format(cmd)
			print err
			
		return p.returncode
	'''		
	envVar = "SAMPLE=%s,BAMFILE=%s,SCRATCH_DIR=%s" % (sample, bamFile, config.scratch_dir)
	cmd = ['qsub', '-N', sample, '-v', envVar, os.path.join('./genotypeSample.sh')]
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()
	'''
	return None

def waitForSamples(sampleSet):
	return None
	
def genotypeSamples(sampleSet, config):
	for sampleIndex in range(len(sampleSet)):
		sample = sampleSet[sampleIndex]
		sampleName = sample["name"]
		tsvFile = os.path.join(config.cache_dir, sampleName + ".tsv")
		if os.path.exists(tsvFile) == False:
			sampleBam = sample["bam"]
			sample['jobId'] = genotypeSample(sampleName, sampleBam, config)

def getCompareResult(resultFile):
	return None
	
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

	genotypeSamples(sampleSet1, config)
	genotypeSamples(sampleSet2, config)
'''
	waitForSamples(sampleSet1)
	waitForSamples(sampleSet2)
	
	compareSamples(sampleSet1, sampleSet2, args.cache_dir, args.results_dir)
	
	if len(sampleSet1) == len(sampleSet2):
		doStraightLineComparison()
	
'''
if __name__ == "__main__":
	sys.exit(main())