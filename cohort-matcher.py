import argparse
import os
import random
import string
import sys
from bammatcher_methods import *

'''
Name: cohort-matcher

Desc: Check two sets of bam files for genotype concordance
Assumptions:
1. BAM files are full URLS and are named as sample.*bam (ie filename up to
   first period is the sample name)
2. 

Steps:
1.  Get list of set1 bam files
2.  Get list of set2 bam files
3.  Read configuration information
4.  Verify configuration and variant caller
5.  For each bam file in set1 and set2:
      If vcf is not cache:
5a.      run variant caller to generate vcf file (by qsubbing call)
6.  Wait for all jobs to finish.
7.  Check that all vcf files are present.
8.  Make matrix of concordance
'''
def parseArguments(argv):
    parser = argparse.ArgumentParser(description="Cohort matching pipeline")

    parser_grp1 = parser.add_argument_group("REQUIRED")
    parser_grp1.add_argument("--bamlist1", action="store", required=True, dest="bamlist1", help="File containing list of bam files in set 1")
    parser_grp1.add_argument("--bamlist2", action="store", required=True, dest="bamlist2", help="File containing list of bam files in set 2")

    parser_grp2 = parser.add_argument_group("CONFIGURATION")
    parser_grp2.add_argument("--config", "-c", required=False, default='', help="Specify configuration file (default = /dir/where/script/is/located/bam-matcher.conf)")
    parser_grp2.add_argument("--generate-config", "-G", required=False, nargs="?", action="store", default=None, const='', dest="generateConfigFile", help="Specify where to generate configuration file template")

    parser_grp3 = parser.add_argument_group("OUTPUT REPORT")
    parser_grp3.add_argument("--output", "-o", required=False, help="Specify output report path (default = /current/dir/bam_matcher.SUBFIX)")
    parser_grp3.add_argument("--short-output","-so", required=False, default = False, action="store_true", help="Short output mode (tab-separated).")
    parser_grp3.add_argument("--html", "-H",  required=False, action="store_true", help="Enable HTML output. HTML file name = report + '.html'")
    parser_grp3.add_argument("--no-report", "-n",  required=False, action="store_true", help="Don't write output to file. Results output to command line only.")
    parser_grp3.add_argument("--scratch-dir", "-s",  required=False, help="Scratch directory for temporary files. If not specified, the report output directory will be used (default = /tmp/[random_string])")

    # Variants VCF file
    parser_grp4 = parser.add_argument_group("VARIANTS")
    parser_grp4.add_argument("--vcf", "-V",  required=False, default=None, help="VCF file containing SNPs to check (default can be specified in config file instead)")

    # caller settings - setting these will override config
    parser_grp5 = parser.add_argument_group("CALLERS AND SETTINGS (will override config values)")
    parser_grp5.add_argument("--caller", "-CL", required=False, choices=('gatk', 'freebayes', 'varscan'), help="Specify which caller to use (default = 'gatk')")
    parser_grp5.add_argument("--dp-threshold", "-DP", required=False, type=int, help="Minimum required depth for comparing variants")
    parser_grp5.add_argument("--number_of_snps", "-N", required=False, type=int, help="Number of SNPs to compare.")
    parser_grp5.add_argument("--fastfreebayes",  "-FF", required=False, action="store_true", help="Use --targets option for Freebayes.")
    parser_grp5.add_argument("--gatk-mem-gb" ,   "-GM", required=False, type=int, help="Specify Java heap size for GATK (GB, int)")
    parser_grp5.add_argument("--gatk-nt" ,       "-GT", required=False, type=int, help="Specify number of threads for GATK UnifiedGenotyper (-nt option)")
    parser_grp5.add_argument("--varscan-mem-gb", "-VM", required=False, type=int, help="Specify Java heap size for VarScan2 (GB, int)")

    # Genome references
    parser_grp6 = parser.add_argument_group("REFERENCES")
    parser_grp6.add_argument("--reference", "-R",  required=False, help="Default reference fasta file. Needs to be indexed with samtools faidx. Overrides config settings.")
    parser_grp6.add_argument("--ref-alternate", "-R2",  required=False, help="Alternate reference fasta file. Needs to be indexed with samtools faidx. Overrides config settings.")
    parser_grp6.add_argument("--chromosome-map", "-M", required=False, help="Required when using alternate reference. Run BAM-matcher with --about-alternate-ref for more details.")
    parser_grp6.add_argument("--about-alternate-ref", "-A", required=False, action="store_true", default="False", help="Print information about using --alternate-ref and --chromosome-map")

    # for batch operations
    parser_grp7 = parser.add_argument_group("BATCH OPERATIONS")
    parser_grp7.add_argument("--do-not-cache", "-NC", required=False, default=False, action="store_true", help="Do not keep variant-calling output for future comparison. By default (False) data is written to /bam/filepath/without/dotbam.GT_compare_data")
    parser_grp7.add_argument("--recalculate", "-RC", required=False, default=False, action="store_true", help="Don't use cached variant calling data, redo variant-calling. Will overwrite cached data unless told not to (-NC)")
    parser_grp7.add_argument("--cache-dir", "-CD", required=False, dest="cache_dir", help="Specify directory for cached data. Overrides configuration")

    # Experimental features
    # parser_grp8 = parser.add_argument_group("EXPERIMENTAL")
    # parser_grp8.add_argument("--allele-freq", "-AF", required=False, default=False, action="store_true", help="Plot variant allele frequency graphs")

    # optional, not in config
    parser.add_argument("--debug", "-d", required=False, action="store_true", help="Debug mode. Temporary files are not removed")
    parser.add_argument("--verbose", "-v", required=False, action="store_true", help="Verbose reporting. Default = False")
    return parser.parse_args(argv)

def parseArguments2(argv):
    parser = argparse.ArgumentParser(description="Cohort matching pipeline")
    
    parser_grp1 = parser.add_argument_group("REQUIRED")
    parser_grp1.add_argument("--bamlist1", action="store", required=False, dest="bamlist1", help="File containing list of bam files in set 1")
    parser_grp1.add_argument("--bamlist2", action="store", required=False, dest="bamlist2", help="File containing list of bam files in set 2")
    
    parser_grp2 = parser.add_argument_group("CONFIGURATION")
    parser_grp2.add_argument("--config", "-c", required=False, default='', help="Specify configuration file (default = /dir/where/script/is/located/bam-matcher.conf)")
    parser_grp2.add_argument("--generate-config", "-G", required=False, nargs="?", action="store", default=None, const='', dest="generateConfigFile", help="Specify where to generate configuration file template")

    parser_grp3 = parser.add_argument_group("OUTPUT REPORT")
    parser_grp3.add_argument("--output", "-o", required=False, help="Specify output report path (default = /current/dir/bam_matcher.SUBFIX)")
    parser_grp3.add_argument("--short-output","-so", required=False, default = False, action="store_true", help="Short output mode (tab-separated).")
    parser_grp3.add_argument("--html", "-H",  required=False, action="store_true", help="Enable HTML output. HTML file name = report + '.html'")
    parser_grp3.add_argument("--no-report", "-n",  required=False, action="store_true", help="Don't write output to file. Results output to command line only.")
    parser_grp3.add_argument("--scratch-dir", "-s",  required=False, help="Scratch directory for temporary files. If not specified, the report output directory will be used (default = /tmp/[random_string])")

    # Variants VCF file
    parser_grp4 = parser.add_argument_group("VARIANTS")
    parser_grp4.add_argument("--vcf", "-V",  required=False, default=None, help="VCF file containing SNPs to check (default can be specified in config file instead)")

    # caller settings - setting these will override config
    parser_grp5 = parser.add_argument_group("CALLERS AND SETTINGS (will override config values)")
    parser_grp5.add_argument("--caller", "-CL", required=False, choices=('gatk', 'freebayes', 'varscan'), help="Specify which caller to use (default = 'gatk')")
    parser_grp5.add_argument("--dp-threshold", "-DP", required=False, type=int, help="Minimum required depth for comparing variants")
    parser_grp5.add_argument("--number_of_snps", "-N", required=False, type=int, help="Number of SNPs to compare.")
    parser_grp5.add_argument("--fastfreebayes",  "-FF", required=False, action="store_true", help="Use --targets option for Freebayes.")
    parser_grp5.add_argument("--gatk-mem-gb" ,   "-GM", required=False, type=int, help="Specify Java heap size for GATK (GB, int)")
    parser_grp5.add_argument("--gatk-nt" ,       "-GT", required=False, type=int, help="Specify number of threads for GATK UnifiedGenotyper (-nt option)")
    parser_grp5.add_argument("--varscan-mem-gb", "-VM", required=False, type=int, help="Specify Java heap size for VarScan2 (GB, int)")

    # Genome references
    parser_grp6 = parser.add_argument_group("REFERENCES")
    parser_grp6.add_argument("--reference", "-R",  required=False, help="Default reference fasta file. Needs to be indexed with samtools faidx. Overrides config settings.")
    parser_grp6.add_argument("--ref-alternate", "-R2",  required=False, help="Alternate reference fasta file. Needs to be indexed with samtools faidx. Overrides config settings.")
    parser_grp6.add_argument("--chromosome-map", "-M", required=False, help="Required when using alternate reference. Run BAM-matcher with --about-alternate-ref for more details.")
    parser_grp6.add_argument("--about-alternate-ref", "-A", required=False, action="store_true", default="False", help="Print information about using --alternate-ref and --chromosome-map")

    # for batch operations
    parser_grp7 = parser.add_argument_group("BATCH OPERATIONS")
    parser_grp7.add_argument("--do-not-cache", "-NC", required=False, default=False, action="store_true", help="Do not keep variant-calling output for future comparison. By default (False) data is written to /bam/filepath/without/dotbam.GT_compare_data")
    parser_grp7.add_argument("--recalculate", "-RC", required=False, default=False, action="store_true", help="Don't use cached variant calling data, redo variant-calling. Will overwrite cached data unless told not to (-NC)")
    parser_grp7.add_argument("--cache-dir", "-CD", required=False, dest="cache_dir", help="Specify directory for cached data. Overrides configuration")

    # Experimental features
    # parser_grp8 = parser.add_argument_group("EXPERIMENTAL")
    # parser_grp8.add_argument("--allele-freq", "-AF", required=False, default=False, action="store_true", help="Plot variant allele frequency graphs")
    
    # optional, not in config
    parser.add_argument("--debug", "-d", required=False, action="store_true", help="Debug mode. Temporary files are not removed")
    parser.add_argument("--verbose", "-v", required=False, action="store_true", help="Verbose reporting. Default = False")
    return parser.parse_args(argv)

def checkArguments(args):
    if os.access(args.bamlist1, os.R_OK) == False:
        print "Unable to read bamlist1"
        return False
    if os.access(args.bamlist2, os.R_OK) == False:
        print "Unable to read bamlist2"
        return False
    return True

def generateConfig(config_template_output):
    if config_template_output == '':
        config_template_output = 'bam-matcher.conf.template'
    
    config_template_output = os.path.abspath(config_template_output)    
    print """
====================================
Generating configuration file template
====================================

Config template will be written to %s

""" % config_template_output

    # make sure that it doesn't overwrite anything!
    if os.path.isfile(config_template_output):
        print "%s\nThe specified path ('%s') for config template exists already." % (FILE_ERROR, config_template_output)
        print "Write to another file."
        return 1
    fout = open(config_template_output, "w")
    fout.write(CONFIG_TEMPLATE_STR)
    fout.close()
    return 0

def readConfiguration(args):
    if args.config == '':
        SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
        config_file = os.path.join(SCRIPT_DIR, "bam-matcher.conf")
    else:
        config_file = os.path.abspath(os.path.expanduser(args.config))
        
    settings = dict()
    settings['CONFIG_FILE'] = config_file
    
    if args.verbose:
        print "Using config file: %s" % config_file
        
    # Test if the config file exists
    if check_file_read(config_file, "config", CONFIG_ERROR) == False: exit(1)

    # Load config and check if all necessary bits are available
    config = ConfigParser.ConfigParser()
    try:
        config.read(config_file)
    except ConfigParser.Error as e:
        print "%s\nUnspecified configuration file error. Please check configuration file.\n\nPython error msg:\n%s" % (CONFIG_ERROR, e)
        exit(1)

    # are all the sections there?
    config_sections = config.sections()
    REQUIRED_CONFIG_SECTIONS = ["GenomeReference", "VariantCallerParameters",
                                "ScriptOptions", "VariantCallers", "Miscellaneous",
                                "BatchOperations"]
    for sect in REQUIRED_CONFIG_SECTIONS:
        if sect not in config_sections:
            print """%s\n    Missing required section in config file: %s\n    """ % (CONFIG_ERROR, sect)
            exit(1)
            
    settings['CACHE_DIR'] = fetch_config_value(config, "BatchOperations", "CACHE_DIR")
    if args.cache_dir:
        settings['CACHE_DIR'] = args.cache_dir 
    
    #-------------------------------------------------------------------------------
    # setting variables using the config file
    settings['CALLER']         = fetch_config_value(config, "VariantCallers", "caller")
    settings['GATK']           = fetch_config_value(config, "VariantCallers", "GATK")
    settings['FREEBAYES']      = fetch_config_value(config, "VariantCallers", "freebayes")
    settings['SAMTOOLS']       = fetch_config_value(config, "VariantCallers", "samtools")
    settings['VARSCAN']        = fetch_config_value(config, "VariantCallers", "varscan")
    settings['JAVA']           = fetch_config_value(config, "VariantCallers", "java")
    settings['DP_THRESH'] = fetch_config_value(config, "ScriptOptions", "DP_threshold")
    settings['NUMBER_OF_SNPS'] = fetch_config_value(config, "ScriptOptions", "number_of_SNPs")
    settings['FAST_FREEBAYES'] = fetch_config_value(config, "ScriptOptions", "fast_freebayes")
    settings['VCF_FILE']       = fetch_config_value(config, "ScriptOptions", "VCF_file")
    settings['REFERENCE']      = fetch_config_value(config, "GenomeReference", "REFERENCE")
    settings['REF_ALTERNATE']  = fetch_config_value(config, "GenomeReference", "REF_ALTERNATE")
    settings['CHROM_MAP']      = fetch_config_value(config, "GenomeReference", "CHROM_MAP")
    settings['GATK_MEM']       = fetch_config_value(config, "VariantCallerParameters", "GATK_MEM")
    settings['GATK_NT']        = fetch_config_value(config, "VariantCallerParameters", "GATK_nt")
    settings['VARSCAN_MEM']    = fetch_config_value(config, "VariantCallerParameters", "VARSCAN_MEM")
    settings['SCRATCH_DIR']    = fetch_config_value(config, "Miscellaneous", "SCRATCH_DIR")
    
    BATCH_RECALCULATE = False
    BATCH_USE_CACHED  = True
    BATCH_WRITE_CACHE = True
       
    # if argument, use argument, else use config
    if args.caller: 
        settings['CALLER'] = args.caller
    # if not set in argument or config, default to Freebayes
    if settings['CALLER'] == '':
        settings['CALLER'] = "freebayes"
        print "No default caller was specified in the configuration file nor at runtime.\nWill default to Freebayes.\n"
    elif settings['CALLER'] not in ["gatk", "freebayes", "varscan"]:
        print "Incorrect caller specified.\nThe only values accepted for the caller parameter are: 'gatk', 'freebayes', and 'varscan'"
        exit(1)

    if args.vcf:
        settings['VCF_FILE'] = os.path.abspath(os.path.expanduser(args.vcf))
    if settings['VCF_FILE'] == "":
        print "No variants file (VCF) has been specified.\nUse --vcf/-V at command line or VCF_FILE in the configuration file."
        exit(1)    
    settings['VCF_FILE'] = os.path.abspath(os.path.expanduser(settings['VCF_FILE']))
    if check_file_read(settings['VCF_FILE'], "variants VCF", CONFIG_ERROR) == False: exit(1)
   
    if args.dp_threshold != None:
        settings['DP_THRESH'] = args.dp_threshold 
    if settings['DP_THRESH'] == "":
        print "DP_threshold value was not specified in the config file or arguments.\nDefault value (15) will be used instead."
        settings['DP_THRESH'] = 15
    try:
        settings['DP_THRESH'] = int(settings['DP_THRESH'])
    except ValueError, e:
        print "DP_threshold value ('%s') in config file is not a valid integer." % settings['DP_THRESH'] 
        exit(1)

    if args.number_of_snps != None:
        settings['NUMBER_OF_SNPS'] = args.number_of_snps
    # get from config file
    else:
        if settings['NUMBER_OF_SNPS'] == "":
            settings['NUMBER_OF_SNPS'] = 0
        else:
            try:
                settings['NUMBER_OF_SNPS'] = int(settings['NUMBER_OF_SNPS'])
            except ValueError, e:
                print "number_of_SNPs value ('%s') in config file is not a valid integer." % settings['NUMBER_OF_SNPS']
                sys.exit(1)    
    if settings['NUMBER_OF_SNPS'] <= 200 and settings['NUMBER_OF_SNPS'] > 0:
        print "Using fewer than 200 SNPs is not recommended, may not be sufficient to correctly discriminate between samples."
    
    # Fast Freebayes
    # only check if actually using Freebayes for variant calling
    if settings['CALLER'] == "freebayes":
        # get from command line
        if args.fastfreebayes:
            settings['FAST_FREEBAYES'] = args.fastfreebayes
        # get from config file
        else:
            if settings['FAST_FREEBAYES'] == "":
                print "fast_freebayes was not set in the configuration file.\nDefault value (False) will be used instead."
                settings['FAST_FREEBAYES'] = False
            else:
                if settings['FAST_FREEBAYES'] in ["False", "false", "F", "FALSE"]:
                    settings['FAST_FREEBAYES'] = False
                elif settings['FAST_FREEBAYES'] in ["True", "true", "T", "TRUE"]:
                    settings['FAST_FREEBAYES'] = True
                else:
                    print "Invalid value ('%s') was specified for fast_freebayes in the configuration file.\nUse 'False' or 'True'" % settings['FAST_FREEBAYES']
                    sys.exit(1)
    
    if settings['REFERENCE'] == "" and args.reference == None:
        print """%s
                 No genome reference has been specified anywhere.
                 Need to do this in either the configuration file or at run time (--reference/-R).
    
                 ALTERNATE_REF (in config) or --alternate_ref/-A should only be used if there is
                 already a default genome reference speficied.
              """ % CONFIG_ERROR
        exit(1)
    
    # overriding config REFERENCE
    if args.reference != None:
        settings['REFERENCE'] = os.path.abspath(os.path.expanduser(args.reference))
        print """%s
                 --reference/-R argument overrides config setting (REFERENCE).
                 Default reference = %s
              """ % (WARNING_MSG, settings['REFERENCE'])

    # overriding config REF_ALTERNATE
    if args.ref_alternate != None:
        settings['REF_ALTERNATE'] = os.path.abspath(os.path.expanduser(args.ref_alternate))
        print """%s
                 --ref-alternate/-R2 argument overrides config setting (REF_ALTERNATE).
                 Alternate reference = %s
              """ % (WARNING_MSG, settings['REF_ALTERNATE'])

    # overriding config CHROM_MAP
    if args.chromosome_map != None:
        settings['CHROM_MAP'] = os.path.abspath(os.path.expanduser(args.chromosome_map))
        print """%s
                 --chromosome-map/-M argument overrides config setting (CHROM_MAP).
                 Chromosome map = %s
              """ % (WARNING_MSG, settings['CHROM_MAP'])

    # check that if ref_alternate is used, then chromosome_map is also supplied
    if settings['REF_ALTERNATE']:
        if not settings['CHROM_MAP']:
            print """%s
                     When using an alternate genome reference (--ref-alternate), chromosome map
                     (--chromosome-map/-M or CHROM_MAP in config) must also be supplied.
                     For more details, run BAM-matcher with --about-alternate-ref/-A.
                  """ % (CONFIG_ERROR)
            exit(1)
    else:
        # check REFERENCE
        if check_file_read(settings['REFERENCE'], "default reference", FILE_ERROR) == False: exit(1)
        if not check_fasta_index(settings['REFERENCE']): exit(1)

    return settings

def printConfigurationSettings(settings):
    print "\n\nCONFIG SETTINGS"
    print "VCF file:                   %s" % settings['VCF_FILE']
    print "DP_threshold:               %d" % settings['DP_THRESH']
    print "number_of_SNPs:             %d (if 0, all variants in VCF file will be used)" % settings['NUMBER_OF_SNPS']
    print "Caller:                     %s" % settings['CALLER']
    if settings['CALLER'] == "freebayes":
        print "path to caller:            ", settings['FREEBAYES']
    if settings['CALLER'] == "freebayes":
        print "fast_freebayes:            ", settings['FAST_FREEBAYES']

    print "default genome reference:   %s" % settings['REFERENCE'] 
    print "alternate genome reference: %s" % settings['REF_ALTERNATE']

    print "cache directory:            %s" % settings['CACHE_DIR']    

def checkVariantCaller(args, settings):
    if args.verbose:
        print "\n---------------\nChecking caller\n---------------"
        
    CALLER = settings['CALLER']
    JAVA = settings['JAVA']
    SCRATCH_DIR = settings['SCRATCH_DIR']
    GATK = settings['GATK']
    FREEBAYES = settings['FREEBAYES']
    VARSCAN = settings['VARSCAN']
    
    if CALLER == "freebayes" and not JAVA:
        print "%s\nJava command was not specified.\nDo this in the configuration file" % CONFIG_ERROR
        sys.exit(1)
    caller_check_log = os.path.join(SCRATCH_DIR, "caller_check.log")
    if CALLER == "gatk":
        check_caller(CALLER, GATK, JAVA, args.verbose, logfile=caller_check_log)
    elif CALLER == "freebayes":
        check_caller(CALLER, FREEBAYES, JAVA, args.verbose, logfile=caller_check_log)
    elif CALLER == "varscan":
        check_caller(CALLER, VARSCAN,   JAVA, args.verbose, logfile=caller_check_log, SAMTL=SAMTOOLS)
    if args.verbose:
        print "Caller settings seem okay.\n"

def generateVcfs(args, settings):
    if args.verbose:
        print "Using cache dir: %s" % settings['CACHE_DIR']
    
    jobScript = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'genotypeSample.sh');
    
    for bamlist in [args.bamlist1, args.bamlist2]:
        with open(bamlist, 'r') as f:
            for line in f:
                line = line.rstrip('\r\n')
    
                if len(line) == 0 or line.startswith('#'):
                    continue
                  
                sampleBamFile = os.path.basename(line)
                sample = sampleBamFile.split('.')[0]
                sampleVcf = "%s/%s.vcf" % (settings['CACHE_DIR'], sample)
                if os.path.exists(sampleVcf) == False:
                    print "qsub -N %s -v CONFIG_FILE=%s,BAMFILE=%s %s" % (sample, settings['CONFIG_FILE'], line, jobScript)
    return None

def downloadFile(bamFile, destDir):
    destPathName = "%s/%s" % (destDir, os.path.basename(bamFile))
    if os.path.exists(destPathName):
        return destPathName
    
    if bamFile.startswith("s3://"):
        cmd = ['aws', 's3' ,'cp', bamFile, destPathName]
        print cmd
        retVal = subprocess.call(cmd)
        if retVal != 0:
            return None
        if os.path.exists(destPathName):
            return destPathName
        else:
            return None

def checkReference(bamFile, settings):
    # TODO:
    return
        
def generateIntervalFile(sampleName, settings):
    intervalListFile = os.path.join(settings['SCRATCH_DIR'], sampleName  + ".intervals")
    if os.path.exists(intervalListFile) == True and os.path.getsize(intervalListFile) > 0:
        return intervalListFile
    convert_vcf_to_intervals(settings['VCF_FILE'], intervalListFile, 0, settings['NUMBER_OF_SNPS'], settings['CALLER'])
    if os.path.getsize(intervalListFile) == 0:
        print """%s
                No intervals were extracted from variants list.
                Genotype calling have no targets and will either fail or generate an empty VCF file.
                
                Check:
                1. input VCF file (whose genomic position format should match the DEFAULT genome reference fasta),
                2. default genome reference fasta
                3. alternate genome reference fasta if it is being used
                4. the chromosome map file if it is being used.
                
                Input VCF file:             %s
                Default genome reference:   %s
                Alternate genome reference: %s
                chromosome map:             %s
                """ % (CONFIG_ERROR, settings['VCF_FILE'], settings['REFERENCE'], settings['REF_ALTERNATE'], settings['CHROM_MAP'])
        exit(1)    
    return intervalListFile

def genotypeSample(configFile, remoteBamFile):
    argv = ['-c', configFile]
    args = parseArguments2(argv)
    settings = readConfiguration(args)
    printConfigurationSettings(settings)
    bamFile = downloadFile(remoteBamFile, settings['SCRATCH_DIR'])
    if bamFile == None:
       return

    remoteBaiFile1 = remoteBamFile.replace(".bam", ".bai")
    remoteBaiFile2 = remoteBamFile + ".bai"
    baiFile1 = os.path.join(settings['SCRATCH_DIR'], os.path.basename(remoteBaiFile1))
    baiFile2 = os.path.join(settings['SCRATCH_DIR'], os.path.basename(remoteBaiFile2))
    if os.path.exists(baiFile1) == False and os.path.exists(baiFile2) == False:
        baiFile = downloadFile(remoteBaiFile1, settings['SCRATCH_DIR'])
        if baiFile == None:
            baiFile = downloadFile(remoteBaiFile2, settings['SCRATCH_DIR'])
            if baiFile == None:
                return
        
    sampleBamFile = os.path.basename(bamFile)
    sampleName = sampleBamFile.split('.')[0]
    
    checkReference(bamFile, settings)
    intervalListFile = generateIntervalFile(sampleName, settings)

    out_vcf = os.path.join(settings['CACHE_DIR'], sampleName + ".vcf")
    caller_log_file = os.path.join(settings['SCRATCH_DIR'], "%s_caller.log" % sampleName)
    caller_log = open(caller_log_file, "w")
    
    if settings['CALLER'] == "freebayes":
        fout = open(out_vcf, "w")

        # fast-Freebayes, single intervals file
        if settings['FAST_FREEBAYES']:
            varcall_cmd = [settings['FREEBAYES'], "--fasta-reference", settings['REFERENCE'], "--targets",
                           intervalListFile, "--no-indels", "--min-coverage", str(settings['DP_THRESH'])]
            varcall_cmd += ["--report-all-haplotype-alleles", "--report-monomorphic", bamFile]
            print "Freebayes variant-calling command:\n" + " ".join(varcall_cmd)
 
            varcall_proc = subprocess.Popen(varcall_cmd, stdout=fout, stderr=caller_log)
            varcall_proc.communicate()
            varcall_proc_returncode = varcall_proc.returncode
            fout.close()
            caller_log.close()

            # check calling was successful
            if varcall_proc_returncode != 0:
                print_caller_failure_message(" ".join(varcall_cmd), caller_log_file)
                print "\n Variant calling failed."
                exit(1)
            else:
                print "\nVariant calling successful."
        else:
           # -----------------------
            # slow Freebayes,calling each site separately
            write_header = True
            fin = open(intervalListFile, "r")
            for line in fin:
                bits = line.strip("\n").split("\t")
                bits[1] = int(bits[1])
                bits[2] = int(bits[2]) + 1
                region_str = "%s:%d-%d" % (bits[0], bits[1], bits[2])
                varcall_cmd = [settings['FREEBAYES'], "--fasta-reference", settings['REFERENCE'], "--region",
                               region_str, "--no-indels", "--min-coverage", str(settings['DP_THRESH'])]
                varcall_cmd += ["--report-all-haplotype-alleles", "--report-monomorphic", bamFile]
                print "Freebayes variant-calling command:\n" + " ".join(varcall_cmd)

                varcall_proc_returncode = subprocess.call(varcall_cmd)

            # if a call failed
            if varcall_proc_returncode != 0:
                print "\n Variant calling failed."
                exit(1)
            else:
                print "\nVariant calling successful."

    # Clean up
    os.remove(bamFile)
    os.remove(baiFile)
    os.remove(intervalListFile)
    
def main(argv):
    ''' 
    If we are running with no command-line parameters, its possible we were called from qsub.
    If we are qsubbed, then a few environment variables will be set, so let's check them.
    '''
    if len(argv) == 0 and os.environ.get('CONFIG_FILE') and os.environ.get('BAMFILE'):
        genotypeSample(os.environ.get('CONFIG_FILE'), os.environ.get('BAMFILE'))
        return
    
    args = parseArguments(argv)
    if args.generateConfigFile != None:
        return generateConfig(args.generateConfigFile)
    if args.about_alternate_ref == True:
        print ABOUT_ALTERNATE_REF_MSG
        return
    if checkArguments(args) == False:
        return
    settings = readConfiguration(args)
    checkVariantCaller(args, settings)
    generateVcfs(args, settings)
    
if __name__ == '__main__':
    main(sys.argv[1:])