import unittest
import logging
from mock import patch, MagicMock, mock
from tempfile import NamedTemporaryFile
import os
import filecmp
from cohort_matcher import checkConfig, compareGenotypes, is_same_gt, main, parseArguments, \
     readSamples, vcfToIntervals, genotypeSamples, genotypeSample, \
     compareSamples, get_tsv_variants
import cohort_matcher

class TestCohortMatcher(unittest.TestCase):
    @patch('os.path.isdir')
    @patch('os.path.exists')
    def test_checkConfig(self, mock_exists, mock_isdir):
        # Set up test case
        config = MagicMock(log_level=logging.INFO, reference='', reference2=None,
                           vcf2=None, chromosome_map=None)
        # Set up supporting mocks
        # Test
        retVal = checkConfig(config)
        # Check results
        self.assertTrue(retVal)

    def test_compareGenotypes(self):
        # Set up test case
        # chr1:1 hom same
        # chr1:2 het same
        # chr1:3 hom diff
        # chr1:4 het diff
        # chr1:5 hom vs het diff
        # chr1:6 het vs hom diff
        # chr1:7 hom vs het subset
        # chr1:8 het vs hom subset
        var_list = {'chr1\t1': {'GT': 'C/C'}, 'chr1\t2': {'GT': 'A/G'},
                    'chr1\t3': {'GT': 'A/A'}, 'chr1\t4': {'GT': 'C/T'},
                    'chr1\t5': {'GT': 'A/A'}, 'chr1\t6': {'GT': 'A/T'},
                    'chr1\t7': {'GT': 'A/A'}, 'chr1\t8': {'GT': 'C/G'}}
        
        var_list2 = {'chr1\t1': {'GT': 'C/C'}, 'chr1\t2': {'GT': 'A/G'},
                     'chr1\t3': {'GT': 'G/G'}, 'chr1\t4': {'GT': 'A/G'},
                     'chr1\t5': {'GT': 'C/T'}, 'chr1\t6': {'GT': 'G/G'},
                     'chr1\t7': {'GT': 'A/T'}, 'chr1\t8': {'GT': 'C/C'}}
        intersection = ['chr1\t1', 'chr1\t2', 'chr1\t3', 'chr1\t4', 'chr1\t5',
                        'chr1\t6', 'chr1\t7', 'chr1\t8']
        # Set up supporting mocks
        # Test
        results = compareGenotypes(var_list, var_list2, intersection, None, None)
        # Check retults
        self.assertEqual(results['total_compared'], 8)
        self.assertEqual(results['ct_common'], 2)
        self.assertEqual(results['comm_hom_ct'], 1)
        self.assertEqual(results['comm_het_ct'], 1)
        self.assertEqual(results['ct_diff'], 6)
        self.assertEqual(results['diff_hom_ct'], 1)
        self.assertEqual(results['diff_het_ct'], 1)
        self.assertEqual(results['diff_hom_het_ct'], 1)
        self.assertEqual(results['diff_het_hom_ct'], 1)
        self.assertEqual(results['diff_1sub2_ct'], 1)
        self.assertEqual(results['diff_2sub1_ct'], 1)

    def test_compareGenotypes_altChroms(self):
		# Set up test case
        var_list = {'chr1\t1': {'GT': 'C/C'}, 'chr1\t2': {'GT': 'A/G'},
                    'chr1\t3': {'GT': 'A/A'}, 'chr1\t4': {'GT': 'C/T'},
                    'chr1\t5': {'GT': 'A/A'}, 'chr1\t6': {'GT': 'A/T'},
                    'chr1\t7': {'GT': 'A/A'}, 'chr1\t8': {'GT': 'C/G'}}

        var_list2 = {'1\t1': {'GT': 'C/C'}, '1\t2': {'GT': 'A/G'},
                     '1\t3': {'GT': 'G/G'}, '1\t4': {'GT': 'A/G'},
                     '1\t5': {'GT': 'C/T'}, '1\t6': {'GT': 'G/G'},
                     '1\t7': {'GT': 'A/T'}, '1\t8': {'GT': 'C/C'}}
        intersection = ['chr1\t1', 'chr1\t2', 'chr1\t3', 'chr1\t4', 'chr1\t5',
                        'chr1\t6', 'chr1\t7', 'chr1\t8']
        def_to_alt = {'chr1': '1'}
        alt_chroms = ['1']
		# Set up supporting mocks
		# Test
        results = compareGenotypes(var_list, var_list2, intersection, alt_chroms, def_to_alt)
		# Check results
        self.assertEqual(results['total_compared'], 8)
        self.assertEqual(results['ct_common'], 2)
        self.assertEqual(results['comm_hom_ct'], 1)
        self.assertEqual(results['comm_het_ct'], 1)
        self.assertEqual(results['ct_diff'], 6)
        self.assertEqual(results['diff_hom_ct'], 1)
        self.assertEqual(results['diff_het_ct'], 1)
        self.assertEqual(results['diff_hom_het_ct'], 1)
        self.assertEqual(results['diff_het_hom_ct'], 1)
        self.assertEqual(results['diff_1sub2_ct'], 1)
        self.assertEqual(results['diff_2sub1_ct'], 1)

    @patch('cohort_matcher.get_tsv_variants')
    @patch('cohort_matcher.getIntersectingVariants')
    @patch('cohort_matcher.compareGenotypes')
    @patch('cohort_matcher.writeSampleComparisonReport')
    @patch('cohort_matcher.writeSimilarityMatrix')
    @patch('os.path.exists')
    def test_compareSamples(self, mock_exists, mock_writeSimilarityMatrix,
			                mock_writeSampleComparisonReport, mock_compareGenotypes,
							mock_getIntersectingVariants, mock_get_tsv_variants):
        # Set up test parameters
        sampleSet1 = [{'name': 'A004AX275-001', 'bam': 'sample1.bam'}]
        sampleSet2 = [{'name': 'A004AX474-001', 'bam': 'sample2.bam'}]
        config = MagicMock(name="config", dp_threshold=10, chromosome_map=None,
                           cache_dir='test_data', scratch_dir='/scratch')
        # Set up supporting mocks
        mock_exists.return_value = False
        # Test
        compareSamples(sampleSet1, sampleSet2, config)
        # Check results
        # Make sure all samples were sample sheet were compared against each other
        mock_writeSampleComparisonReport.assert_called()

    def test_downloadBAMFile(self):
        self.skipTest("not yet implemented")

    @patch('os.path.exists')
    @patch('cohort_matcher.downloadBAMFile')
    @patch('cohort_matcher.get_chrom_names_from_BAM')
    @patch('cohort_matcher.get_chrom_names_from_REF')
    @patch('cohort_matcher.get_chrom_names_from_VCF')
    @patch('subprocess.Popen')
    @patch('cohort_matcher.VCFtoTSV')
    @patch('os.remove')
    def test_genotypeSample(self, mock_osremove, mock_VCFtoTSV, mock_Popen, 
                            mock_getChromNamesFromVCF, mock_getChromNamesFromREF, 
                            mock_getChromNamesFromBAM,  mock_downloadBAMFile, 
                            mock_pathexists):
        # Set up test parameters
        sample = 'asdf'
        bamFile = 's3://some/path/to/asdf.bam'
        # Set up supporting mocks
        mock_pathexists.side_effect = [False, True, False]
        mock_downloadBAMFile.return_value = "/tmp/asdf.bam", "/tmp/asdf.bam.bai"
        mock_getChromNamesFromBAM.return_value = ['chr1', 'chr2', 'chr3']
        mock_getChromNamesFromREF.return_value = ['chr1', 'chr2', 'chr3']
        mock_getChromNamesFromVCF.return_value = ['chr1', 'chr2', 'chr3']
        p = MagicMock(name="caller_cmd", returncode=0)
        p.communicate.return_value = None, None
        mock_Popen.return_value = p
        # Test
        genotypeSample(sample, bamFile, "reference", "vcf", "intervalsFile",
                       MagicMock(name="config", cache_dir="./cache", caller="freebayes"))
        # Check results
        mock_Popen.assert_called_once()
        mock_VCFtoTSV.assert_called_once()
        mock_osremove.assert_any_call('/tmp/asdf.bam')
        mock_osremove.assert_any_call('/tmp/asdf.bam.bai')
        
    @patch('multiprocessing.cpu_count')
    @patch('multiprocessing.Pool')
    @patch('os.path.join')
    @patch('os.path.exists')
    def test_genotypeSamples(self, mock_exists, mock_pathjoin, mock_pool, mock_cpucount):
        # Set up test parameters
        sampleSet = [{'name': 'sample1', 'bam': 'sample1.bam'},
                     {'name': 'sample2', 'bam': 'sample2.bam'}]
        # Set up supporting mocks
        mock_cpucount.return_value = 1
        pool = MagicMock(name="pool")
        mock_pool.return_value = pool

        mock_exists.side_effect = [False, True]
        # Test
        genotypeSamples(sampleSet, MagicMock(name='reference'), MagicMock(name='vcf'),
                        MagicMock(name='intervalsFile'), MagicMock(name='config'))
        # Check results
        pool.apply_async.assert_called_once()

    def test_get_chrom_names_from_BAM(self):
        self.skipTest('nyi')
    
    def test_get_chrom_names_from_REF(self):
        self.skipTest('nyi')
        
    def test_get_chrom_names_from_VCF(self):
        self.skipTest('nyi')

    def test_get_tsv_variants(self):
        # Set up test parameters
        tsvFile = 'test_data/A004AX275-001.tsv'
        dp_threshold = 15
        # Set up supporting mocks
        # Test
        variants = get_tsv_variants(tsvFile, dp_threshold)
        # Check results
        self.assertEqual(len(variants), 5899)
        self.assertEqual(variants["chr1\t881627"], {'ALT': '.', 'GT': 'G/G', 'REF': 'G', 'DP': 50})

    def test_getIntersectingVariants(self):
        # Set up test case
        var_list = ['1\t1', '1\t2', '2\t4']
        var_list2 = ['chr1\t1', 'chr1\t3', 'chr2\t4']
        def_to_alt = {'1': 'chr1', '2': 'chr2'}
        alt_to_def = {'chr1': '1', 'chr2': '2'}
        # Set up supporting mocks
        # Test
        actual = cohort_matcher.getIntersectingVariants(var_list, var_list2,
                                                        def_to_alt, alt_to_def)
        # Check results
        expected = ['1\t1', '2\t4']
        self.assertEqual(actual, expected)
            
    def test_is_same_gt(self):
        self.assertTrue(is_same_gt('A/A', 'A/A'))
        self.assertTrue(is_same_gt('C/T', 'C/T'))
        self.assertFalse(is_same_gt('A/A', 'G/G'))
        self.assertTrue(is_same_gt('A/C', 'C/A'))
        
    @patch('cohort_matcher.parseArguments')
    @patch('cohort_matcher.checkConfig')
    @patch('cohort_matcher.readSamples')
    @patch('cohort_matcher.vcfToIntervals')
    @patch('cohort_matcher.genotypeSamples')
    @patch('cohort_matcher.compareSamples')
    @patch('cohort_matcher.plotResults')
    def test_main(self, mock_plotResults, mock_compareSamples,
                  mock_genotypeSamples, mock_vcfToIntervals,
                  mock_readSamples, mock_checkConfig,
                  mock_parseArguments):
        # Set up test case
        argv = []
        # Set up supporting mocks
        # Test
        with patch('logging.basicConfig') as mock_basicConfig:
            retval = main(argv)
        # Check results
        self.assertEqual(retval, 0)

    @patch('os.path.isfile')
    def test_readSamples(self, mock_isfile):
        # Set up test parameters
        sampleSheetFile = MagicMock(name="sampleSheetFile")
        # Set up supporting mocks
        mock_isfile.return_value = True
        sampleSheetText = ['sample1\ts3://some/path/to/sample1.bam\n',
                           'sample2\ts3://some/path/to/sample2.bam\n']
        mock_open = mock.mock_open(read_data=''.join(sampleSheetText))
        mock_open.return_value.__iter__ = lambda self: iter(self.readline, '')
        # Test
        with patch('__builtin__.open', mock_open):
            actual_sampleSheet = readSamples(sampleSheetFile)
        # Check results
        self.assertEqual(len(actual_sampleSheet), 2)
        self.assertEqual(actual_sampleSheet[0]['name'], 'sample1')
        self.assertEqual(actual_sampleSheet[1]['name'], 'sample2')
        self.assertEqual(actual_sampleSheet[0]['bam'], 's3://some/path/to/sample1.bam')
        self.assertEqual(actual_sampleSheet[1]['bam'], 's3://some/path/to/sample2.bam')
        
    def test_parseArguments(self):
        # Set up test case
        args = ['-S1', 'set1', '-S2', 'set2', '-V', 'somevcf', '-R', 'someref']
        # Set up supporting mocks
        # Test
        actual_args = parseArguments(args)
        # Check results
        self.assertEqual(actual_args.set1, 'set1')
        self.assertEqual(actual_args.set2, 'set2')
        self.assertEqual(actual_args.vcf, 'somevcf')
        self.assertEqual(actual_args.reference, 'someref')

    def test_vcfToIntervals(self):
        # Set up test parameters
        vcfFile = "hg19.exome.highAF.1511.vcf"
        bedFile = NamedTemporaryFile()
        bedFile.close()
        # Set up supporting mocks
        # Test
        vcfToIntervals(vcfFile, bedFile.name)
        # Check results
        self.assertTrue(filecmp.cmp(bedFile.name, "test_data/hg19.exome.highAF.1511.bed"))
        # Clean up
        os.remove(bedFile.name)
    
    def test_VCFtoTSV(self):
        self.skipTest('nyi')
