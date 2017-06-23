import unittest
import logging
from mock import patch, MagicMock, mock
from tempfile import NamedTemporaryFile
import os
import filecmp
from cohort_matcher import checkConfig, main, parseArguments, readSamples, \
     vcfToIntervals, genotypeSamples, genotypeSample, compareSamples, get_tsv_variants

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

    def test_compareSamples(self):
        # Set up test parameters
        sampleSet1 = [{'name': 'A004AX275-001', 'bam': 'sample1.bam'}]
        sampleSet2 = [{'name': 'A004AX474-001', 'bam': 'sample2.bam'}]
        config = MagicMock(name="config", dp_threshold=10, chromosome_map=None,
                           cache_dir='test_data', scratch_dir='/scratch')
        # Set up supporting mocks
        #mock_get_tsv_variants.return_value = [{"chr1\t123": 1},
        #                                      {"chr1\t1234": 1}]
        # Test
        compareSamples(sampleSet1, sampleSet2, config)
        # Check results

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
    def test_genotypeSample(self, mock_osremove, mock_VCFtoTSV, mock_Popen, mock_getChromNamesFromVCF, 
                            mock_getChromNamesFromREF, mock_getChromNamesFromBAM, 
                            mock_downloadBAMFile, mock_pathexists):
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
        genotypeSample(sample, bamFile, "reference", "vcf", "intervalsFile", MagicMock(name="config", cache_dir="./cache", caller="freebayes"))
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
        tsvFile = 'test_data/sample1.tsv'
        dp_threshold = 15
        # Set up supporting mocks
        # Test
        variants = get_tsv_variants(tsvFile, dp_threshold)
        # Check results
        self.assertEqual(len(variants), 5)
        self.assertEqual(variants["chr1\t881627"], 1)
        
    @patch('cohort_matcher.parseArguments')
    @patch('cohort_matcher.checkConfig')
    @patch('cohort_matcher.readSamples')
    @patch('cohort_matcher.vcfToIntervals')
    @patch('cohort_matcher.genotypeSamples')
    @patch('cohort_matcher.compareSamples')
    def test_main(self, mock_compareSamples, mock_genotypeSamples,
                  mock_vcfToIntervals, mock_readSamples,
                  mock_checkConfig, mock_parseArguments):
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