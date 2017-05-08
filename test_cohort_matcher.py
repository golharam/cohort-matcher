import unittest
from mock import patch

from cohort_matcher import main, parseArguments

class TestCohortMatcher(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_checkConfig(self):
        self.skipTest("not yet implemented")

    def test_compareSamples(self):
        self.skipTest("not yet implemented")

    def test_genotypeSamples(self):
        self.skipTest("not yet implemented")

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
        retval = main(argv)
        # Check results
        self.assertEqual(retval, 0)

    def test_readSamples(self):
        self.skipTest("not yet implemented")

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
        self.skipTest("not yet implemented")
