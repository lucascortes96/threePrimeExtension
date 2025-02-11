import unittest
from unittest.mock import MagicMock
import pandas as pd
import pyBigWig

# Assuming getIntervals function is defined in findBigWigCoverage.py
from findBigWigCoverage import getIntervals

class TestGetIntervals(unittest.TestCase):
    def setUp(self):
        # Create a mock BigWig object
        self.bw = MagicMock(spec=pyBigWig.pyBigWig)
        
        # Mock intervals method
        self.bw.intervals = MagicMock(side_effect=self.mock_intervals)
        
        # Create test DataFrames
        self.df1 = pd.DataFrame({
            'Chromosome': ['chr1', 'chr1'],
            'Start': [110, 200],
            'End': [120, 210],
            'Strand1': ['forward', 'reverse'],
            'Other': ['data1', 'data2'],
            'MoreData': ['data3', 'data4'],
            'Strand': ['forward', 'reverse'],
            'Extra': ['data5', 'data6'],
            'MoreExtra': ['data7', 'data8'],
            'StartPos': [110, 200],
            'EndPos': [120, 210]
        })
        
        self.df2 = pd.DataFrame({
            'Chromosome': ['chr2', 'chr2'],
            'Start': [300, 400],
            'End': [310, 410],
            'Strand1': ['forward', 'reverse'],
            'Other': ['data9', 'data10'],
            'MoreData': ['data11', 'data12'],
            'Strand': ['forward', 'reverse'],
            'Extra': ['data13', 'data14'],
            'MoreExtra': ['data15', 'data16'],
            'StartPos': [310, 400],
            'EndPos': [320, 410]
        })

    def mock_intervals(self, chrom, start, end):
        # Mock interval data
        if chrom == 'chr1':
            if start == 110 and end == 120:
                return [(110, 120, 15)]
            elif start == 90 and end == 100:
                return [(90, 100, 5)]
        elif chrom == 'chr2':
            if start == 310 and end == 320:
                return [(310, 320, 20)]
            elif start == 390 and end == 400:
                return [(390, 400, 8)]
        return []

    def test_get_intervals_df1(self):
        result = getIntervals(self.bw, self.df1, 'chr1')
        expected = pd.DataFrame({
            'Chromosome': ['chr1', 'chr1'],
            'Start': [110, 200],
            'End': [120, 210],
            'Strand1': ['forward', 'reverse'],
            'Other': ['data1', 'data2'],
            'MoreData': ['data3', 'data4'],
            'Strand': ['forward', 'reverse'],
            'Extra': ['data5', 'data6'],
            'MoreExtra': ['data7', 'data8'],
            'StartPos': [110, 200],
            'EndPos': [120, 210],
            'longreadSupport': [True, False]
        })
        pd.testing.assert_frame_equal(result, expected)

    #def test_get_intervals_df2(self):
        result = getIntervals(self.bw, self.df2, 'chr2')
        expected = pd.DataFrame({
            'Chromosome': ['chr2', 'chr2'],
            'Start': [310, 400],
            'End': [320, 410],
            'Strand1': ['forward', 'reverse'],
            'Other': ['data9', 'data10'],
            'MoreData': ['data11', 'data12'],
            'Strand': ['forward', 'reverse'],
            'Extra': ['data13', 'data14'],
            'MoreExtra': ['data15', 'data16'],
            'StartPos': [310, 400],
            'EndPos': [320, 410],
            'longreadSupport': [True, False]
        })
        pd.testing.assert_frame_equal(result, expected)

    def test_has_value_greater_than_10(self):
        # Test cases for has_value_greater_than_10 logic
        interval_data_1 = [(120, 130, 15)]
        interval_data_2 = [(110, 120, 5)]
        interval_data_3 = [(310, 320, 20)]
        interval_data_4 = [(390, 400, 8)]
            
        self.assertTrue(any(value > 10 for _, _, value in interval_data_1))
        self.assertFalse(any(value > 10 for _, _, value in interval_data_2))
        self.assertTrue(any(value > 10 for _, _, value in interval_data_3))
        self.assertFalse(any(value > 10 for _, _, value in interval_data_4))

        test = self.bw.intervals('chr1', 110, 120)
        self.assertTrue(any(value > 10 for _, _, value in test))

if __name__ == '__main__':
    unittest.main()