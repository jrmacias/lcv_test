"""
Document
"""
import unittest
from pipeline_illumina import fastqc


class TestFastqc(unittest.TestCase):
    """_summary_

    Args:
        unittest (_type_): _description_
    """

    def test_fastqc(self):
        """_summary_
        """

        output = ""
        sample = ""
        threads = ""
        reads1 = ""
        reads2 = ""

        fastqc(output, sample, threads, reads1, reads2)


if __name__ == '__main__':
    unittest.main()
