import neteval.processing_functions as pf
import unittest
import pandas as pd

class Test(unittest.TestCase):
    def setUp(self):
        # Set up test environment before running tests
        pass
    
    def tearDown(self):
        # Clean up test environment after running tests
        pass
    
    def test_prefix_extraction(self):
        # Test function 1
        self.assertEqual("Q59482", pf.extract_id_with_prefix("uniprot:Q59482|asdfa:asdfa", [("uniprot:", "|")]))
        self.assertEqual("Q59482", pf.extract_id_with_prefix("uniprot:Q59482", [("uniprot:", "|")]))
        self.assertEqual("Q59482", pf.extract_id_with_prefix("asdfa:asdfa|uniprot:Q59482", [("uniprot:", "|")]))
        self.assertEqual("Q59482", pf.extract_id_with_prefix("uniprot:Q59482", [("uniprot:")]))
        self.assertEqual("asdfa", pf.extract_id_with_prefix("kb:Q59482|dip:asdfa", [("uniprot:", "|"), ("dip:", "|")]))
        self.assertTrue(pd.isna(pf.extract_id_with_prefix("Q59482", [("uniprot:", "|")])))
        self.assertEqual("DIP-asdfa", pf.extract_id_with_prefix("kb:Q59482|DIP-asdfa", [("DIP-", "|"), ("kb:", "|")]))


if __name__ == '__main__':
    unittest.main()
