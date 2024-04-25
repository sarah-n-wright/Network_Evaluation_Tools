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
    
    def test_parsing(self):
        pass
    
    def test_clean_score(self):
        self.assertEqual(0.7, pf.clean_score("0.7")) #string to float
        self.assertEqual(0.7, pf.clean_score(0.7)) #float to float
        self.assertEqual(7, pf.clean_score("7")) #int string to float
        self.assertIsInstance(pf.clean_score("asdgasga"), type(pd.NA)) #no number in the string, should return NA
        self.assertIsInstance(pf.clean_score(dict()), type(pd.NA)) # invalid type
        self.assertEqual(8, pf.clean_score("asfasga8agasd")) #TODO, should this one work???
        self.assertEqual(8, pf.clean_score("score:8")) # separator, ending number
        self.assertEqual(8, pf.clean_score("score|8")) # different separator
        self.assertEqual(8, pf.clean_score("8=score")) #     leading number, separator
        self.assertEqual(8, pf.clean_score("8|score|9")) # take just the first number
        self.assertEqual(0.4, pf.clean_score("4e-1")) # scientific notation
        self.assertEqual(0.4, pf.clean_score("4E-1")) # scientific notation capitalized
    
    def test_binarize_complexes(self):
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
