import neteval.gene_mapper as gmap
import unittest
import pandas as pd
import os


class Test(unittest.TestCase):
    def setUp(self):
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
    
    def tearDown(self):
        file_list = []
        for f in file_list:
            if os.path.exists(self.dir_path + "/data/"+f):
                os.remove(self.dir_path + "/data/"+f)
                
    def test_mygene_api(self):
        pass
    
    def test_uniprot_api(self):
        pass
    
    def test_hgnc_api(self):
        pass
    
    def test_ensembl_api(self):
        pass
    
    def test_update_nodes_general(self):
        pass
    
    def test_convert_nodes_general(self):
        pass
    
    def test_update_uniprot(self):
        pass
    
    def test_update_ensembl(self):
        pass
    
    def test_update_symbol(self):
        pass
    
    def test_symbol_to_entrez(self):
        pass
    
    def test_ensembl_to_entrez(self):
        pass
    
    def test_uniprot_to_entrez(self):
        pass
    
    def test_refseq_to_entrex(self):
        pass
    
    def test_ensembl_to_symbol(self):
        pass
    
    def test_uniprot_to_symbol(self):
        pass
    
    def test_refseq_to_symbol(self):
        pass
    
    def test_entrez_to_symbol(self):
        pass
    
    def test_outdated_uniprot(self):
        pass
    
    def test_previous_versus_alias(self):
        pass
    

if __name__ == '__main__':
    unittest.main()
