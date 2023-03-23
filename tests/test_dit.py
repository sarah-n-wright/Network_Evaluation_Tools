import neteval.data_import_export_tools as dit
import unittest
import pandas as pd
import networkx as nx
import os


class Test(unittest.TestCase):
    def setUp(self):
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        G_df = pd.read_csv(self.dir_path + "/data/example_graph.txt", sep="\t")
        self.G = nx.from_pandas_edgelist(G_df, source="A", target="B")
        self.shuff_G = nx.from_pandas_edgelist(G_df, source="A", target="B")        
        # Set up test environment before running tests
    
    def tearDown(self):
        file_list = []
        for f in file_list:
            if os.path.exists(self.dir_path + "/data/"+f):
                os.remove(self.dir_path + "/data/"+f)
        # delete files generated
    
    def test_ndex_download(self):
        pass
    
    def test_ndex_upload(self):
        pass
    
    def test_load_genesets(self):
        genesets = dit.load_node_sets(self.dir_path + '/data/example_genesets.txt', id_type='Entrez')
        self.assertEqual(len(genesets), 26)
    
    def test_load_network_to_networkx(self):
        #TODO
        pass
        #G_test = dit.load_edgelist_to_networkx(self.dir_path + "/data/example_graph.txt", node_prefix="")
        #self.assertTrue(nx.is_isomorphic(G_test, self.G), "Function does not load equivalent to manual loading")

    def test_load_network_to_dataframe(self):
        pass
        
    def test_write_plain_edgelist(self):
        pass
    
    def test_write_edgelist_with_attributes(self):
        pass


if __name__ == '__main__':
    unittest.main()
