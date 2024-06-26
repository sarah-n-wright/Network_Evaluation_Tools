import neteval.data_import_export_tools as dit
from neteval.shuffle_networks import shuffle_network, write_shuffled_network, parse_arguments
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
        #self.genesets = dit.load_node_sets(self.dir_path + '/data/example_genesets.txt', id_type='Entrez')

        # Set up test environment before running tests
    
    def tearDown(self):
        file_list = ["example_graph_shuffled.txt"]
        for f in file_list:
            if os.path.exists(self.dir_path + "/data/"+f):
                os.remove(self.dir_path + "/data/"+f)
        # delete files generated
    
    def test_shuffle_original_unmodified(self):
        _ = shuffle_network(self.shuff_G, n_swaps=1)
        self.assertTrue(nx.is_isomorphic(self.shuff_G, self.G), "Original graph is being modified by shuffling function")

    def test_shuffle_nswaps(self):
        G_shuff1 = shuffle_network(self.shuff_G, n_swaps=0.1)
        G_shuff2 = shuffle_network(self.shuff_G, n_swaps=0.8)
        G_shuff3 = shuffle_network(self.shuff_G, n_swaps=2)
        shared_edges1 = len(set(self.shuff_G.edges()).intersection(set(G_shuff1.edges())))
        shared_edges2 = len(set(self.shuff_G.edges()).intersection(set(G_shuff2.edges())))
        shared_edges3 = len(set(self.shuff_G.edges()).intersection(set(G_shuff3.edges())))
        self.assertGreater(shared_edges1, shared_edges2)
        self.assertGreater(shared_edges2, shared_edges3)
    
    def test_write_shuffled_network(self):
        write_shuffled_network(self.G, datafile=self.dir_path + "/data/example_graph.txt", outpath=self.dir_path + "/data/")
        self.assertTrue(os.path.exists(self.dir_path + "/data/example_graph_shuffled.txt"), "File does not exist in desired location.")
        G_test = nx.from_pandas_edgelist(pd.read_csv(self.dir_path + "/data/example_graph_shuffled.txt", sep="\t"), source="Node_A", target="Node_B")
        self.assertTrue(nx.is_isomorphic(G_test, self.G), "Written network is incorrectly modified from original")
        
    def test_parsing_defaults(self):
        test_string = ["-o", "out.path", "datafile.txt"]
        args = parse_arguments(test_string)
        self.assertEqual(args.o, "out.path")
        self.assertEqual(args.datafile, "datafile.txt")
        self.assertEqual(args.nSwaps, 1)
        self.assertTrue(args.testMode)
        self.assertFalse(args.verbose)
    
    def test_parsing_non_default(self):
        test_string = ["--testMode", "0", "--nSwaps", "4", "-o", "out.path", "--verbose", "1", "datafile.txt"]
        args = parse_arguments(test_string)
        self.assertEqual(args.o, "out.path")
        self.assertEqual(args.datafile, "datafile.txt")
        self.assertEqual(args.nSwaps, 4)
        self.assertFalse(args.testMode) 
        self.assertTrue(args.verbose)
        


if __name__ == '__main__':
    unittest.main()
