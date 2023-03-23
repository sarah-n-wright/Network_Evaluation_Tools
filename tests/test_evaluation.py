from neteval.shuffle_networks import load_network, write_network, shuffle_network, parse_arguments
from neteval.data_import_tools import load_node_sets
from neteval.network_evaluation_functions import calculate_p, small_network_AUPRC_wrapper, construct_prop_kernel, calculate_network_performance_score, calculate_network_performance_gain
from neteval.network_propagation import closed_form_network_propagation, calculate_alpha, normalize_network
import unittest
import pandas as pd
import networkx as nx
import os
import numpy as np
import warnings
warnings.filterwarnings(action='ignore', category=DeprecationWarning)

class Test(unittest.TestCase):
    def setUp(self):
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        G_df = pd.read_csv(self.dir_path + "/data/example_graph.txt", sep="\t")
        self.G = nx.from_pandas_edgelist(G_df, source="A", target="B")
        self.shuff_G = nx.from_pandas_edgelist(G_df, source="A", target="B")
        self.genesets = load_node_sets('/cellar/users/snwright/Git/Network_Evaluation_Tools/tests/data/example_genesets.txt', id_type='Entrez')
        
        # Set up test environment before running tests
    
    def tearDown(self):
        file_list = ["example_graph_shuffled.txt"]
        for f in file_list:
            if os.path.exists(self.dir_path + "/test_data/"+f):
                os.remove(self.dir_path + "/test_data/"+f)
        # delete files generated
        
    def test_load_genesets(self):
        genesets = load_node_sets('/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/DisGeNET_genesets_entrez.txt', id_type='Entrez')
        self.assertEqual(len(genesets), 446)
    
    def test_calculate_p(self):
        genesets_p = calculate_p(self.G, self.genesets, id_type='Entrez')
        self.assertEqual(len(genesets_p), len(self.genesets))
    
    def test_calculate_alpha(self):
        m=-0.02935302
        b=0.74842057
        alpha = calculate_alpha(self.G, m=m, b=b)
        self.assertEqual(alpha, 0.656)
    
    def test_closed_form_propagation(self):
        binary = pd.DataFrame(data=np.identity(len(self.G.nodes())), index=self.G.nodes(), columns=self.G.nodes())
        net_matrix = closed_form_network_propagation(network=self.G, binary_matrix=binary, network_alpha=0.5, verbose=True)
        self.assertEqual(binary.shape, net_matrix.shape)
        
    def test_closed_form_propagation_symmetric(self):
        binary = pd.DataFrame(data=np.identity(len(self.G.nodes())), index=self.G.nodes(), columns=self.G.nodes())
        x = closed_form_network_propagation(network=self.G, binary_matrix=binary, network_alpha=0.5, verbose=False,
                                            symmetric_norm=True)
        
    def test_normalize_network(self):
        array_norm = normalize_network(self.G, symmetric_norm=False)
        array_sym_norm = normalize_network(self.G, symmetric_norm=True)
        self.assertEqual(array_norm.shape, array_sym_norm.shape)

    def test_load_network(self):
        G_test = load_network(self.dir_path + "/data/example_graph.txt", node_prefix="")
        self.assertTrue(nx.is_isomorphic(G_test, self.G), "Function does not load equivalent to manual loading")
    
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
        
    def test_small_network_AUPRC_wrapper(self):
        genesets_p = calculate_p(self.G, self.genesets, id_type='Entrez')
        m=-0.02935302
        b=0.74842057
        alpha = calculate_alpha(self.G, m=m, b=b)
        kernel = construct_prop_kernel(self.G, alpha=alpha, verbose=True)
        AUPRC_values = small_network_AUPRC_wrapper(kernel, self.genesets, genesets_p, n=30, cores=4, verbose=True)
        null_AUPRCs = []
        print("loop")
        for i in range(10):
            shuffNet = shuffle_network(self.G, n_swaps=1)
            shuffNet_kernel = construct_prop_kernel(shuffNet, alpha=alpha, verbose=False)
            shuffNet_AUPRCs = small_network_AUPRC_wrapper(shuffNet_kernel, self.genesets, genesets_p, n=30, cores=4, verbose=False)
            null_AUPRCs.append(shuffNet_AUPRCs)
            print('shuffNet', repr(i+1), 'AUPRCs calculated')
        null_AUPRCs_table = pd.concat(null_AUPRCs, axis=1)
        null_AUPRCs_table.columns = ['shuffNet'+repr(i+1) for i in range(len(null_AUPRCs))]
        network_performance = calculate_network_performance_score(AUPRC_values, null_AUPRCs_table, verbose=True)
        network_perf_gain = calculate_network_performance_gain(AUPRC_values, null_AUPRCs_table, verbose=True)
        network_perf_gain.name = 'Test Network'
        network_performance.name = 'Test Network'
    
    def test_calculate_network_performance_score(self):
        pass
    
    def test_calculate_network_performance_gain(self):
        pass
        
    
    def test_write_network(self):
        write_network(self.G, datafile=self.dir_path + "/data/example_graph.txt", outpath=self.dir_path + "/data/")
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
    
    def test_parsing_non_default(self):
        test_string = ["--testMode", "0", "--nSwaps", "4", "-o", "out.path", "datafile.txt"]
        args = parse_arguments(test_string)
        self.assertEqual(args.o, "out.path")
        self.assertEqual(args.datafile, "datafile.txt")
        self.assertEqual(args.nSwaps, 4)
        self.assertFalse(args.testMode) 
        
if __name__=="__main__":
    unittest.main()
