from neteval.shuffle_networks import shuffle_network
from neteval.data_import_export_tools import load_node_sets
from neteval.network_evaluation_functions import calculate_p, small_network_AUPRC_wrapper, construct_prop_kernel, calculate_network_performance_score, calculate_network_performance_gain, calculate_precision_recall
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
        file_list = []
        for f in file_list:
            if os.path.exists(self.dir_path + "/data/"+f):
                os.remove(self.dir_path + "/data/"+f)
        # delete files generated
        
    def test_calculate_precision_recall(self):
        pass
    
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
        
if __name__=="__main__":
    unittest.main()
