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
        G_df_attr = pd.read_csv(self.dir_path + "/data/example_graph_attributes.txt", sep="\t")
        self.G_attr = nx.from_pandas_edgelist(G_df_attr, source="A", target="B", edge_attr=True)
        self.shuff_G = nx.from_pandas_edgelist(G_df, source="A", target="B")
        self.test_uuid = 'f7a218c0-2376-11ea-bb65-0ac135e8bacf' # bioplex 3 shared interactions
        self.export_file = self.dir_path + "/data/export_test_network.txt"
        self.pcnet_file = self.dir_path + "/data/export_test_pcnet.txt"
        self.db_names = {'m': 'M', 'n': 'N', 'p': 'P', 'o': 'O', 'q': 'Q', 'r': 'R'}
        # Set up test environment before running tests
    
    def tearDown(self):
        file_list = ['test_edgelist.txt', 'test_edgelist_attr.txt', 'test_genesets.txt']
        for f in file_list:
            if os.path.exists(self.dir_path + "/data/"+f):
                os.remove(self.dir_path + "/data/"+f)
    
    def test_ndex_download(self):
        net = dit.load_public_network_from_ndex(self.test_uuid)
        self.assertEqual(len(net.nodes), 8364)
        self.assertEqual(len(net.edges), 35704)
    
    def test_load_genesets(self):
        genesets = dit.load_node_sets(self.dir_path + '/data/example_genesets.txt', id_type='Entrez')
        self.assertEqual(len(genesets), 26)
        
    def test_write_genesets(self):
        genesets = {'A': {1, 2, 3}, 'B': {4, 5, 6}}
        dit.write_node_sets(genesets, self.dir_path + '/data/test_genesets.txt')
        genesets_test = dit.load_node_sets(self.dir_path + '/data/test_genesets.txt', id_type='Entrez')
        self.assertEqual(genesets, genesets_test)
    
    def test_load_network_to_networkx_default(self):
        G_test = dit.load_edgelist_to_networkx(self.dir_path + "/data/example_graph.txt")
        self.assertTrue(nx.is_isomorphic(G_test, self.G), "Function does not load equivalent to manual loading")
        self.assertIsInstance(list(G_test.nodes)[0], int, "Nodes are not integers")
        
    def test_load_network_to_networkx_with_attributes(self):
        G_test = dit.load_edgelist_to_networkx(self.dir_path + "/data/example_graph_attributes.txt", keep_attributes=True)
        self.assertTrue(nx.is_isomorphic(G_test, self.G_attr), "Function does not load equivalent to manual loading")
        
    def test_write_plain_edgelist(self):
        dit.write_networkx_to_file(self.G, self.dir_path + "/data/test_edgelist.txt", source='A', target='B')
        G_test = dit.load_edgelist_to_networkx(self.dir_path + "/data/test_edgelist.txt")
        self.assertTrue(nx.is_isomorphic(G_test, self.G), "Changes occured in writing and reading edgelist")
    
    def test_write_edgelist_with_attributes(self):
        dit.write_networkx_to_file(self.G_attr, self.dir_path + "/data/test_edgelist_attr.txt", source='A', target='B')
        G_test = dit.load_edgelist_to_networkx(self.dir_path + "/data/test_edgelist_attr.txt")
        self.assertTrue(nx.is_isomorphic(G_test, self.G_attr), "Changes occured in writing and reading edgelist")
    
    def test_create_networkx_for_export_basic(self):
        G_export = dit.create_networkx_for_export(self.export_file, id_type='Entrez')
        self.assertEqual(len(G_export.nodes), 300, 'Number of nodes is not correct')
        self.assertEqual(len(G_export.edges), 299, 'Number of edges is not correct')
        self.assertIn('name', G_export.edges[(3172, 4041)], 'Attributes are not correct')
        
    def test_create_networkx_for_export_with_attributes(self):
        G_export = dit.create_networkx_for_export(self.export_file, id_type='Entrez', add_symbols=True, attributes=['Score'])
        self.assertEqual(len(G_export.nodes), 300, 'Number of nodes is not correct')
        self.assertEqual(len(G_export.edges), 299, 'Number of edges is not correct')
        self.assertIn('name', G_export.edges[(3172, 4041)], 'Edge attributes are not correct')
        self.assertIn('name', G_export.nodes[3172], 'Node attributes are not correct')
    
    def test_create_pcnet_networkx_for_export(self):
        G_export = dit.create_pcnet_networkx_for_export(self.pcnet_file, id_type='Entrez', db_name_dict=self.db_names)
        self.assertEqual(len(G_export.nodes), 300, 'Number of nodes is not correct')
        self.assertEqual(len(G_export.edges), 299, 'Number of edges is not correct')
        self.assertIn('name', G_export.edges[(3172, 4041)], 'Edge attributes are not correct')
        self.assertIn('Supporting_Databases', G_export.edges[(3172, 4041)], 'Edge attributes are not correct')
        self.assertIsInstance(G_export.edges[(3172, 4041)]['Supporting_Databases'], list, 'Supporting databases is not a list')
        self.assertIn('name', G_export.nodes[3172], 'Node attributes are not correct')
    

if __name__ == '__main__':
    unittest.main()
