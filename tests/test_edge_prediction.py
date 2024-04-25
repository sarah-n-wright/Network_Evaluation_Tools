import neteval.edge_prediction as ep
import unittest
import pandas as pd
import os

class Test(unittest.TestCase):
    def setUp(self):
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        self.example_edgelist = self.dir_path + "/data/example_graph.txt"
        self.example_GS = self.dir_path + "/data/example_GS_edges.txt"
        # Set up test environment before running tests
    
    def tearDown(self):
        # Clean up test environment after running tests
        pass
    
    def test_edge_list(self):
        x = ep.load_edge_list(self.example_edgelist, edge_tuples=False)
        self.assertEqual(len(x), 1458)
        x = ep.load_edge_list(self.example_edgelist, edge_tuples=True)
        self.assertIsInstance(x, set)
        
    def test_edge_filter(self):
        test_edges = ep.load_edge_list(self.example_GS, edge_tuples=False)
        x = ep.filter_test_edges(test_edges, self.example_edgelist)
        node_exclude = set([(100, 101), (100, 102), (100, 105), (4, 101), (2, 100)])
        edge_exclude = set([(4, 52), (4, 56), (4, 59), (4, 61), (4, 62), (4, 69)])
        self.assertEqual(len(x-node_exclude), len(x), "Node filtering is not working")
        self.assertEqual(len(x-edge_exclude), len(x), "Edge filtering is not working")
        

    
if __name__ == '__main__':
    unittest.main()
