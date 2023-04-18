import neteval.type_mapper as tmap
import unittest
import pandas as pd
import os


class Test(unittest.TestCase):
    def setUp(self):
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        self.pmid = tmap.PmidClassifications(self.dir_path+"/data/pmid_unassigned.txt", write_inplace=False)
        self.pmid_assigned = tmap.PmidClassifications(self.dir_path+"/data/pmid_assigned.txt", write_inplace=False)
    
    def tearDown(self):
        file_list = []
        for f in file_list:
            if os.path.exists(self.dir_path + "/data/"+f):
                os.remove(self.dir_path + "/data/"+f)

    def test_global_variables(self):
        self.assertIsInstance(tmap.EDGE_TYPES, dict)
        self.assertIsInstance(tmap.SUB_CLASSES, dict)
        edge_dict = tmap.reverse_dict(tmap.EDGE_TYPES)
        sub_dict = tmap.reverse_dict(tmap.SUB_CLASSES)
        # test that all sub-classes are associated with a type
        for subC in tmap.SUB_CLASSES.keys():
            self.assertTrue(subC in edge_dict, "Subclass: %s does not have a type mapping" % subC)


    def test_reverse_dict(self):
        test_dict = {'a':{1, 4, 5}, 'b':{2}, 'c':{3, 7, 8}}
        rev_dict = tmap.reverse_dict(test_dict)
        # test that all values have been reversed
        for k in test_dict.keys():
            for v in test_dict[k]:
                self.assertTrue(v in rev_dict, "Value: %s not in reversed dictionary" % v)
        # test that self edges are added
        for k in test_dict.keys():
            self.assertEqual(k,  rev_dict[k], "Self edge not added for %s" % k)
        # test final dictionary has correct values
        self.assertDictEqual(rev_dict, {1:'a', 2:'b', 3:'c', 4:'a', 5:'a', 7:'c', 8:'c', 'a':'a', 'b':'b', 'c':'c'})
    
    def test_pmid_object_initialization(self):
        pmid = tmap.PmidClassifications(self.dir_path+"/data/pmid_unassigned.txt", write_inplace=False)
        self.assertEqual(pmid.pmid_map.shape[1], 2, "PMID map does not have correct number of columns")
    
    def test_add_pmid_classifications(self):
        # test exception raised if pmid is already in the data
        self.assertRaises(AssertionError, self.pmid.add_pmid, pmid='15546613', classification='genetic', sub_classification='Phenotypic Enhancement')
        # test adding a new pmid
        self.pmid.add_pmid('123', 'genetic', 'Phenotypic Enhancement')
        # check that it is added to pmid list
        self.assertIn('123', self.pmid.pmids, "PMID not added to pmids list")
        # check that classifications added correctly
        self.assertEqual(self.pmid.pmid_map.loc['123', "Classification"], "Genetic", "Classification added does not match")
        self.assertEqual(self.pmid.pmid_map.loc['123', "SubClassification"], "Phenotypic", "Sub-Classification added does not match")
        # test adding an integer pmid
        self.pmid.add_pmid(456, 'genetic', 'Phenotypic Enhancement')
        self.assertIn('456', self.pmid.pmids, "PMID not added to pmids list")

    def test_update_pmid_classification(self):
        # test exception raised if pmid is not present
        self.assertRaises(AssertionError, self.pmid.update_classification, pmid='789', classification='genetic')
        self.assertRaises(AssertionError, self.pmid.update_sub_classification, pmid='789', sub_classification='genetic')
        # update from U-U to X-U
        self.pmid.update_classification('15546613', 'genetic')
        self.assertEqual(self.pmid.get_pmid('15546613').tolist(), ['Genetic', 'Unassigned'], "PMID classification not updated correctly")
        self.assertEqual
        # updated form U-U to U-X
        self.pmid.update_sub_classification('14633992', 'Phenotypic Enhancement')
        self.assertEqual(self.pmid.get_pmid('14633992').tolist(), ['Unassigned', 'Phenotypic'], "PMID subclassification not updated correctly")
        # updated from U-U to X-X
        self.pmid.classify('15829968', 'physical', 'FRET')
        self.assertEqual(self.pmid.get_pmid('15829968').tolist(), ['Physical', 'FRET'], "PMID classification not updated correctly")
        # updated from X to M
        self.pmid.classify('12808131', 'physical', 'PCA')
        self.pmid.classify('12808131', 'genetic', 'PCA')
        self.assertEqual(self.pmid.get_pmid('12808131').tolist(), ['Mixed', 'Mixed'], "PMID classification not updated correctly")
        # updated from M-M to M-M
        self.pmid.classify('15616580', 'physical', 'PCA')
        self.pmid.update_sub_classification('15616580', 'FRET')
        self.assertEqual(self.pmid.get_pmid('15616580').tolist(), ['Physical', 'Mixed'], "PMID classification not updated correctly")
        
    def test_map_with_pmids(self):
        data = pd.read_csv(self.dir_path+"/data/pmids.txt", sep="\t", header=None, names=["PMID"], index_col=None)
        mapped_data = tmap.map_with_pmids(data, self.pmid_assigned, pmid_col="PMID")
        #how to test the results?
        self.assertEqual(mapped_data.shape[1], 3, "Mapped data does not have correct number of columns")
        
            
    def test_assign_types_from_col(self):
        data = pd.read_csv(self.dir_path+"/data/pmid_evidence_test_data.txt", sep="\t", index_col=None)
        # assign both types
        mapped_data = tmap.assign_type_from_col(data, type_col="Type", sub_col="Subtype")
        self.assertIn('EdgeType', mapped_data.columns)
        self.assertIn('EdgeSubType', mapped_data.columns)

    def test_assign_subtype_from_col(self):
        data = pd.read_csv(self.dir_path+"/data/pmid_evidence_test_data.txt", sep="\t", index_col=None)
        mapped_data = tmap.assign_type_from_col(data, sub_col="Subtype")
        self.assertIn('EdgeType', mapped_data.columns)
        self.assertIn('EdgeSubType', mapped_data.columns)
        
    def test_assign_type_only_from_col(self):
        data = pd.read_csv(self.dir_path+"/data/pmid_evidence_test_data.txt", sep="\t", index_col=None)
        mapped_data = tmap.assign_type_from_col(data, type_col="Type")
        self.assertIn('EdgeType', mapped_data.columns)
        self.assertIn('EdgeSubType', mapped_data.columns)
        # check that subtype is added and = unassigned
        self.assertSetEqual(set(mapped_data["EdgeSubType"].values), {"Unassigned"})
        
    def test_missing_types_assign_from_col(self):
        data  = pd.DataFrame({"Subtype": 'Reconstituted Complex', "Type": 'test1'}, index=[1])
        self.assertWarns(UserWarning, tmap.assign_type_from_col, data_to_map=data, type_col="Type", sub_col="Subtype")

if __name__ == '__main__':
    unittest.main()