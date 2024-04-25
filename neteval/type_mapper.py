import pandas as pd
import numpy as np
import warnings

# TODO: Figure out the right rules for updated classifications and sub classifications, particularly when only one is supplied.
# Right now I am getting different results if I use the classify function versus if I use the update_classification and update_sub_classification functions alone.

EDGE_TYPES = {"Genetic":{'Phenotypic', 'Synthetic;Dosage', 'Pos;Neg Genetic', 'genetic'}, 
                "Physical":{'physical', 'Two-hybrid;Reconstituted Complex', 'Protein-RNA','Affinity Capture (MS;Western)','Co-crystal;Protein-peptide', 
                        'Affinity Capture-Luminescence', 'Biochemical Activity', 'FRET','Co-purification', 'Co-localization', 'Far Western', 'Co-fractionation',
                        'PCA', 'Proximity Label-MS'}, 
                "Regulatory":{}, 
                "Literature":{}, 'Mixed':{"Mixed", "mixed"}}
SUB_CLASSES = {'Two-hybrid;Reconstituted Complex': ["Two-hybrid", "Reconstituted Complex"],
                'Protein-RNA':["Affinity Capture-RNA", "Protein-RNA"],
                'Phenotypic':["Phenotypic Enhancement", "Phenotypic Suppression"],
                'Synthetic;Dosage': ['Synthetic Rescue','Synthetic Lethality', 'Synthetic Growth Defect', 'Dosage Lethality', 'Dosage Rescue', 'Dosage Growth Defect'],
                'Affinity Capture (MS;Western)': ["Affinity Capture-MS", "Affinity Capture-Western"],
                'Co-crystal;Protein-peptide':['Co-crystal Structure', 'Protein-peptide'],
                'Affinity Capture-Luminescence': ['Affinity Capture-Luminescence'],
                'Biochemical Activity':['Biochemical Activity'],
                'FRET':['FRET'],
                'Co-purification':['Co-purification'],
                'Co-localization': ['Co-localization'],
                'Far Western':['Far Western'],
                'Co-fractionation':['Co-fractionation'],
                'PCA': ['PCA'],
                'Pos;Neg Genetic': ["Negative Genetic", "Positive Genetic"],
                'Proximity Label-MS' : ["Proximity Label-MS"],
                'Mixed': ['Mixed', 'mixed']}

def reverse_dict(d):
    new_dict = {}
    for k, v in d.items():
        new_dict[k] = k
        for i in v:
            new_dict[i] = k
    return new_dict

class PmidClassifications:
    def __init__(self, filepath, write_inplace=True):
        self.pmid_map = pd.read_csv(filepath, sep="\t", index_col=0)
        self.pmids = [str(x) for x in self.pmid_map.index.tolist()]
        self.pmid_map.index = self.pmid_map.index.astype(str)
        self.filepath = filepath
        self.inplace = write_inplace
        self.classification_dict = reverse_dict(EDGE_TYPES)
        self.sub_classification_dict = reverse_dict(SUB_CLASSES)
        
    def classify(self, pmid, classification, sub_classification):
        if type(pmid) != str:
            pmid = str(int(pmid))
        if pmid not in self.pmids:
            self.add_pmid(pmid, classification, sub_classification)
        else:
            current_classification = self.pmid_map.loc[pmid, 'Classification']
            current_sub_classification = self.pmid_map.loc[pmid, 'SubClassification']
            if current_classification != "Unassigned":
                if current_classification != self.classification_dict[classification]:
                    classification = "Mixed"
                    sub_classification = "Mixed"
                else:
                    if current_sub_classification != "Unassigned":
                        if current_sub_classification != self.sub_classification_dict[sub_classification]:
                            sub_classification = "Mixed"   
            else:
                if current_sub_classification != "Unassigned":
                        if current_sub_classification != self.sub_classification_dict[sub_classification]:
                            sub_classification = "Mixed" 
            self.update_classification(pmid, classification)
            self.update_sub_classification(pmid, sub_classification)
    
    def add_pmid(self, pmid, classification, sub_classification):
        if type(pmid) != str:
            pmid = str(int(pmid))
        assert pmid not in self.pmids, f'{pmid} is already in the list of PMIDs'
        classification = self.classification_dict[classification]
        sub_classification = self.sub_classification_dict[sub_classification]
        add_line = pd.DataFrame([[classification, sub_classification]], index=[pmid], columns=['Classification', 'SubClassification'])
        self.pmid_map = self.pmid_map.append(add_line)
        self.pmids.append(pmid)
    
    def update_classification(self, pmid, classification):
        if classification == "Unassigned":
            return
        if type(pmid) != str:
            pmid = str(int(pmid))
        assert pmid in self.pmids, f'{pmid} is not in the list of PMIDs'
        classification = self.classification_dict[classification]
        # Check if the classification is already assigned
        current_classification = self.pmid_map.loc[pmid, 'Classification']
        if current_classification == "Mixed":
            return
        if (current_classification != 'Unassigned') and (current_classification != classification):
            if classification != "Mixed":
                print(f'{pmid} is already classified as {current_classification} and cannot be reclassified as {classification}. Assigning mixed.')
                classification = 'Mixed'
        self.pmid_map.loc[pmid, 'Classification'] = classification
    
    def update_sub_classification(self, pmid, sub_classification):
        if sub_classification == "Unassigned":
            return
        if type(pmid) != str:
            pmid = str(int(pmid))
        assert pmid in self.pmids, f'{pmid} is not in the list of PMIDs'
        sub_classification = self.sub_classification_dict[sub_classification]
        current_sub_classification = self.pmid_map.loc[pmid, 'SubClassification']
        if current_sub_classification == "Mixed":
            return
        if (current_sub_classification != 'Unassigned') and (sub_classification != current_sub_classification):
            if sub_classification != "Mixed":
                print(f'{pmid} is already classified as {current_sub_classification} and cannot be reclassified as {sub_classification}. Assigning mixed.')
                sub_classification = 'Mixed'
        self.pmid_map.loc[pmid, 'SubClassification'] = sub_classification
    
    def get_pmid(self, pmid):
        if type(pmid) != str:
            pmid = str(int(pmid))
        assert pmid in self.pmids, f'{pmid} is not in the list of PMIDs'
        return self.pmid_map.loc[pmid]
    
    def get_classifications(self, pmid_list):
        pmid_list = [str(int(pm)) if type(pm) != str else pm for pm in pmid_list]
        present_pmids = [pm for pm in pmid_list if pm in self.pmids]
        missing_pmids = [pm for pm in pmid_list if pm not in self.pmids]
        present_df = self.pmid_map.loc[present_pmids]
        missing_df = pd.DataFrame([['Unassigned', 'Unassigned']]*len(missing_pmids), index=missing_pmids, columns=['Classification', 'SubClassification'])
        return pd.concat([present_df, missing_df])
    
    def write(self, new_filepath=None):
        if self.pmid_map.index.name != 'PMID':
            self.pmid_map.index.name = 'PMID'
        if self.inplace:
            self.pmid_map.to_csv(self.filepath, sep="\t", index=True)
        else:
            assert new_filepath is not None, "Must provide a new filepath to write to if not writing in place"
            self.pmid_map.to_csv(new_filepath, sep="\t", index=True)

def map_with_pmids(data_to_map, pmid_obj, pmid_col="PMID"):
    data_to_map[pmid_col] = data_to_map[pmid_col].astype(str)
    type_map = pmid_obj.get_classifications(data_to_map[pmid_col].tolist())
    type_map['PMID'] = type_map.index.astype(str)
    data_to_map = data_to_map.merge(type_map, left_on=pmid_col, right_on='PMID', how='left')
    return data_to_map

def assign_type_from_col(data_to_map, type_col=None, sub_col=None):
    assert (type_col is not None) or (sub_col is not None), "Must provide a type OR sub-type column to assign types from"
    edge_dict = reverse_dict(EDGE_TYPES)
    sub_dict = reverse_dict(SUB_CLASSES)
    if (type_col is not None):
        missing_types = [x for x in data_to_map[type_col].unique() if x not in edge_dict]
        if len(missing_types) > 0:
            warnings.warn("Missing types: " +", ".join(missing_types))
    if (sub_col is not None):
        missing_sub_types = [x for x in data_to_map[sub_col].unique() if x not in sub_dict]
        if len(missing_sub_types) > 0:
            warnings.warn("Missing sub-types: "+ ", ".join(missing_sub_types))
    if (type_col is not None) and (sub_col is not None):
        data_to_map['EdgeType'] = data_to_map[type_col].apply(lambda x: edge_dict[x] if x in edge_dict else 'Unassigned')
        data_to_map['EdgeSubType'] = data_to_map[sub_col].apply(lambda x: sub_dict[x] if x in sub_dict else 'Unassigned')
    elif (type_col is not None):
        data_to_map['EdgeType'] = data_to_map[type_col].apply(lambda x: edge_dict[x] if x in edge_dict else 'Unassigned')
        data_to_map['EdgeSubType'] = 'Unassigned'
    elif (sub_col is not None):
        data_to_map['EdgeSubType'] = data_to_map[sub_col].apply(lambda x: sub_dict[x] if x in sub_dict else 'Unassigned')
        data_to_map['EdgeType'] = data_to_map['EdgeSubType'].apply(lambda x: edge_dict[x] if x in edge_dict else 'Unassigned')
    return data_to_map


def assign_type():
    pass

def import_type_classifications():
    pass

def clean_type_data():
    pass

def process_type_data():
    pass

def type_mapper_wrapper(no_type):
    if no_type:
        assign_type()
    else:
        import_type_classifications()
        clean_type_data()
        process_type_data()
    pass

if __name__=="__main__":
    data = pd.read_csv("/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/BIOGRID/BIOGRID-ORGANISM-Homo_sapiens-4.4.217.tab3.txt", 
                    sep="\t", nrows=100)
    data.head()
    
    path = "/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/pmid_classification.txt"
    pmid = PmidClassifications(path)
    data = pd.read_csv("/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/test_biogrid_evidence.txt", sep="\t", index_col=None)
    for i in range(len(data)):
        pmid.classify(data.loc[i, 'PMID'], data.loc[i, 'Experimental System Type'], data.loc[i, 'Experimental System'])
    pmid.write()
    