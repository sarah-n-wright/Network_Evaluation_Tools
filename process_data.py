import pandas as pd

def update_nodes(nodes, id_type, keep="present"):
    """ Takes a set of node identifiers and updates them to the latest version of the same identifier type.

    Args:
        nodes (set): The set of nodes to be updated
        id_type (str): The type of identifier being used and updated
        keep (str): Which nodes should be kept? "updated" for only those that have been updated, "present" for all those present in the dataset, 
        whether updated or not, "all" to keep all nodes (even if they are not present in the mapping data)

    Returns:
        dict: Mapping between input nodes (key) and updated identifiers (values)
    """
    # must return 1:1
    updated_node_map = {}
    if id_type == "Uniprot":
        pass
    elif id_type == "Symbol":
        pass
    elif id_type == "Ensembl":
        pass
    elif id_type == "Entrez":
        pass
    return updated_node_map
    
def convert_nodes(nodes, initial_id, target_id):
    """ Converts nodes between two different identifier types

    Args:
        nodes (set): Set of nodes to be converted
        initial_id (str): Identifier type of input nodes
        target_id (str): Identifier type to be converted to
        
    Returns:
        dict: mapping between input nodes and new identifier
        set: nodes that were not able to be mapped to new identifiers.
    """
    pass


class NetworkData:
    def __init__(self, datafile, node_a, node_b, identifiers, target_id_type, net_name, score=None, species=None, species_code=None, sep="\t"):

        available_id_types = ["Symbol", "Entrez", "Uniprot", "Ensembl"]
        assert species is None or ((species is not None) and (species_code is not None)), "If species column is given, a corresponding value representing the species must be given"
        assert target_id_type in available_id_types, "Target identifier must be in " + "; ".join(available_id_types)
        
        self.species = species
        self.score = score
        self.species_code = species_code
        self.node_a = node_a
        self.node_b = node_b
        self.cols = [col for col in [node_a, node_b, score, species] if col is not None]
        self.target_id_type = target_id_type
        self.net_name = net_name
        
        # check the identifiers given
        if type(identifiers) == str:
            assert identifiers in available_id_types, "Identifiers must be in " + "; ".join(available_id_types)
            self.mixed_identifiers = False
            self.indentifiers = identifiers
        elif type(identifiers) == list:
            for id_type in identifiers:
                assert id_type in available_id_types, "Identifiers must be in " + "; ".join(available_id_types)
            self.mixed_identifiers = True
            self.identifiers = identifiers
        
        # Load the data
        self.raw_data = pd.read_csv(datafile, sep=sep)
        self.data = pd.read_csv(datafile, sep=sep)
        
        # Check that the columns are in the data
        assert self.node_a in self.raw_data.columns, node_a + " is not present as a column in the data"
        assert self.node_b in self.raw_data.columns, node_b + " is not present as a column in the data"
        
        if self.score is not None:
            assert self.score in self.raw_data.columns, score + " is not present as a column in the data"
            self.score_subset = pd.DataFrame()
        
        if self.species is not None:
            assert self.species in self.raw_data.columns, species+ " is not present as a column in the data"
            assert species_code in self.raw_data[self.species].values, "The `species_code` is not present in the `species` column"
            
        # create dictinary for tracking the stats
        self.stats = {"edges":{"raw":len(self.raw_data)}, "nodes":{}}
        
    def clean_data(self):
        self.data = self.data.loc[:, self.cols]
        self.sort_node_pairs()
        self.stats["nodes"]["raw"] = self.get_unique_nodes()
        # keep just the desired species data, then remove the column
        if self.species is not None:
            self.data = self.data.loc[self.data[self.species] == self.species_code]
            self.data.drop(columns = self.species, inplace=True)
            self.stats["edges"]["species"] = len(self.data)
            self.cols.remove(self.species)
        # drop duplicates, in the case of dscored edges keep the highest scoring if there is a duplicate
        if self.score is not None:
            self.data.sort_values(by=self.score, ascending=False, inplace=True)
        self.data.drop_duplicates(inplace=True)
        self.stats["edges"]["de-duped1"] = len(self.data)       
            
    def sort_node_pairs(self):
        node_array = self.data.iloc[:, (self.node_a, self.node_b)].to_numpy()
        node_array.sort()
        sorted_data = pd.DataFrame(node_array)
        sorted_data.columns = [self.node_a, self.node_b]
        if self.data.shape[1] > 2:
            self.data = pd.concat([sorted_data.reset_index(drop=True), self.data.drop(columns=[self.node_a, self.node_b]).reset_index(drop=True)], axis=1)
        else:
            self.data = sorted_data
    
    def get_unique_nodes(self):
        n_a = self.data[self.node_a].unique
        n_b = self.data[self.node_b].unique
        all_nodes = set(n_a).union(set(n_b))
        return all_nodes
    
    def convert_nodes(self):
        """Updates and converts node identifiers, and then converts all edges to new identifiers. When the input data contains multiple input
        types, they are converted in order. 
        """
        gene_map = {}
        all_nodes = self.get_unique_nodes()
        if self.mixed_identifiers:
            unmapped = all_nodes
            for id_type in self.identifiers:
                id_map, unmapped = self.get_node_conversion_map(unmapped, id_type, self.target_id_type)
                gene_map = {**gene_map, **id_map}
        else:
            node_map = self.get_node_conversion_map(all_nodes, id_type, self.target_id_type)
        unmapped_nodes = all_nodes.difference(set(node_map.keys()))
        print(len(unmapped_nodes), "nodes are unmapped")
        print("UNMAPPED NODES")
        print("\n".join(unmapped_nodes))
        print("END UNMAPPED NODES")
        self.stats["nodes"]["unmapped"] = len(unmapped_nodes)
        self.stats["nodes"]["mapped"] = len(all_nodes) - len(unmapped_nodes)
        self.convert_edges(node_map)
    
    def get_node_conversion_map(self, nodes, initial_id, target_id):
        """_summary_

        Args:
            nodes (_type_): _description_
            initial_id (_type_): _description_
            target_id (_type_): _description_

        Returns:
            _type_: _description_
        """
        # this must return a 1:1 mapping
        if initial_id == target_id:
            print("No identifier conversion needed, but IDs will be checked to match latest version")
            node_map = update_nodes(nodes, initial_id)
            return node_map, nodes.difference(set(node_map.keys()))
        else:
            updated_id_map = update_nodes(nodes, initial_id)
            node_map = convert_nodes(list(updated_id_map.values()), initial_id, target_id)
            
            updated_df = pd.DataFrame() # to allow conversion of successive mappings
            return #TODO
    
    def convert_edges(self, node_map):
        pass
        #TODO edge and node stats
    
    def subset_on_score(self, percentile):
        if self.score is not None:
            pass
        # TODO edge and node stats
        
    def write_network_data(self, outpath, percentile=""):
        self.data.to_csv(outpath + self.net_name + "txt", sep="\t")
        if self.score is not None:
            self.score_subset.to_csv(outpath + self.net_name + "_"+ str(percentile) + "txt", sep="\t")
        

if __name__=="__main__":
    df = pd.DataFrame({"node1":[1,2,3,4,5], "node2":[6,3,1,2,6], "species": [9606, 2, 9606, 9606, 9606],
                       "score":[6,3,4,2,1]})
    df.to_csv("testfile.tsv", sep="\t")
    nd = NetworkData("testfile.tsv", "node1", "node2", "Entrez", "Symbol" , score="score", species="species", species_code=9606)
    nd.clean_data()
    
