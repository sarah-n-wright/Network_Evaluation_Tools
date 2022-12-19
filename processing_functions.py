import pandas as pd
import numpy as np
import query_uniprot as uni
import query_hgnc as hgnc
import query_ensembl as ensg
import mygene
import csv

from datetime import datetime

def print_time(message=""):
    now = datetime.now()
    print(now.strftime("%H:%M:%S"), message)
    
    
def print_delta(a, b, message=""):
    diff = b-a
    print(str(diff), message)
    
    
class Timer:
    def __init__(self):
        self.start_times = {}
        self.finish_times = {}
        self.elapsed_times = {}
        self.tasks = []
        self.current_task_stack = []
        self.indents = {}
        
    def start(self, taskstr):
        if taskstr in self.start_times.keys():
            taskstr=taskstr + "1"
            i=1
            while taskstr in self.start_times.keys():
                i += 1
                taskstr = taskstr[0:-1] + str(i) 
        self.current_task_stack.append(taskstr)
        self.indents[taskstr] = len(self.current_task_stack) - 1
        self.tasks.append(taskstr)
        self.start_times[taskstr] = datetime.now()
   
        
    def end(self, taskstr):
        if taskstr in self.finish_times:
            matching_tasks = [task for task in self.start_times.keys() if taskstr in task]
            taskstr = matching_tasks[-1]
        self.current_task_stack.remove(taskstr)
        self.finish_times[taskstr] = datetime.now()
        self.elapsed_times[taskstr] = str(self.finish_times[taskstr] - self.start_times[taskstr])
        
    def print_all_times(self):
        try:
            for task in self.tasks:
                if task not in self.elapsed_times:
                    self.end(task)
                if self.indents[task] > 0:
                    print("".join(["|", "---"*self.indents[task], ">"]),self.elapsed_times[task], task)
                else:
                    print(self.elapsed_times[task], task)
        except:
            print(self.elapsed_times)
        

def update_nodes(nodes, id_type, keep="present", timer=None):
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
    if timer is None:
        timer = Timer()
    timer.start("Update Nodes")
    updated_node_map = {}
    if id_type == "Uniprot":
        results, failed = uni.perform_uniprot_query(ids = nodes, from_db="UniProtKB_AC-ID", to_db="Uniprot")
        print("DUPLICATED UNIPROT MAPPINGS")
        print(results.loc[results.duplicated()])
        #remove duplicates
        results = results.drop_duplicates(subset=["from"])
    elif id_type == "Symbol":
        results, failed = hgnc.perform_hgnc_query(nodes, "Symbol", "Symbol")
        results = pd.DataFrame.from_dict(results, orient="index", columns = ["to"])
        results["from"] = results.index.values
    elif id_type == "Ensembl":
        results, failed = ensg.get_latest_ensembl_id(nodes)
    elif id_type == "Entrez":
        results, failed = get_mygene(nodes, 'entrezgene')
        results.index = results.index.astype(str)
        results["to"] = results["to"].astype(str)
        results["from"] = results["from"].astype(str)

    # process the final data
    if keep == "updated":
        results = results.loc[results["from"] != results["to"]]
    elif keep == "all":
        results = pd.concat([results, pd.DataFrame({"from": failed, "to": np.nan})], axis=0)
    # convert to a dictionary
    results = results.set_index('from')
    updated_node_map = results["to"].to_dict()
    timer.end("Update Nodes")
    return updated_node_map, failed
    
def convert_node_ids(nodes, initial_id, target_id, timer=None):
    """ Converts nodes between two different identifier types

    Args:
        nodes (set): Set of nodes to be converted
        initial_id (str): Identifier type of input nodes
        target_id (str): Identifier type to be converted to
        
    Returns:
        dict: mapping between input nodes and new identifier
        set: nodes that were not able to be mapped to new identifiers.
    """
    # TODO can any of these be looped together?
    # TODO for multiple Ids will need to split and do each separately. 
    mygene_fields = {"Symbol": "symbol", "Entrez": 'entrezgene', "Uniprot": "uniprot", "Ensembl": "ensembl.gene"}
    if timer is None:
        timer = Timer()
    timer.start("Convert node IDs")
    if (initial_id == "Symbol") and (target_id == 'Entrez'):
        # we will use mygeneinfo to do the conversion...
        converted_df, missing = query_mygene(nodes, "symbol", "entrezgene")
        converted_node_map = converted_df.dropna(subset=["_id"])["_id"].to_dict()
        if len(missing) > 0:
            missing_map, still_missing = hgnc.query_other_id(missing, "Entrez")
            converted_node_map = {**converted_node_map, **missing_map}
        else:
            still_missing=missing
    elif (initial_id == "Entrez") and (target_id == "Symbol"):
        converted_df, still_missing = get_mygene(nodes, "symbol")
        converted_df["from"] = converted_df["from"].astype(str)
        converted_df.index = converted_df.index.astype(str)
        converted_node_map = converted_df["to"].to_dict()
    elif (initial_id == "Uniprot"):
        converted_df, still_missing = uni.perform_uniprot_query(ids = nodes, from_db="UniProtKB_AC-ID", to_db=target_id)
        converted_df['from'] = converted_df['from'].astype(str)
        converted_df.index = converted_df["from"]
        converted_df['to'] = converted_df['to'].astype(str)
        converted_node_map = converted_df['to'].to_dict()
        # secondary check for missing ids. 
        if len(still_missing) > 0:
            secondary, still_missing = query_mygene(still_missing, scopes='uniprot', fields=mygene_fields[target_id])
            secondary = secondary.dropna(subset=[mygene_fields[target_id]])
            secondary[mygene_fields[target_id]] = secondary[mygene_fields[target_id]].astype(str)
            # add to the node map
            converted_node_map = {**converted_node_map, **secondary[mygene_fields[target_id]].to_dict()}
        # third check for missing ids (convert first to symbol via uiport and then to Entrez
        if len(still_missing) > 0:
            missing_df, still_missing = uni.perform_uniprot_query(ids = set(still_missing), from_db="UniProtKB_AC-ID", to_db='Symbol')
            tertiary, still_missing = query_mygene(missing_df["to"].values, scopes='symbol', fields=mygene_fields[target_id])
            missing_df.index = missing_df["to"]
            missing_dict= missing_df["from"].to_dict()
            tertiary = tertiary.dropna(subset=[mygene_fields[target_id]])
            tertiary_dict = {}
            for node in tertiary.index.values:
                tertiary_dict[missing_dict[node]] = tertiary.loc[node, mygene_fields[target_id]]
            converted_node_map = {**converted_node_map, **tertiary_dict}
        
        
        
    elif initial_id == "Ensembl":
        converted_df, still_missing = query_mygene(nodes, scopes=mygene_fields["Ensembl"], fields=mygene_fields[target_id])
        converted_df = converted_df.dropna(subset=[mygene_fields[target_id]])
        converted_node_map = converted_df[mygene_fields[target_id]].to_dict()
        
    timer.end("Convert node IDs")
    return converted_node_map, still_missing

def query_mygene(gene_list, scopes, fields):
    mg = mygene.MyGeneInfo()
    results_df = mg.querymany(qterms=gene_list, scopes=scopes, fields=fields, species='human', 
                                returnall=True, verbose=False, as_dataframe=True, entrezonly=True)
    mapped = results_df["out"]
    dups = results_df["dup"]
    missing = results_df["missing"]
    unmapped = []
    if len(dups) > 0:
        unmapped += list(results_df["dup"]["query"].values)
    if len(missing) > 0:
        unmapped += list(results_df["missing"]["query"].values)
    return mapped, unmapped

def get_mygene(gene_list, target_id):
    mg = mygene.MyGeneInfo()
    results = mg.getgenes(gene_list, as_dataframe=True, fields=target_id)
    failed = list(results.loc[results[target_id].isna()].index.values)
    results = results.dropna(subset=[target_id])
    results["from"] = results.index.values
    results = results.loc[:, ("from", target_id)]
    results.columns = ["from", "to"]
    
    return results, failed


class NetworkData:
    def __init__(self, datafile, node_a, node_b, identifiers, target_id_type, net_name, score=None, species=None, species_code=None, sep="\t", header=0, test_mode=False):
        self.T = Timer()
        self.T.start("Total")
        available_id_types = ["Symbol", "Entrez", "Uniprot", "Ensembl"]
        assert species is None or ((species is not None) and (species_code is not None)), "If species column is given, a corresponding value representing the species must be given"
        assert target_id_type in available_id_types, "Target identifier must be in " + "; ".join(available_id_types)
        
        self.species = species
        self.two_species_cols = True if type(species) == list else False
            
        self.score = score
        self.node_a = node_a
        self.node_b = node_b
        all_cols = [node_a, node_b, score] + species if self.two_species_cols else [node_a, node_b, score, species]
        self.cols = [col for col in all_cols if col is not None]
        
        self.species_code = species_code
        
        self.target_id_type = target_id_type
        self.net_name = net_name
        
        # check the identifiers given
        if (type(identifiers) == str) or (len(identifiers) == 1):
            if type(identifiers) == list:
                identifiers = identifiers[0]
            assert identifiers in available_id_types, "Identifiers must be in " + "; ".join(available_id_types)
            self.mixed_identifiers = False
            self.identifiers = identifiers
        elif (type(identifiers) == list) and (len(identifiers) > 1):
            for id_type in identifiers:
                assert id_type in available_id_types, "Identifiers must be in " + "; ".join(available_id_types)
            self.mixed_identifiers = True
            self.identifiers = identifiers
        
        # Load the data
        self.T.start("Load data")
        if test_mode:
            self.raw_data = pd.read_csv(datafile, sep=sep, header=header, engine='python', nrows=100, quoting=csv.QUOTE_ALL, quotechar='"')
            self.data = pd.read_csv(datafile, sep=sep, header=header, engine='python', nrows=100, quoting=csv.QUOTE_ALL, quotechar='"')
        else:
            self.raw_data = pd.read_csv(datafile, sep=sep, header=header, engine='python', quoting=csv.QUOTE_ALL, quotechar='"')
            self.data = pd.read_csv(datafile, sep=sep, header=header, engine='python', quoting=csv.QUOTE_ALL, quotechar='"')
        self.T.end("Load data")
        # Check that the columns are in the data
        print(self.raw_data.head())
        assert self.node_a in self.raw_data.columns,str(node_a) + " is not present as a column in the data"
        assert self.node_b in self.raw_data.columns, str(node_b) + " is not present as a column in the data"
        
        if self.score is not None:
            assert self.score in self.raw_data.columns, str(score) + " is not present as a column in the data"
            self.score_subset = pd.DataFrame()
        
        if self.species is not None:
            if self.two_species_cols:
                assert self.species[0] in self.raw_data.columns, str(self.species[0])+ " is not present as a column in the data"
                assert self.species[1] in self.raw_data.columns, str(self.species[1])+ " is not present as a column in the data"

            else:
                assert self.species in self.raw_data.columns, str(species)+ " is not present as a column in the data"

            #print(self.raw_data[self.species].values)
            print(self.species_code, type(self.species_code))
            assert species_code in self.raw_data[self.species].values, "The `species_code` "+ str(species_code) + " is not present in the `species` column"
            
        # create dictinary for tracking the stats
        self.stats = {"edges":{"raw":len(self.raw_data)}, "nodes":{}}
    
    #@profile
    def clean_data(self):
        self.T.start("Clean data")
        self.data = self.data.loc[:, self.cols]
        self.data = self.data.dropna(subset=[self.node_a, self.node_b])
        self.sort_node_pairs()
        self.stats["nodes"]["raw"] = len(self.get_unique_nodes())
        # keep just the desired species data, then remove the column
        if self.species is not None:
            if self.two_species_cols:
                self.data = self.data.loc[((self.data[self.species[0]] == self.species_code) &
                                        (self.data[self.species[1]] == self.species_code))]
                for col in self.species:
                    self.cols.remove(col)
            else:
                self.data = self.data.loc[self.data[self.species] == self.species_code]
                self.cols.remove(self.species)
            self.data.drop(columns = self.species, inplace=True)
            self.stats["edges"]["species"] = len(self.data)
            
        # drop duplicates, in the case of dscored edges keep the highest scoring if there is a duplicate
        if self.score is not None:
            if len(self.data[self.score].unique()) <= 1:
                self.data.drop(columns=[self.score])
                self.cols.remove(self.score)
                self.score = None
            else:
                self.data.sort_values(by=self.score, ascending=False, inplace=True)
        self.data.drop_duplicates(inplace=True)
        self.stats["edges"]["de-duped1"] = len(self.data)
        
        # convert identifiers to strings
        self.data[self.node_a] = self.data[self.node_a].astype(str)
        self.data[self.node_b] = self.data[self.node_b].astype(str)
        self.T.end("Clean data")
                
                
            
    def sort_node_pairs(self):
        self.T.start("Sort nodes")
        node_array = self.data.loc[:, (self.node_a, self.node_b)].to_numpy()
        try:
            node_array.sort()
        except TypeError:
            print("TYPE ERROR")
            print(self.data.head())
            print(node_array)
            raise
        sorted_data = pd.DataFrame(node_array)
        sorted_data.columns = [self.node_a, self.node_b]
        if self.data.shape[1] > 2:
            self.data = pd.concat([sorted_data.reset_index(drop=True), self.data.drop(columns=[self.node_a, self.node_b]).reset_index(drop=True)], axis=1)
        else:
            self.data = sorted_data
        self.T.end("Sort nodes")
    
    def get_unique_nodes(self, score_subset=False):
        self.T.start("Unique nodes")
        if not score_subset:
            n_a = self.data[self.node_a].unique()
            n_b = self.data[self.node_b].unique()
        else:
            n_a = self.score_subset[self.node_a].unique()
            n_b = self.score_subset[self.node_b].unique()
        all_nodes = set(n_a).union(set(n_b))
        self.T.end("Unique nodes")
        return all_nodes
    
    #@profile
    def convert_nodes(self):
        """Updates and converts node identifiers, and then converts all edges to new identifiers. When the input data contains multiple input
        types, they are converted in order. 
        """
        self.T.start("Convert nodes")
        gene_map = {}
        all_nodes = self.get_unique_nodes()
        if self.mixed_identifiers:
            unmapped = all_nodes
            for id_type in self.identifiers:
                if len(unmapped) > 0:
                    id_map = self.get_node_conversion_map(unmapped, id_type, self.target_id_type)
                    gene_map = {**gene_map, **id_map}
                    unmapped = all_nodes.difference(set(gene_map.keys()))
            node_map = gene_map
        else:
            node_map = self.get_node_conversion_map(all_nodes, self.identifiers, self.target_id_type)
        unmapped_nodes = all_nodes.difference(set(node_map.keys()))
        print(len(unmapped_nodes), "nodes are unmapped")
        print("UNMAPPED NODES")
        print("\n".join(unmapped_nodes))
        print("END UNMAPPED NODES")
        self.stats["nodes"]["unmapped"] = len(unmapped_nodes)
        self.stats["nodes"]["mapped"] = len(all_nodes) - len(unmapped_nodes)
        self.convert_edges(node_map, unmapped_nodes)
        print("Converted")
        self.T.end("Convert nodes")
    
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
        self.T.start("Get node conversion map")
        if initial_id == target_id:
            print("No identifier conversion needed, but IDs will be checked to match latest version")
            node_map, _ = update_nodes(nodes, initial_id, timer=self.T)
        else:
            # don't keep ids that are not able to be updated/queried in their own identifier
            updated_id_map, _ = update_nodes(nodes, initial_id, timer=self.T)
            converted_map, _ = convert_node_ids(set(updated_id_map.values()), initial_id, target_id, timer=self.T)
            node_map = {}
            for node in updated_id_map.keys():
                if node in converted_map.keys(): # to catch nodes that were updated but not mapped to target identifier
                    node_map[node] = converted_map[updated_id_map[node]]
            #match updated_id_map and node_map to remove intermediate updates
            #updated_df = pd.DataFrame() # to allow conversion of successive mappings
        self.T.end("Get node conversion map")
        return node_map
    
    #@profile
    def convert_edges_old(self, node_map, unmapped_nodes):
        # first remove edges with unmapped nodes
        self.T.start("Convert edges")
        self.T.start("Remove unmapped")
        self.data = self.data.loc[((~self.data[self.node_a].isin(unmapped_nodes)) & (~self.data[self.node_b].isin(unmapped_nodes)))]
        self.T.end("Remove unmapped")
        self.stats['edges']['mapped'] = len(self.data)
        # convert node a and node b
        self.T.start("Replace")
        self.data.replace(to_replace={self.node_a:node_map, self.node_b: node_map}, inplace=True)
        self.T.end("Replace")
        self.T.end("Convert edges")
        
    def convert_edges(self, node_map, unmapped_nodes):
        # did I need the deep copy??
        self.T.start("Convert edges")
        self.T.start("Remove unmapped")
        self.data = self.data.loc[((~self.data[self.node_a].isin(unmapped_nodes)) & (~self.data[self.node_b].isin(unmapped_nodes)))]
        self.T.end("Remove unmapped")
        self.stats['edges']['mapped'] = len(self.data)
        node_names = [self.node_a, self.node_b]
        self.data.dropna(subset=node_names, inplace=True)
        self.T.start("Map")
        self.data[node_names[0]] = self.data[node_names[0]].map(node_map, na_action="ignore")
        self.data[node_names[1]] = self.data[node_names[1]].map(node_map, na_action="ignore")
        self.T.end("Map")
        self.T.end("Convert edges")

    def subset_on_score(self, score_col, percentile):
        if self.score is not None:
            cutoff = np.percentile(self.data[self.score].values, percentile)
            # are all scores high=better?
            self.score_subset = self.data.loc[self.data[self.score] > cutoff]
            self.stats['edges']['score_subset'] = len(self.score_subset)
            self.stats['nodes']['score_subset'] = len(self.get_unique_nodes(score_subset=True))
        
        
    def write_network_data(self, outpath, percentile=95):
        # rename columns
        self.T.start("Write data")
        final_names= {self.node_a: self.target_id_type + "_A", self.node_b: self.target_id_type + "_B"}
        if self.score is not None:
            final_names[self.score] = "Score"
        if self.score is not None:
            self.subset_on_score("Score", percentile)
            
        self.data.rename(columns = final_names, inplace=True)
        self.data.to_csv(outpath + self.net_name + ".txt", sep="\t", index=False)
        
        if self.score is not None:
            self.score_subset.rename(columns = final_names, inplace=True)
            self.score_subset.to_csv(outpath + self.net_name + "_"+ str(percentile) + ".txt", sep="\t", index=False)

        self.T.end("Write data")
        self.T.end("Total")
            
    def write_stats(self, outpath):
        stats_df = pd.DataFrame.from_dict(self.stats)
        stats_df.to_csv(outpath + self.net_name + ".stats", sep="\t", index=True)
        self.T.print_all_times()
        


#if __name__=="__main__":

if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/HI-II-14.tsv"
    nd = NetworkData(datafile, node_a=0, node_b=1, 
                    target_id_type='Entrez',  identifiers='Ensembl',
                    header=None, net_name="HI_test", test_mode=True, sep="\t")
    #nd.data = pd.concat([nd.data, pd.DataFrame({"Gene_A": 'A6NKP2', "Gene_B": "P48506", "Weight": 1}, index=[101])])
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)

if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/GeneMania_2021_COMBINED.DEFAULT_NETWORKS.BP_COMBINING.txt"
    nd = NetworkData(datafile, node_a="Gene_A", node_b="Gene_B", 
                    target_id_type='Entrez',  identifiers=['Uniprot', 'Ensembl'], score="Weight",
                    header=0, net_name="genemania_test", test_mode=True, sep="\t")
    #nd.data = pd.concat([nd.data, pd.DataFrame({"Gene_A": ['A6NKP2', "ENSG00000260342"], "Gene_B": ["Q96PG1", "ENSG00000188897"], "Weight": [1,1]}, index=[100,101])])
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)

if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/BioPlex_293T_Network_10K_Dec_2019.tsv"
    nd = NetworkData(datafile, node_a="GeneA", node_b="GeneB", 
                    target_id_type='Symbol',  identifiers='Entrez', 
                    header=0, net_name="bioplex_test", test_mode=True, sep="\t")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)


if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/BIOGRID/BIOGRID-ORGANISM-Homo_sapiens-4.4.213.tab3.txt"
    nd = NetworkData(datafile, node_a="Entrez Gene Interactor A", node_b="Entrez Gene Interactor B", species="Organism ID Interactor A", 
                    target_id_type='Symbol',  identifiers='Entrez', species_code=9606, 
                    score="Score", header=0, net_name="BIOGRID_entrex_test", test_mode=True, sep="\t")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)