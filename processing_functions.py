import pandas as pd
import numpy as np
import query_uniprot as uni
import query_hgnc as hgnc
import query_ensembl as ensg
import mygene
import csv
import re
from itertools import combinations

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
        query_nodes = [node for node in nodes if (("_HUMAN" in node) or (re.fullmatch("[a-zA-Z0-9\.-]+", node) is not None))]
        results, failed = uni.perform_uniprot_query(ids = query_nodes, from_db="UniProtKB_AC-ID", to_db="Uniprot")
        #print("DUPLICATED UNIPROT MAPPINGS")
        #print(results.loc[results.duplicated()])
        #remove duplicates
        failed = list(failed) + [node for node in nodes if node not in query_nodes]
        results = results.drop_duplicates(subset=["from"])
    elif id_type == "Symbol":
        exclude_prefix_suffix = ["CHEBI:", "_HUMAN"]
        query_nodes = [node for node in list(nodes) if re.search("|".join(exclude_prefix_suffix), node) is None]
        results, failed = hgnc.perform_hgnc_query(query_nodes, "Symbol", "Symbol")
        results = pd.DataFrame.from_dict(results, orient="index", columns = ["to"])
        results["from"] = results.index.values
        failed = list(failed) + [node for node in list(nodes) if re.search("|".join(exclude_prefix_suffix), node) is not None]
    elif id_type in ["Ensembl", "EnsemblProtein"]:
        results, failed = ensg.get_latest_ensembl_id(nodes)
    elif id_type == "Entrez":
        results, failed = get_mygene(nodes, 'entrezgene')
        results.index = results.index.astype(str)
        results["to"] = results["to"].astype(str)
        results["from"] = results["from"].astype(str)
    elif id_type == "DIP":
        dip_ids = [n for n in nodes if "DIP-" in n]
        results = pd.DataFrame({"to":dip_ids, "from":dip_ids})
        failed = [n for n in nodes if "DIP-" not in n]
    elif id_type == "Refseq":
        refseq_ids = [n for n in nodes if "_" in n]
        results = pd.DataFrame({"to":refseq_ids, "from":refseq_ids})
        failed = [n for n in nodes if "_" not in n]

    # process the final data
    if keep == "updated":
        results = results.loc[results["from"] != results["to"]]
    elif keep == "all":
        results = pd.concat([results, pd.DataFrame({"from": failed, "to": np.nan})], axis=0)
    # convert to a dictionary
    if len(results) > 0:
        results = results.set_index('from')
        updated_node_map = results["to"].to_dict()
    else:
        updated_node_map = {}
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
    mygene_fields = {"Symbol": "symbol", "Entrez": 'entrezgene', "Uniprot": "uniprot", "Ensembl": "ensembl.gene",
                    "Refseq":"refseq", "EnsemblProtein":"ensembl.protein"}
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
    elif (initial_id == "Uniprot") or (initial_id == "DIP"):
        if (initial_id == "DIP"):
            dip_df, missing_dip = uni.perform_uniprot_query(ids = nodes, from_db="DIP", to_db="Uniprot")
            dip_df['from'] = dip_df['from'].astype(str)
            dip_df['to'] = dip_df['to'].astype(str)
            if (target_id == "Uniprot"):
                return
            else:
                nodes = dip_df["to"].unique()
        
        converted_df, still_missing = uni.perform_uniprot_query(ids = nodes, from_db="UniProtKB_AC-ID", to_db=target_id)
        converted_df['from'] = converted_df['from'].astype(str)
        converted_df.index = converted_df["from"]
        converted_df['to'] = converted_df['to'].astype(str)
        converted_node_map = converted_df['to'].to_dict()
        # secondary check for missing ids. 
        if len(still_missing) > 0:
            secondary, still_missing = query_mygene(still_missing, scopes='uniprot', fields=mygene_fields[target_id])
            if mygene_fields[target_id] in secondary.columns:  #otherwise none were found
                
                secondary = secondary.dropna(subset=[mygene_fields[target_id]])
                secondary[mygene_fields[target_id]] = secondary[mygene_fields[target_id]].astype(str)
                # add to the node map
                converted_node_map = {**converted_node_map, **secondary[mygene_fields[target_id]].to_dict()}
        # third check for missing ids (convert first to symbol via uniprot and then to Entrez
        if len(still_missing) > 0:
            missing_df, still_missing = uni.perform_uniprot_query(ids = set(still_missing), from_db="UniProtKB_AC-ID", to_db='Symbol')
            if len(missing_df) > 0:
                tertiary, still_missing = query_mygene(missing_df["to"].values, scopes='symbol', fields=mygene_fields[target_id])
                still_missing = list(missing_df[missing_df["to"].isin(still_missing)]["from"])
                missing_df.index = missing_df["from"]
                missing_dict= missing_df["to"].to_dict()
                if mygene_fields[target_id] in tertiary.columns:
                    tertiary = tertiary.dropna(subset=[mygene_fields[target_id]])
                    tertiary["input"] = tertiary.index.values
                    tertiary = tertiary.drop_duplicates(subset=[mygene_fields[target_id], "input"])
                    tertiary_dict = {}
                    for node in missing_dict:
                        if missing_dict[node] in tertiary.index.values:
                            tertiary_dict[node] = tertiary.loc[missing_dict[node], mygene_fields[target_id]]
                    converted_node_map = {**converted_node_map, **tertiary_dict}
        if (initial_id == "DIP"):
            uniprot_df = pd.DataFrame.from_dict(converted_node_map, orient="index", columns=["target"])
            full_df = dip_df.join(uniprot_df, on="to", how="left")
            full_df.dropna(inplace=True)
            full_df.index = full_df["from"]
            converted_node_map = full_df["target"].to_dict()
            still_missing = missing_dip + list(set(dip_df["from"]).difference(set(full_df["from"])))
        
    elif initial_id in ["Ensembl", "Refseq", "EnsemblProtein"]:
        converted_df, still_missing = query_mygene(nodes, scopes=mygene_fields[initial_id], fields=mygene_fields[target_id])
        converted_df = converted_df.dropna(subset=[mygene_fields[target_id]])
        converted_node_map = converted_df[mygene_fields[target_id]].to_dict()
        
    timer.end("Convert node IDs")
    return converted_node_map, still_missing

def query_mygene(gene_list, scopes, fields, retries=10):
    mg = mygene.MyGeneInfo()
    for retry in range(retries):
        try:
            results_df = mg.querymany(qterms=gene_list, scopes=scopes, fields=fields, species='human', 
                                returnall=True, verbose=False, as_dataframe=True, entrezonly=True)
            break
        except Exception as e:
            if retry < retries - 1:
                #print("PF", gene_list)
                print(f"Retrying mg.querymany: {e}")
            else:
                print("Max retries reach for mg.querymany")
                raise e
    mapped = results_df["out"]
    dups = results_df["dup"]
    missing = results_df["missing"]
    unmapped = []
    if len(dups) > 0:
        unmapped += list(results_df["dup"]["query"].values)
    if len(missing) > 0:
        unmapped += list(results_df["missing"]["query"].values)
    return mapped, unmapped

def get_mygene(gene_list, target_id, retries=10):
    mg = mygene.MyGeneInfo()
    for retry in range(retries):
        try:
            results = mg.getgenes(gene_list, as_dataframe=True, fields=target_id)
            break
        except Exception as e:
            if retry < retries - 1:
                print(f"Retrying mg.getgenes: {e}")
            else:
                print("Max retries reach for mg.getgenes")
                raise e
    failed = list(results.loc[results[target_id].isna()].index.values)
    results = results.dropna(subset=[target_id])
    results["from"] = results.index.values
    results = results.loc[:, ("from", target_id)]
    results.columns = ["from", "to"]
    
    return results, failed

    
def extract_id_with_prefix(id_str, pref_sep):
    add_back_prefixes = ["DIP-", "ESNG", "ESNP"]
    if isinstance(pref_sep[0], str):
        # Handle case where pref_sep contains a single string
        prefix = pref_sep[0]
        separator = ""
        if prefix in id_str:  # otherwise try next prefix if available
            id_str = id_str.split(prefix)[1]
            if separator and separator in id_str:
                id_str = id_str.split(separator)[0]
            id_str = prefix + id_str if prefix in add_back_prefixes else id_str
            return id_str
    else:
        # Handle case where pref_sep contains a list of tuples
        for prefix, separator in pref_sep:
            if prefix in id_str:  # otherwise try next prefix if available
                id_str = id_str.split(prefix)[1]
                if separator and separator in id_str:
                    id_str = id_str.split(separator)[0]
                id_str = prefix + id_str if prefix in add_back_prefixes else id_str
                return id_str
    # if no prefixes found return NA
    return pd.NA


def clean_score(score):
    try:
        score = float(score)
        return score
    except ValueError:
        if isinstance(score, str):
            if re.search(r"^\D+$", score):  # no numeric score can be found
                return pd.NA
            else: #look for the first numeric occurence
                m = re.search(r"(\d+\.?\d*(?:e[+-]?\d+)?)", score).group(1)
                return float(m)
        else:
            return pd.NA
    


class NetworkData:
    def __init__(self, datafile, node_a, node_b, identifiers, target_id_type, net_name, score=None, species=None, species_code=None, sep="\t", header=0, test_mode=False, prefixes=None):
        test_size=10000
        self.T = Timer()
        self.T.start("Total")
        available_id_types = ["Symbol", "Entrez", "Uniprot", "Ensembl", "DIP", "Refseq", "EnsemblProtein"]
        assert species is None or ((species is not None) and (species_code is not None)), "If species column is given, a corresponding value representing the species must be given"
        assert target_id_type in available_id_types, "Target identifier must be in " + "; ".join(available_id_types)
        
        self.species = species
        self.two_species_cols = True if type(species) == list else False
            
        self.score = score
        if (self.score is not None) and (type(self.score) == str):
            if "[p]" in self.score:
                self.score = self.score.split("[p]")[0]
                self.score = int(self.score) if self.score.isnumeric() else self.score
                self.score_is_p = True
            else:
                self.score_is_p = False
        else:
            self.score_is_p = False
            
        self.node_a = node_a
        self.node_b = node_b
        self.complexes = True if self.node_b is None else False
        
        all_cols = [self.node_a, self.node_b, self.score] + self.species if self.two_species_cols else [self.node_a, self.node_b, self.score, self.species]
        self.cols = [col for col in all_cols if col is not None]
        
        self.species_code = species_code
        
        self.target_id_type = target_id_type
        self.net_name = net_name
        self.prefix = prefixes
        
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
        chars = ['"', "@"]
        for i, quotechar in enumerate(chars):
            try:
                if test_mode:
                    self.raw_data = pd.read_csv(datafile, sep=sep, header=header, engine='python', nrows=test_size, quoting=csv.QUOTE_ALL, quotechar=quotechar)
                    self.data = pd.read_csv(datafile, sep=sep, header=header, engine='python', nrows=test_size, quoting=csv.QUOTE_ALL, quotechar=quotechar)
                    if (header is not None) and any(["Unnamed" in col for col in self.raw_data.columns]):
                # Possible malformed file, first column being incorrectly used as an index? Force it to not assign index
                        self.raw_data = pd.read_csv(datafile, sep=sep, index_col=False, header=header, engine='python', nrows=test_size, quoting=csv.QUOTE_ALL, quotechar=quotechar)
                        self.data = pd.read_csv(datafile, sep=sep, index_col=False, header=header, engine='python', nrows=test_size, quoting=csv.QUOTE_ALL, quotechar=quotechar)
                
                else:
                    self.raw_data = pd.read_csv(datafile, sep=sep, header=header, engine='python', quoting=csv.QUOTE_ALL, quotechar=quotechar)
                    self.data = pd.read_csv(datafile, sep=sep, header=header, engine='python', quoting=csv.QUOTE_ALL, quotechar=quotechar)
                    if (header is not None) and any(["Unnamed" in col for col in self.raw_data.columns]):
                        self.raw_data = pd.read_csv(datafile, sep=sep, index_col=False, header=header, engine='python', quoting=csv.QUOTE_ALL, quotechar=quotechar)
                        self.data = pd.read_csv(datafile, sep=sep, index_col=False, header=header, engine='python', quoting=csv.QUOTE_ALL, quotechar=quotechar)
                break
            except:
                print("WARNING: Data not able to be loaded with quotechar", quotechar)
                if i == len(chars) - 1:
                    raise pd.errors.ParserError("Data could not be loaded with available quote characters")
                
        # Intialize an index column:
        self.data[net_name+"_ID"] = self.data.index.values
        self.cols.append(net_name+"_ID")
        self.T.end("Load data")
        # Check that the columns are in the data
        #print(self.raw_data.head())
        assert self.node_a in self.raw_data.columns,str(node_a) + " is not present as a column in the data"
        assert (self.node_b in self.raw_data.columns, str(node_b)) or self.node_b is None, + " is not present as a column in the data"
        
        if self.score is not None:
            assert self.score in self.raw_data.columns, str(score) + " is not present as a column in the data"
            self.score_subset = pd.DataFrame()
        
        #print(self.raw_data.columns)
        #print(self.raw_data.iloc[1])
        if self.species is not None:
            if self.two_species_cols:
                assert self.species[0] in self.raw_data.columns, str(self.species[0])+ " is not present as a column in the data"
                assert self.species[1] in self.raw_data.columns, str(self.species[1])+ " is not present as a column in the data"
                #print(self.species_code, type(self.species_code))
                #print(self.raw_data[self.species[0]].values)
                assert species_code in self.raw_data[self.species[0]].values, "The `species_code` "+ str(species_code) + " is not present in the `species[0]` column"
                assert species_code in self.raw_data[self.species[1]].values, "The `species_code` "+ str(species_code) + " is not present in the `species[1]` column"
                
            else:
                assert self.species in self.raw_data.columns, str(species)+ " is not present as a column in the data"
                #print(self.species_code, type(self.species_code))
                assert species_code in self.raw_data[self.species].values, "The `species_code` "+ str(species_code) + " is not present in the `species` column"

            
        # create dictinary for tracking the stats
        self.stats = {"edges":{"raw":len(self.raw_data)}, "nodes":{}}
    
    #@profile
    def clean_data(self):
        self.T.start("Clean data")
        self.data = self.data.loc[:, self.cols]
        if self.complexes:
            self.binarize_complexes()
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
            elif self.score_is_p:
                # set pseudo count to 1/100 smallest non-zero p-value.
                lowest_score = self.data[self.data[self.score] > 0][self.score].min()
                self.data[self.score] = self.data[self.score].apply(lambda x: -1 * np.log10(x + lowest_score/100))
                self.data.sort_values(by=self.score, ascending=False, inplace=True)
            else:
                try:
                    self.data[self.score] = self.data[self.score].astype(float)
                except ValueError:
                    self.data[self.score] = self.data[self.score].apply(lambda x: clean_score(x))
                self.data.sort_values(by=self.score, ascending=False, inplace=True)
        self.extract_from_prefixes()        
        self.stats["edges"]["Node prefix"] = len(self.data)
        self.data.drop_duplicates(inplace=True)
        self.stats["edges"]["de-duped1"] = len(self.data)
        
        # convert identifiers to strings
        self.data[self.node_a] = self.data[self.node_a].astype(str)
        self.data[self.node_b] = self.data[self.node_b].astype(str)
        self.T.end("Clean data")
        
    def binarize_complexes(self):
        results = []
        separator = None
        idx=0
        while separator is None:
            comp = self.data.iloc[idx][self.node_a]
            try:
                if "_HUMAN" in comp:
                    separator = re.search("[:,\|]", comp).group()
                else:
                    separator = re.search("[:,\|_]", comp).group()
            except AttributeError:
                idx += 1
        # Iterate over the rows of the input data frame to create a mapping of complex string to pairwise ids
        for comp in self.data[self.node_a]:
            # Get the list of IDs and the values for the other columns
            complex_id_list = comp.split(separator)
            if len(complex_id_list) > 1:
                # Iterate over all pairwise combinations of IDs
                for id1, id2 in combinations(complex_id_list, 2):
                    # Append the pairwise combination to the results list
                    results.append((comp, id1, id2))
            else:
                print("SELF EDGE:", comp)
        results_df = pd.DataFrame(results, columns=[self.node_a, "NodeA", "NodeB"])
        self.data = pd.merge(self.data, results_df, on=self.node_a, how="right")
        self.cols.remove(self.node_a)
        self.data.drop(columns=self.node_a, inplace=True)
        self.node_a = "NodeA"
        self.node_b = "NodeB"
        self.cols = self.cols + ["NodeA", "NodeB"]
        

                
    def extract_from_prefixes(self):
        if (self.prefix is not None) and (not self.mixed_identifiers):
            new_nodes_a = self.data[self.node_a].apply(lambda x: extract_id_with_prefix(x, self.prefix))
            new_nodes_b = self.data[self.node_b].apply(lambda x: extract_id_with_prefix(x, self.prefix))
            self.data[self.node_a] = new_nodes_a
            self.data[self.node_b] = new_nodes_b
            self.data.dropna(subset=[self.node_a,self.node_b])  
            
    def sort_node_pairs(self):
        self.T.start("Sort nodes")
        node_array = self.data.loc[:, (self.node_a, self.node_b)].to_numpy()
        try:
            node_array.sort()
        except TypeError as e:
            print("TYPE ERROR")
            print(self.data.head())
            print(node_array)
            raise e
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
        #add_nodes = ['CBWD2_HUMAN', 'BUB1_HUMAN', 'MYO3A_HUMAN', 'NUD10_HUMAN', 'NHLC1_HUMAN', 'SRP72_HUMAN', 'NOX4_HUMAN', 'CACO1_HUMAN', 'HDGF_HUMAN', 'ZC3HE_HUMAN', 'ABI2_HUMAN', 'UBN1_HUMAN', 'PI42C_HUMAN', 'DUOX2_HUMAN', 'Myosin 5A', 'SC24C_HUMAN', 'DHE4_HUMAN', 'THIO_HUMAN', 'CIP4_HUMAN', 'PLK4_HUMAN', 'TRI27_HUMAN', 'NOT7', 'Fibrocystin;PARD3', 'KCNA2_HUMAN', 'MTCH1_HUMAN', 'ZSWM7_HUMAN', 'FBLN3_HUMAN', 'PAHX_HUMAN', 'CAH2_HUMAN', 'CTDS2_HUMAN', 'UBP53_HUMAN', 'GRP75_HUMAN', 'TM9S1_HUMAN', 'NMT2_HUMAN', 'STXB1_HUMAN', 'NELL2_HUMAN', 'MMP8_HUMAN', 'U1SBP_HUMAN', 'RPP30_HUMAN', 'WFS1_HUMAN', 'AGO1_HUMAN', 'ERN1_HUMAN', 'CFTR_HUMAN', 'ACOD_HUMAN', 'QARS', 'UB2V1_HUMAN;Ubiquitin conjugating enzyme E2 Kua-UEV', 'CENPB_HUMAN', 'PP2AB_HUMAN', 'FLT3_HUMAN', 'DYL1_HUMAN', 'STX4_HUMAN', 'PNO1_HUMAN', 'HXA2_HUMAN', 'HIST1H2AC', 'NAPSA_HUMAN', 'CDK1_HUMAN', 'PSMF1_HUMAN', 'OLIG3_HUMAN', 'Clathrin adaptor complex AP2, MU subunit', 'PADC1_HUMAN', 'Chorionic gonadotropin beta polypeptide 5;Chorionic gonoadotropin, beta chain', 'FA86C_HUMAN', 'Intercellular adhesion molecule 5', 'GLRX2_HUMAN', 'SHIP1_HUMAN', 'DDX58_HUMAN', 'Transcription initiation factor IIB', 'SPYA_HUMAN', 'FOSB_HUMAN', 'HNF6_HUMAN', 'RD23B_HUMAN', 'CO8A1_HUMAN', 'CX04A_HUMAN', 'TPM1_HUMAN', 'EFHC1_HUMAN', 'COIA1_HUMAN', 'FGF receptor 3', 'Placental ribonuclease inhibitor']
        #all_nodes = all_nodes.union(set(add_nodes))
        if self.mixed_identifiers:
            unmapped = all_nodes
            print("# UNMAPPED", len(unmapped))
            for i, id_type in enumerate(self.identifiers):
                if len(unmapped) > 0:
                    if self.prefix is not None:
                        clean_nodes = {}
                        for node in unmapped:
                            node_no_pref = extract_id_with_prefix(node,[self.prefix[i]] )
                            if not pd.isna(node_no_pref):
                                clean_nodes[node] = node_no_pref
                        clean_map = self.get_node_conversion_map(set(clean_nodes.values()), id_type, self.target_id_type)
                        id_map = {node: clean_map[clean_nodes[node]] for node in clean_nodes if clean_nodes[node] in clean_map}
                    else:
                        id_map = self.get_node_conversion_map(unmapped, id_type, self.target_id_type)
                        #print(len(id_map))
                    # catch ids not in all_nodes
                    gene_map = {**gene_map, **id_map}
                    #print("# UNMAPPED - ID_MAP = ", len(unmapped) - len(id_map))
                    unmapped = all_nodes.difference(set(gene_map.keys()))
                    modified = [node for node in id_map.keys() if node not in all_nodes]
                    if len(modified) > 0:
                        modified_map = {}
                        for node in modified:
                            match = [x for x in unmapped if node in x]
                            if len(match) > 0:
                                for matched_node in match:
                                    unmapped.remove(matched_node)
                                    modified_map[matched_node] = id_map[node]
                        gene_map = {**gene_map, **modified_map}
                    #print("ACTUAL UNMAPPED", len(unmapped))
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
            node_map = {}
            if len(updated_id_map) > 0:
                converted_map, _ = convert_node_ids(set(updated_id_map.values()), initial_id, target_id, timer=self.T)
                for node in updated_id_map.keys():
                    if updated_id_map[node] in converted_map.keys(): # to catch nodes that were updated but not mapped to target identifier
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
            # assign minimum score to NA values 
            replace_na = min(0, self.data[self.score].min())
            self.data.loc[self.data[self.score].isna(), self.score] = replace_na
            cutoff = np.percentile(self.data[self.score].values, percentile)
            # are all scores high=better?
            self.score_subset = self.data.loc[self.data[self.score] > cutoff]
            self.stats['edges']['score_subset'] = len(self.score_subset)
            self.stats['nodes']['score_subset'] = len(self.get_unique_nodes(score_subset=True))
        
        
    def write_network_data(self, outpath, percentile=95):
        # rename columns
        self.T.start("Write data")
        final_names= {self.node_a: self.target_id_type + "_A", self.node_b: self.target_id_type + "_B", self.net_name+"_ID": "ID"}
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
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/PROPER_v1.csv"
    nd = NetworkData(datafile, node_a="Gene1", node_b="Gene2", 
                    target_id_type='Entrez',  identifiers='Symbol', species= "Cell line specificity", 
                    species_code="shared", score="BH-corrected p-value[p]",
                    header=0, net_name="proper", test_mode=True, sep=",")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)


if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/PathwayCommons/PathwayCommons12_uniprot_test.txt"
    nd = NetworkData(datafile, node_a=0, node_b=2, 
                    target_id_type='Entrez',  identifiers=['Uniprot', 'Symbol'], 
                    header=None, net_name="hprd_test", test_mode=True, sep="\t")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)
    
    
if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/ConsensusPathDB_human_PPI_v35.tsv"
    nd = NetworkData(datafile, node_a="interaction_participants", node_b=None, score=None,
                    target_id_type='Entrez',  identifiers='Uniprot',
                    header=1, net_name="consensusDB", test_mode=True, sep="\t")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)
    
if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/iRefIndex/9606.mitab.06-11-2021.txt.zip"
    nd = NetworkData(datafile, node_a="#uidA", node_b="uidB", score=None, species=["taxa", "taxb"], species_code="taxid:9606(Homo sapiens)",
                    target_id_type='Entrez',  identifiers=['Uniprot', 'Refseq'], prefixes=[("uniprotkb:"), ("refseq:")],
                    header=0, net_name="iref_test", test_mode=True, sep="\t")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)

if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/IntAct/intact_test.txt"
    nd = NetworkData(datafile, node_a=0, node_b=1, score=None,
                    target_id_type='Entrez',  identifiers='Uniprot', prefixes=[("uniprotkb:")],
                    header=None, net_name="intact_b_test", test_mode=True, sep="\t")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)

if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/InBio_Map_core_2016_09_12/core.psimitab"
    nd = NetworkData(datafile, node_a=0, node_b=1, score=14,
                    target_id_type='Entrez',  identifiers='Uniprot', prefixes=[("uniprotkb:")],
                    header=None, net_name="InBio_test", test_mode=True, sep="\t")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)


if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/HI-II-14.tsv"
    nd = NetworkData(datafile, node_a=0, node_b=1, 
                    target_id_type='Entrez',  identifiers='Ensembl', 
                    header=None, net_name="HI_test", test_mode=True, sep="\t")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)


if False:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/DIP_Hsapi20170205.txt"
    nd = NetworkData(datafile, node_a="ID interactor A", node_b="ID interactor B", 
                    target_id_type='Entrez',  identifiers=["DIP", "Uniprot"],
                    header=0, net_name="dip_test", test_mode=True, sep="\t", 
                    prefixes=[("DIP-", "|"), ("uniprotkb:", "|")], species_code="taxid:9606(Homo sapiens)",
                    species=["Taxid interactor A", "Taxid interactor B"])
    #nd.data = pd.concat([nd.data, pd.DataFrame({"Gene_A": 'A6NKP2', "Gene_B": "P48506", "Weight": 1}, index=[101])])
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)


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
                    target_id_type='Entrez',  identifiers='Entrez', species_code=9606, 
                    score="Score", header=0, net_name="BIOGRID_entrex_test", test_mode=True, sep="\t")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)


if True:
    datafile = "/cellar/users/snwright/Data/Network_Analysis/Network_Data_Raw/PathwayCommons/PathwayCommons.8.bind.BINARY_SIF.hgnc.txt.sif"
    nd = NetworkData(datafile, node_a=0, node_b=2, species=None, 
                    target_id_type='Entrez',  identifiers='Symbol', 
                    header=None, net_name="bind_test", test_mode=True, sep="\t")
    outpath = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
    nd.clean_data()
    nd.convert_nodes()
    final_nodes = nd.get_unique_nodes()
    nd.write_network_data(outpath)
    nd.write_stats(outpath)