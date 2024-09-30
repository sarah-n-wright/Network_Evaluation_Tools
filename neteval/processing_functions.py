import pandas as pd
import numpy as np
from neteval.gene_mapper import *
from neteval.Timer import Timer
import csv
import re
from itertools import combinations
import os

    
def extract_id_with_prefix(id_str, pref_sep):
    """Extracts the gene identifier from a string containing a prefix and separator
    
    Args:
        id_str (str): String containing the gene identifier
        pref_sep (str): Prefix and separator for the gene identifier
        
    Returns:
        str: Extracted gene identifier
    
    """
    add_back_prefixes = ["DIP-", "ESNG", "ESNP"] # prefixes that are actually part of gene identifiers
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
            if prefix == '':
                if separator and separator in id_str:
                    id_str = id_str.split(separator)[0]
                    return id_str
            elif prefix in id_str:  # otherwise try next prefix if available
                id_str = id_str.split(prefix)[1]
                if separator and separator in id_str:
                    id_str = id_str.split(separator)[0]
                id_str = prefix + id_str if prefix in add_back_prefixes else id_str
                return id_str
    # if no prefixes found return NA
    return pd.NA


def clean_score(score):
    """Extracts a numeric score from a string with a score
    
    Args:
        score (str): String containing a score
        
    Returns:
        float: Extracted score or NA if no score is found
    
    """
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
    except TypeError:
        return pd.NA


class NetworkData:
    """Class for processing and converting network data based on configuration parameters
    
    Args:
        datafile (str): File path to the data file
        node_a (str): Column name of the first node
        node_b (str): Column name of the second node
        identifiers (str or list): Identifier type(s) for the input data
        target_id_type (str): Identifier type to convert to
        net_name (str): Name to give the network
        score (str): Column name for column containing scores for the interactions. Default None.
        species (str): Column name(s) containing species information for a given interaction. Default None.
        species_code (str): Value with species columns indicating interactions to keep. Default None.
        sep (str): Separator for the input data file. Default '\t'.
        header (int): Row number of the header. Default 0
        test_mode (bool): True if in test mode (process first 10000 interactions only). Default False.
        prefixes (str): Prefix and separator for the gene identifiers. Default None.
        
    """
    def __init__(self, datafile, node_a, node_b, identifiers, target_id_type, net_name, score=None, species=None, species_code=None, sep="\t", header=0, test_mode=False, prefixes=None):
        test_size=10000
        self.T = Timer()
        self.T.start("Total")
        self.datafile = datafile
        assert os.path.exists(datafile), "File not found, check path"
        available_id_types = ["Symbol", "Entrez", "Uniprot", "Ensembl", "DIP", "Refseq", "EnsemblProtein"]
        assert species is None or ((species is not None) and (species_code is not None)), "If species column is given, a corresponding value representing the species must be given"
        assert target_id_type in available_id_types, "Target identifier must be in " + "; ".join(available_id_types)
        
        self.species = species
        self.two_species_cols = True if type(species) == list else False
            
        self.score = score
        if (self.score is not None) and (type(self.score) == str):
            if "[p]" in self.score:
                # check if scores are p-values (as the best scores are the smallest values in this case)
                self.score = self.score.split("[p]")[0]
                self.score = int(self.score) if self.score.isnumeric() else self.score
                self.score_is_p = True
            else:
                self.score_is_p = False
        else:
            self.score_is_p = False
            
        self.node_a = node_a
        self.node_b = node_b
        # if only one node column is given, the data is assumed to be a list of complexes
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
        
        self.import_data(datafile, header, sep, net_name, test_mode, test_size)
            
        self.T.end("Load data")
        # Check that the columns are in the data
        assert self.node_a in self.raw_data.columns,str(node_a) + " is not present as a column in the data"
        assert (self.node_b in self.raw_data.columns, str(node_b)) or self.node_b is None, + " is not present as a column in the data"
        
        if self.score is not None:
            assert self.score in self.raw_data.columns, str(score) + " is not present as a column in the data"
            self.score_subset = pd.DataFrame()

        print(self.raw_data.head())
        if self.species is not None:
            if self.two_species_cols:
                assert self.species[0] in self.raw_data.columns, str(self.species[0])+ " is not present as a column in the data"
                assert self.species[1] in self.raw_data.columns, str(self.species[1])+ " is not present as a column in the data"
                #print(self.species_code, type(self.species_code))
                #print(self.raw_data[self.species[0]].values)
                assert species_code in self.raw_data[self.species[0]].values, "The `species_code` "+ str(species_code) + " is not present in the " + str(self.species[0]) +  " column"
                assert species_code in self.raw_data[self.species[1]].values, "The `species_code` "+ str(species_code) + " is not present in the " + str(self.species[1]) +  " column"
                
            else:
                assert self.species in self.raw_data.columns, str(species)+ " is not present as a column in the data"
                #print(self.species_code, type(self.species_code))
                assert species_code in self.raw_data[self.species].values, "The `species_code` "+ str(species_code) + " is not present in the `species` column"

        # create dictinary for tracking the stats
        self.stats = {"edges":{"raw":len(self.raw_data)}, "nodes":{}}
        
    def import_data(self, datafile, header, sep, net_name, test_mode=True, test_size=10000):
        """Imports the data from a file and loads it to the class
        
        Args:
            datafile (str): File path to the data file
            header (int): Row number of the header
            sep (str): Separator for the input data file
            net_name (str): Name to give the network
            test_mode (bool): True if in test mode (process first 10000 interactions only)
            test_size (int): Number of interactions to process in test mode. Default 10000.
            
        Returns:
            None
        """
        chars = ['"', "@", ""]
        quoting = csv.QUOTE_ALL
        # Try to load the data with different quote characters
        for i, quotechar in enumerate(chars):
            quoting = csv.QUOTE_NONE if i == (len(chars) - 1) else csv.QUOTE_ALL
            index=None
            nrows = test_size if test_mode else None
            try:
                test_data = pd.read_csv(datafile, sep=sep, header=header, engine='python', nrows=10, quoting=quoting, quotechar=quotechar)
                if (header is not None) and any(["Unnamed" in col for col in test_data.columns]):
                    index=False
                try:
                    self.raw_data = pd.read_csv(datafile, sep=sep, header=header, engine='python', nrows=nrows, 
                                                quoting=quoting, quotechar=quotechar, usecols=self.cols, index_col=index)   
                    self.data = self.raw_data.copy()
                except pd.errors.ParserError as err:
                    if "NULL" in str(err):
                        self.raw_data = pd.read_csv(datafile, sep=sep, header=header, engine='c', nrows=nrows, 
                                                quoting=quoting, quotechar=quotechar, usecols=self.cols, index_col=index)   
                        self.data = self.raw_data.copy()
                    else:
                        raise err
                except ValueError as e:
                    print(test_data.cols)
                    raise e
                break
            except:
                print("WARNING: Data not able to be loaded with quotechar", quotechar)
                if i == len(chars) - 1:
                    raise pd.errors.ParserError("Data could not be loaded with available quote characters")
                
        # Intialize an index column:
        self.data[net_name+"_ID"] = self.data.index.values
        self.cols.append(net_name+"_ID")

    def clean_data(self):
        """Cleans the data by removing columns, filtering by species, binarizing complexes, and sorting node pairs
        
        Args:
            None
            
        Returns:
            None
        """
        self.T.start("Clean data")
        self.data = self.data.loc[:, self.cols]
        # binarize complexes
        if self.complexes:
            self.binarize_complexes()
        self.data = self.data.dropna(subset=[self.node_a, self.node_b])
        # sort node pairs for easier identification of duplicates
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
            
        # drop duplicates, in the case of discordant edges keep the highest scoring if there is a duplicate
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
                
        # extract gene identifiers from prefixes
        self.extract_from_prefixes()        
        self.stats["edges"]["Node prefix"] = len(self.data)
        self.data.drop_duplicates(inplace=True, subset=[self.node_a, self.node_b])
        self.stats["edges"]["de-duped1"] = len(self.data)
        
        # convert identifiers to strings
        self.data[self.node_a] = self.data[self.node_a].astype(str)
        self.data[self.node_b] = self.data[self.node_b].astype(str)
        self.T.end("Clean data")
        
    def binarize_complexes(self):
        """Binarizes the complex data by creating pairwise combinations of the complex members
        
        Args:
            None
            
        Returns:
            None
        """
        results = []
        if any(["_HUMAN" in node for node in self.data[self.node_a].values]):
            separators = "[:,\|\.]"
        else:
            separators = "[:,\|_]"
        # Iterate over the rows of the input data frame to create a mapping of complex string to pairwise ids
        for comp in self.data[self.node_a]:
            # Get the list of IDs and the values for the other columns
            complex_id_list = re.split(separators, comp)
            # complex_id_list = comp.split(separator)
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
        """Extracts gene identifiers from prefixes
        
        Args:
            None
            
        Returns:
            None
        """
        if (self.prefix is not None) and (not self.mixed_identifiers):
            new_nodes_a = self.data[self.node_a].apply(lambda x: extract_id_with_prefix(x, self.prefix))
            new_nodes_b = self.data[self.node_b].apply(lambda x: extract_id_with_prefix(x, self.prefix))
            self.data[self.node_a] = new_nodes_a
            self.data[self.node_b] = new_nodes_b
            self.data.dropna(subset=[self.node_a,self.node_b])  
            
    def sort_node_pairs(self):
        """Sorts the node pairs in the data
        
        Args:
            None
            
        Returns:
            None
        """
            
        self.T.start("Sort nodes")
        node_array = self.data.loc[:, (self.node_a, self.node_b)].to_numpy()
        try:
            node_array.sort()
        except TypeError as e:
            if "not supported between instances" in str(e):
                try:
                    self.data[self.node_a] = self.data[self.node_a].astype(str)
                    self.data[self.node_b] = self.data[self.node_b].astype(str)
                    node_array = self.data.loc[:, (self.node_a, self.node_b)].to_numpy()
                    node_array.sort()
                except Exception as e:
                    raise e
            else:
                print("TYPE ERROR")
                print(self.data.head())
                print(node_array)
                raise e
        except ValueError as e:
                print(self.data.head())
                print(node_array)
                # save to the directory of datafile
                np.savetxt(X=node_array, fname=os.path.dirname(self.datafile) + "/node_array.txt", fmt='%s')
                raise e
        sorted_data = pd.DataFrame(node_array)
        sorted_data.columns = [self.node_a, self.node_b]
        if self.data.shape[1] > 2:
            self.data = pd.concat([sorted_data.reset_index(drop=True), self.data.drop(columns=[self.node_a, self.node_b]).reset_index(drop=True)], axis=1)
        else:
            self.data = sorted_data
        self.T.end("Sort nodes")
    
    def get_unique_nodes(self, score_subset=False):
        """Gets the unique nodes in the data
        
        Args:
            score_subset (bool): True if the unique nodes are for the score subset. Default False.
            
        Returns:
            set: Unique nodes in the data
            
        """
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
    
    def convert_nodes(self):
        """Wrapper for updating and converting node identifiers, and then converting all edges to new identifiers. When the input data contains multiple input
        types, they are converted in order. 
        
        Args:
            None
            
        Returns:
            None
        """
        self.T.start("Convert nodes")
        gene_map = {}
        all_nodes = self.get_unique_nodes()
        # create a node map
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

                    # catch ids not in all_nodes
                    gene_map = {**gene_map, **id_map}

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
        # use the node map to convert node identifiers
        self.convert_edges(node_map, unmapped_nodes)
        print("Converted")
        self.T.end("Convert nodes")
        return unmapped
    
    def get_node_conversion_map(self, nodes, initial_id, target_id):
        """Gets a mapping of node identifiers from the initial identifier to the target identifier. First updates the identifiers to the latest version, then converts the identifiers.

        Args:
            nodes (set): Set of node identifiers
            initial_id (str or list): Initial identifier type
            target_id (str): Target identifier type
            
        Returns:
            dict: Mapping of node identifiers from the initial identifier to the target identifier
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
        self.T.end("Get node conversion map")
        return node_map
        
    def convert_edges(self, node_map, unmapped_nodes):
        """Converts the node identifiers in the data
        
        Args:
            node_map (dict): Mapping of node identifiers from the initial identifier to the target identifier
            unmapped_nodes (set): Set of unmapped node identifiers
            
        Returns:
            None
        """
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
        """Subsets the data based on a score column and a percentile
        
        Args:
            score_col (str): Column name for the score
            percentile (int): Percentile to use for the subset
            
        Returns:
            None
            
        """
        if self.score is not None:
            # assign minimum score to NA values 
            replace_na = min(0, self.data[self.score].min())
            self.data.loc[self.data[self.score].isna(), self.score] = replace_na
            cutoff = np.percentile(self.data[self.score].values, percentile)
            self.score_subset = self.data.loc[self.data[self.score] > cutoff]
            self.stats['edges']['score_subset'] = len(self.score_subset)
            self.stats['nodes']['score_subset'] = len(self.get_unique_nodes(score_subset=True))
    
    def remove_duplicates(self):
        """Removes duplicate edges"""
        self.sort_node_pairs()
        self.data = self.data.drop_duplicates(subset=[self.node_a, self.node_b])
        self.stats["edges"]["de-duped_final"] = len(self.data)
        
    def remove_self_edges(self):
        """Removes self edges"""
        self.data = self.data.loc[self.data[self.node_a] != self.data[self.node_b]]
        self.stats["edges"]["removed_self_edges"] = len(self.data)
        
    def write_network_data(self, outpath, percentile=90):
        """Writes the processed network data to a file
        
        Args:
            outpath (str): Output path for the network data
            percentile (int): Percentile to use for the score subset. Default 90.
            
        Returns:
            None
        """
        # rename columns
        self.T.start("Write data")
        final_names= {self.node_a: self.target_id_type + "_A", self.node_b: self.target_id_type + "_B", self.net_name+"_ID": "ID"}
        if self.score is not None:
            final_names[self.score] = "Score"
        if self.score is not None:
            self.subset_on_score("Score", percentile)
            score_subset_nodes = pd.DataFrame({"Unique_Nodes": list(self.get_unique_nodes(score_subset=True))})
        pd.DataFrame({"Unique_Nodes": list(self.get_unique_nodes())}).to_csv(outpath + self.net_name + ".nodelist", sep="\t", index=False) 
        self.data.rename(columns = final_names, inplace=True)
        self.data.to_csv(outpath + self.net_name + "_net.txt", sep="\t", index=False)
        
        if self.score is not None:
            self.score_subset.rename(columns = final_names, inplace=True)
            score_subset_nodes.to_csv(outpath + self.net_name + "_"+ str(percentile) + ".nodelist", sep="\t", index=False)    
            self.score_subset.to_csv(outpath + self.net_name + "_"+ str(percentile) + "_net.txt", sep="\t", index=False)

        self.T.end("Write data")
        self.T.end("Total")
            
    def write_stats(self, outpath, unmapped=None):
        """Writes the processing stats to a file
        
        Args:
            outpath (str): Output path for the stats file
            
        Returns:
            None
        """
        stats_df = pd.DataFrame.from_dict(self.stats)
        stats_df.to_csv(outpath + self.net_name + ".stats", sep="\t", index=True)
        if unmapped is not None:
            with open(outpath + self.net_name + "_unmapped_nodes.txt", "w") as f:
                f.write("\n".join([str(x) for x in unmapped])+"\n")
        self.T.print_all_times()
    