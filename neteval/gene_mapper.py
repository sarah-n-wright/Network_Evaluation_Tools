import pandas as pd
import numpy as np
import neteval.query_uniprot as uni
import neteval.query_hgnc as hgnc
import neteval.query_ensembl as ensg
import mygene
import csv
import re
from itertools import combinations
import os

from datetime import datetime


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

    
