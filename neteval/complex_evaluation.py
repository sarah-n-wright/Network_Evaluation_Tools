import itertools
from goatools.semantic import semantic_similarity
from collections import defaultdict
from neteval.gsea_functions import *
import obonet as obo
from goatools.obo_parser import GODag
from neteval.data_import_export_tools import *

# Load the GO DAG (Directed Acyclic Graph) for semantic similarity calculations
# go_dag = GODag("go-basic.obo")

import networkx as nx
import pandas as pd
from itertools import combinations
from sklearn.metrics import jaccard_score
import numpy as np
from tqdm import tqdm
import os
import argparse
# set working directory

def load_go_data(go_data_path=None):
    assert go_data_path is not None, 'Please provide a path to the GO data file (gene2go).'
    assert os.path.exists(go_data_path), f'GO data file not found at {go_data_path}.'
    data = pd.read_csv('go_data_path', sep='\t', names=['TaxID', 'GeneID', 'GO_term', 'Evidence', 'Qualifier', 'GO_name', 'Pubmed', 'Category'])
    geneid2go = data.groupby('GeneID')['GO_term'].apply(set).to_dict()
    godag = GODag('go.obo')
    branches = {'BP': data[data.Category=='Process'].GO_term.unique(), 'MF': data[data.Category=='Function'].GO_term.unique(), 'CC': data[data.Category=='Component'].GO_term.unique()}
    return geneid2go, godag, branches


def calculate_modularity(G, clusters):
    """
    Calculate the modularity of each cluster in the graph G.

    Parameters:
    G (networkx.Graph): The graph representing the protein interactions.
    clusters (dict): A dictionary with cluster IDs as keys and lists of proteins (nodes) as values.

    Returns:
    dict: A dictionary with cluster IDs as keys and modularity scores as values.
    """
    modularity_dict = {}
    for cluster_id, nodes in clusters.items():
        if len(nodes) < len(G.nodes):
            modularity = nx.algorithms.community.quality.modularity(G=G, communities=[set(nodes), set(G.nodes) - set(nodes)])
            modularity_dict[cluster_id] = modularity
    return modularity_dict

def calculate_clustering_coefficient(G, clusters):
    """
    Calculate the average clustering coefficient of each cluster in the graph G.

    Parameters:
    G (networkx.Graph): The graph representing the protein interactions.
    clusters (dict): A dictionary with cluster IDs as keys and lists of proteins (nodes) as values.

    Returns:
    dict: A dictionary with cluster IDs as keys and clustering coefficient scores as values.
    """
    clustering_coefficient_dict = {}
    for cluster_id, nodes in clusters.items():
        if len(nodes) < len(G.nodes):
            subgraph = G.subgraph(nodes)
            clustering_coefficient = nx.average_clustering(subgraph)
            clustering_coefficient_dict[cluster_id] = clustering_coefficient
    return clustering_coefficient_dict


def jaccard_index(set1, set2):
    """
    Calculate the Jaccard Index between two sets.

    Parameters:
    set1 (set): First set of elements (proteins).
    set2 (set): Second set of elements (proteins).

    Returns:
    float: Jaccard index between the two sets.
    """
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union

def calculate_precision_recall(clusters, benchmark_complexes):
    """
    Calculate the Jaccard index, precision, and recall for each cluster against benchmark complexes.

    Parameters:
    clusters (dict): A dictionary with cluster IDs as keys and lists of proteins (nodes) as values.
    benchmark_complexes (dict): A dictionary with complex IDs as keys and lists of proteins (nodes) as values.

    Returns:
    dict: A dictionary with cluster IDs as keys and a tuple (jaccard_index, precision, recall) as values.
    """
    evaluation_dict = {}
    for cluster_id, cluster_nodes in clusters.items():
        best_jaccard = 0
        best_precision = 0
        best_recall = 0
        cluster_set = set(cluster_nodes)
        for complex_id, complex_nodes in benchmark_complexes.items():
            complex_set = set(complex_nodes)
            j_index = jaccard_index(cluster_set, complex_set)
            precision = len(cluster_set.intersection(complex_set)) / len(cluster_set)
            recall = len(cluster_set.intersection(complex_set)) / len(complex_set)
            if j_index > best_jaccard:
                best_jaccard = j_index
                best_precision = precision
                best_recall = recall
        evaluation_dict[cluster_id] = (best_jaccard, best_precision, best_recall)
    return evaluation_dict


def calculate_go_similarity(go_terms1, go_terms2, go_dag, already_calculated={}):
    """
    Calculate the semantic similarity between two sets of GO terms using Wang et al.'s method.

    Parameters:
    go_terms1 (list): List of GO terms for the first protein.
    go_terms2 (list): List of GO terms for the second protein.
    go_dag (GODag): The GO DAG loaded from an OBO file.
    method (str): The method used for semantic similarity ('Wang' is the default).

    Returns:
    float: The average semantic similarity score between the two sets of GO terms.
    """
    scores = []
    for term1 in go_terms1:
        for term2 in go_terms2:
            if (term1, term2) in already_calculated:
                score = already_calculated[(term1, term2)]
            elif (term2, term1) in already_calculated:
                score = already_calculated[(term2, term1)]
            else:
                score = semantic_similarity(term1, term2, go_dag)
                scores.append(score)
                already_calculated[(term1, term2)] = score
    #scores = [score for score in scores if score is not None]
    return sum(scores) / len(scores) if scores else 0.0, already_calculated

def calculate_homogeneity(clusters, go_annotations, go_dag, go_branches):
    """
    Calculate the homogeneity of each cluster based on GO semantic similarity.

    Parameters:
    clusters (dict): A dictionary with cluster IDs as keys and lists of proteins (nodes) as values.
    go_annotations (pd.DataFrame): DataFrame with columns ['Protein', 'GO_term', 'Category']

    Returns:
    dict: A dictionary with cluster IDs as keys and homogeneity scores as values.
    """
    homogeneity_dict = {}
    already_calculated_go = defaultdict(float)
    already_calculated_genes = defaultdict(dict)
    
    for cluster_id, genes in tqdm(clusters.items()):
        print(len(already_calculated_go), len(already_calculated_genes))
        genes = [g for g in genes if g in go_annotations]
        gene_pairs = list(itertools.combinations(genes, 2))
        
        mf_scores = []
        bp_scores = []
        cc_scores = []
        
        for gene1, gene2 in gene_pairs:
            if (gene1, gene2) in already_calculated_genes:
                mf_scores.append(already_calculated_genes[(gene1, gene2)]['mf'])
                bp_scores.append(already_calculated_genes[(gene1, gene2)]['bp'])
                cc_scores.append(already_calculated_genes[(gene1, gene2)]['cc'])
                continue
            elif (gene2, gene1) in already_calculated_genes:
                mf_scores.append(already_calculated_genes[(gene2, gene1)]['mf'])
                bp_scores.append(already_calculated_genes[(gene2, gene1)]['bp'])
                cc_scores.append(already_calculated_genes[(gene2, gene1)]['cc'])
                continue
            
            go_terms1 = go_annotations[gene1]
            go_terms2 = go_annotations[gene2]
            mf_terms1 = [t for t in go_terms1 if t in go_dag and (go_dag[t].namespace == 'molecular_function')]
            mf_terms2 = [t for t in go_terms2 if t in go_dag and (go_dag[t].namespace == 'molecular_function')]
            bp_terms1 = [t for t in go_terms1 if t in go_dag and (go_dag[t].namespace == 'biological_process')]
            bp_terms2 = [t for t in go_terms2 if t in go_dag and (go_dag[t].namespace == 'biological_process')]
            cc_terms1 = [t for t in go_terms1 if t in go_dag and (go_dag[t].namespace == 'cellular_component')]
            cc_terms2 = [t for t in go_terms2 if t in go_dag and (go_dag[t].namespace == 'cellular_component')]
            already_calculated_genes[(gene1, gene2)] = {'mf':0, 'bp':0, 'cc':0}
            if mf_terms1 and mf_terms2:
                score, already_calculated_go = calculate_go_similarity(mf_terms1, mf_terms2, go_dag, already_calculated_go)
                mf_scores.append(score)
                already_calculated_genes[(gene1, gene2)]['mf'] = score
            if bp_terms1 and bp_terms2:
                score, already_calculated_go = calculate_go_similarity(bp_terms1, bp_terms2, go_dag, already_calculated_go)
                bp_scores.append(score)
                already_calculated_genes[(gene1, gene2)]['bp'] = score
            if cc_terms1 and cc_terms2:
                score, already_calculated_go = calculate_go_similarity(cc_terms1, cc_terms2, go_dag, already_calculated_go)
                cc_scores.append(score)
                already_calculated_genes[(gene1, gene2)]['cc'] = score

        
        # Average the scores across all protein pairs in the cluster
        mf_score_avg = sum(mf_scores) / len(mf_scores) if mf_scores else 0.0
        bp_score_avg = sum(bp_scores) / len(bp_scores) if bp_scores else 0.0
        cc_score_avg = sum(cc_scores) / len(cc_scores) if cc_scores else 0.0
        
        # Calculate the final GO score for the cluster
        go_score = (mf_score_avg + bp_score_avg + cc_score_avg) / 3
        homogeneity_dict[cluster_id] = go_score
    
    return homogeneity_dict

def load_clusters(datadir, network_pref, max_res):
    # Load the clusters
    clusters = {}
    with open(os.path.join(datadir, 'Reviewer_Updates', 'complex_pred', f'{network_pref}.{max_res}.nodes')) as f:
        for line in f:
            clusterid, size, nodes, persistence = line.strip().split('\t')
            nodes = [int(x) for x in nodes.split(' ') if x.isnumeric()]
            clusters[clusterid] = nodes
    return clusters

if __name__=='__main__':
    #Use add parser to read int he network file path and max res
    parser = argparse.ArgumentParser(description='Evaluate a network against a benchmark.')
    parser.add_argument('--datadir', type=str, help='Path to the parent directory containing the network and cluster files')
    parser.add_argument('--network_pref', type=str, help='Network prefix.')
    parser.add_argument('--outpref', type=str, help='Path to write outputs')
    parser.add_argument('--max_res', type=int, help='Maximum number of clusters to extract.')
    parser.add_argument('--max_cluster_size', type=int, required=False, default = 30000, help='Maximum size of clusters to extract.')
    parser.add_argument('--go_data', type=str, required=False, help='Path to the GO data file (gene2go).')
    args = parser.parse_args()
    clusters = load_clusters(args.datadir, args.network_pref, args.max_res)
    G = load_edgelist_to_networkx(os.path.join(args.datadir, 'Processed_Data', 'v2_fixed', args.network_pref+'_net.txt'), id_type="Entrez", testmode=False, timer=None, delimiter="\t", node_cols=[0,1],
                            keep_attributes=False, verbose=False)
    
    mods = calculate_modularity(G, clusters)
    coeffs = calculate_clustering_coefficient(G, clusters)
    # load the clusters   
    gd, godag, go_branches = load_go_data(args.go_data)
    sub_clusters = {k: v for k, v in clusters.items() if len(v) < args.max_cluster_size}
    homogs = calculate_homogeneity(sub_clusters, gd, godag, go_branches)
    
    output_df = pd.DataFrame({'Modularity': mods, 'ClusteringCoeff': coeffs, 'Homogeneity': homogs})
    output_df['Cluster_size'] = output_df.index.map(lambda x: len(clusters[x]))
    
    output_df.to_csv(args.outpref +'.evaluation.tsv', sep='\t')

