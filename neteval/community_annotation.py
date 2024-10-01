import pandas as pd
import itertools
import numpy as np
from goatools.obo_parser import GODag
from goatools.semantic import semantic_similarity
from collections import defaultdict
import networkx as nx


def load_corum_data(corum_file, outdir):
    corum_df = pd.read_csv(corum_file, sep='\t', usecols=['ComplexID', 'ComplexName', 'subunits(Entrez IDs)', 'GO ID', 'GO description'])
    corum_df['NumSubunits'] = corum_df['subunits(Entrez IDs)'].apply(lambda x: len(x.split(';')))
    complex_info = corum_df.loc[:,['ComplexID', 'ComplexName', 'GO ID', 'GO description', 'NumSubunits']]
    complex_info.to_csv(os.path.join(outdir,'corum_complex_info.tsv'), sep='\t')
    all_complexes = []
    for index, row in corum_df.iterrows():
        subunits = row['subunits(Entrez IDs)'].split(';')
        subunits = [int(x) for x in subunits if x.isnumeric()]
        all_complexes.append(pd.DataFrame({'Complex': row['ComplexID'], 'Subunit': subunits}))
    complex_df = pd.concat(all_complexes).reset_index(drop=True)
    complex_df.Subunit = complex_df.Subunit.astype(int)
    complex_df['value'] = 1
    complex_mat = complex_df.pivot_table(index='Subunit', columns='Complex', values='value')
    complex_mat = complex_mat.fillna(0)
    complex_mat.index.name=None
    complex_mat.columns.name=None
    complex_mat.to_csv(os.path.join(outdir, 'corum_processed.tsv'), sep='\t')
    return complex_mat


def process_network_clusters(prefix, outdir):
    test_df = pd.read_csv(os.path.join(outdir, f'{prefix}.nodes'), usecols=[0,2], sep='\t', header=None, names=['ID', 'Genes'])
    all_clusters = []
    for index, row in test_df.iterrows():
        subunits = row['Genes'].split(' ')
        subunits = [int(x) for x in subunits if x.isnumeric()]
        all_clusters.append(pd.DataFrame({'Complex': row['ID'], 'Subunit': subunits}))
    test_df = pd.concat(all_clusters).reset_index(drop=True)
    test_df.Subunit = test_df.Subunit.astype(int)
    test_df['value'] = 1
    test_mat = test_df.pivot_table(index='Subunit', columns='Complex', values='value')
    test_mat = test_mat.fillna(0)
    test_mat.index.name=None
    test_mat.columns.name=None
    return test_mat
    
    
def standardize_mats(test_mat, complex_mat):
    not_in_test = [x for x in complex_mat.index if x not in test_mat.index]
    not_in_corum = [x for x in test_mat.index if x not in complex_mat.index]
    complex_filled = pd.concat([complex_mat, pd.DataFrame({comp: 0 for comp in complex_mat.columns}, index=not_in_corum)])
    test_filled = pd.concat([test_mat, pd.DataFrame({comp: 0 for comp in test_mat.columns}, index=not_in_test)])
    complex_filled.sort_index(inplace=True)
    test_filled.sort_index(inplace=True)
    return complex_filled, test_filled


def calculate_Jaccard(test_filled, complex_filled):
    intersection = np.dot(test_filled.T, complex_filled)
    sum1 = np.repeat([np.sum(test_filled, axis=0).values], [complex_filled.shape[1]], axis=0).T
    sum2 = np.repeat([np.sum(complex_filled, axis=0).values], [test_filled.shape[1]], axis=0)
    union = sum1 + sum2 - intersection
    J = intersection / union
    return J


def get_best_matches(J, test_ids, complex_ids):
    best_match_to_test = {test_ids[i]: complex_ids[x]for i, x in enumerate(np.argmax(J, axis=1))}
    best_similarity_test = {test_ids[i]: x for i, x in enumerate(np.max(J, axis=1))}

    # Find the best match for each complex in dfB against dfA
    best_match_to_complex = {complex_ids[i]: test_ids[x] for i, x in enumerate(np.argmax(J, axis=0))}
    best_similarity_complex = {complex_ids[i]: x for i, x in enumerate(np.max(J, axis=0))}
    test_results = pd.DataFrame({'CORUM_match':best_match_to_test, 'J': best_similarity_test})
    complex_results = pd.DataFrame({'Test_match': best_match_to_complex, 'J': best_similarity_complex})
    return test_results, complex_results


def map_complex_ids(res_df, outdir ):
    complex_info = pd.read_csv(os.path.join(outdir,'corum_complex_info.tsv'), sep='\t', index_col=0)
    if 'CORUM_match' in res_df.columns:
        res_df['TestID'] = res_df.index
        res_df = res_df.merge(complex_info.loc[:, ('ComplexID', 'ComplexName', 'NumSubunits')], left_on='CORUM_match', right_on='ComplexID', how='left')
    else:
        res_df = res_df.merge(complex_info.loc[:, ('ComplexID', 'ComplexName', 'NumSubunits')], left_index=True, right_on='ComplexID', how='left')
    return res_df

def write_results(test_df, complex_df, prefix, outdir):
    test_df.to_csv(os.path.join(outdir, prefix+ '.test_to_corum.tsv'), sep='\t')
    complex_df.to_csv(os.path.join(outdir, prefix+ '.corum_to_test.tsv'), sep='\t')
    
    
def corum_analysis_wrapper(corum_file, prefix, outdir):
    # load the necessary data
    complex_mat = load_corum_data(corum_file, outdir=outdir)
    test_mat = process_network_clusters(prefix, outdir)
    # calculate all Jaccard indexes
    complex_filled, test_filled = standardize_mats(test_mat, complex_mat)
    J = calculate_Jaccard(test_filled, complex_filled)
    # identify and get information on best matches
    test_res, complex_res = get_best_matches(J, test_filled.columns, complex_filled.columns)
    complex_out = map_complex_ids(complex_res, outdir)
    test_out = map_complex_ids(test_res, outdir)
    #save the results
    write_results(test_out, complex_out, prefix, outdir)
    
    
def load_go_data():
    data = pd.read_csv('gene2go', sep='\t', names=['TaxID', 'GeneID', 'GO_term', 'Evidence', 'Qualifier', 'GO_name', 'Pubmed', 'Category'])
    geneid2go = data.groupby('GeneID')['GO_term'].apply(set).to_dict()
    godag = GODag('go.obo')
    branches = {'BP': data[data.Category=='Process'].GO_term.unique(), 'MF': data[data.Category=='Function'].GO_term.unique(), 'CC': data[data.Category=='Component'].GO_term.unique()}
    return geneid2go, godag, branches


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
        genes = [g for g in genes if g in go_annotations]
        gene_pairs = list(itertools.combinations(genes, 2))
        scores = {'mf' : [], 'bp': [],'cc': []}
        
        for gene1, gene2 in gene_pairs:
            if (gene1, gene2) in already_calculated_genes:
                for branch in ['mf', 'bp', 'cc']:
                    scores[branch].append(already_calculated_genes[(gene1, gene2)][branch])
                continue
            elif (gene2, gene1) in already_calculated_genes:
                for branch in ['mf', 'bp', 'cc']:
                    scores[branch].append(already_calculated_genes[(gene2, gene1)][branch])
                continue
            terms1 = {'go' :go_annotations[gene1]}
            terms2 = {'go' :go_annotations[gene2]}
            branch_names = {'mf':'molecular_function', 'bp':'biological_process', 'cc':'cellular_component'}
            for branch in ['mf', 'bp', 'cc']:
                terms1[branch] = [t for t in terms1['go'] if t in go_dag and (go_dag[t].namespace == branch_names[branch])]
                terms2[branch] = [t for t in terms2['go'] if t in go_dag and (go_dag[t].namespace == branch_names[branch])]
      
            already_calculated_genes[(gene1, gene2)] = {'mf':0, 'bp':0, 'cc':0}
            for branch in ['mf', 'bp', 'cc']:
                if terms1[branch] and terms2[branch]:
                    score, already_calculated_go = calculate_go_similarity(terms1[branch], terms2[branch], 
                                                                           go_dag, already_calculated_go)
                    scores[branch].append(score)
                    already_calculated_genes[(gene1, gene2)][branch] = score

        avg_scores = {}
        for branch in ['mf', 'bp', 'cc']:
            avg_scores[branch] = sum(scores[branch]) / len(scores[branch]) if scores[branch] else 0.0
        
        # Calculate the final GO score for the cluster
        go_score = (avg_scores['mf'] + avg_scores['bp'] + avg_scores['cc']) / 3
        homogeneity_dict[cluster_id] = go_score
    
    return homogeneity_dict


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