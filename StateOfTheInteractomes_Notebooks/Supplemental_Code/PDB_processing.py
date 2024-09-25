import pandas as pd
import os
from tqdm import tqdm
import numpy as np
import glob
import itertools
import sys
from neteval.gene_mapper import *

def get_pdb_files(datadir, download_id):
    return glob.glob(os.path.join(datadir, f'*{download_id}*'))


def process_pdb_data(file_list):
    pdb_files = []
    for f in file_list:
        pdb = pd.read_csv(f, header=1)
        pdb_files.append(pdb)
    pdb = pd.concat(pdb_files)
    pdb['Entry ID'] = pdb['Entry ID'].fillna(method='ffill')
    pdb_filt = pdb[pdb['Source Organism'] == 'Homo sapiens']
    if 'Entity Polymer Type' in pdb_filt.columns:
        pdb_filt = pdb_filt[pdb_filt['Entity Polymer Type'] == 'Protein']
    return pdb_filt

def extract_pdb_complexes(pdb_df):
    entries = pdb_df['Entry ID'].unique()
    complexes = {}
    for ent in tqdm(entries):
        proteins = pdb_df[pdb_df['Entry ID'] == ent]['Accession Code(s)'].dropna().unique()
        pairs = list(itertools.combinations(proteins, 2))
        pairs = [sorted(p) for p in pairs]
        protA = [x[0] for x in pairs]
        protB = [x[1] for x in pairs]
        if len(proteins) > 1:
            complexes[ent] = {'ProtA':protA, 'ProtB':protB}
    return make_edgelist(complexes)
            
def make_edgelist(data):
    # Initialize an empty list to store edges
    edges = []

    # Iterate over the dictionary to create edges
    for key, value in data.items():
        protA_list = value['ProtA']
        protB_list = value['ProtB']
        for protA, protB in zip(protA_list, protB_list):
            edges.append({'Source': protA, 'Target': protB, 'Key': key})

    # Convert the list of edges to a DataFrame
    edgelist_df = pd.DataFrame(edges)

    # Display the DataFrame
    return edgelist_df

def convert_uniprot_ids(pdb_net):
    print(pdb_net.head())
    unique_prots = set(pdb_net.Source).union(set(pdb_net.Target))
    updated, failed = update_nodes(unique_prots, 'Uniprot', keep="present", timer=None)
    converted, failed = convert_node_ids(list(updated.values()), 'Uniprot', 'Entrez')
    node_map = {x:converted[updated[x]] for x in updated.keys() if (x in updated) and (updated[x] in converted)}
    pdb_net['GeneA'] = pdb_net.Source.apply(lambda x: int(node_map[x]) if x in node_map else np.nan)
    pdb_net['GeneB'] = pdb_net.Target.apply(lambda x: int(node_map[x]) if x in node_map else np.nan)
    pdb_net['Edge'] = pdb_net.apply(lambda x: sorted([x.GeneA, x.GeneB]), axis=1)
    pdb_net['GeneA'] = pdb_net.Edge.apply(lambda x: x[0])
    pdb_net['GeneB'] = pdb_net.Edge.apply(lambda x: x[1])
    pdb_net = pdb_net.dropna()
    return pdb_net.loc[:, ('GeneA', 'GeneB')].astype(int).drop_duplicates().reset_index(drop=True)

def extract_individual_proteins(pdb_df):
    prot_ids = pdb_df['Accession Code(s)'].dropna().unique()
    updated_ids, failed_uni = update_nodes(prot_ids, 'Uniprot', keep="present", timer=None)
    converted_ids, failed_ent = convert_node_ids(list(updated_ids.values()), 'Uniprot', 'Entrez')
    node_map_ids = {x:converted_ids[updated_ids[x]] for x in updated_ids.keys() if (x in updated_ids) and (updated_ids[x] in converted_ids)}
    pdb_prots = set(node_map_ids.values())
    return pdb_prots
    

if __name__ == '__main__':
    datadir = sys.argv[-3]
    download_id = sys.argv[-2]
    file_type = sys.argv[-1]
    pdb_files = get_pdb_files(datadir, download_id)
    print(pdb_files)
    pdb_df = process_pdb_data(pdb_files)
    if file_type == 'complex':
        complex_df = extract_pdb_complexes(pdb_df)
        pdb_edgelist = convert_uniprot_ids(complex_df)
        pdb_edgelist.to_csv(os.path.join(datadir, 'pdb_edgelist.tsv'), sep='\t', index=False)
    elif file_type == 'protein':
        pdb_prots = extract_individual_proteins(pdb_df)
        with open(os.path.join(datadir, 'pdb_individual_prot_list.txt'), 'w') as f:
            f.write('\n'.join(list(pdb_prots))+'\n')