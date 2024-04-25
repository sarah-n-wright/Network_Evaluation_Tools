###############################################
# ---------- Data Import Functions ---------- #
###############################################

import pandas as pd
import networkx as nx
import ndex2
import ndex2.client
from getpass import getpass
from neteval.gene_mapper import convert_node_ids

def load_public_network_from_ndex(uuid, remove_self_edges=False, verbose=True):
    interactome_uuid=uuid # for PCNet
    # interactome_uuid='275bd84e-3d18-11e8-a935-0ac135e8bacf' # for STRING high confidence
    ndex_server='public.ndexbio.org'
    ndex_user=None
    ndex_password=None
    G_int = ndex2.create_nice_cx_from_server(
            ndex_server, 
            username=ndex_user, 
            password=ndex_password, 
            uuid=interactome_uuid
        ).to_networkx()
    nodes = list(G_int.nodes)
    # pcnet appears to have some self edges... should remove them. 
    if remove_self_edges:
        G_int.remove_edges_from(nx.selfloop_edges(G_int))
    # print out interactome num nodes and edges for diagnostic purposes
    if verbose:
        print('number of nodes:')
        print(len(G_int.nodes))
        print('\nnumber of edges:')
        print(len(G_int.edges))
    return G_int


def load_private_network_from_ndex(uuid, username=None, password=None, remove_self_edges=False, verbose=True):
    if username is None:
        username = getpass("NDEx Username")
    if password is None:
        password = getpass("NDEx Password")
    interactome_uuid=uuid # for PCNet
    # interactome_uuid='275bd84e-3d18-11e8-a935-0ac135e8bacf' # for STRING high confidence
    ndex_server='public.ndexbio.org'
    G_int = ndex2.create_nice_cx_from_server(
            ndex_server, 
            username=username, 
            password=password, 
            uuid=interactome_uuid
        ).to_networkx()
    if remove_self_edges:
        G_int.remove_edges_from(nx.selfloop_edges(G_int))
    # print out interactome num nodes and edges for diagnostic purposes
    if verbose:
        print('number of nodes:')
        print(len(G_int.nodes))
        print('\nnumber of edges:')
        print(len(G_int.edges))
    return G_int
    

def create_networkx_for_export(edge_file, verbose=True, id_type='Entrez', add_symbols=False, attributes=None):
    if isinstance(edge_file, str):
        G_df = pd.read_csv(edge_file, sep='\t')
    elif isinstance(edge_file, pd.DataFrame):
        G_df= edge_file
    else:
        raise ValueError("edge_file must be a string (filename to load to dataframe) or a pandas DataFrame")
    if id_type == 'Entrez':
        G_df['Gene_A'] = G_df['Entrez_A'].astype(int)
        G_df['Gene_B'] = G_df['Entrez_B'].astype(int)
    else:
        G_df.columns[0:2] = ['Gene_A', 'Gene_B']
    if attributes is None:
        G_df = G_df.iloc[:, 0:2]
    else:
        G_df = G_df.loc[:, ['Gene_A', 'Gene_B', 'Entrez_A', 'Entrez_B'] + attributes]
    if add_symbols:
        node_list = list(set(G_df['Gene_A']).union(set(G_df['Gene_B'])))
        node_map, missing = convert_node_ids(node_list, 'Entrez', 'Symbol')
        node_map_int = {int(k):v for k,v in node_map.items()}
        G_df['name']=G_df.apply(lambda x: node_map_int[x['Gene_A']] + '_' + node_map_int[x['Gene_B']], axis=1)
        
    else:
        G_df['name'] = G_df['Gene_A'].astype(str) + '_' + G_df['Gene_B'].astype(str)
    ## check for duplicate names
    name_counts = pd.DataFrame.from_dict(node_map_int, orient='index').value_counts()
    dups = name_counts[name_counts>1]
    if len(dups) > 0:
        for dup in dups.index:
            dup_ids = sorted([k for k in node_map_int.keys() if node_map_int[k] == dup[0]])
            for i, idx in enumerate(dup_ids):
                if i > 0:
                    node_map_int[idx] = node_map_int[idx] + '_' + str(idx)
        
    G_df['GeneID_Edge'] = G_df['Gene_A'].astype(str) + '_' + G_df['Gene_B'].astype(str)
    # create networkx representation
    G_export = nx.from_pandas_edgelist(G_df, source='Gene_A', target='Gene_B', create_using=nx.Graph, edge_attr=attributes+['name', 'Entrez_A', 'Entrez_B'])
    
    if add_symbols:
        nx.set_node_attributes(G_export, node_map_int, 'name')
        #nx.set_node_attributes(G_export, {k:k for k in node_map_int.keys()}, 'GeneID')
        nx.set_node_attributes(G_export, {k:'ncbigene:'+str(k) for k in node_map_int.keys()}, 'represents')
    nx.set_edge_attributes(G_export, 'interacts-with', 'interaction')
        
    return G_export
    
def create_pcnet_networkx_for_export(edge_file, verbose=True, id_type='Entrez', db_name_dict=None):
    G_df = pd.read_csv(edge_file, sep='\t')
    if db_name_dict is not None:
        db_combos = G_df.dbs.unique()
        db_map = {}
        for db_group in db_combos:
            dbs = db_group.split(',')
            db_map[db_group] = ','.join([db_name_dict[db] for db in dbs])
        G_df['dbs'] = G_df['dbs'].map(db_map)
    G_df.rename(columns={'dbs': 'Supporting_Databases', 'num_dbs': 'Number_of_Supporting_Databases'}, inplace=True)
    G_df['Supporting_Databases'] = G_df['Supporting_Databases'].apply(lambda x: x.split(','))
    G_export = create_networkx_for_export(G_df, verbose=verbose, id_type=id_type, add_symbols=True, attributes=['Supporting_Databases', 'Number_of_Supporting_Databases'])
    return G_export


def export_networkx_to_ndex(G, G_name, username=None, password=None, contexts={'ncbigene':'https://identifiers.org/ncbigene/'}, properties=None):
    G.name = G_name
    if username is None:
        username = getpass("NDEx Username")
    if password is None:
        password = getpass("NDEx Password")
    my_ndex=ndex2.client.Ndex2("http://public.ndexbio.org", username, password)
    cx = ndex2.create_nice_cx_from_networkx(G)
    cx.set_context(contexts)
    #cx.set_network_attributes({'prov:generated_by': 'test network', 'Notes': 'test notes'})
    url = cx.upload_to(client=my_ndex)
    network_id = url.split('/')[-1]
    if properties is not None:
        my_ndex.set_network_properties(network_id, properties)
    print(url)
    return url
    

def load_edgelist_to_networkx(datafile, id_type="Entrez", testmode=False, timer=None, delimiter="\t", node_cols=[0,1],
                            keep_attributes=False, verbose=False):
    # Assumes node columns are columns 0,1
    # Do I really want it to return a graph when there is an error?
    valid_id_types = ["Entrez", "Symbol"]
    assert id_type in valid_id_types, "id_type must be one of:" + ", ".join(valid_id_types)
    if timer is not None:
        timer.start("Load network")
    if testmode:
        net_df = pd.read_csv(datafile, sep=delimiter, index_col=None, nrows=1000)
    else:    
        net_df = pd.read_csv(datafile, sep=delimiter, index_col=None)
    #has_edge_attributes = True if len(net_df.columns) > 2 else None
    source_col, target_col = [net_df.columns[i] for i in node_cols]
    if id_type == "Entrez":
        net_df[source_col] = net_df[source_col].astype(int)
        net_df[target_col] = net_df[target_col].astype(int)
    try:
        if keep_attributes:
            G = nx.from_pandas_edgelist(net_df, source=source_col, target=target_col, edge_attr=keep_attributes)
        else:
            G = nx.from_pandas_edgelist(net_df, source=source_col, target=target_col, edge_attr=None)
    except KeyError: 
        print("FILE LOAD ERROR:", datafile)
        G = nx.Graph()
    except pd.errors.EmptyDataError:
        print("EMPTY DATA ERROR:", datafile)
        G = nx.Graph()
    G.graph['file'] = datafile
    if verbose:
        print('Network File Loaded:', datafile)
        print("# Nodes:", len(G.nodes))
        print("# Edges:", len(G.edges))
    if timer is not None:
        timer.end("Load network")
    return G


def write_networkx_to_file(G, outfilepath, timer=None, verbose=False, source="Node_A", target="Node_B"):
    if timer is not None:
        timer.start("Write Network")
    net_df = nx.to_pandas_edgelist(G, source=source, target=source)
    net_df.to_csv(outfilepath, index=False, sep="\t")
    if timer is not None:
        timer.end("Write Network")
    return


# Construct dictionary of node sets from input text file to perform AUPRC analysis on for network of interest
# File format: Each line is a delimited list with the first item in the list is the name of the node set
# All other nodes in the list follow the node set name
def load_node_sets(node_set_file, delimiter='\t', verbose=False, id_type="Entrez"):
    """_summary_

    Args:
        node_set_file (_type_): _description_
        delimiter (str, optional): _description_. Defaults to '\t'.
        verbose (bool, optional): _description_. Defaults to False.
        id_type (str, optional): _description_. Defaults to "Symbol".

    Returns:
        _type_: _description_
    """
    f = open(node_set_file)
    node_set_lines = f.read().splitlines()
    node_set_lines_split = [line.split(delimiter) for line in node_set_lines]
    f.close()
    node_sets = {node_set[0]:set(node_set[1:]) for node_set in node_set_lines_split}
    if id_type == "Entrez":
        for set_id in node_sets:
            node_sets[set_id] = {int(node) for node in list(node_sets[set_id]) if node.isnumeric()}
    if verbose:
        print('Node cohorts loaded:', node_set_file)
    return node_sets


def write_node_sets(genesets, outfile, sep="\t"):
    """ Write gene sets to file.
    Args:
        genesets (dict): dictionary of gene sets
        outfile (str): path to output file
        sep (str): separator for output file
    
    Returns:
        None
    """
    out_strings = [sep.join([set_id] + [str(g) for g in genesets[set_id]]) for set_id in genesets]
    with open(outfile,'w') as tfile:
        tfile.write('\n'.join(out_strings))
