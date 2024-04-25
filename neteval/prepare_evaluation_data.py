import neteval.data_import_export_tools as dit
from neteval.gene_mapper import update_nodes, convert_node_ids
import pandas as pd
import argparse
import os


def check_min_genes_per_network(gene_set, network_file_list, datadir, min_genes=20):
    """Check that a gene set has at least min_genes in each network in the network_file_list.

    Args:
        gene_set (set): set of gene ids
        network_file_list (str): path to file containing list of network prefixes
        datadir (str): path to directory containing network files
        min_genes (int): minimum number of genes required in each network
    
    Returns:
        bool: True if gene set has at least min_genes in each network, False otherwise

    """
    with open(network_file_list, 'r') as f:
        prefs = f.readlines()
    node_files = [datadir + p.strip()+".nodelist" for p in prefs]
    # iterate over networks
    for net_file in node_files:
        net_genes = set(pd.read_csv(net_file).Unique_Nodes.astype(str).values)
        intersect = net_genes.intersection(gene_set)
        if len(intersect) < min_genes:
            print("Failed on", net_file, len(intersect))
            return False  # not all networks have at least min genes in the gene set
    return True


def check_all_genesets_against_network(gene_sets, network_file, min_genes=20):
    """Check that all gene sets have at least min_genes in the network specified by network_file.
    
    Args:
        gene_sets (dict): dictionary of gene sets
        network_file (str): path to network file
        min_genes (int): minimum number of genes required in each network
    
    Returns:
        dict: dictionary of gene sets that have at least min_genes in the network specified by network_file
    """
    net_genes = pd.read_csv(network_file)
    net_genes['Unique_Nodes'] = net_genes['Unique_Nodes'].astype(int)
    remove_sets = []
    for set_id in gene_sets:
        result = check_single_geneset(gene_sets[set_id], net_genes, min_genes)
        if result is False:
            remove_sets.append(set_id)
    for set_id in remove_sets:
        gene_sets.pop(set_id)
    print("Gene sets removed by", network_file, ":", len(remove_sets))
    return gene_sets
    
def check_single_geneset(set_nodes, net_genes, min_genes=20):
    """ Check that a single gene set has at least min_genes in the network specified by net_genes.
    Args:
        set_nodes (set): set of gene ids
        net_genes (pd.DataFrame): dataframe containing network genes
        min_genes (int): minimum number of genes required in each network
    
    Returns:   
        bool: True if gene set has at least min_genes in the network specified by net_genes, False otherwise
    """
    intersection_df = net_genes[net_genes.Unique_Nodes.isin(set_nodes)]
    return len(intersection_df) >= min_genes
    
        
def filter_gene_sets(genesets, datadir, network_file_list=None, min_genes=1, max_genes=20000):
    """ Filter gene sets to remove those that do not have at least min_genes in each network specified by network_file_list.
    Args:
        genesets (dict): dictionary of gene sets
        datadir (str): path to directory containing network files
        network_file_list (str): path to file containing list of network prefixes
        min_genes (int): minimum number of genes required in each network
        max_genes (int): maximum number of genes allowed in each network
        
    Returns:
        dict: dictionary of gene sets that have at least min_genes in each network specified by network_file_list
    
    """
    print(network_file_list)
    if (network_file_list is not None) and (network_file_list != ""):
        genesets = filter_gene_sets_by_network(genesets, network_file_list, datadir, min_genes=min_genes)
    out_genesets = {}
    for set_id in genesets:
        if len(genesets[set_id]) >= min_genes and len(genesets[set_id]) <= max_genes:
            out_genesets[set_id] = genesets[set_id]
    return genesets

def filter_gene_sets_by_network(genesets, network_file_list, datadir, min_genes=20):
    """ Filter gene sets to remove those that do not have at least min_genes in each network specified by network_file_list.
    Args:
        genesets (dict): dictionary of gene sets
        network_file_list (str): path to file containing list of network prefixes
        datadir (str): path to directory containing network files
        min_genes (int): minimum number of genes required in each network
        
    Returns:
        dict: dictionary of gene sets that have at least min_genes in each network specified by network_file_list
    """
    with open(network_file_list, 'r') as f:
        prefs = f.readlines()
    print(prefs)
    print(datadir) 
    #print(genesets['Reperfusion Injury'])
    node_files = [datadir + p.strip()+".nodelist" for p in prefs]
    keep_gene_sets = {set_id:genesets[set_id] for set_id in genesets}
    for network in node_files:
        keep_gene_sets = check_all_genesets_against_network(keep_gene_sets, network, min_genes)
    print(len(keep_gene_sets), "/", len(genesets), "retained after filtering.")
    return keep_gene_sets


def convert_genesets(genesets, initial_id, target_id):
    """ Convert gene sets from one id type to another.
    Args:
        genesets (dict): dictionary of gene sets
        initial_id (str): initial id type
        target_id (str): target id type
    Returns:
        dict: dictionary of gene sets with ids converted from initial_id to target_id
    """
    all_nodes = set()
    for set_id in genesets:
        all_nodes = all_nodes.union(genesets[set_id])
    updated_nodes, unmapped1 = update_nodes(all_nodes, id_type=initial_id)
    converted_nodes, unmapped2 = convert_node_ids(updated_nodes.values(),initial_id=initial_id, target_id=target_id)
    node_map = {g: converted_nodes[updated_nodes[g]] for g in updated_nodes if updated_nodes[g] in converted_nodes}
    new_gene_sets = {}
    for set_id in genesets:
        new_gene_sets[set_id] = {node_map[g] for g in genesets[set_id] if g in node_map} 
    print("Unmapped genes:", len(unmapped1 + list(unmapped2)), "/", len(all_nodes))
    return new_gene_sets 
    
def write_gene_sets(genesets, outfile, sep="\t"):
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
        
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Process Evaluation Data')
    parser.add_argument('-s', metavar='set_source_file', required=True, type=str)
    parser.add_argument('-i', metavar='id_type', required=True, type=str)
    parser.add_argument('-n', metavar='network_list', required=False, default='', type=str)
    parser.add_argument('-o', metavar='output', required=True, type=str)
    parser.add_argument('-C', action='store_true')
    parser.add_argument('-F', action='store_true')
    parser.add_argument('-d', metavar='datadir', required=False, default='', type=str)
    parser.add_argument('-m', metavar='min_genes', required=False, default=10, type=int)
    parser.add_argument('-M', metavar='max_genes', required=False, default=500, type=int)
    
    args = parser.parse_args()
    assert (args.C or args.F), "Must specify at least one of -C or -F"
    
    genesets = dit.load_node_sets(args.s, id_type=args.i)
    if args.C:
        genesets = convert_genesets(genesets, initial_id=args.i, target_id='Entrez')
    if args.i != 'Entrez':
        for set_id in genesets:
            genesets[set_id] = {int(node) for node in list(genesets[set_id])}
    if args.F:
        print(args.n)
        # get file path of current file
        os.path.dirname(os.path.abspath(__file__)) + '/../Data/'
        genesets = filter_gene_sets(genesets, args.d, network_file_list=args.n, min_genes=args.m, max_genes=args.M)
    write_gene_sets(genesets, args.o)
