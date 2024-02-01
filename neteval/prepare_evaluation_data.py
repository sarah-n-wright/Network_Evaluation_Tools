import neteval.data_import_export_tools as dit
from neteval.gene_mapper import update_nodes, convert_node_ids
import pandas as pd
import csv
from tqdm import tqdm
import requests
from getpass import getpass
import argparse
import re


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
        datadir='/cellar/users/snwright/Data'
        genesets = filter_gene_sets(genesets, args.d, network_file_list=args.n, min_genes=args.m, max_genes=args.M)
    write_gene_sets(genesets, args.o)
        
# def parse_geneset_arguments(args):
#     if args.s == "disgen":
#         disgen_args = args.A.split(";")
#         disgen_args = {arg.split("=")[0]:arg.split("=")[1] for arg in disgen_args}
#         disgen_args['min'] = int(disgen_args['min'])
#         disgen_args['max'] = int(disgen_args['max'])
#         disgen_args['types'] = disgen_args['types'].split(",")
#         args.A = disgen_args
#     if args.s != "disgen":
#         raise NotImplementedError("Only disgenet sets are currently supported")
#     args.m = int(args.m)
#     return args

# if __name__=="__main__":
#     parser = argparse.ArgumentParser(description='Process Evaluation Data')
#     parser.add_argument('-e', metavar='email', required=False, default=None)
#     parser.add_argument('-p', metavar="password", required=False, default=None)
#     parser.add_argument('-s', metavar="set_type", required=True)
#     parser.add_argument('-A', metavar='disgen_args', required=False, default="source=BEFREE;min=10;max=500;types=disease;file=/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/updated_gda_2023_BEFREE.tsv")
#     parser.add_argument('-o', metavar='outpath', required=True)
#     parser.add_argument('-n', metavar='network_list', required=False, default=None)
#     parser.add_argument('-d', metavar='network_dir', required=False, default=None)
#     parser.add_argument('-m', metavar='net_min', required=False, default=20)
#     parser.add_argument('-u', metavar='update', default=False)
#     parser.add_argument('-S', metavar='set_file', required=False, default=None)

#     args = parser.parse_args()
#     args = parse_geneset_arguments(args)   
    
#     # TODO this file needs a lot of cleanup and testing, particularly to add capability for other genesets
    
#     # Update the genesets
#     if args.u:
#         if args.set_type == "disgen":
#             get_disgenet_sets(args.disgen_args['file'], args.email, args.password, 
#                     outfile=args.outpath + 'Disgenet_' + args.disgen_args['source'] + '_gda.tsv', source=args.disgen_args['source'],
#                     types=args.disgen_args['types'], min_genes=args.disgen_args['min'] , max_genes=args.disgen_args['max'])
    
#         # Process the genesets
#             process_disgenet_data(args.outpath + 'Disgenet_' + args.disgen_args['source'] + '_gda.tsv', 
#                                 args.outpath + 'Disgenet_' + args.disgen_args['source'] + '_genesets.tsv', id_type='gene_symbol',
#                                 min_genes=args.disgen_args['min'] )
            
#     if args.n is not None:
#         gene_sets = dit.load_node_sets(args.S, id_type="Entrez")
#         keep_gene_sets = filter_gene_sets(gene_sets, network_file_list= args.n, datadir= args.d, min_genes= args.m)
#         write_gene_sets(keep_gene_sets, args.o + 'Disgenet_' + args.A['source'] + '_filtered_genesets.tsv')
    
    #password = getpass("User Password:")
    #process_disgenet_data("/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/updated_gda_2023_BEFREE.tsv",
    #                        outfile='/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/DisGeNET_genesets_BEFREE_2023.txt')
    #process_disgenet_data("/cellar/users/snwright/Data/Transfer/Disgenet_GWASDB_gda.tsv",
    #                        outfile='/cellar/users/snwright/Data/Transfer/Disgenet_GWASDB_genesets.tsv', id_type="gene_symbol")
    #get_disgenet_sets('/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/HPO_diseases.tsv', "snwright@ucsd.edu", "Ballon44!", 
    #                outfile='/cellar/users/snwright/Data/Transfer/Disgenet_HPO_gda.tsv', source="HPO",
    #                types=["disease", "phenotype"], min_genes=10, max_genes=300)
    #process_disgenet_data("/cellar/users/snwright/Data/Transfer/Disgenet_HPO_gda.tsv",
    #                        outfile='/cellar/users/snwright/Data/Transfer/Disgenet_HPO_genesets.tsv', id_type="gene_symbol", min_genes=10)
    #get_intersecting_genes('Git/Network_Evaluation_Tools/Data/v1_net_prefixes.txt', 
                            #datadir='/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/')
    # process_gwas_catalog_data('/cellar/users/snwright/Data/Network_Analysis/Reference_Data/gwas_cat_Jan13_2023.txt', "", pval_th=5e-8, include_intergenic=False, min_genes=20, max_genes=300)
    
    # datadir='/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/'
    # gitdir='/cellar/users/snwright/Git/Network_Evaluation_Tools/'
    # genesets = dit.load_node_sets(gitdir+ 'Data/DisGeNET_genesets_BEFREE_2023.txt')
    # filtered_genesets = filter_gene_sets(genesets, gitdir+'Data/v1_net_prefixes.txt', min_genes=20,datadir=datadir)
    # write_gene_sets(filtered_genesets, gitdir+'Data/DisGeNET_genesets_BEFREE_2023_v1_filtered_genesets_20.txt')
    
    # genesets = dit.load_node_sets(gitdir+ 'Data/DisGeNET_genesets_BEFREE_2023.txt')
    # filtered_genesets = filter_gene_sets(genesets, gitdir+'Data/v2_net_prefixes.txt', min_genes=20,datadir=datadir)
    # write_gene_sets(filtered_genesets, gitdir+'Data/DisGeNET_genesets_BEFREE_2023_v2_filtered_genesets_20.txt')
    
    # genesets = dit.load_node_sets(gitdir+ 'Data/DisGeNET_genesets_ALL_2023.txt')
    # filtered_genesets = filter_gene_sets(genesets, gitdir+'Data/v1_net_prefixes.txt', min_genes=20,datadir=datadir)
    # write_gene_sets(filtered_genesets, gitdir+'Data/DisGeNET_genesets_ALL_2023_v1_filtered_genesets_20.txt')
    
    # genesets = dit.load_node_sets(gitdir+ 'Data/DisGeNET_genesets_ALL_2023.txt')
    # filtered_genesets = filter_gene_sets(genesets, gitdir+'Data/v2_net_prefixes.txt', min_genes=20,datadir=datadir)
    # write_gene_sets(filtered_genesets, gitdir+'Data/DisGeNET_genesets_ALL_2023_v2_filtered_genesets_20.txt')
    # new_genesets = convert_genesets(genesets, initial_id='Symbol', target_id='Entrez')
    # write_gene_sets(new_genesets, 'Git/Network_Evaluation_Tools/Data/Oncogeneic_Components_genesets_entrez.txt')
