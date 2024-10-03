#!/usr/bin/env python3
import networkx as nx
import pandas as pd
from collections import defaultdict
import argparse
import sys

def parse_arguments(args):
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('prefix_file', metavar='p', help='')
    parser.add_argument('-d', metavar='datapath', required=True)
    parser.add_argument('-o', metavar='outpath', required=True)
    parser.add_argument('-n',  type=int, required=True, nargs='+')
    parser.add_argument('-m', metavar='method', default='parsimonious')
    parser.add_argument('--name', type=str, required=True)
    parser.add_argument('--nodepref', default='Entrez_')
    args = parser.parse_args(args)
    return args
    

def get_node_edge_files(datadir, prefix_file):
    with open(prefix_file, 'r') as f:
        network_list = f.readlines()
    node_files = {net.split("\n")[0]: datadir + net.split("\n")[0]+".nodelist" for net in network_list}
    edge_files = {net.split("\n")[0]: datadir + net.split("\n")[0]+"_net.txt" for net in network_list}
    nets = [net.split("\n")[0] for net in network_list]
    return node_files, edge_files, nets


def get_unique_nodes_edges(node_files, edge_files, node_pref):
    nodes = defaultdict(int)
    edges = defaultdict(int)
    edge_data = defaultdict(list)
    for net in node_files.keys():
        print(net)
        nodes, edges, edge_data = parse_file_network(node_files[net], edge_files[net], nodes, edges, edge_data, net, node_pref)
        print(net, "finished.")
    return nodes, edges, edge_data


def parse_file_network(node_file, edge_file, node_dict, edge_dict, edge_data_dict=None, net_name=None, node_pref='Entrez_'):
    n = pd.read_csv(node_file)
    g = pd.read_csv(edge_file, sep="\t")
    e = g.apply(lambda x: (x[node_pref +"A"], x[node_pref+"B"]), axis=1)
    node_dict = count_nodes(node_dict, n.Unique_Nodes.values)
    edge_dict, edge_data_dict = count_edges(edge_dict, e, edge_data_dict = edge_data_dict, net_name=net_name )
    return node_dict, edge_dict, edge_data_dict

def count_nodes(node_dict, nodes):
    for node in nodes:
        node_dict[node] += 1
    return node_dict

def count_edges(edge_dict, edges, edge_data_dict=None, net_name=None):
    for edge in edges:
        edge = list(edge)
        edge.sort()
        edge = tuple(edge)
        edge_dict[edge] += 1
        if edge_data_dict is not None:
            edge_data_dict[edge].append(net_name)
    return edge_dict, edge_data_dict
    

def create_network_subset(edges, edge_data, outpath ,min_dbs=2, node_pref="Entrez_", name=""):
    print("Creating network subset:", min_dbs)
    keep_edges = [e for e in edges if edges[e] >= min_dbs]
    G = nx.Graph()
    use_edges = [(e[0], e[1], {"num_dbs":edges[e], "dbs":",".join(edge_data[e])}) for e in keep_edges]
    G.add_edges_from(use_edges)
    G_df = nx.to_pandas_edgelist(G, source=node_pref +"A", target=node_pref +"B")
    if node_pref == "Entrez_":
        G_df["Entrez_A"] = G_df["Entrez_A"].astype(int)
        G_df["Entrez_B"] = G_df["Entrez_B"].astype(int)
    G_df.to_csv(outpath +name+"_composite_min"+str(min_dbs)+"_net.txt", sep="\t", index=False)
    node_df = pd.DataFrame({"Unique_Nodes": list(G.nodes())})
    node_df["Unique_Nodes"] = node_df["Unique_Nodes"].astype(int)
    node_df.to_csv(outpath + name+ "_composite_min"+str(min_dbs) + ".nodelist", sep="\t", index=False)  
    return G


if __name__=="__main__":
    args = parse_arguments(sys.argv[1:])
    # args = parse_arguments(['-d', '/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/',
    #                         '-o', '/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/',
    #                         '-n', '2', '--name', 'test', '--nodepref', 'Entrez_',
    #                         '/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/v2_net_prefixes.txt'
    #                         ])
    print(args)
    node_f, edge_f, nets = get_node_edge_files(args.d, args.prefix_file)
    print(node_f)
    assert args.m in ['parsimonious', 'ordered'], "Method must be one of 'parsimonious' or 'ordered'."
    if args.m == 'parsimonious':
        n_dict, e_dict, data_dict = get_unique_nodes_edges(node_f, edge_f, node_pref=args.nodepref)
        for n in args.n:
            create_network_subset(e_dict, data_dict, args.o, min_dbs=n, node_pref=args.nodepref, name=args.name)
    elif args.m == 'ordered':
        n_dict = defaultdict(int)
        e_dict = defaultdict(int)
        data_dict = defaultdict(list)
        n_dict, e_dict, data_dict = get_unique_nodes_edges({net:node_f[net] for net in nets[:(args.n[0]-1)]}, {net:edge_f[net] for net in nets[:(args.n[0]-1)]}, node_pref=args.nodepref)
        assert len(args.n) == 1, "Only one minimum number of databases can be specified for ordered method."
        for i in range(args.n[0]-1, len(node_f)):
            n_dict, e_dict, data_dict = parse_file_network(node_f[nets[i]], edge_f[nets[i]], n_dict, e_dict, data_dict, nets[i], args.nodepref)
            #top_nets = nets[:i]
            #print(i, top_nets)
            #n_dict, e_dict, data_dict = get_unique_nodes_edges({net:node_f[net] for net in top_nets}, {net:edge_f[net] for net in top_nets}, node_pref=args.nodepref)
            create_network_subset(e_dict, data_dict, args.o, min_dbs=args.n[0], node_pref=args.nodepref, name=args.name+str(i+1)+'_')
    
    
