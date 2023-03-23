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
    parser.add_argument('--name', type=str, required=True)
    parser.add_argument('--nodepref', default='Entrez_')
    args = parser.parse_args(args)
    return args
    

def get_node_edge_files(datadir, prefix_file):
    with open(prefix_file, 'r') as f:
        network_list = f.readlines()
    node_files = {net.split("\n")[0]: datadir + net.split("\n")[0]+".nodelist" for net in network_list}
    edge_files = {net.split("\n")[0]: datadir + net.split("\n")[0]+"_net.txt" for net in network_list}
    return node_files, edge_files


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
    node_f, edge_f = get_node_edge_files(args.d, args.prefix_file)
    n_dict, e_dict, data_dict = get_unique_nodes_edges(node_f, edge_f, node_pref=args.nodepref)
    for n in args.n:
        create_network_subset(e_dict, data_dict, args.o, min_dbs=n, node_pref=args.nodepref, name=args.name)
    
    