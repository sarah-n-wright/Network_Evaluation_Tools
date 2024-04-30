import argparse
import networkx as nx
from neteval.Timer import Timer
import os, sys
import warnings
import neteval.data_import_export_tools as dit

def shuffle_network(G, n_swaps, timer=None):
    if timer is not None:
        timer.start("Shuffle Network")
    edge_len=len(G.edges())
    G_shuff = G.copy()
    try:
        nx.double_edge_swap(G_shuff, nswap=edge_len*n_swaps, max_tries=edge_len*n_swaps*10, seed=None)
    except nx.NetworkXAlgorithmError:
        warning_string = 'Maximum number of swap attempts ('+str(edge_len*10)+') exceeded before desired swaps achieved ('+str(edge_len*n_swaps)+') for file' + G.graph['file'] +'.'
        warnings.warn(warning_string)
    except nx.NetworkXError:
        print("NETWORK ERROR:", G.graph['file'])
    shared_edges = len(set(G.edges()).intersection(set(G_shuff.edges())))
    try:
        print('Edge Similarity:', shared_edges/float(edge_len), G.graph['file'])
    except KeyError:
        print('Edge Similarity:', shared_edges/float(edge_len))
    if timer is not None:
        timer.end("Shuffle Network")
    return G_shuff

def write_shuffled_network(G, datafile, outpath, timer=None):
    outfile = outpath + os.path.split(datafile)[1].split(".txt")[0] + "_shuffled.txt"
    dit.write_networkx_to_file(G, outfilepath=outfile, timer=timer)
    
def load_shuffled_network(datafile, outpath, timer=None):
    outfile = outpath + os.path.split(datafile)[1].split(".txt")[0] + "_shuffled.txt"
    G = dit.load_edgelist_to_networkx(outfile, timer=timer)
    return G



def parse_arguments(args):
    parser = argparse.ArgumentParser(description='Shuffle the edges of a network')
    parser.add_argument('datafile', metavar='d', help='String with the full file path storing the network data')
    parser.add_argument('-o', metavar='outpath', required=True)
    parser.add_argument('--nSwaps', default='1', type=float)
    parser.add_argument('--testMode', default='1', type=int)
    parser.add_argument('--verbose', default='0', type=int)
    args = parser.parse_args(args)
    args.testMode = bool(args.testMode)
    args.verbose = bool(args.verbose)
    return args


if __name__=="__main__":
    args = parse_arguments(sys.argv[1:])
    
    T = Timer()
    T.start("Total")    
    #print(args.A, args.B)
    if args.verbose:
        print(args)
        print("Analysis of", args.datafile)  
    G = dit.load_edgelist_to_networkx(args.datafile, testmode=args.testMode, timer=T)
    if args.verbose:
        print("Data Loaded")
    if len(G.edges) > 0:
        G_shuff = shuffle_network(G, args.nSwaps, timer=T)
        if args.verbose:
            print("Network Shuffled")
    # get output file name from the datafile, then append to outpath
        write_shuffled_network(G_shuff, args.datafile, args.o, timer=T)
    else:
        print("NO EDGES:", args.datafile)
    T.end("Total")
    if args.verbose:
        T.print_all_times()
        print("Complete.")    
    
