###################################################################
# Command line script to analyze network on node sets of interest #
###################################################################

import neteval.network_evaluation_functions as nef
import neteval.data_import_export_tools as dit
import neteval.shuffle_networks as shuf
import neteval.network_propagation as prop
import argparse
import os
import pandas as pd
import numpy as np


def restricted_float(x):
    """Check that numerical value is a float in the range (0.0, 1.0) exclusive. Value can also be None"""
    if x is not None:
        x = float(x)
        if x <= 0.0 or x >= 1.0:
            raise argparse.ArgumentTypeError("%r not in range (0.0, 1.0) exclusive"%(x,))
    return x


def positive_int(x):
    """Check that numerical value is a positive integer. Value must be >0."""
    x = int(x)
    if x <= 0:
        raise argparse.ArgumentTypeError("%s must be a positive integer" % x)
    return x


def valid_infile(in_file):
    """Check that input file path exists and is readable"""
    if not os.path.isfile(in_file):
        raise argparse.ArgumentTypeError("{0} is not a valid input file path".format(in_file))	
    if os.access(in_file, os.R_OK):
        return in_file
    else:
        raise argparse.ArgumentTypeError("{0} is not a readable input file".format(in_file))


def valid_outfile(out_file):
    """Check that output file path is writable. Note: This uses '/' character for splitting pathnames on Linux and Mac OSX. 
    The character may need to be changed to '\' for Windows executions"""
    outdir = '/'.join(out_file.split('/')[:-1])
    if not os.path.isdir(outdir):
        try: 
            # create the directory
            os.makedirs(outdir)
        except OSError:
            raise argparse.ArgumentTypeError("{0} is not a valid output directory".format(outdir))
    if os.access(outdir, os.W_OK):
        return out_file
    else:
        raise argparse.ArgumentTypeError("{0} is not a writable output directory".format(outdir))

if __name__ == "__main__":
    # Network Evaluation Setup Variables
    parser = argparse.ArgumentParser(description='Analyze network performance on ability to aggregate sets of nodes in network space.')
    parser.add_argument("network_path", type=valid_infile, 
        help='Path to file of network to be evaluated. File must be 2-column edge list where each line is a gene interaction separated by a common delimiter.')
    parser.add_argument("node_sets_file", type=valid_infile, 
        help='Path to file of node sets. Each line is a list, separated by a common delimiter. The first item in each line will be the name of the node set.')
    parser.add_argument("actual_AUPRCs_save_path", type=valid_outfile, 
        help='CSV file path of network evaluation result scores (AUPRCs). This script minimally returns these values to save. Must have a writable directory.')		
    parser.add_argument('-v', '--verbose', default=False, action="store_true", required=False,
        help='Verbosity flag for reporting on patient similarity network construction steps.')	
    parser.add_argument('-netd', '--net_file_delim', type=str, default='\t', required=False,
        help='Delimiter used in network file between columns. Default is tab white space.')
    parser.add_argument('-setd', '--set_file_delim', type=str, default='\t', required=False,
        help='Delimiter used in node set file to delimit lists. Default is tab white space.')	
    parser.add_argument("-p", "--sample_p", type=restricted_float, default=None, required=False,
        help='Sub-sampling percentage for node sets of interest. Default is None. Each gene set''s p is automatically determined by the network in this case.')
    parser.add_argument("-a", "--alpha", type=restricted_float, default=None, required=False,
        help='Propagation constant to use in the propagation of node sub-samples over given network. Overrides alpha calculation model if given.')
    parser.add_argument("-n", "--sub_sample_iter", type=positive_int, default=30, required=False,
        help='Number of times to perform sub-sampling during performance recovery (AUPRC) calculation for each node set. Default is 30.')
    parser.add_argument('-c', '--cores', type=positive_int, default=1, required=False,
        help='Number of cores to be utilized by machine for performance calculation step. NOTE: Each core must have enough memory to store at least network-sized square matrix and given node sets to perform calculations.')	
    parser.add_argument('-bg', '--background', type=str, default='network', choices=['genesets', 'network'], required=False,
        help='Establishes the background gene set to calculate AUPRC over. Default is to use all genes in the network, can change to use only genes from the union of all gene sets tested (i.e. disease genes only).')	
    parser.add_argument('-ming', '--min_genes', type=int, default=None, required=False,
        help='Enforce a minimum number of nodes from the network in each geneset. Any genesets without sufficient seed genes will be removed from the analysis')
    parser.add_argument('-w', '--weighted',default=False, action="store_true", required=False,
        help='Should weighted network propagation be performed?')
    # Network performance score calculations (with null networks)
    parser.add_argument("-i", "--null_iter", type=positive_int, default=30, required=False,
        help='Number of times to perform degree-preserved shuffling of network to construct performance value null distribution. Default is 30. If this value is >0, --null_AUPRCs_save_path will be required')
    parser.add_argument('-nno', '--null_network_outdir', type=valid_outfile, default=None, required=False,
        help='File directory to save null networks after generation.')
    parser.add_argument('-nsp', '--null_AUPRCs_save_path', type=valid_outfile, default=None, required=False,
        help='CSV file path of where to save null network evaluation results. Used in the calculation of network performance score and perfomance gain scores')
    parser.add_argument('-psp', '--performance_save_path', type=valid_outfile, default=None, required=False,
        help='CSV file path of where to save network evaluation results as z-scores.')
    parser.add_argument('-gsp', '--performance_gain_save_path', type=valid_outfile, default=None, required=False,
        help='CSV file path of where to save network evaluation results as gain in AUPRC over median null AUPRCs.')
    
    args = parser.parse_args()

    ####################################
    ##### Network Evaluation Setup #####
    ####################################
    if args.null_iter > 0:
    # A file path must be given to either save the null networks or the null network performance
        if (args.null_AUPRCs_save_path is None) and (args.null_network_outdir is None):
            parser.error('Save path required for null network edge lists or null network evaluation results.')
    # Load Network
    network = dit.load_edgelist_to_networkx(args.network_path, verbose=args.verbose)
    network_size = len(network.nodes())

    # Load Gene sets
    genesets = dit.load_node_sets(args.node_sets_file, verbose=args.verbose, id_type="Entrez")

    # Calculate gene set sub-sample rate with network (if not set)
    if args.sample_p is None:
        genesets_p, mean_coverage = nef.calculate_p(network, genesets, id_type="Entrez")
    else:
        _, mean_coverage = nef.calculate_p(network, genesets, id_type="Entrez")
        genesets_p = {geneset:args.sample_p for geneset in genesets}
    # Calculate network kernel (also determine propagation constant if not set)
    if args.alpha is None:
        alpha = prop.calculate_alpha(np.log10(len(network.edges)), mean_coverage, np.mean(list(genesets_p.values())))
    else:
        alpha = args.alpha
    kernel = nef.construct_prop_kernel(network, alpha=alpha, verbose=True)
    # Change background gene list if needed
    if args.background == 'genesets':
        background_node_set = set()
        for geneset in genesets:
            background_node_set = background_node_set.union(genesets[geneset])
        background_nodes = list(background_node_set.intersection(set(kernel.index)))
    else:
        background_nodes = list(kernel.index)

    ############################################
    ##### Network Performance Calculations #####
    ############################################

    # Calculate AUPRC for each gene set on actual network (large networks are >=10k nodes)
    if network_size < 10000:
        #TODO take background variable again
        actual_AUPRC_values, actual_FDR_values = nef.small_network_AUPRC_wrapper(kernel, genesets, genesets_p, n=args.sub_sample_iter, cores=args.cores, verbose=True, min_nodes=args.min_genes)
    else:
        #TODO take backgrond variable again
        actual_AUPRC_values, actual_FDR_values = nef.large_network_AUPRC_wrapper(kernel, genesets, genesets_p, n=args.sub_sample_iter, cores=args.cores, verbose=True, min_nodes=args.min_genes)

    # Save the actual network's AUPRC values
    actual_AUPRC_values.to_csv(args.actual_AUPRCs_save_path)

    #################################################
    ##### Null Network Performance Calculations #####
    #################################################

    # If number of null networks > 0:
    if args.null_iter > 0:
        null_AUPRCs = []
        null_FDRs = []
        for i in range(args.null_iter):
            print('Null iteration:', i+1)
            # Construct null networks and calculate AUPRCs for each gene set on each null network
            if args.null_network_outdir is not None:
                # check if samed netowrk exists
                try:
                    shuffNet = shuf.load_shuffled_network(datafile=args.network_path, outpath=args.null_network_outdir+'shuffNet_'+repr(i+1)+'_')
                except FileNotFoundError:
                    shuffNet = shuf.shuffle_network(network, n_swaps=1)
                    shuf.write_shuffled_network(shuffNet, datafile=args.network_path, outpath=args.null_network_outdir+'shuffNet_'+repr(i+1)+'_')
                    print('Shuffled Network', i+1, 'written to file')
            else:
                shuffNet = shuf.shuffle_network(network, n_swaps=1)
            # Save null network if null network output directory is given
                if args.null_network_outdir is not None:
                    shuf.write_shuffled_network(shuffNet, datafile=args.network_path, outpath=args.null_network_outdir+'shuffNet_'+repr(i+1)+'_')
                    if args.verbose:
                        print('Shuffled Network', i+1, 'written to file')
            # Construct null network kernel
            shuffNet_kernel = nef.construct_prop_kernel(shuffNet, alpha=alpha, verbose=False)
            # Calculate null network AUPRCs
            if network_size < 10000:
                #TODO make background an input again
                shuffNet_AUPRCs, shuffNet_FDRs = nef.small_network_AUPRC_wrapper(shuffNet_kernel, genesets, genesets_p, n=args.sub_sample_iter, cores=args.cores, verbose=True, min_nodes=args.min_genes)
            else:
                #TODO make background an input again
                shuffNet_AUPRCs, shuffNet_FDRs = nef.large_network_AUPRC_wrapper(shuffNet_kernel, genesets, genesets_p, n=args.sub_sample_iter, cores=args.cores, verbose=True, min_nodes=args.min_genes)
            null_AUPRCs.append(shuffNet_AUPRCs)
            null_FDRs.append(shuffNet_FDRs)
        # Construct table of null AUPRCs
        print('Calculating AUPRCs')
        null_AUPRCs_table = pd.concat(null_AUPRCs, axis=1)
        null_AUPRCs_table.columns = ['shuffNet'+repr(i+1) for i in range(len(null_AUPRCs))]
        null_FDRs_table = pd.concat(null_AUPRCs, axis=1)
        null_FDRs_table.columns = ['shuffNet'+repr(i+1) for i in range(len(null_AUPRCs))]
        if args.verbose:
            print('All null network gene set AUPRCs calculated')
        # Save null network AUPRCs if save path is given
        if args.null_AUPRCs_save_path is not None:
            null_AUPRCs_table.to_csv(args.null_AUPRCs_save_path)
        # Calculate performance score for each gene set's AUPRC if performance score save path is given
        if args.performance_save_path is not None:
            network_performance = nef.calculate_network_performance_score(actual_AUPRC_values, null_AUPRCs_table, verbose=args.verbose)			
            network_performance.to_csv(args.performance_save_path)
        # Calculate network performance gain over median null AUPRC if AUPRC performance gain save path is given
        if args.performance_gain_save_path is not None:
            network_perf_gain = nef.calculate_network_performance_gain(actual_AUPRC_values, null_AUPRCs_table, verbose=args.verbose)			
            network_perf_gain.to_csv(args.performance_gain_save_path)
    
