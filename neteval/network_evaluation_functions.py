#####################################################################
# ---------- Node Set-Based Network Evaluation Functions ---------- #
#####################################################################
from multiprocessing import Pool
import neteval.data_import_export_tools as dit
import neteval.network_propagation as prop
import neteval.shuffle_networks as shuf
import numpy as np
import os
import random
import scipy.stats as stats
import sklearn.metrics as metrics
import pandas as pd
import time
import warnings
warnings.filterwarnings(action='ignore', category=DeprecationWarning)


def precision_at_k(precision, ranks, k):
    """Calculate precision at k given a list of precision values and a dictionary of ranks.
    
    Args:
        precision (list): List of precision values.
        ranks (dict): Dictionary of all set genes, and their predicted rank in the network {gene: rank}
        k (int): Rank value to calculate precision at.
        
    Returns:
        float: Precision at k.
    """
    rank_values = list(ranks.values())
    if k in rank_values:
        return precision[rank_values.index(k)+1]
    elif k < min(rank_values):
        return 0
    elif k > max(rank_values):
        warnings.warn("k is greater than the maximum rank value. Precision will be calculated at the maximum rank value.")
        return precision[-1]
    else: # extrapolate the precision at k from last true prediction
        for i, rank in enumerate(rank_values):
            if rank >= k:
                k_idx = i
                break
    p_idx = precision[k_idx]
    tp = p_idx*rank_values[k_idx-1]
    p_at_k = tp/k
    return p_at_k


def calculate_p(network, nodesets,  size_coef=0.00933849, coverage_coef=-0.00126822, intercept=0.4369, id_type='Symbol', minp=0.1, maxp=0.8):
    """Calculate optimal sub-sampling proportion for test/train based on the optimization performed across all networks. Network coverage is 
    calculated as the number of nodes from a given set that are present in the network. 
    
    Args:
        network (NetworkX object): Network to calculate sub-sampling proportion for.
        nodesets (dict): Dictionary of {geneset name:list of genes}.
        size_coef (float): Coefficient for network size in p calculation.
        coverage_coef (float): Coefficient for coverage in p calculation.
        intercept (float): Intercept for p calculation.
        id_type (str): Type of gene identifier used in nodesets.
        minp (float): Minimum p (subsampling proportion) value.
        maxp (float): Maximum p (subsampling proportion) value.
        
    Returns:
        dict: Dictionary of {geneset name: optimal sub-sampling proportion}.
    
    """
    coverages = []
    if id_type == 'Entrez':
        network_nodes = [gene for gene in network.nodes()]
    else:
        network_nodes = [str(gene) for gene in network.nodes()]
    nodesets_p = {}
    for nodeset in nodesets:
        nodesets_coverage = len([node for node in nodesets[nodeset] if node in network_nodes])
        if nodesets_coverage == 0:
            nodesets_p[nodeset] = 0
        else:
            coverages.append(nodesets_coverage)
            nodesets_p[nodeset] = fit_p(net_size=np.log10(len(network_nodes)), coverage=nodesets_coverage, size_coef=size_coef, coverage_coef=coverage_coef, intercept=intercept, minp=minp, maxp=maxp)
    return nodesets_p, np.mean(np.array(coverages))


def fit_p(net_size, coverage, size_coef=0.00933849, coverage_coef=-0.00126822, intercept=0.4369, minp=0.1, maxp=0.8):
    """Fit a p value based on network size and coverage.
    
    Args:
        net_size (float): Log10 of network size.
        coverage (int): Number of nodes in the gene set that are in the network.
        size_coef (float): Coefficient for network size in p calculation.
        coverage_coef (float): Coefficient for coverage in p calculation.
        intercept (float): Intercept for p calculation.
        minp (float): Minimum p (subsampling proportion) value.
        maxp (float): Maximum p (subsampling proportion) value.
    
    Returns:
        float: Optimal sub-sampling proportion.
    
    """
    p =  intercept + net_size*size_coef + coverage * coverage_coef
    if p < minp:
        return minp
    elif p > maxp:
        return maxp
    return p


def construct_prop_kernel(network, alpha=0.5, m=-0.02935302, b=0.74842057, verbose=False, 
                        save_path=None, weighted=False):
    """Construct influence matrix of each network node propagated across network to use as kernel in AUPRC analysis.
    No propagation constant or alpha model required, can be calculated.
    
    Args:
        network (NetworkX object): Network to propagate.
        alpha (float): Propagation constant.
        m (float): Coefficient for alpha calculation.
        b (float): Intercept for alpha calculation.
        verbose (bool): Whether to print progress to stdout.
        save_path (str): Path to save network kernel.
        weighted (bool): Whether to use weighted network propagation. NotImplemented.
    
    Returns:
        pandas.DataFrame: Propagated network kernel.
    """
    # initialize network kernel
    network_Fo = pd.DataFrame(data=np.identity(len(network.nodes())), index=network.nodes(), columns=network.nodes())
    if weighted:
        raise NotImplementedError("Weighted network propagation not yet implemented.")
    if alpha is None:
        alpha_val = prop.calculate_alpha(np.log10(len(network.edges)), m=m, b=b)
    else:
        alpha_val = alpha
    network_Fn = prop.closed_form_network_propagation(network, network_Fo, alpha_val, verbose=verbose)
    network_Fn = network_Fn.loc[network_Fn.columns]
    if verbose:
        print('Propagated network kernel constructed')
    if save_path is not None:
        if save_path.endswith('.hdf'):
            network_Fn.to_hdf(save_path, key='Kernel', mode='w')
        else:
            network_Fn.to_csv(save_path)
    return network_Fn


def global_var_initializer(global_net_kernel):
    """Global variable initialization function for small network AUPRC calculations"""
    global kernel
    kernel = global_net_kernel


def calculate_small_network_AUPRC(params):
    """Calculate AUPRC of a single node set's recovery for small networks (<250k edges). This method is faster 
    for smaller networks, but still has a relatively large memory footprint. The parallel setup for this situation 
    requires passing the network kernel to each individual thread
    
    Args:
        params (list): List of parameters for AUPRC calculation. [node_set_name, node_set, p, n, bg, verbose]
    
    Returns:
        list: [node_set_name, AUPRC, recall_at_FDR]
    
    """
    node_set_name, node_set, p, n, bg, verbose, exclude_genes = params[0], params[1], params[2], params[3], params[4], params[5], params[6]
    runtime = time.time()
    intersect = [nodes for nodes in node_set if nodes in kernel.index]
    AUPRCs = []
    FDRs = []
    sample_size = int(round(p*len(intersect)))
    for i in range(n):																				# Number of times to run the sampling
        sample = random.sample(intersect, sample_size)													# get node set sample
        intersect_non_sample = [node for node in intersect if node not in sample]					   	# nodes in intersect not in sample
        bg_non_sample = [node for node in bg if node not in sample]							 			# nodes in background gene list not in sample
        bg_sample_sum = kernel.loc[sample][bg_non_sample].sum().sort_values(ascending=False)				# summed prop value for all nodes in background
        bg_sample_sum = bg_sample_sum[~bg_sample_sum.index.isin(exclude_genes)]
        y_actual = pd.Series(0, index=bg_sample_sum.index, dtype=int)									# nodes sorted by mean prop value
        y_actual.loc[intersect_non_sample]+=1															# which nodes in sorted list are in intersect_non_sample
        intersect_non_sample_sorted = y_actual[y_actual==1].index
        # intersect_non_sample sorted
        P_totals = {node:float(y_actual.loc[:node].shape[0]) for node in intersect_non_sample_sorted}
        geneset, auprc, fdr = calculate_precision_recall(node_set_name+'(rep '+str(i)+")", 
                        intersect_non_sample_sorted, P_totals, verbose=verbose)
        AUPRCs.append(auprc)
        FDRs.append(fdr)

    if verbose:
        print('AUPRC Analysis for given node set', '('+repr(len(intersect))+' nodes in network) complete:', round(time.time()-runtime, 2), 'seconds.')
    return [node_set_name, np.mean(AUPRCs), np.mean(FDRs)]


def calculate_precision_recall(geneset, intersect_non_sample_sorted, P_totals, verbose=False):
    """Calculate precision-recall curve for a given node set's recovery
    
    Args:
        geneset (str): Name of node set to calculate AUPRC for.
        intersect_non_sample_sorted (list): List of nodes in intersect_non_sample sorted by summed prop value.
        P_totals (dict): Dictionary of {node: total number of nodes in intersect_non_sample up to that node}.
        verbose (bool): Whether to print progress to stdout.
    
    Returns:
        tuple: (geneset, AUPRC, recall_at_FDR)
    
    """
    #geneset, intersect_non_sample_sorted, P_totals, verbose = params[0], params[1], params[2], params[3]
    runtime = time.time()
    TP, FN = 0, len(intersect_non_sample_sorted)	# initialize true positives and false negatives
    precision, recall = [1], [0]					# initialize precision and recall curves
    fdr = [1]
    for node in intersect_non_sample_sorted:		# Step down sorted nodes by summed prop value by nodes that are in intersect_non_sample
        TP += 1.0										# Calculate true positives found at this point in list
        FN -= 1.0										# Calculate false negatives found at this point in list
        precision.append(TP/float(P_totals[node]))		# Calculate precision ( TP / TP+FP ) and add point to curve
        recall.append(TP/float(TP+FN))		
        fdr.append(1-TP/float(P_totals[node]))# Calculate recall ( TP / TP+FN ) and add point to curve
    try:
        AUPRC = metrics.auc(recall, precision)
        recall_at_fdr = find_recall_at_fdr(fdr, recall, threshold=0.25)# Calculate Area Under Precision-Recall Curve (AUPRC)
    except ValueError:
        print("Geneset:", geneset)
        print("Intersect Length", len(intersect_non_sample_sorted))
        print("Precision:", precision)
        print("Recall:", recall)
        AUPRC = np.nan
        recall_at_fdr = np.nan
    if verbose:
        print('AUPRC Analysis for given node set:', geneset, 'complete:', round(time.time()-runtime, 2), 'seconds.')
    return geneset, AUPRC, recall_at_fdr


def find_recall_at_fdr(fdr_vals, recall, threshold=0.25):
    """Calculate recall at a given FDR threshold
    
    Args:
        fdr_vals (list): List of false discovery rates.
        recall (list): List of recall values.
        threshold (float): FDR threshold to calculate recall at.
        
    Returns:
        float: Recall at FDR threshold.
    
    """
    results = pd.DataFrame({'fdr':fdr_vals, 'recall':recall})
    low_fdr = results[results['fdr'] < threshold]
    if low_fdr.shape[0] > 0:
        return low_fdr['recall'].max()
    else:
        min_fdr = results['fdr'].min()
        return results[results['fdr'] == min_fdr]['recall'].max()


def calculate_large_network_AUPRC(params):
    """Calculate AUPRC of a single node set's recovery for large networks (>=250k edges). This method is slower than the 
    small network case, as the parallel setup would require too much memory for large networks.
    
    Args:
        params (list): List of parameters for AUPRC calculation. [geneset, intersect_non_sample_sorted, P_totals, verbose]
    
    Returns:
        list: [geneset, AUPRC, recall_at_FDR]
    """
    geneset, intersect_non_sample_sorted, P_totals, verbose = params[0], params[1], params[2], params[3]
    return calculate_precision_recall(geneset, intersect_non_sample_sorted, P_totals, verbose=verbose)


def small_network_AUPRC_wrapper(net_kernel, genesets, genesets_p, n=30, cores=1, bg=None, verbose=True, min_nodes=None, full_genesets=None):
    """Top level function to calculate AUPRC of multiple node sets' recovery for small networks (<10k edges)
        Args:
            net_kernel (pandas.DataFrame): Network kernel to use for AUPRC analysis
            genesets (dict): Dictionary of node sets to calculate AUPRC for
            genesets_p (dict): Dictionary sub-sampling parameters for each node set
            n (int): Number of sub-sampling iterations to perform
            cores (int): Number of cores to use for parallelization
            bg (list): List of nodes to use as background for AUPRC analysis
            verbose (bool): Whether to print progress to stdout
    """
    # Construct params list
    if bg is None:
        bg_intersect = list(net_kernel.index)
    else:
        bg_intersect = list(set(bg).intersection(set(net_kernel.index)))
    if min_nodes is not None:
        print("NUMBER INPUT GENESETS: ", len(genesets)	)
        genesets = {geneset:genesets[geneset] for geneset in genesets if len(set(genesets[geneset]).intersection(set(net_kernel.index))) >= min_nodes}
        print("NUMBER FILTERED GENESETS: ", len(genesets))
    # Parameters:
    # geneset (str): Name of node set to calculate AUPRC for
    # genesets (list/set): Nodes in the geneset
    # genesets_p (float): Proportion of nodes in geneset to sample
    # n (int): Number of sub-sampling iterations to perform
    # bg_intersect (list): List of nodes to use as background for AUPRC analysis
    if full_genesets is not None:
        AUPRC_Analysis_params = [[geneset, genesets[geneset], genesets_p[geneset], n, bg_intersect, verbose, [n for n in full_genesets[geneset] if n not in genesets[geneset]]] for geneset in genesets]
    else:
        AUPRC_Analysis_params = [[geneset, genesets[geneset], genesets_p[geneset], n, bg_intersect, verbose, []] for geneset in genesets]
    # Determine parallel calculation status
    if cores == 1:
        # Set network kernel
        global_var_initializer(net_kernel)
        # Calculate AUPRC values for all gene sets
        AUPRC_results = []
        for params_list in AUPRC_Analysis_params:
            AUPRC_results.append(calculate_small_network_AUPRC(params_list))
    else:
        # Initialize worker pool
        pool = Pool(cores, global_var_initializer, [net_kernel])
        # Run the AUPRC analysis for each geneset
        AUPRC_results = pool.map(calculate_small_network_AUPRC, AUPRC_Analysis_params)
        # Close worker pool
        pool.close()
    # Construct AUPRC results
    geneset_AUPRCs = {result[0]:result[1] for result in AUPRC_results}
    geneset_FDRs = {result[0]:result[2] for result in AUPRC_results}		
    AUPRCs_table = pd.Series(geneset_AUPRCs, name='AUPRC')
    FDRs_table = pd.Series(geneset_FDRs, name='Recall at FDR')
    return AUPRCs_table, FDRs_table


def large_network_AUPRC_wrapper(net_kernel, genesets, genesets_p, n=30, cores=1, bg=None, verbose=True, min_nodes=None, full_genesets=None):
    """ Wrapper to calculate AUPRC of multiple node sets' recovery for large networks (>=250k edges)
    
    Args:
        net_kernel (pandas.DataFrame): Network kernel to use for AUPRC analysis
        genesets (dict): Dictionary of node sets to calculate AUPRC for
        genesets_p (dict): Dictionary sub-sampling parameters for each node set
        n (int): Number of sub-sampling iterations to perform
        cores (int): Number of cores to use for parallelization
        bg (list): List of nodes to use as background for AUPRC analysis
        verbose (bool): Whether to print progress to stdout
        min_nodes (int): Minimum number of nodes in a gene set to be considered for AUPRC analysis
        
    Returns:
        tuple: (AUPRCs_table, FDRs_table)
    """
    starttime = time.time()
    # Construct binary gene set sub-sample matrix
    if min_nodes is not None:
        print("NUMBER INPUT GENESETS: ", len(genesets)	)
        genesets = {geneset:genesets[geneset] for geneset in genesets if len(set(genesets[geneset]).intersection(set(net_kernel.index))) >= min_nodes}
        print("NUMBER FILTERED GENESETS: ", len(genesets))
    geneset_list = list(genesets.keys())
    m, c = len(geneset_list), net_kernel.shape[0]
    subsample_mat = np.zeros((n*m, c))
    y_actual_mat = np.zeros((n*m, c))
    # Each block of length n rows is a sub-sampled binary vector of the corresponding gene set
    for i in range(m):
        geneset = geneset_list[i]
        # Get indices of gene set genes in kernel
        intersect = [gene for gene in genesets[geneset] if gene in net_kernel.index]
        index_dict = dict((gene, idx) for idx, gene in enumerate(net_kernel.index))
        intersect_idx = [index_dict[gene] for gene in intersect]
        # Generate n sub-samples
        for j in range(n):
            # Sub-sample gene set indices
            sample_size = int(round(genesets_p[geneset]*len(intersect)))
            sample_idx = random.sample(intersect_idx, sample_size)
            non_sample_idx = [idx for idx in intersect_idx if idx not in sample_idx]
            # Set sub-sampled list to 1
            row = (i*n)+j
            subsample_mat[row, sample_idx] = 1
            y_actual_mat[row, non_sample_idx] = 1
    if verbose:
        print('Binary gene set sub-sample matrix constructed')
    # Propagate sub-samples
    prop_subsamples = np.dot(subsample_mat, net_kernel)
    if verbose:
        print('Binary gene set sub-sample matrix propagated')
    # Construct parameter list to be passed
    AUPRC_Analysis_params = []
    for i in range(len(geneset_list)):
        if full_genesets is not None:
            exclude_genes = [n for n in full_genesets[geneset_list[i]] if n not in genesets[geneset_list[i]]]
        else:
            exclude_genes = []
        for j in range(n):
            row = (i*n)+j
            prop_result_full = pd.DataFrame(np.array((subsample_mat[row], y_actual_mat[row], prop_subsamples[row])), 
                                                    index=['Sub-Sample', 'Non-Sample', 'Prop Score'], columns=net_kernel.columns).T
            # Set background gene sets from a predefined gene set or all network genes
            if bg is None:
                prop_result = prop_result_full.sort_values(by=['Sub-Sample', 'Prop Score', 'Non-Sample'],
                                                            ascending=[False, False, False]).iloc[int(sum(subsample_mat[row])):]['Non-Sample']
            else:
                prop_result = prop_result_full.loc[bg].dropna().sort_values(by=['Sub-Sample', 'Prop Score', 'Non-Sample'],
                                                                    ascending=[False, False, False])
                print("ROW:", subsample_mat[row])
                print("SUM", sum(subsample_mat[row]))
                print("INDEXER:", int(sum(subsample_mat[row])))			
                prop_result = prop_result.iloc[int(sum(subsample_mat[row])):]['Non-Sample']
            prop_result = prop_result[~prop_result.index.isin(exclude_genes)]
            intersect_non_sample_sorted = prop_result[prop_result==1].index
            P_totals = {node:float(prop_result.loc[:node].shape[0]) for node in intersect_non_sample_sorted}
            AUPRC_Analysis_params.append([geneset_list[i], intersect_non_sample_sorted, P_totals, verbose])
    # Determine parallel calculation status
    if cores == 1:
        # Calculate AUPRC values for all gene sets
        AUPRC_results = []
        for params_list in AUPRC_Analysis_params:
            AUPRC_results.append(calculate_large_network_AUPRC(params_list))
    else:
        # Initialize worker pool
        pool = Pool(cores)
        # Run the AUPRC analysis for each geneset
        AUPRC_results = pool.map(calculate_large_network_AUPRC, AUPRC_Analysis_params)
        # Close worker pool
        pool.close()		  
    # Construct AUPRC results
    geneset_AUPRCs = pd.DataFrame(AUPRC_results, columns=['Gene Set', 'AUPRCs', "Recall_at_FDR"]).set_index('Gene Set', drop=True)
    geneset_AUPRCs_merged = {geneset:geneset_AUPRCs.loc[geneset]['AUPRCs'].mean() for geneset in geneset_list}
    geneset_FDRs_merged = {geneset:geneset_AUPRCs.loc[geneset]['Recall_at_FDR'].mean() for geneset in geneset_list}
    AUPRCs_table = pd.Series(geneset_AUPRCs_merged, name='AUPRC')
    FDRs_table = pd.Series(geneset_FDRs_merged, name='Recall_at_FDR')
    return AUPRCs_table, FDRs_table


def calculate_network_performance_score(actual_net_AUPRCs, shuff_net_AUPRCs, verbose=True, save_path=None):
    """ Calculate robust z-score metric for a network on given node sets given results directory of AUPRC calculations.
    Requires the AUPRCs calculated for the actual network in a pandas Series. Also requires the AUPRCs calculated for 
    the same gene sets on the shuffled networks in a pandas DataFrame
    
    Args:
        actual_net_AUPRCs (pandas.Series): AUPRC values for the actual network
        shuff_net_AUPRCs (pandas.DataFrame): AUPRC values for the shuffled networks
        verbose (bool): Whether to print progress to stdout
        save_path (str): Path to save network performance score
        
    Returns:
        pandas.Series: Network performance score
    """
    # Align data (only calculate for gene sets with full data on both actual networks and all shuffled networks)
    genesets = sorted(list(set(actual_net_AUPRCs.index).intersection(set(shuff_net_AUPRCs.index))), key=lambda s: s.lower())
    actual_net_AUPRCs = actual_net_AUPRCs.loc[genesets]
    shuff_net_AUPRCs = shuff_net_AUPRCs.loc[genesets]
    # Compute robust z-score for composite network performances
    k = 1/stats.norm.ppf(0.75)	# Mean absolute deviation scaling factor to make median absolute deviation behave similarly to the standard deviation of a normal distribution
    AUPRC_null_median = shuff_net_AUPRCs.median(axis=1)
    AUPRC_null_MAD = abs(shuff_net_AUPRCs.subtract(AUPRC_null_median, axis=0)).median(axis=1)
    AUPRC_null_MAD_scaled = k*AUPRC_null_MAD
    AUPRC_ZNorm = (actual_net_AUPRCs - AUPRC_null_median).divide(AUPRC_null_MAD_scaled)
    if save_path is not None:
        AUPRC_ZNorm.to_csv(save_path)
    if verbose:
        print('AUPRC values z-normalized')
    return AUPRC_ZNorm


def calculate_network_performance_gain(actual_net_AUPRCs, shuff_net_AUPRCs, verbose=True, save_path=None):
    """Calculate relative gain of actual network AUPRC over median random network AUPRC performance for each gene set
    Requires the AUPRCs calculated for the actual network in a pandas Series. Also requires the AUPRCs calculated for
    the same gene sets on the shuffled networks in a pandas DataFrame
    
    Args:
        actual_net_AUPRCs (pandas.Series): AUPRC values for the actual network
        shuff_net_AUPRCs (pandas.DataFrame): AUPRC values for the shuffled networks
        verbose (bool): Whether to print progress to stdout
        save_path (str): Path to save relative gain
        
    Returns:
        pandas.Series: Relative gain of actual network AUPRC over median random network AUPRC performance
    """
    # Align data (only calculate for gene sets with full data on both actual networks and all shuffled networks)
    genesets = sorted(list(set(actual_net_AUPRCs.index).intersection(set(shuff_net_AUPRCs.index))), key=lambda s: s.lower())
    actual_net_AUPRCs = actual_net_AUPRCs.loc[genesets]
    shuff_net_AUPRCs = shuff_net_AUPRCs.loc[genesets]	
    # Compute relative gain
    AUPRC_null_median = shuff_net_AUPRCs.median(axis=1)
    AUPRC_gain = (actual_net_AUPRCs - AUPRC_null_median).divide(AUPRC_null_median)
    if save_path is not None:
        AUPRC_gain.to_csv(save_path)
    if verbose:
        print('AUPRC relative performance gain calculated')
    return AUPRC_gain
