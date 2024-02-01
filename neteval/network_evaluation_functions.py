#####################################################################
# ---------- Node Set-Based Network Evaluation Functions ---------- #
#####################################################################
from multiprocessing import Pool
import neteval.data_import_export_tools as dit
import neteval.network_propagation as prop
import neteval.shuffle_networks as shuf
from neteval.Timer import Timer
import networkx as nx
import numpy as np
import os
import random
import scipy.stats as stats
import sklearn.metrics as metrics
import pandas as pd
import time
import warnings
from collections import defaultdict
warnings.filterwarnings(action='ignore', category=DeprecationWarning)

class NetworkAnalyzer:
    def __init__(self, network, genesets, outdir, out_pref, network_path, alpha=0.5, sample_proportion=0.2, 
                    num_samples=50, null_iterations=200, weighted=False, background=None,
                    verbose=False, k=20, min_genes=0):
        self.timer = Timer()
        self.timer.start('Overall')
        self.timer.start('Initialization')
        self.num_samples = num_samples
        self.null_iterations = null_iterations
        self.alpha = alpha
        if isinstance(sample_proportion, float):
            self.sample_proportions = {geneset:sample_proportion for geneset in genesets}
        else:
            self.sample_proportions = sample_proportion
        self.genesets = genesets
        self.network_path = network_path
        self.network = network
        self.out_pref = out_pref
        self.weighted = weighted
        self.verbose = verbose
        self.k = k
        self.min_genes = min_genes
        ## Get the network propagation kernel
        self.timer.start('Kernel Construction')
        self.kernel = construct_prop_kernel(self.network, alpha=self.alpha, 
                                            verbose=self.verbose, weighted=self.weighted)
        self.timer.end('Kernel Construction')
        ## Construct the background set as the background of all genes in at least one geneset
        if background is None:
            self.background = list(network.nodes())
        else:
            self.background = self.create_geneset_background()
        # check that all output paths exist
        self.outdir = outdir
        if not os.path.exists(os.path.dirname(self.outdir)):
            os.makedirs(os.path.dirname(self.outdir))
            if verbose:
                print("Created directory:", os.path.dirname(self.outdir))
        for folder in ['AURPCs', 'AUPRCs/Null_AUPRCs', 'Performance', 'Performance_Gain', 'Null_Networks']:
            if not os.path.exists(os.path.join(self.outdir, folder)):
                os.makedirs(os.path.join(self.outdir, folder))
                if verbose:
                    print("Created directory:", os.path.join(self.outdir, folder))    
        self.timer.start("Filter GeneSets")
        self.filter_genesets()
        self.timer.end("Filter GeneSets")
        
        #self.out_pref = os.path.join(outdir, "_".join([self.net_pref, str(self.alpha), str(self.sample_proportion), str(self.num_samples), str(self.null_iterations)]))

        self.timer.end('Initialization')
    
    def create_geneset_background(self):
        background_node_set = set()
        for geneset in self.genesets:
            background_node_set = background_node_set.union(self.genesets[geneset])
        background_nodes = list(background_node_set.intersection(set(self.kernel.index)))
        return background_nodes
    
    def filter_genesets(self):
        if self.min_genes is not None:
            if self.min_genes > 0:
                if self.verbose:
                    print("NUMBER INPUT GENESETS: ", len(self.genesets)	)
                filtered_genesets = {geneset:self.genesets[geneset] for geneset in self.genesets if len(set(self.genesets[geneset]).intersection(set(self.kernel.index))) >= self.min_genes}
                if self.verbose:
                    print("NUMBER FILTERED GENESETS: ", len(filtered_genesets))
                self.genesets = filtered_genesets
                self.sample_proportions = {geneset:self.sample_proportions[geneset] for geneset in self.sample_proportions if geneset in self.genesets}
        
    #     def _parellelize_analysis(self, node_set):
    #     intersect = [nodes for nodes in self.genesets[node_set] if nodes in self.kernel.index]
    #     subsample_result = self.perform_geneset_subsample_propagations(intersect)
    #     null_result = self.perform_null_geneset_propagations(intersect, node_set)
    #     return subsample_result, null_result
    
    # def perform_analysis(self):
    #     self.timer.start('Analysis')
    #     subsample_results = {}
    #     null_results = {}
    #     batch = 0
    #     # create a multiprocessing pool with a max of 10 processes
    #     pool = Pool(processes=5)
    #     geneset_names = list(self.genesets.keys())
    #     for i, output in enumerate(pool.map(self._parellelize_analysis, geneset_names)):
    #         subsample_results[geneset_names[i]], null_results[geneset_names[i]] = output[0], output[1]
    
    def perform_analysis(self):
        self.timer.start('Analysis')
        subsample_results = {}
        batch = 0
        self.timer.start("Batch " + str(batch*100) + "-" + str((batch+1)*100))
        for i, node_set in enumerate(self.genesets):
            intersect = [nodes for nodes in self.genesets[node_set] if nodes in self.kernel.index]
            if len(intersect) < self.min_genes:
                continue
            if (i+1) % 100 == 0:
                self.timer.end("Batch " + str(batch*100) + "-" + str((batch+1)*100))
                batch += 1
                self.timer.start("Batch " + str(batch*100) + "-" + str((batch+1)*100))
            subsample_results[node_set]= self.perform_geneset_subsample_propagations(self.kernel, intersect, self.sample_proportions[node_set])
        self.timer.end("Batch " + str(batch*100) + "-" + str((batch+1)*100))
        self.subsample_results = self.create_multiIndex_results(subsample_results)
        self.timer.end('Analysis')
        
        
    def perform_shuffling_analysis(self):
        self.timer.start('Shuffled Propagations')
        all_null_results = {}
        for i in range(self.null_iterations):
            self.timer.start('Shuffled Propagation '+str(i))
            if self.verbose:
                print('Null Iteration:', i, 'of', self.null_iterations)
            all_null_results[i] = self.perform_shuffled_propagations(i)
            self.timer.end('Shuffled Propagation '+str(i))
        ## TODO check that this is the right format. 
        self.null_results = pd.concat(all_null_results, axis=0).swaplevel(0,1).sort_index() 
        self.timer.end('Shuffled Propagations')
        
    def perform_shuffled_propagations(self, i):
        # get shuffled network
        self.timer.start('Load shuffled network')
        try:
            shuffNet = shuf.load_shuffled_network(datafile=self.network_path, outpath=os.path.join(self.outdir, 'Null_Networks', 'shuffNet_' +repr(i+1)+'_'))
        except FileNotFoundError:
            warnings.warn("Shuffled network not found. Creating new shuffled network.")
            shuffNet = shuf.shuffle_network(self.network, n_swaps=1)
            shuf.write_shuffled_network(shuffNet, datafile=self.network_path, outpath=os.path.join(self.outdir, 'Null_Networks', 'shuffNet_' +repr(i+1)+'_'))
            if self.verbose:
                print('Shuffled Network', i+1, 'written to file')
        self.timer.end('Load shuffled network')
        # get shuffled kernel
        self.timer.start('Construct shuffled kernel')
        shuff_kernel = construct_prop_kernel(shuffNet, alpha=self.alpha, 
                                            verbose=self.verbose, weighted=self.weighted)
        self.timer.end('Construct shuffled kernel')
        self.timer.start('Propagation Analysis')
        null_results = {}
        for j, node_set in enumerate(self.genesets):
            intersect = [nodes for nodes in self.genesets[node_set] if nodes in shuff_kernel.index]
            if len(intersect) < self.min_genes:
                continue
            if (j+1) % 100 == 0:
                if self.verbose:
                    print( j, "of", len(self.genesets), 'completed.')
            null_results[node_set] = self.perform_geneset_subsample_propagations(shuff_kernel, intersect, self.sample_proportions[node_set])
        self.null_results = self.create_multiIndex_results(null_results)
        self.timer.end('Propagation Analysis')

        return self.create_multiIndex_results(null_results)
    
    # def perform_analysis_old(self):
    #     self.timer.start('Analysis')
    #     subsample_results = {}
    #     null_results = {}
    #     batch = 0
    #     self.timer.start("Batch " + str(batch*100) + "-" + str((batch+1)*100))
    #     for i, node_set in enumerate(self.genesets):
    #         if self.verbose:
    #             print("Performing analysis for node set:", node_set, i, "of", len(self.genesets))
    #         intersect = [nodes for nodes in self.genesets[node_set] if nodes in self.kernel.index]
    #         if len(intersect) < self.min_genes:
    #             continue
    #         subsample_results[node_set]= self.perform_geneset_subsample_propagations(intersect, self.sample_proportions[node_set])
    #         null_results[node_set] = self.perform_null_geneset_propagations(intersect, self.sample_proportions[node_set])
    #         if (i+1) % 100 == 0:
    #             self.timer.end("Batch " + str(batch*100) + "-" + str((batch+1)*100))
    #             batch += 1
    #             self.timer.start("Batch " + str(batch*100) + "-" + str((batch+1)*100))
    #     self.timer.end("Batch " + str(batch*100) + "-" + str((batch+1)*100))
    #     self.subsample_results = self.create_multiIndex_results(subsample_results)
    #     #self.subsample_results_stds = pd.DataFrame.from_dict(subsample_results['stds'], orient='index')
    #     self.null_results = self.create_multiIndex_results(null_results)
    #     #self.null_results_stds = pd.DataFrame.from_dict(null_results['stds'], orient='index')
    #     self.timer.end('Analysis')
        

    
    def create_multiIndex_results(self, results_dict):
        rekeyed_dict= {(outerKey, innerKey): values for outerKey, innerDict in results_dict.items() for innerKey, values in innerDict.items()}
        return pd.DataFrame.from_dict(rekeyed_dict, orient='index')
    
    def perform_geneset_subsample_propagations(self, kernel, node_set, sample_p):
        sample_size = int(round(sample_p*len(node_set)))
        results_dict = {}
        for i in range(self.num_samples):
            sample = random.sample(node_set, sample_size)
            results_dict[i] = self.perform_single_propagation(kernel, node_set, sample)
        summary_results = pd.DataFrame.from_dict(results_dict, orient='index').describe(percentiles=[]).loc[['mean', 'std', '50%']]
        return summary_results.to_dict(orient='index')
    
    def perform_single_propagation(self, kernel, node_set, sample):
        test_set = [node for node in node_set if node not in sample]
        bg_non_sample = [node for node in self.background if node not in sample]							 			# nodes in background gene list not in sample
        propagation_results = kernel.loc[sample][bg_non_sample].sum().sort_values(ascending=False)
        y_actual = pd.Series(0, index=propagation_results.index, dtype=int)									# nodes sorted by mean prop value
        y_actual.loc[test_set]+=1
        sorted_test_set = y_actual[y_actual==1].index # get the test_set ordered by propagation score
        test_ranks = {node:float(y_actual.loc[:node].shape[0]) for node in sorted_test_set}
        metrics = calculate_precision_recall_metrics(sorted_test_set, test_ranks, verbose=self.verbose, k=self.k)
        return metrics
    
    # def perform_null_geneset_propagations(self, geneset, sample_p):
    #     sample_size = int(round(sample_p*len(geneset)))
    #     results_dict = {}
    #     for i in range(self.null_iterations):
    #         null_set = self.get_degree_matched_nodeset(geneset)
    #         null_sample = random.sample(null_set, sample_size)
    #         results_dict[i] = self.perform_single_propagation(geneset, null_sample)
    #     results_df = pd.DataFrame.from_dict(results_dict, orient='index')
    #     summary_results = results_df.describe(percentiles=[]).loc[['mean', 'std', '50%']]
    #     null_medians = summary_results.loc['50%']
    #     null_MAD = abs(results_df.subtract(null_medians)).median(axis=0)
    #     k = 1/stats.norm.ppf(0.75)	# Mean absolute deviation scaling factor to make median absolute deviation behave similarly to the standard deviation of a normal distribution
    #     null_MAD_scaled = k*null_MAD
    #     return {**summary_results.to_dict(orient='index'), **{'scaled_MAD':null_MAD_scaled.to_dict()}}
    #     #return summary_results.loc['mean'].to_dict(), summary_results.loc['std'].to_dict()

            
    # def get_degree_matched_nodeset(self, nodeset):
    #     null_set = []
    #     for node in nodeset:
    #         random_node = self.BinnedNodes.random_degree_matched_node(self.network.degree(node))
    #         while random_node in null_set: # in case a duplicate node is selected
    #             random_node = self.BinnedNodes.random_degree_matched_node(self.network.degree(node))
    #         null_set.append(random_node)
    #     return null_set
    
    def calculate_performance_metrics(self):
        assert hasattr(self, 'subsample_results'), "Subsample results have not been calculated. Please run perform_analysis()"
        self.timer.start('Performance Calculation')
        self.t_results = self.perform_t_analysis()
        self.robust_z_results = self.perform_robust_z_analysis()
        self.timer.end('Performance Calculation')
        
    def write_results(self):
        self.timer.start("Write outputs")
        # TODO do I want to switch the format of this? take a look at the two and decide which is more intuitive. 
        # get AUPRCs metics
        pd.DataFrame(self.get_all_metrics('AUPRC')).to_csv(self.out_pref + '_AUPRCs.tsv', sep="\t")
        pd.DataFrame(self.get_all_metrics('recall_at_fdr')).to_csv(self.out_pref+'_Recall_at_FDR.tsv', sep="\t")
        pd.DataFrame(self.get_all_metrics('precision_at_'+ str(self.k))).to_csv(self.out_pref+'_Precision_at'+ str(self.k)+'.tsv', sep="\t")
    
    def get_all_metrics(self, measure):
        means, medians = self.subsample_results.xs('mean', level=1)[measure], self.subsample_results.xs('50%', level=1)[measure]
        t, pt, gaint = self.t_results.xs('stat', level=1)[measure], self.t_results.xs('pval', level=1)[measure], self.t_results.xs('gain', level=1)[measure]
        z, pz, gainz = self.robust_z_results.xs('stat', level=1)[measure], self.robust_z_results.xs('pval', level=1)[measure], self.robust_z_results.xs('gain', level=1)[measure]
        return {'mean':means, 'median':medians, 't':t, 'pval_t':pt, 'gain_t':gaint, 'z':z, 'pval_z':pz, 'gain_z':gainz}
        
        
    def perform_t_analysis(self):
        subsample_ses = self.subsample_results.xs('std', level=1) / np.sqrt(self.num_samples)
        null_ses = self.null_results.xs('std', level=1) / np.sqrt(self.null_iterations)
        sqrt_se1_se0 = np.sqrt((subsample_ses ** 2) + (null_ses ** 2))
        observed_means = self.subsample_results.xs('mean', level=1)
        null_means = self.null_results.xs('mean', level=1)
        t_scores = (observed_means - null_means) / sqrt_se1_se0
        t_sigs = t_scores.applymap(lambda x: stats.t.sf(np.abs(x), self.num_samples + self.null_iterations - 2))
        gain = (observed_means - null_means) / null_means
        t_results = pd.concat({'stat': t_scores, 'pval':t_sigs, 'gain':gain}, axis=0).swaplevel(0,1).sort_index()   
        return t_results
    
    def perform_robust_z_analysis(self):
        mean_minus_median = self.subsample_results.xs('mean', level=1) - self.null_results.xs('50%', level=1)
        z_norm = mean_minus_median / self.null_results.xs('scaled_MAD', level=1)
        z_gain = mean_minus_median / self.null_results.xs('50%', level=1)
        z_sig = z_norm.applymap(lambda x: stats.norm.sf(np.abs(x)))
        z_results = pd.concat({'stat': z_norm, 'pval':z_sig, 'gain':z_gain}, axis=0).swaplevel(0,1).sort_index()
        return z_results

    # class DegreeBinnedNodes:
    #     def __init__(self, network, min_bin_size):
    #         degree_dict = defaultdict(list)
    #         for node in network.nodes():
    #             degree_dict[network.degree(node)].append(node)
    #         self.bins = [0]
    #         self.bin_sizes = []
    #         bin_size = 0
    #         for degree in sorted(degree_dict.keys()):
    #             if bin_size < min_bin_size:
    #                 bin_size += len(degree_dict[degree])
    #             else:
    #                 self.bins.append(degree)
    #                 self.bin_sizes.append(bin_size)
    #                 bin_size = 0
    #         if bin_size < min_bin_size:
    #             # there are not sufficient nodes in final bin so remove it, and add the final bin
    #             self.bins.pop()
    #             self.bin_sizes[-1] += bin_size
    #         self.bins.append(degree+1)

    #         self.bin_dict = defaultdict(list)
    #         for i, left_edge in enumerate(self.bins[:-1]):
    #             for degree in range(left_edge, self.bins[i+1]):
    #                 self.bin_dict[left_edge] += degree_dict[degree]
    #         # map all degrees to their bin
    #         self.degree_to_bin = {deg: left_edge for i, left_edge in enumerate(self.bins[:-1]) for deg in range(left_edge, self.bins[i+1])}
            
    #     def random_degree_matched_node(self, degree):
    #         return random.choice(self.bin_dict[self.degree_to_bin[degree]])
    
    # def finish_analysis(self):
    #     self.timer.end('Overall')
    #     self.timer.print_all_times()
            

def calculate_precision_recall_metrics(test_set_sorted, test_ranks, verbose=False, k=None):
    #geneset, intersect_non_sample_sorted, P_totals, verbose = params[0], params[1], params[2], params[3]
    runtime = time.time()
    TP, FN = 0, len(test_set_sorted)	# initialize true positives and false negatives
    precision, recall = [1], [0]					# initialize precision and recall curves
    fdr = [1]
    results = {}
    for node in test_set_sorted:		# Step down sorted nodes by summed prop value by nodes that are in intersect_non_sample
        TP += 1.0										# Calculate true positives found at this point in list
        FN -= 1.0										# Calculate false negatives found at this point in list
        precision.append(TP/float(test_ranks[node]))		# Calculate precision ( TP / TP+FP ) and add point to curve
        recall.append(TP/float(TP+FN))		
        fdr.append(1-TP/float(test_ranks[node]))# Calculate recall ( TP / TP+FN ) and add point to curve
    results['AUPRC'] = metrics.auc(recall, precision)
    results['recall_at_fdr'] = find_recall_at_fdr(fdr, recall, threshold=0.25)# Calculate Area Under Precision-Recall Curve (AUPRC)
    if k is not None:
        results['precision_at_' + str(k)] = precision_at_k(precision, test_ranks, k)
    if verbose:
        print('AUPRC Analysis for given node set complete:', round(time.time()-runtime, 2), 'seconds.')
    return results


def precision_at_k(precision, ranks, k):
    rank_values = list(ranks.values())
    if k in rank_values:
        return precision[rank_values.index(k)+1]
    elif k < min(rank_values):
        return 0
    elif k > max(rank_values):
        warnings.warn("k is greater than the maximum rank value. Precision will be calculated at the maximum rank value.")
        return precision[-1]
    else: # impute the precision at k from last true prediction
        for i, rank in enumerate(rank_values):
            if rank >= k:
                k_idx = i
                break
    p_idx = precision[k_idx]
    tp = p_idx*rank_values[k_idx-1]
    p_at_k = tp/k
    return p_at_k



##### OLD FUNCTIONS #####

# Shuffle network in degree-preserving manner
# Input: network - networkx formatted network
# For large networks this can be slow: may need to be sped up to prevent bottlenecking
def old_shuffle_network(network, max_tries_n=10, verbose=False):
    # Shuffle Network
    shuff_time = time.time()
    edge_len=len(network.edges())
    shuff_net=network.copy()
    try:
        nx.double_edge_swap(shuff_net, nswap=edge_len, max_tries=edge_len*max_tries_n)
    except:
        if verbose:
            print('Note: Maximum number of swap attempts ('+repr(edge_len*max_tries_n)+') exceeded before desired swaps achieved ('+repr(edge_len)+').')
    if verbose:
        # Evaluate Network Similarity
        shared_edges = len(set(network.edges()).intersection(set(shuff_net.edges())))
        print('Network shuffled:', time.time()-shuff_time, 'seconds. Edge similarity:', shared_edges/float(edge_len))
    return shuff_net

# Calculate optimal sub-sampling proportion for test/train
# Input: NetworkX object and dictionary of {geneset name:list of genes}
def calculate_p(network, nodesets,  size_coef=0.00933849, coverage_coef=-0.00126822, intercept=0.4369, id_type='Symbol', minp=0.1, maxp=0.8):
    coverages = []
    if id_type == 'Entrez':
        network_nodes = [gene for gene in network.nodes()]
    else:
        network_nodes = [str(gene) for gene in network.nodes()]
    nodesets_p = {}
    for nodeset in nodesets:
        nodesets_coverage = len([node for node in nodesets[nodeset] if node in network_nodes])
        if nodesets_coverage == 0:
            # TODO this should no longer be needed when all gene sets have been filtered against all networks.
            nodesets_p[nodeset] = 0
        else:
            coverages.append(nodesets_coverage)
            nodesets_p[nodeset] = fit_p(net_size=np.log10(len(network_nodes)), coverage=nodesets_coverage, size_coef=size_coef, coverage_coef=coverage_coef, intercept=intercept, minp=minp, maxp=maxp)
    return nodesets_p, np.mean(np.array(coverages))

def fit_p(net_size, coverage, size_coef=0.00933849, coverage_coef=-0.00126822, intercept=0.4369, minp=0.1, maxp=0.8):
    p =  intercept + net_size*size_coef + coverage * coverage_coef
    if p < minp:
        return minp
    elif p > maxp:
        return maxp
    return p


# Construct influence matrix of each network node propagated across network to use as kernel in AUPRC analysis
# Input: NetowkrkX object. No propagation constant or alpha model required, can be calculated
def construct_prop_kernel(network, alpha=0.5, m=-0.02935302, b=0.74842057, verbose=False, 
                        save_path=None, weighted=False):
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


# Global variable initialization function for small network AUPRC calculations
def global_var_initializer(global_net_kernel):
    global kernel
    kernel = global_net_kernel


def calculate_network_AUPRC(params):
    runtime = time.time()
    
    
def network_AUPRC_wrapper(params):
    runtime = time.time()


# Calculate AUPRC of a single node set's recovery for small networks (<250k edges)
# This method is faster for smaller networks, but still has a relatively large memory footprint
# The parallel setup for this situation requires passing the network kernel to each individual thread
def calculate_small_network_AUPRC(params):
    node_set_name, node_set, p, n, bg, verbose = params[0], params[1], params[2], params[3], params[4], params[5]
    runtime = time.time()
    # print("Kernel_Index:",  kernel.index)
    # print("Node set:", node_set)
    intersect = [nodes for nodes in node_set if nodes in kernel.index]
    AUPRCs = []
    FDRs = []
    sample_size = int(round(p*len(intersect)))
    for i in range(n):																				# Number of times to run the sampling
        sample = random.sample(intersect, sample_size)													# get node set sample
        intersect_non_sample = [node for node in intersect if node not in sample]					   	# nodes in intersect not in sample
        bg_non_sample = [node for node in bg if node not in sample]							 			# nodes in background gene list not in sample
        bg_sample_sum = kernel.loc[sample][bg_non_sample].sum().sort_values(ascending=False)				# summed prop value for all nodes in background
        y_actual = pd.Series(0, index=bg_sample_sum.index, dtype=int)									# nodes sorted by mean prop value
        y_actual.loc[intersect_non_sample]+=1															# which nodes in sorted list are in intersect_non_sample
        intersect_non_sample_sorted = y_actual[y_actual==1].index
        # intersect_non_sample sorted
        P_totals = {node:float(y_actual.loc[:node].shape[0]) for node in intersect_non_sample_sorted}
        geneset, auprc, fdr = calculate_precision_recall(node_set_name+'(rep '+str(i)+")", 
                        intersect_non_sample_sorted, P_totals, verbose=verbose)
        AUPRCs.append(auprc)
        FDRs.append(fdr)
        # TP, FN = 0, len(intersect_non_sample_sorted)													# initialize precision and recall curves
        # precision, recall = [1], [0]																	# initialize true positives and false negatives
        # for node in intersect_non_sample_sorted:														# Slide down sorted nodes by summed prop value by nodes that are in intersect_non_sample
        # 	TP += 1.0									   													# Calculate true positives found at this point in list
        # 	FN -= 1.0																					   	# Calculate false negatives found at this point in list
        # 	precision.append(TP/float(y_actual.loc[:node].shape[0]))										 	# Calculate precision ( TP / TP+FP ) and add point to curve
        # 	recall.append(TP/float(TP+FN))
        # try:# Calculate recall ( TP / TP+FN ) and add point to curve
        # 	AUPRCs.append(metrics.auc(recall, precision))
        # except ValueError as e:
        # 	print("Recall:", recall)
        # 	print("Preceision:", precision)
        # 	print("Intersect_non_sample_sorted:", intersect_non_sample_sorted)
        # 	print("Sample:", sample)
        # 	raise(e)# Calculate Area Under Precision-Recall Curve (AUPRC)
    if verbose:
        print('AUPRC Analysis for given node set', '('+repr(len(intersect))+' nodes in network) complete:', round(time.time()-runtime, 2), 'seconds.')
    return [node_set_name, np.mean(AUPRCs), np.mean(FDRs)]


def calculate_precision_recall(geneset, intersect_non_sample_sorted, P_totals, verbose=False):
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
    results = pd.DataFrame({'fdr':fdr_vals, 'recall':recall})
    low_fdr = results[results['fdr'] < threshold]
    if low_fdr.shape[0] > 0:
        return low_fdr['recall'].max()
    else:
        min_fdr = results['fdr'].min()
        print("!!MIN FDR is:", min_fdr)
        return results[results['fdr'] == min_fdr]['recall'].max()


# Caclulate AUPRC of a single node set's recovery for large networks (>=250k edges)
# This method is slower than the small network case, as well as forces the memory footprint to be too large
# The parallel setup for this situation requries 
def calculate_large_network_AUPRC(params):
    geneset, intersect_non_sample_sorted, P_totals, verbose = params[0], params[1], params[2], params[3]
    return calculate_precision_recall(geneset, intersect_non_sample_sorted, P_totals, verbose=verbose)
    # runtime = time.time()
    # TP, FN = 0, len(intersect_non_sample_sorted)	# initialize true positives and false negatives
    # precision, recall = [1], [0]					# initialize precision and recall curves
    # for node in intersect_non_sample_sorted:		# Step down sorted nodes by summed prop value by nodes that are in intersect_non_sample
    # 	TP += 1.0										# Calculate true positives found at this point in list
    # 	FN -= 1.0										# Calculate false negatives found at this point in list
    # 	precision.append(TP/float(P_totals[node]))		# Calculate precision ( TP / TP+FP ) and add point to curve
    # 	recall.append(TP/float(TP+FN))					# Calculate recall ( TP / TP+FN ) and add point to curve
    # AUPRC = metrics.auc(recall, precision)				# Calculate Area Under Precision-Recall Curve (AUPRC)
    # if verbose:
    # 	print('AUPRC Analysis for given node set:', geneset, 'complete:', round(time.time()-runtime, 2), 'seconds.')
    # return [geneset, AUPRC]

# Wrapper to calculate AUPRC of multiple node sets' recovery for small networks (<250k edges)
def small_network_AUPRC_wrapper(net_kernel, genesets, genesets_p, n=30, cores=1, bg=None, verbose=True, min_nodes=None):
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
    AUPRC_Analysis_params = [[geneset, genesets[geneset], genesets_p[geneset], n, bg_intersect, verbose] for geneset in genesets]
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

# Wrapper to calculate AUPRC of multiple node sets' recovery for large networks (>=250k edges)
def large_network_AUPRC_wrapper(net_kernel, genesets, genesets_p, n=30, cores=1, bg=None, verbose=True, min_nodes=None):
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

# Wrapper to calculate AUPRCs of multiple node sets given network and node set files
def AUPRC_Analysis_single(network_file, genesets_file, shuffle=False, kernel_file=None, prop_constant=None, 
                        subsample_iter=30, cores=1, geneset_background=False, save_path=None, verbose=True):
    starttime = time.time()
    # Load network
    network = dit.load_edgelist_to_networkx(network_file, verbose=verbose)
    # Shuffle network?
    if shuffle:
        network = shuf.shuffle_network(network, verbose=verbose)
    # Get network size
    net_nodes = network.nodes()
    net_size = len(net_nodes)
    if verbose:
        print('Network size:', net_size, 'Nodes')
    # Calculate or load network propagation kernel
    if kernel_file is None:
        # Determine propagation constant
        if prop_constant is None:
            alpha = prop.calculate_alpha(network)
        else:
            alpha = prop_constant
        # Calculate network propagation kernel
        net_kernel = construct_prop_kernel(network, alpha=alpha, verbose=verbose)
    else:
        # Load network propagation kernel
        if kernel_file.endswith('.hdf'):
            net_kernel = pd.read_hdf(kernel_file)
        else:
            net_kernel = pd.read_csv(kernel_file)
    # Load node sets to recover
    genesets = dit.load_node_sets(genesets_file, verbose=verbose)
    # Calculate sub-sample rate for each node set given network
    genesets_p = calculate_p(network, genesets)
    non_represented_gene_sets = [gset for gset in genesets_p if genesets_p[gset]==0]
    for gset in non_represented_gene_sets:
        print("WARNING:", gset, "removed due to insufficient coverage.")
        genesets.pop(gset)
        genesets_p.pop(gset)
    # Set background of genes to recover as all network nodes or union of all gene sets' genes
    if geneset_background:
        background_gene_set = set()
        for geneset in genesets:
            background_gene_set = background_gene_set.union(genesets[geneset])
        background_genes = list(background_gene_set.intersection(set(net_nodes)))
    else:
        background_genes = list(net_nodes)
    # if network is small:
    if net_size < 10000:
        AUPRC_table = small_network_AUPRC_wrapper(net_kernel, genesets, genesets_p, n=subsample_iter, cores=cores, bg=background_genes, verbose=verbose)
    # if network is large:
    elif (net_size >= 10000) & (net_size < 15000):
        AUPRC_table = large_network_AUPRC_wrapper(net_kernel, genesets, genesets_p, n=subsample_iter, cores=cores, bg=background_genes, verbose=verbose)
    # if network is large:
    else:
        #TODO why was cores=1 set for large networks?
        AUPRC_table = large_network_AUPRC_wrapper(net_kernel, genesets, genesets_p, n=subsample_iter, cores=cores, bg=background_genes, verbose=verbose)
    if verbose:
        print('AUPRC values calculated', time.time()-starttime, 'seconds')
    # Save table
    if save_path is not None:
        AUPRC_table.to_csv(save_path)
    if verbose:
        print('AUPRC table saved:', save_path)
    return AUPRC_table

# The function will take all files containing the filename marker given to shuff_net_AUPRCs_fn and construct a single null AUPRCs table from them (in wd)
# shuff_net_AUPRCs_fn is a generic filename marker (assumes all shuff_net_AUPRCs files have the same file name structure)
def get_null_AUPRCs_table(wd, shuff_net_AUPRCs_fn, geneset_list=None):
    shuff_net_AUPRCs = [pd.read_csv(wd+fn, index_col=0, header=-1) for fn in os.listdir(wd) if shuff_net_AUPRCs_fn in fn]
    shuff_net_AUPRCs = pd.concat(shuff_net_AUPRCs, axis=1)
    if geneset_list is None:
        return shuff_net_AUPRCs
    else:
        return shuff_net_AUPRCs.loc[geneset_list].dropna(axis=1)

# Calculate robust z-score metric for a network on given node sets given results directory of AUPRC calculations
# Requires the AUPRCs calculated for the actual network in a pandas Series
# Also requires the AUPRCs calculated for the same gene sets on the shuffled networks in a pandas DataFrame
def calculate_network_performance_score(actual_net_AUPRCs, shuff_net_AUPRCs, verbose=True, save_path=None):
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

# Calculate relative gain of actual network AUPRC over median random network AUPRC performance for each gene set
# Requires the AUPRCs calculated for the actual network in a pandas Series
# Also requires the AUPRCs calculated for the same gene sets on the shuffled networks in a pandas DataFrame
def calculate_network_performance_gain(actual_net_AUPRCs, shuff_net_AUPRCs, verbose=True, save_path=None):
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


if __name__=='__main__':
    if False:
        ranks = {'a':1, 'b':4, 'c': 9}
        precision = [1, 1, 0.5, 1/3]
        true_precisions = [1, 1, 1/2, 1/3, 2/4, 2/5, 2/6, 2/7, 2/8, 3/9, 3/10]
        precision_at_k(precision, ranks, 3)
    if False:
        network_path = '/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/dip.PC_net.txt'
        node_sets_file = '/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/disgen_all_genesets_aim2.txt'
        outdir = '/cellar/users/snwright/Data/Network_Analysis/Evaluation/dev'
        network = dit.load_edgelist_to_networkx(network_path, verbose=False, keep_attributes=False)
        network_size = len(network.nodes())
        alpha = prop.calculate_alpha(network)
        genesets = dit.load_node_sets(node_sets_file, verbose=False, id_type="Entrez")
        # get just the first 10 genesets
        genesets = {k:genesets[k] for k in list(genesets.keys())[:10]}
        sample_p = calculate_p(network, genesets)
        na = NetworkAnalyzer(network, genesets, outdir, 'dip_test', alpha=alpha, sample_proportion=sample_p, 
                            num_samples=5, null_iterations=5, k=5)
        na.perform_analysis()
        na.calculate_performance_metrics()