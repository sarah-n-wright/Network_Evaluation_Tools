###############################################
# ---------- Data Import Functions ---------- #
###############################################

import pandas as pd
import networkx as nx
import time
import os
import matplotlib.pyplot as plt

def load_network_from_ndex(uuid):
    raise NotImplementedError


def export_network_to_ndex():
    raise NotImplementedError


def load_network_to_dataframe():
    raise NotImplementedError

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


def write_networkx_to_file(G, outfilepath, timer=None, verbose=False):
    if timer is not None:
        timer.start("Write Network")
    net_df = nx.to_pandas_edgelist(G, source="Node_A", target="Node_B")
    net_df.to_csv(outfilepath, index=False, sep="\t")
    if timer is not None:
        timer.end("Write Network")
    return


# Construct dictionary of node sets from input text file to perform AUPRC analysis on for network of interest
# File format: Each line is a delimited list with the first item in the list is the name of the node set
# All other nodes in the list follow the node set name
def load_node_sets(node_set_file, delimiter='\t', verbose=False, id_type="Symbol"):
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
			node_sets[set_id] = {int(node) for node in list(node_sets[set_id])}
	if verbose:
		print('Node cohorts loaded:', node_set_file)
	return node_sets


######################################
##	OLD FUNCTIONS - TO BE REMOVED 	##
######################################
# Load network from file as unweighted network
# Can set delimiter, but default delimiter is tab
# Only will read edges as first two columns, all other columns will be ignored
def load_network_file_DELETE(network_file_path, delimiter='\t', verbose=False, id_type='Symbol', keep_attributes=False):
	"""_summary_

	Args:
		network_file_path (_type_): _description_
		delimiter (str, optional): _description_. Defaults to '\t'.
		verbose (bool, optional): _description_. Defaults to False.
		id_type (str, optional): _description_. Defaults to 'Symbol'.
		keep_attributes (bool, optional): _description_. Defaults to False.

	Returns:
		_type_: _description_
	"""
	net_df = pd.read_csv(network_file_path, sep=delimiter)
	source_col, target_col = net_df.columns[0:2]
	if id_type == "Entrez":
		net_df[source_col] = net_df[source_col].astype(int)
		net_df[target_col] = net_df[target_col].astype(int)
	if keep_attributes:
		network = nx.from_pandas_edgelist(net_df, source=source_col, target=target_col, edge_attr=keep_attributes)
	else:
		network = nx.from_pandas_edgelist(net_df, source=source_col, target=target_col)
	if verbose:
		print('Network File Loaded:', network_file_path)
	return network


# Filter extended sif file where all edges are weighted by a specific quantile
# Return the filtered network edge list and save it to a file if desired (for import by load_network_file)
def filter_weighted_network_sif_DELETE(network_file_path, nodeA_col=0, nodeB_col=1, score_col=2, q=0.9, delimiter='\t', verbose=False, save_path=None):
    data = pd.read_csv(network_file_path, sep=delimiter, header=None, low_memory=False)
    # Filter edges by score quantile
    q_score = data[score_col].quantile(q)
    if verbose:
        print(str(round(q*100,2))+'%', 'score:', q_score)
    data_filt = data[data[score_col]>q_score][data.columns[[nodeA_col, nodeB_col, score_col]]]
    data_filt.columns = ['nodeA', 'nodeB', 'edgeScore']
    if verbose:
        print(data_filt.shape[0], '/', data.shape[0], 'edges retained')
    if save_path is not None:
        data_filt.to_csv(save_path, sep='\t', header=False, index=False)
    return data_filt


def plot_changes_to_dataset_DELETE(input_raw, input_raw_v2, edgelist_filt, edgelist_filt_v2, input_human=None, input_human_v2=None):
    plt.rcParams['font.size'] = '14'
    stats = pd.DataFrame({"v1":[len(input_raw)],
                        "v2":[len(input_raw_v2)]}, index=["Input"])
    if input_human is not None:
        stats = pd.concat([stats, pd.DataFrame({"v1":len(input_human), "v2":len(input_human_v2)}, index=["Human only"])])
    else:
        print("No human filtering done")
    stats = pd.concat([stats, pd.DataFrame({"v1":len(edgelist_filt), "v2":len(edgelist_filt_v2)}, index=["Filtered"])])
    # get node stats
    nodes_v1 = set(edgelist_filt["symbol_n1"].values).union(set(edgelist_filt["symbol_n2"].values))
    nodes_v2 = set(edgelist_filt_v2["symbol_n1"].values).union(set(edgelist_filt_v2["symbol_n2"].values))
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,5))
    stats.plot.bar(ax=ax1, fontsize=14)
    ax1.set_ylabel("Number of edges", fontsize=14 )
    ax2.bar(["v1", "v2", "change", "+ V2", "- V2"], [len(nodes_v1), len(nodes_v2), len(nodes_v2)-len(nodes_v1),
                                                        len(nodes_v2.difference(nodes_v1)), -1* len(nodes_v1.difference(nodes_v2))])
    ax2.set
    ax2.set_ylabel("Number of nodes")
    

    
# Get full paths to all networks in directory with a given file name structure:
# e.g. If filename = 'BIND_Symbol.sif', then network_name='BIND', suffix='_Symbol', ext='.sif
def get_networks_DELETE(wd, suffix=None, file_ext='.sif'):
	network_files = {}
	for fn in os.listdir(wd):
		if suffix==None:
			if fn.endswith(file_ext):
				network_files[fn.split(file_ext)[0]]=wd+fn			
		else:
			if fn.endswith(file_ext) and fn.split(file_ext)[0].endswith(suffix):
				network_files[fn.split(suffix)[0]]=wd+fn
	return network_files

# Companion function with get_networks(), loads all of the network files found in a directory
# Uses the load_network_file() function to load each network, also only imports first two columns, no edge data
# Constructs a dictionary of useful network items for each network in the directory:
#  - Actual networkx object representation of network
#  - List of nodes by name for each network
#  - List of edges by node name for each network
def load_networks_DELETE(network_file_map, delimiter='\t', verbose=False):
	# Initialize dictionaries
	networks, network_edges, network_nodes = {}, {}, {}
	# Loading network and network properties
	for network_name in network_file_map:
		loadtime = time.time()
		# Load network
		network = load_network_file(network_file_map[network_name], verbose=verbose)
		networks[network_name]=network
		# Construct network node list
		network_nodes[network_name] = network.nodes()
		# Construct network edge list
		network_edges[network_name] = network.edges()
	if verbose:
		print('All given network files loaded')
	# Return data structure
	return networks, network_edges, network_nodes

# Convert and save MAF from Broad Firehose
# Can produce 2 types of filetypes: 'matrix' or 'list', matrix is a full samples-by-genes binary csv, 'list' is a sparse representaiton of 'matrix'
# This is a conversion tool, so the result must be saved (most tools will require a path to a processed MAF file and load it separately)
# Gene naming can be 'Symbol' or 'Entrez'
def process_TCGA_MAF_DELETE(maf_file, save_path, filetype='matrix', gene_naming='Symbol', verbose=False):
	loadtime = time.time()
	# Load MAF File
	TCGA_MAF = pd.read_csv(maf_file,sep='\t',low_memory=False)
	# Get all patient somatic mutation (sm) pairs from MAF file
	if gene_naming=='Entrez':
		TCGA_sm = TCGA_MAF.groupby(['Tumor_Sample_Barcode', 'Entrez_Gene_Id']).size()
	else:
		TCGA_sm = TCGA_MAF.groupby(['Tumor_Sample_Barcode', 'Hugo_Symbol']).size()
	# Turn somatic mutation data into binary matrix
	TCGA_sm_mat = TCGA_sm.unstack().fillna(0)
	TCGA_sm_mat = (TCGA_sm_mat>0).astype(int)
	# Trim TCGA barcodes
	TCGA_sm_mat.index = [pat[:12] for pat in TCGA_sm_mat.index]
	# Filter samples with duplicate IDs
	non_dup_IDs = list(TCGA_sm_mat.index.value_counts().index[TCGA_sm_mat.index.value_counts()==1])
	dup_IDs = list(TCGA_sm_mat.index.value_counts().index[TCGA_sm_mat.index.value_counts()>1])
	# Save file as binary matrix or sparse list
	if filetype=='list':
		# Now try to construct two-column/sparse representation of binary sm data
		# Get list of all patient somatic mutations
		index_list = list(TCGA_sm.index)
		# Filter list of patient somatic mutations of duplicate patient barcodes
		index_list_filt = [i for i in index_list if not any([True if barcode in i[0] else False for barcode in dup_IDs])]
		# Save patient somatic mutations list to file
		f = open(save_path, 'w')
		for sm in index_list_filt:
			f.write(sm[0][:12]+'\t'+sm[1]+'\n')
		f.close()
		if verbose:
			print('Binary somatic mutations list saved')
	else:
		# Save non-duplicate patients' binary TCGA somatic mutation matrix to csv
		TCGA_sm_mat_filt = TCGA_sm_mat.loc[non_dup_IDs]
		# Remove all genes that have no more mutations after patient filtering
		nonempty_cols = [col for col in TCGA_sm_mat_filt.columns if not all(TCGA_sm_mat_filt[col]==0)]
		TCGA_sm_mat_filt2 = TCGA_sm_mat_filt[nonempty_cols]
		# Remove columns with bad names like '0'
		named_cols = [col for col in TCGA_sm_mat_filt.columns if col!='0']
		TCGA_sm_mat_filt3 = TCGA_sm_mat_filt2[nonempty_cols]
		TCGA_sm_mat_filt3.to_csv(save_path)
		if verbose:
			print('Binary somatic mutation matrix saved')
	if verbose:
		print('MAF file processed:', maf_file, round(time.time()-loadtime, 2), 'seconds.')
	return

# Load binary mutation data with 2 file types (filetype= 'matrix' or 'list')
# filetype=='matrix' is a csv or tsv style matrix with row and column headers, rows are samples/patients, columns are genes
# filetype=='list' is a 2 columns text file separated by the delimiter where 1st column is sample/patient, 2nd column is one gene mutated in that patient
# Line example in 'list' file: 'Patient ID','Gene Mutated'
def load_binary_mutation_data(filename, filetype='matrix', delimiter=',', verbose=False):
	if filetype=='list':
		f = open(filename)
		binary_mat_lines = f.read().splitlines()
		binary_mat_data = [(line.split('\t')[0], line.split('\t')[1]) for line in binary_mat_lines]
		binary_mat_index = pd.MultiIndex.from_tuples(binary_mat_data, names=['Tumor_Sample_Barcode', 'Hugo_Symbol'])
		binary_mat_2col = pd.DataFrame(1, index=binary_mat_index, columns=[0])[0]
		binary_mat = binary_mat_2col.unstack().fillna(0)
	else:
		binary_mat = pd.read_csv(filename, delimiter=delimiter, index_col=0).astype(int)
	if verbose:
		print('Binary Mutation Matrix Loaded:', filename)
	return binary_mat

# Concatinate multiple mutation matrices together
# All file type structures and delimiters must be the same (see load_binary_mutation_matrix()) across all files
def concat_binary_mutation_matrices_DELETE(filename_list, filetype='matrix', delimiter=',', verbose=False, save_path=None):
	binary_mat_list = [load_binary_mutation_data(fn, filetype=filetype, delimiter=delimiter, verbose=verbose) for fn in filename_list]  
	binary_mat_concat = pd.concat(binary_mat_list).fillna(0)
	if verbose:
		print('All binary mutation matrices loaded and concatenated')
	if save_path==None:
		return binary_mat_concat
	else:
		binary_mat_concat.to_csv(save_path)
		return binary_mat_concat
