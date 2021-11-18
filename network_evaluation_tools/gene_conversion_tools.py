################################################################
# ---------- Network Gene Name Conversion Functions ---------- #
################################################################
import requests
import re
import time
import pandas as pd
import numpy as np
from tqdm import tqdm
import json

# Determine if id to be input is a valid gene name (does not contain parentheses or quotations or whitespace)
def exclude_id(name, bad_prefixes=None):
	excluded_id_regex = re.compile('[(),\'\"\s\/\|\.<>]+')
	# Remove genes that may also have prefixes that we do not want (e.g. CHEBI)
	if bad_prefixes:
		for prefix in bad_prefixes:
			if name.startswith(prefix):
				return True
	return excluded_id_regex.search(name)

# Remove the naming system prefix, if there is one
def get_identifier_without_prefix(string):
	elements = string.split(':')
	length = len(elements)
	if length == 2:
		return str(elements[1])
	elif length > 2:
		return None
	else:
		return string

# Construct string for bach query to MyGene.Info v3.0.0 API
def query_constructor(gene_list, exclude_prefixes=None, print_invalid_genes=False):
	# Find genes that are valid and return only gene identifiers
	valid_query_genes = [get_identifier_without_prefix(gene) for gene in gene_list if exclude_id(gene, exclude_prefixes)==None]
	# Find all genes that have invalid names
	invalid_query_genes = [gene for gene in gene_list if exclude_id(gene, exclude_prefixes)!=None]
	print(len(valid_query_genes), "Valid Query Genes")
	if print_invalid_genes:
		print(len(invalid_query_genes), "Invalid Query Genes:")
		print(invalid_query_genes)
	else:
		print(len(invalid_query_genes), "Invalid Query Genes")
	query_string = ' '.join(valid_query_genes) # Build string of names to input into MyGene.Info
	return query_string, valid_query_genes, invalid_query_genes

# Function for posting batch query to MyGene.info v3 API (2021)
def query_batch(query_string, tax_id='9606', scopes="symbol,entrezgene,alias,uniprot", fields="symbol,entrezgene"):
    query_split = query_string.split(' ')
    query_n = len(query_split)
    query_time = time.time()
    if query_n <=1000:
        q = ','.join(query_split)
        headers = {'content-type': 'application/x-www-form-urlencoded'}
        req = 'q='+q+"&scopes="+scopes+'&fields='+fields+"&taxid="+tax_id
        res = requests.post('http://mygene.info/v3/query', data=req, headers=headers)
        results = json.loads(res.text)
    else:
        # If the query is too long, we will need to break it up into chunks of 1000 query genes (MyGene.info cap)
        if query_n % 1000 == 0:
            chunks = round(query_n / 1000)
        else:
            chunks = round(np.floor(query_n / 1000) + 1)
        query_chunks = []
        for i in range(chunks):
            start_i, end_i = i*1000, (i+1)*1000
            query_chunks.append(','.join(query_split[start_i:end_i]))
        results = []
        for chunk in query_chunks:
            headers = {'content-type': 'application/x-www-form-urlencoded'}
            req = 'q='+chunk+"&scopes="+scopes+'&fields='+fields+"&taxid="+tax_id
            res = requests.post('http://mygene.info/v3/query', data=req, headers=headers)
            results = results + json.loads(res.text)         
    print(len(results), 'Matched query results')
    print('Batch query complete:', round(time.time()-query_time,2), 'seconds')
    return results

# Construct matched queries maps
def construct_query_map_table(query_result, query_genes, display_unmatched_queries=False):
	if type(query_genes) == str:
		query_genes = query_genes.split(' ')
	construction_time = time.time()
	# Construct DataFrame of matched queries (only keep the results for each query where both symbol and entrez id were mapped)
	matched_data, matched_genes=[], []
	for match in query_result:
		if match.get('entrezgene') and match.get('symbol'):
			matched_data.append([match.get('query'), match.get('_score'), match.get('symbol'), str(match.get('entrezgene'))])
			matched_genes.append(match.get('query'))
	# Add all other partial mappings or non-mappings to the list
	partial_match_genes = [gene for gene in query_genes if gene not in matched_genes]
	partial_match_results = []
	for match in query_result:
		if match.get('query') in partial_match_genes:
			partial_match_results.append(match)
			if match.get('entrezgene'): # If there if an entrez gene, we want that that in string form, otherwise we want None
				matched_data.append([match.get('query'), match.get('_score'), match.get('symbol'), str(match.get('entrezgene'))])
			else:
				matched_data.append([match.get('query'), match.get('_score'), match.get('symbol'), match.get('entrezgene')])
	print('Queries without full matching results found:', len(partial_match_results))
	if display_unmatched_queries:
		for entry in partial_match_results:
			print(entry)
	# Convert matched data list into data frame table
	match_table = pd.DataFrame(data=matched_data, columns=['Query','Score','Symbol','EntrezID'])
	match_table = match_table.set_index('Query')
	# Some genes will be matched in duplicates (due to alias mapping, generally the highest scoring matches will be correct)
	# Therefore we remove duplicate mappings to create 1-to-1 mappings for query to genes.
	duplicate_matched_genes = []
	for gene in query_genes:
		if type(match_table.loc[gene])==pd.DataFrame:
			duplicate_matched_genes.append(gene)
	print()
	print(len(duplicate_matched_genes), "Queries with mutliple matches found")
	# Construct mapping table of genes with only one full result
	single_match_genes = [gene for gene in query_genes if gene not in duplicate_matched_genes]

	match_table_single = match_table.loc[single_match_genes]
	# Keep matches of queries matched only once if there are duplicate matches for genes
	if len(duplicate_matched_genes) > 0:
		# Keep maximum scored matches of queries matched more than once
		max_score_matches=[]
		for gene in duplicate_matched_genes:
			matched_duplicates = match_table.loc[gene]
			max_score = max(matched_duplicates['Score'])
			max_score_matches.append(matched_duplicates[matched_duplicates['Score']==max_score])
		match_table_duplicate_max = pd.concat(max_score_matches)
		# Construct Query maps for symbol and entrez
		match_table_trim = pd.concat([match_table_single, match_table_duplicate_max])
	else:
		match_table_trim = match_table_single.copy(deep=True)
	# Construct query map dictionaries
	query_to_symbol = match_table_trim['Symbol'].to_dict()
	query_to_entrez = match_table_trim['EntrezID'].to_dict()
	print()
	print('Query mapping table/dictionary construction complete:', round(time.time()-construction_time,2), 'seconds')
	return match_table_trim, query_to_symbol, query_to_entrez

# Filter edgelist to remove all genes that contain invalid query names
# This function is only required if there are any invalid genes found by query_constructor()
def filter_query_edgelist(query_edgelist, invalid_genes):
	edgelist_filt = []
	count=0
	for edge in query_edgelist:
		if edge[0] in invalid_genes or edge[1] in invalid_genes:
			count+=1
		else:
			edgelist_filt.append(edge)
	print(count, '/', len(query_edgelist), 'edges with invalid nodes removed')
	return edgelist_filt

# Convert network edge lists
# Third column is for weights if desired to pass weights forward
def convert_edgelist(query_edgelist, query_to_symbol, weighted=False, fill_na=False):
    edge_df = pd.DataFrame(query_edgelist, columns=["n1", "n2", "weight"])
    gene_map = pd.DataFrame.from_dict(query_to_symbol, orient='index')
    gene_map.columns=['symbol']
    edge_df = edge_df.merge(gene_map, left_on="n1", right_index=True, how="left", suffixes=("", "_n1"))
    edge_df = edge_df.merge(gene_map, left_on="n2", right_index=True, how="left", suffixes=("", "_n2"))
    edge_df.columns = ["n1", "n2", "weight", "symbol_n1", "symbol_n2"]
    if fill_na:
        edge_df.loc[edge_df.symbol_n1.isna(), "symbol_n1"] = edge_df.loc[edge_df.symbol_n1.isna(), ("n1")]
        edge_df.loc[edge_df.symbol_n1.isna(), "symbol_n2"] = edge_df.loc[edge_df.symbol_n1.isna(), ("n2")]
    if weighted:
        edge_df = edge_df.loc[:, ("symbol_n1", "symbol_n2", "weight")]
    else:
        edge_df = edge_df.loc[:, ("symbol_n1", "symbol_n2")]
    return edge_df

# Sometimes each node needs to be converted by its best match if there are multiple names per node
# This function uses the match_table constructed earlier to convert genes to either symbol or entrez format only
def convert_custom_namelist(names, field, match_table):
	# Keep only mappings defined for field of interest
	if field=='symbol':
		# Return match table values that have matched symbol
		conversion = match_table.ix[names][~(match_table.ix[names]['Symbol'].isnull())]
		if conversion.shape[0]==0:
			return None
		else:
			# Return conversion with max score or None if no conversion
			max_score = conversion['Score'].max()
			return conversion[conversion['Score']==max_score].ix[0]['Symbol']
	elif field=='entrez':
		# Return match table values that have matched symbol
		conversion = match_table.ix[names][~(match_table.ix[names]['EntrezID'].isnull())]
		if conversion.shape[0]==0:
			return None
		else:
			# Return conversion with max score or None if no conversion
			max_score = conversion['Score'].max()
			return conversion[conversion['Score']==max_score].ix[0]['EntrezID']

# Filter converted edge lists
def filter_converted_edgelist(edgelist, remove_self_edges=True, weighted=False):
	filter_time = time.time()
	print(len(edgelist),'input edges')
	# Remove self-edges
	if remove_self_edges:
		edgelist_filt1 = edgelist.loc[(edgelist.symbol_n1 != edgelist.symbol_n2)]
		print(len(edgelist)-len(edgelist_filt1), 'self-edges removed')
	else:
		edgelist_filt1 = edgelist
		print('Self-edges not removed')
	if weighted:
		# Remove edges where one or both nodes are "None"
		edgelist_filt2 = edgelist_filt1.dropna()
		print(len(edgelist_filt1)-len(edgelist_filt2), 'edges with un-mapped genes removed')
		# Remove duplicates by keeping the max score
		edgelist_filt3_scoremap = {}
		edgelist_filt2 = edgelist_filt2.sort_values(by="weight", ascending=False)
		edgelist_filt3 = edgelist_filt2.drop_duplicates(subset=['symbol_n1', 'symbol_n2'], keep='first')
		# Convert dictionary of scores to list
		print(len(edgelist_filt2)-len(edgelist_filt3), 'duplicate edges removed')
	else:
		# Remove edges where one or both nodes are "None"
		edgelist_filt2 = edgelist_filt1.dropna()
		print(len(edgelist_filt1)-len(edgelist_filt2), 'edges with un-mapped genes removed')
		# Remove duplicate edges
		edgelist_filt3 = edgelist_filt2.drop_duplicates()
		print(edgelist_filt2.shape[0]-len(edgelist_filt3), 'duplicate edges removed')
	print('Edge list filtered:',round(time.time()-filter_time,2),'seconds')
	print(len(edgelist_filt3), 'Edges remaining')
	return edgelist_filt3

# Write edgelist to file
def write_edgelist(edgelist, output_file, delimiter='\t', binary=True):
    write_time=time.time()
    if type(edgelist) == list:
        f = open(output_file,'w')
        for edge in edgelist:
            if binary:
                f.write(delimiter.join([edge[0], edge[1]])+'\n')
            else:
                f.write(delimiter.join([str(val) for val in edge])+'\n')
        f.close()
    else:
        if binary:
            edgelist.loc[:, ("symbol_n1", "symbol_n2")].to_csv(output_file, sep=delimiter, index=False, header=False)
        else:
            edgelist.to_csv(output_file, sep=delimiter, index=False, header=False)
    print('Edge list saved:', round(time.time()-write_time,2),'seconds')