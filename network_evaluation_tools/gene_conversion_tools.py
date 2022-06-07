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
import urllib.parse
import urllib.request

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
        for chunk in tqdm(query_chunks):
            headers = {'content-type': 'application/x-www-form-urlencoded'}
            req = 'q='+chunk+"&scopes="+scopes+'&fields='+fields+"&taxid="+tax_id
            res = requests.post('http://mygene.info/v3/query', data=req, headers=headers)
            results = results + json.loads(res.text)         
    print(len(results), 'Matched query results')
    print('Batch query complete:', round(time.time()-query_time,2), 'seconds')
    return results


def query_uniprot(ids2map, source_fmt='ACC+ID', target_fmt="GENENAME", output_fmt='tab', return_as_dict=True):
    """
    See https://www.uniprot.org/help/api_idmapping for data codes
    """
    url = 'https://www.uniprot.org/uploadlists/'
    if hasattr(ids2map, 'pop'):
        ids2map = ' '.join(ids2map)
    if not hasattr(target_fmt, 'pop'):
        target_fmt = [target_fmt]
    for i, target in enumerate(target_fmt):
        params = {'from': source_fmt,
                  'to': target,
                  'format': output_fmt,
                  'query': ids2map,
                 'taxon': '9606'}

        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
            response = f.read()
        temp_file = open("uniprot_temp.txt", "w")
        temp_file.write(response.decode('utf-8'))
        temp_file.close()
        response_data = pd.read_csv("uniprot_temp.txt", sep="\t", dtype=str)
        response_data.columns = ["query", target]
        if i == 0:
            output = response_data
        else:
            output = output.merge(response_data, on="query", how="outer")
        column_map = {"P_ENTREZGENEID": 'entrezgene', "GENENAME":'symbol'}
        output.rename(axis="columns", mapper=column_map, inplace=True)
        output = output.assign(_score=0)
    unmapped = [gene for gene in ids2map.split(' ') if gene not in output["query"].values]
    if return_as_dict:
        output = output.to_dict(orient="records")
    return output, unmapped


def query_against_dataset(data_path, queries, query_priority=["Approved symbol", "Previous symbols", "Alias symbols"], 
                          to_id=["Approved symbol", "NCBI Gene ID"], split_cols = ["Previous symbols", "Alias symbols"],
                          filter_kws={"Status":["Approved"]}, return_as_dict=False):
    clean = lambda x: "NaN" if x=="NaN" or x=='' else str(int(x))
    queries = set(queries)
    ref_data = pd.read_csv(data_path, sep="\t", converters={"NCBI Gene ID": clean})
    # perform desired filtering
    if filter_kws is not None:
        for filter_col in filter_kws:
            ref_data = ref_data.loc[ref_data[filter_col].isin(filter_kws[filter_col])]
    # remove entries that do not have the desired outputs ------
    #ref_data.dropna(subset = to_id, inplace=True)
    ref_data = ref_data.loc[:, set(query_priority+to_id)]
    results = []
    for search_col in query_priority:
        if search_col in split_cols:
            search_data = ref_data.dropna(subset=[search_col])
            search_data = search_data.assign(temp=[x.split(",") for x in search_data[search_col]])
            search_data = search_data.explode(column="temp")
            search_data["temp"] = search_data["temp"].apply(lambda x: x.strip())
            search_data.drop(columns=[search_col], inplace=True)
            search_data.rename(columns={"temp":search_col}, inplace=True)
            search_results = search_data.loc[search_data[search_col].isin(queries), set([search_col]+to_id)]
            search_results = search_results.assign(Query=search_results[search_col])
            results.append(search_results)
            queries = queries.difference(set(search_results["Query"].values))
        else:
            search_results = ref_data.loc[ref_data[search_col].isin(queries), set([search_col]+to_id)]
            search_results = search_results.assign(Query=search_results[search_col])
            results.append(search_results)
            queries = queries.difference(set(search_results["Query"].values))
        if len(queries) == 0:
            break    
    output = pd.concat(results)
    output = output.loc[:, ["Query"] + to_id]
    if return_as_dict:
        output.columns = ["query", "symbol", "entrezgene"]
        output["entrezgene"] = output["entrezgene"]
        output["_score"] = 0
        output = output.to_dict(orient="records")
        for missing_gene in queries:
            output.append({"query":missing_gene, "symbol": "NaN", "entrezgene": "NaN", "_score": 0})
    return output, queries


def construct_query_map_table(query_result, query_genes, display_unmatched_queries=False):
    if not isinstance(query_result, pd.DataFrame):
        query_result = pd.DataFrame.from_dict(query_result)
    query_result = query_result.replace("NaN", np.nan)
    query_result = query_result.fillna(value=np.nan)
    if type(query_genes) == str:
        query_genes = query_genes.split(' ')
    # Duplicated        
    duplicated_data = query_result[query_result.duplicated(subset="query")]
    print("Number of genes with multiple matches:", len(duplicated_data))
    if len(duplicated_data) > 0:
        unique_data = query_result.sort_values(by=["_score", "entrezgene", "symbol"])
        unique_data.drop_duplicates(subset="query", keep="first", inplace=True)
    else:
        unique_data = query_result
    # UNMATCHED
    unmatched_genes = query_result[((query_result["symbol"].isna()) & (query_result["entrezgene"].isna()))]["query"].unique().tolist()
    print("Number of unmatched genes:", len(unmatched_genes))
    if display_unmatched_queries:
        for gene in unmatched_genes:
            print(gene)

    # FULLY MATCHED
    print("Number of fully matched genes:", len(unique_data.dropna(subset=["symbol", "entrezgene"])))
    ## No entrezID
    print("Number of partially matched genes:", len(unique_data[(unique_data.symbol.isna()) | (unique_data.entrezgene.isna())])) 
    match_table_trim = unique_data.loc[:, ("query", "symbol", "_score", "entrezgene")]
    match_table_trim.columns = ['Query','Symbol',"Score", 'EntrezID']
    match_table_trim = match_table_trim.set_index('Query')
    query_to_symbol = match_table_trim.dropna(subset=["Symbol"])['Symbol'].to_dict()
    query_to_entrez = match_table_trim.dropna(subset=["EntrezID"])['EntrezID'].to_dict()
    return match_table_trim, query_to_symbol, query_to_entrez


# Construct matched queries maps
def construct_query_map_table_old(query_result, query_genes, display_unmatched_queries=False):
	if type(query_genes) == str:
		query_genes = query_genes.split(' ')
	construction_time = time.time()
	# Construct DataFrame of matched queries (only keep the results for each query where both symbol and entrez id were mapped)
	matched_data, matched_genes=[], []
	for match in query_result:
		if match.get('entrezgene') and match.get('entrezgene') != "NaN" and match.get('symbol') and match.get('symbol') != "NaN":
			matched_data.append([match.get('query'), match.get('_score'), match.get('symbol'), str(match.get('entrezgene'))])
			matched_genes.append(match.get('query'))
	# Add all other partial mappings or non-mappings to the list
	partial_match_genes = [gene for gene in query_genes if gene not in matched_genes]
	partial_match_results = []
	for match in query_result:
		if match.get('query') in partial_match_genes:
			partial_match_results.append(match)
			if match.get('entrezgene') and match.get('entrezgene') != "NaN": # If there if an entrez gene, we want that that in string form, otherwise we want None
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
	query_to_symbol = match_table_trim.loc[match_table_trim.Symbol != "NaN"]['Symbol'].to_dict()
	query_to_entrez = match_table_trim.loc[match_table_trim.EntrezID != "NaN"]['EntrezID'].to_dict()
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
def convert_edgelist(query_edgelist, query_to_symbol, node_cols=[0,1]):
    if not isinstance(query_edgelist, pd.DataFrame):
        edge_df = pd.DataFrame(query_edgelist)
    else:
        edge_df = query_edgelist.copy(deep=True)
    node_names = [edge_df.columns[i] for i in node_cols]
    edge_df.dropna(subset=node_names, inplace=True)
    edge_df[node_names[0]] = edge_df[node_names[0]].astype(str)
    edge_df[node_names[1]] = edge_df[node_names[1]].astype(str) 
    edge_df[node_names[0]] = edge_df[node_names[0]].map(query_to_symbol, na_action="ignore")
    edge_df[node_names[1]] = edge_df[node_names[1]].map(query_to_symbol, na_action="ignore")
    name_map = {node_names[0]:"symbol_n1", node_names[1]:"symbol_n2"}
    edge_df.rename(columns=name_map, inplace=True)
    return edge_df
    

# Third column is for weights if desired to pass weights forward
def convert_edgelist_old(query_edgelist, query_to_symbol, weighted=False, fill_na=False):
    if not isinstance(query_edgelist, pd.DataFrame):
        edge_df = pd.DataFrame(query_edgelist)
    else:
        edge_df = query_edgelist
    if edge_df.shape[1] == 3:
        edge_df.columns=["n1", "n2", "weight"]
        final_cols = ["n1", "n2", "weight", "symbol_n1", "symbol_n2"]
    else:
        edge_df.columns=["n1", "n2"]
        final_cols = ["n1", "n2", "symbol_n1", "symbol_n2"]
    edge_df.drop_na(subset=["n1", "n2"], inplace=True) 
    edge_df["n1"] = edge_df["n1"].astype(str)
    edge_df["n2"] = edge_df["n2"].astype(str)    
    gene_map = pd.DataFrame.from_dict(query_to_symbol, orient='index')
    gene_map.columns=['symbol']
    edge_df = edge_df.merge(gene_map, left_on="n1", right_index=True, how="left", suffixes=("", "_n1"))
    edge_df = edge_df.merge(gene_map, left_on="n2", right_index=True, how="left", suffixes=("", "_n2"))
    edge_df.columns = final_cols
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
def filter_converted_edgelist_old(edgelist, remove_self_edges=True, weighted=False, node_cols=['symbol_n1', 'symbol_n2']):
    filter_time = time.time()
    print(len(edgelist),'input edges')
    # Remove self-edges
    if remove_self_edges:
        edgelist_filt1 = edgelist.loc[(edgelist[node_cols[0]] != edgelist[node_cols[1]])]
        print(len(edgelist)-len(edgelist_filt1), 'self-edges removed')
    else:
        edgelist_filt1 = edgelist
        print('Self-edges not removed')
    if weighted:
        # Remove edges where one or both nodes are "None"
        edgelist_filt2 = edgelist_filt1.dropna(subset=node_cols)
        print(len(edgelist_filt1)-len(edgelist_filt2), 'edges with un-mapped genes removed')
        # Remove duplicates by keeping the max score
        # Create a sorted pair identifier to account for node 1-node 2 ordering
        # edgelist_filt3_scoremap = {}
        edgelist_filt2 = edgelist_filt2.assign(symbol_pair=["--".join(sorted(pair)) for pair in zip(edgelist_filt2[node_cols[0]], edgelist_filt2[node_cols[1]])])
        edgelist_filt2 = edgelist_filt2.sort_values(by=['symbol_pair', 'weight'], ascending=False)
        edgelist_filt3 = edgelist_filt2.drop_duplicates(subset=['symbol_pair'], keep='first')
        print(len(edgelist_filt2)-len(edgelist_filt3), 'duplicate edges removed')
    else:
        # Remove edges where one or both nodes are "None"
        edgelist_filt2 = edgelist_filt1.dropna(subset=node_cols)
        print(len(edgelist_filt1)-len(edgelist_filt2), 'edges with un-mapped genes removed')
        # Remove duplicate edges
        edgelist_filt2 = edgelist_filt2.assign(symbol_pair=["--".join(sorted(pair)) for pair in zip(edgelist_filt2[node_cols[0]], edgelist_filt2[node_cols[1]])])
        edgelist_filt3 = edgelist_filt2.drop_duplicates(subset=['symbol_pair'])
        print(edgelist_filt2.shape[0]-len(edgelist_filt3), 'duplicate edges removed')
    print('Edge list filtered:',round(time.time()-filter_time,2),'seconds')
    print(len(edgelist_filt3), 'Edges remaining')
    return edgelist_filt3


def filter_converted_edgelist(edgelist, remove_self_edges=True, node_cols=['symbol_n1', 'symbol_n2'],
                             weight_col=None):
    filter_time = time.time()
    node_indices = [edgelist.columns.get_loc(col) for col in node_cols]
    print(len(edgelist),'input edges')
    # Remove self-edges
    if remove_self_edges:
        edgelist_filt1 = edgelist.loc[(edgelist[node_cols[0]] != edgelist[node_cols[1]])]
        print(len(edgelist)-len(edgelist_filt1), 'self-edges removed')
    else:
        edgelist_filt1 = edgelist
        print('Self-edges not removed')
    # remove edges with unmapped genes
    edgelist_filt2 = edgelist_filt1.dropna(subset=node_cols)
    print(len(edgelist_filt1)-len(edgelist_filt2), 'edges with un-mapped genes removed')
    edgelist_filt2_sorted = sort_node_pairs(edgelist_filt2, node_indices)
    if weight_col is not None:
        edgelist_filt2_sorted = edgelist_filt2_sorted.sort_values(by=weight_col, ascending=False)
    edgelist_filt3 = edgelist_filt2_sorted.drop_duplicates(subset=node_cols, keep="first")
    print(len(edgelist_filt2_sorted)-len(edgelist_filt3), 'duplicate edges removed')    
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
    
    
def sort_node_pairs(data, node_col_indices=[0,1]):
    node_names = [data.columns[i] for i in node_col_indices]
    node_array = data.iloc[:, node_col_indices].to_numpy()
    node_array.sort()
    sorted_data = pd.DataFrame(node_array)
    sorted_data.columns = node_names
    if data.shape[1] > 2:
        sorted_data = pd.concat([sorted_data.reset_index(drop=True), data.drop(columns=node_names).reset_index(drop=True)], axis=1)
    return sorted_data
