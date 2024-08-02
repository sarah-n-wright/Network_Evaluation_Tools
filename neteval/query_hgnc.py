import pandas as pd
import httplib2 as http
import json
import mygene

from urllib.parse import urlparse
http.RETIRES=10

def query_mygene(gene_list, scopes, fields, retries=10):
    mg = mygene.MyGeneInfo()
    for retry in range(retries):
        try:
            results_df = mg.querymany(qterms=gene_list, scopes=scopes, fields=fields, species='human', 
                                returnall=True, verbose=False, as_dataframe=True, entrezonly=True)
            break
        except Exception as e:
            if retry < retries - 1:
                #print("QH", gene_list)
                print(f"Retrying mg.querymany: {e}")
            else:
                print("Max retries reach for mg.querymany")
                raise e
    mapped = results_df["out"]
    dups = results_df["dup"]
    missing = results_df["missing"]
    unmapped = []
    if len(dups) > 0:
        unmapped += list(results_df["dup"]["query"].values)
    if len(missing) > 0:
        unmapped += list(results_df["missing"]["query"].values)
    return mapped, unmapped


def search_approved_symbols(ids):
    headers = {'Accept': 'application/json'}
    uri = 'https://rest.genenames.org'
    path = '/search/symbol/*+AND+status:Approved'
    target = urlparse(uri+path)
    method = 'GET'
    body = ''
    h = http.Http()
    print("Checking approved symbols")
    response, content = h.request(target.geturl(), method, body, headers)
    if response['status'] == '200':
        print("Response received")
# assume that content is a json reply
# parse content with the json module 
        data = json.loads(content)
        approved_df = pd.DataFrame.from_dict(data['response']['docs'])
        approved_ids = approved_df.loc[approved_df.symbol.isin(ids), "symbol"]
        approved_map = {sym:sym for sym in approved_ids}
        missing = set(ids).difference(set(approved_ids))
    else:
        print('Error detected: ' + response['status'])
    return approved_map, missing, approved_df


def query_previous_symbols(ids, approved_df=pd.DataFrame()):
    headers = {'Accept': 'application/json'}
    uri = 'https://rest.genenames.org'
    previous_map = {}
    print("Checking previous symbols")
    #print("Number of ids to check", len(ids))
    #print(ids)
    #raise TimeoutError
    for i, symbol in enumerate(ids):
        #print(i, "Previous")
        path = '/search/prev_symbol/' + symbol
        target = urlparse(uri+path)
        method = 'GET'
        body = ''
        h = http.Http()
        response, content = h.request(target.geturl(), method, body, headers)
        if response['status'] == '200':
# assume that content is a json reply
# parse content with the json module 
            data = json.loads(content)
            for entry in data['response']['docs']:
                if entry['symbol'] in approved_df.symbol.values:
                    previous_map[symbol] = entry['symbol']
        else:
            print('Error detected: ' + response['status'], symbol)
    missing = set(ids).difference(set(previous_map.keys()))
    return previous_map, missing


def query_alias_symbols(ids, approved_df=pd.DataFrame()):
    headers = {'Accept': 'application/json'}
    uri = 'https://rest.genenames.org'
    alias_map = {}
    print("Searching aliases")
    for symbol in ids:
        path = '/search/alias_symbol/' + symbol
        target = urlparse(uri+path)
        method = 'GET'
        body = ''
        h = http.Http()
        response, content = h.request(target.geturl(), method, body, headers)
        if response['status'] == '200':
# assume that content is a json reply
# parse content with the json module 
            data = json.loads(content)
            for entry in data['response']['docs']:
                if entry['symbol'] in approved_df.symbol.values:
                    alias_map[symbol] = entry['symbol']
        else:
            print('Error detected: ' + response['status'], symbol)
    missing = set(ids).difference(set(alias_map.keys()))
    return alias_map, missing


def query_other_id(ids, target_id):
    headers = {'Accept': 'application/json'}
    uri = 'https://rest.genenames.org'
    if target_id == "Entrez":
        field = 'entrez_id'
    target_map = {}
    print("Searching", target_id)
    for symbol in ids:
        path = '/fetch/symbol/' + symbol
        target = urlparse(uri+path)
        method = 'GET'
        body = ''
        h = http.Http()
        response, content = h.request(target.geturl(), method, body, headers)
        if response['status'] == '200':
# assume that content is a json reply
# parse content with the json module 
            data = json.loads(content)
            for entry in data['response']['docs']:
                if 'status' in entry:
                    if entry['status'] == "Approved":
                        if field in entry.keys():
                            target_map[symbol] = entry[field]
        else:
            print('Error detected: ' + response['status'], symbol)
    missing = set(ids).difference(set(target_map.keys()))
    return target_map, missing


def search_gene_names(ids, approved_df=pd.DataFrame()):
    try:
        name_df, _ = query_mygene(ids, scopes="name,other_names", fields='symbol')
    except KeyError as e:
        print(ids)
        raise e
    if 'symbol' in name_df.columns:
        name_df = name_df.dropna(subset=['symbol'])
        name_df = name_df.sort_values(by='_score', ascending=False)
        name_df = name_df[~name_df.index.duplicated(keep='first')]
        name_map = name_df['symbol'].to_dict()
        missing = [g for g in ids if g not in name_map.keys()]
    else:
        name_map = {}
        missing = ids
    return name_map, missing


def perform_hgnc_query(ids, from_id, to_id):
    if (from_id == "Symbol") and (to_id == "Symbol"):
        print("Initial Ids", len(ids))
        approved_map, missing, approved_df = search_approved_symbols(ids)
        if len(missing) > 0:
            print("Check names", len(missing))
            name_map, missing = search_gene_names(missing, approved_df)
            if len(missing) > 0:        
                print("Previous Ids", len(missing))
                previous_map, missing = query_previous_symbols(missing, approved_df)
            
                if len(missing) > 0:
                    print("Alias Ids", len(missing))
                    alias_map, missing = query_alias_symbols(missing, approved_df)
                else:
                    alias_map = {}
            else:
                previous_map = {}
                alias_map = {}
        else:
            name_map = {}
            previous_map = {}
            alias_map = {}
        id_map = {**approved_map, **alias_map, **previous_map, **name_map}
        return id_map, missing
    else:
        # use my gene info to retrieve Entrez ids
        
        # then any missing query HGNC on a case by case basis
        
        raise(NotImplementedError, "Only symbol updating supported")
