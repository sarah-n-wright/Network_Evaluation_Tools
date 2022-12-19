import pandas as pd
import httplib2 as http
import json

from urllib.parse import urlparse
http.RETIRES=10

def search_approved_symbols(ids):
    headers = {'Accept': 'application/json'}
    uri = 'http://rest.genenames.org'
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
    return approved_map, missing



def query_previous_symbols(ids):
    headers = {'Accept': 'application/json'}
    uri = 'http://rest.genenames.org'
    previous_map = {}
    print("Checking previous symbols")
    print("Number of ids to check", len(ids))
    print(ids)
    #raise TimeoutError
    for symbol in ids:
        print("Previous")
        path = '/fetch/prev_symbol/' + symbol
        target = urlparse(uri+path)
        method = 'GET'
        body = ''
        h = http.Http()
        response, content = h.request(target.geturl(), method, body, headers)
        if response['status'] == '200':
            print("Response received.")
# assume that content is a json reply
# parse content with the json module 
            data = json.loads(content)
            for entry in data['response']['docs']:
                if entry['status'] == "Approved":
                    previous_map[symbol] = entry['symbol']
        else:
            print('Error detected: ' + response['status'])
    missing = set(ids).difference(set(previous_map.keys()))
    return previous_map, missing


def query_alias_symbols(ids):
    headers = {'Accept': 'application/json'}
    uri = 'http://rest.genenames.org'
    alias_map = {}
    print("Searching aliases")
    for symbol in ids:
        path = '/fetch/alias_symbol/' + symbol
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
                if entry['status'] == "Approved":
                    alias_map[symbol] = entry['symbol']
        else:
            print('Error detected: ' + response['status'])
    missing = set(ids).difference(set(alias_map.keys()))
    return alias_map, missing


def query_other_id(ids, target_id):
    headers = {'Accept': 'application/json'}
    uri = 'http://rest.genenames.org'
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
                if entry['status'] == "Approved":
                    target_map[symbol] = entry[field]
        else:
            print('Error detected: ' + response['status'])
    missing = set(ids).difference(set(target_map.keys()))
    return target_map, missing


def perform_hgnc_query(ids, from_id, to_id):
    if (from_id == "Symbol") and (to_id == "Symbol"):
        approved_map, missing = search_approved_symbols(ids)
        previous_map, missing = query_previous_symbols(missing)
        alias_map, missing = query_alias_symbols(missing)
        id_map = {**approved_map, **alias_map, **previous_map}
        return id_map, missing
    else:
        # use my gene info to retrieve Entrez ids
        
        # then any missing query HGNC on a case by case basis
        
        raise(NotImplementedError, "Only symbol updating supported")
