import requests, sys
import pandas as pd
import json

def get_latest_ensembl_id(ids):
    server = "https://rest.ensembl.org"
    ext = "/archive/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    #r = requests.post(server+ext, headers=headers, data='{ "id" : ["ENSG00000157764", "ENSG00000248378"] }')
    r = requests.post(server+ext, headers=headers, data='{ "id" :' + json.dumps(list(ids)) +'}')
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    # if node is missing nothing is returned
    # When an entry is withdrawn it will get current=False/None. Can keep 
    decoded_df = pd.DataFrame.from_dict(decoded)
    decoded_df["to"] = decoded_df.apply(parse_archive_results, axis=1)
    print(decoded_df)
    results_df = decoded_df.loc[:, ("id", "to")]
    results_df.columns = ["from", "to"]
    missing = [node for node in ids if node not in results_df['from'].values]
    return results_df, missing


def parse_archive_results(entry):
    if entry.is_current == 1:
        return entry.id
    else:
        new_id = entry.latest
        new_id = new_id.split(".")[0]
        return new_id


def ensembl_to_other():
    ensemble_dbs = {'Symbol': "HGNC", 'Entrez': "EntrezGene"}
    server = "https://rest.ensembl.org"
    ext = "/xrefs/id/ENST00000288602?all_levels=0;external_db=HGNC"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r. raise_for_status()
        sys.exit()
    decoded = r.json()
    print(repr(decoded))

def ensembl_to_uniprot():
    pass

def ensembl_to_entrez():
    pass