import requests, sys
import pandas as pd
import json

def get_latest_ensembl_id(ids):
    server = "https://rest.ensembl.org"
    ext = "/archive/id"
    ids = list(ids)
    ids = [x.split('.')[0] for x in ids]
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    # Initialize an empty list to store the results
    results_df_list = []
    # Split the list of IDs into batches of 1000
    for i in range(0, len(ids), 1000):
        print("Query batch", i, "-", min(i+1000, len(ids)))
        batch_ids = ids[i:min(i+1000, len(ids))]
        # Send the API request with the current batch of IDs with retry logic
        retries = 15
        session = requests.Session()
        session.mount('https://', requests.adapters.HTTPAdapter(max_retries=retries))
        r = session.post(server+ext, headers=headers, data='{ "id" :' + json.dumps(list(batch_ids)) +'}')
        #r = requests.post(server+ext, headers=headers, data='{ "id" :' + json.dumps(list(batch_ids)) +'}')
        if not r.ok:
            #print(r.request)
            #print("\n".join(batch_ids))
            r.raise_for_status()
            sys.exit()
        # Process the API response
        try:
            decoded = r.json()
            decoded_df = pd.DataFrame.from_dict(decoded)
            decoded_df["to"] = decoded_df.apply(parse_archive_results, axis=1)
        except:
            print("Error decoding response for batch", i, "-", min(i+1000, len(ids)))
            print(decoded_df.head())
            print(ids)
            raise Exception("Error generating/decoding response for batch")
        batch_results_df = decoded_df.loc[:, ("id", "to")]
        batch_results_df.columns = ["from", "to"]
        results_df_list.append(batch_results_df)
    # Concatenate all the results into a single DataFrame
    results_df = pd.concat(results_df_list)
    # Find the IDs that were missing from the API response
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
    raise NotImplementedError
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
    raise NotImplementedError

def ensembl_to_entrez():
    raise NotImplementedError
