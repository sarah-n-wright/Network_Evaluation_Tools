import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import obonet as obo
from collections import defaultdict
import requests
from requests.adapters import HTTPAdapter, Retry
import re
from datetime import datetime


def get_node_database_counts(file_list, id_type="Entrez"):
    node_dict = defaultdict(int)
    for file in file_list:
        with open(file) as f:
            for line in f:
                if id_type == "Entrez":
                    node_dict[int(line.strip())] += 1 # add as integer
                else:
                    node_dict[line.strip()] += 1 # add as string
    return node_dict

def get_uniprot_annotation_data(fields, outdir, page_size=500, max_retries=5, taxid=9606, verbose=False):
    if taxid != 9606:
        raise(NotImplementedError("Only human data is currently supported"))
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=max_retries, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return_fields = ",".join(fields)
    url = 'https://rest.uniprot.org/uniprotkb/search?fields='+return_fields+'&format=tsv&query=%28%28organism_id%3A9606%29%29&size=' + str(page_size)
    results = {}
    for batch, total in get_batch(url, re_next_link, session, verbose=verbose):
        for line in batch.text.splitlines()[1:]:
            data = line.split("\t")
            if len(data[0]) > 0 :
                geneid = data[0].split(";")
                if geneid[1] == '':
                    results[int(geneid[0])] = {'aa_length':int(data[1]), 'mass':int(data[2])}
    results_df = pd.DataFrame.from_dict(results, orient='index')
    
    file_name = 'uniprot_data_'+"_".join(fields) + "_" + datetime.today().strftime('%Y-%m-%d')
    if outdir[-1] != '/':
        outdir += '/'
    results_df.to_csv(outdir+file_name+".tsv", sep="\t")
    print('Results saved to: '+outdir+file_name+".tsv")

def get_next_link(headers, re_next_link):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url, re_next_link, session, verbose=False):
    count = 0
    while batch_url:
        if verbose:
            print("Batch: {}".format(count))
        count += 1
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers, re_next_link)

