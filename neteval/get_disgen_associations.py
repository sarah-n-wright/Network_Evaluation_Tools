import pandas as pd
from tqdm import tqdm
import requests
import argparse
from datetime import datetime

def query_disgenet_disease(email, password, source, disease_code):
    """Query disgenet for a genes associated with a specific disease code and database source
    
    Args:
        email (str): email for disgenet
        password (str): password for disgenet
        source (str): database source to query
        disease_code (str): disease code to query
        
    Returns:
        pd.DataFrame: dataframe of gene-disease associations
    """
    query = '/gda/disease/'+disease_code
    result = query_disgenet(email, password, source=source, query=query)
    if result is not None:
        result['diseaseId'] = disease_code
        return result
    else:
        print("No results for disease:", disease_code)
        
def get_latest_disgenet_disease_list(email, password, source, min_genes=5, outfile=None):
    """Get the latest DisGeNET disease list for a specific source
    
    Args:
        email (str): email for disgenet
        password (str): password for disgenet
        source (str): database source to query
        min_genes (int): minimum number of genes per gene set
        outfile (str): output file
    
    Returns:
        pd.DataFrame: dataframe of diseases
    """
    assert source in ['CURATED', 'INFERRED', 'ANIMAL_MODELS', 'ALL', 'BEFREE', 'CGI', 'CLINGEN', 'CLINVAR', 'CTD_human', 
                    'CTD_mouse', 'CTD_rat', 'GENOMICS_ENGLAND', 'GWASCAT', 'GWASDB', 'HPO', 'LHGDN', 'MGD', 'ORPHANET', 
                    'PSYGENET', 'RGD', 'UNIPROT']
    query = '/disease/source/'+source
    result = query_disgenet(email, password, query)
    result = result[result.num_genes >= min_genes]
    if outfile is not None:
        result.to_csv(outfile, sep="\t", index=False)
    return result

def query_disgenet(email, password, query, source=None):
    """ Wrapper function for API query with DisGeNET
    
    Args:
        email (str): email for disgenet.org
        password (str): password for disgenet.org
        query (str): query string for API
        source (str): source to restrict query to specific database
        
    Returns:
        pd.DataFrame: dataframe of results
    """
    # credit disgenet website
    auth_params = {"email":email,"password":password}
    api_host = "https://www.disgenet.org/api"
    api_key = None
    s = requests.Session()
    try:
        r = s.post(api_host+'/auth/', data=auth_params)
        if(r.status_code == 200):
            #Lets store the api key in a new variable and use it again in new requests
            json_response = r.json()
            api_key = json_response.get("token")
            #print(api_key + "This is your user API key.") #Comment this line if you don't want your API key to show up in the terminal
        else:
            print(r.status_code)
            print(r.text)
    except requests.exceptions.RequestException as req_ex:
        print(req_ex)
        print("Something went wrong with the request.")

    if api_key:
        #Add the api key to the requests headers of the requests Session object in order to use the restricted endpoints.
        s.headers.update({"Authorization": "Bearer %s" % api_key}) 
        #Lets get all the diseases associated to a gene eg. APP (EntrezID 351) and restricted by a source.
        if source is not None:
            gda_response = s.get(api_host+query, params={'source':source})
        else:
            gda_response = s.get(api_host+query)
        try:
            result = pd.DataFrame.from_records(gda_response.json())
            return result
        except ValueError:
            print("ERROR:", gda_response.json())
    if s:
        s.close()


def get_disgenet_associations(diseasefile, email, password, outfile, source, types=["disease"], min_genes=20, max_genes=300):
    """Given a list of disgenet diseases, download the gene associations
    
    Args:
        diseasefile (str): file path for list of diseases
        email (str): email for disgenet.org
        password (str): password for disgenet.org
        outfile (str): output file for associations
        source (str): database source to query
        types (list): list of disease types to query
        min_genes (int): minimum number of genes per gene set
        max_genes (int): maximum number of genes per gene set
        
    Returns:
        None
    
    """
    # first get a list of diseases that have between 20 and 300 genes
    disease_stats = pd.read_csv(diseasefile, sep="\t")
    # filter by disease type
    if "diseaseType" in disease_stats.columns:
        disease_stats = disease_stats[disease_stats.diseaseType.isin(types)]
        keep_diseases = disease_stats[(disease_stats.NofGenes >= min_genes) & (disease_stats.NofGenes <= max_genes)]
        disease_ids = keep_diseases.diseaseId.values
    else:
        disease_stats = disease_stats[disease_stats["type"].isin(types)]
        keep_diseases = disease_stats[(disease_stats.num_genes >= min_genes) & (disease_stats.num_genes <= max_genes)]
        disease_ids = keep_diseases.diseaseid.values
    disease_gene_dfs = []
    for disease in tqdm(disease_ids):
        disease_gene_dfs.append(query_disgenet_disease(email, password, disease_code=disease, source=source))
    gda_df = pd.concat(disease_gene_dfs)
    gda_df.to_csv(outfile, sep="\t")
    


def create_disgenet_genesets(datafile, outfile, sep="\t", id_type="geneid", min_genes=20):
    """Create gene set file from DisGeNET data
    
    Args:
        datafile (str): file path for DisGeNET association data
        outfile (str): output file for gene sets
        sep (str): separator for data file
        id_type (str): column name for gene ids
        min_genes (int): minimum number of genes per gene set
        
    Returns:
        None
    """
    data = pd.read_csv(datafile, sep=sep, index_col=0)
    disease_counts = data.diseaseId.value_counts()
    keep_ids = list(disease_counts[disease_counts>=min_genes].index)
    data_keep = data[data.diseaseId.isin(keep_ids)]
    disease_sets = {}
    id_map = data_keep.loc[:, ('diseaseId', 'disease_name')].drop_duplicates().set_index("diseaseId")['disease_name'].to_dict()
    for set_id in keep_ids:
        disease_sets[id_map[set_id]] = sep.join(data_keep[data_keep.diseaseId == set_id][id_type].astype(str).values) 
    out_strings = [sep.join([set_id, disease_sets[set_id]]) for set_id in disease_sets]
    with open(outfile,'w') as tfile:
        tfile.write('\n'.join(out_strings))
        

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Get and clean latest DisGeNET data')
    parser.add_argument('-o', metavar='out_directiory', required=True, type=str)
    parser.add_argument('-m', metavar='min_genes', default=5, type=int)
    parser.add_argument('-M', metavar='max_genes', default=500, type=int)
    parser.add_argument('-u', action='store_true')
    parser.add_argument('-E', metavar='email', required=True, type=str)
    parser.add_argument('-P', metavar='password', required=True, type=str)
    parser.add_argument('-S', metavar='source', required=False, type=str, default='BEFREE')
    parser.add_argument('-d', metavar='disgenfile', required=False, type=str, default=None)    
    args = parser.parse_args()
    current_date = datetime.now().strftime("%b%d_%Y")
    f"gwas_catalog_{current_date}.txt"
    if args.d is None:
        args.u = True
    if args.u:
        get_latest_disgenet_disease_list(args.E, args.P, args.S, outfile=args.o + 'disease_list_'+args.S+current_date+'.tsv')
    
    get_disgenet_associations(args.o + 'disease_list_'+args.S+current_date+'.tsv', args.E, args.P, 
                    args.o+'gda_'+args.S+current_date+'.tsv', source=args.S,
                    types=["disease"], min_genes=args.m, max_genes=args.M)
        
    create_disgenet_genesets(args.o+'gda_'+args.S+current_date+'.tsv', args.o + 'disease_list_'+args.S+current_date+'.tsv.genesets', 
                            min_genes=args.m)
    
    # email='snwright@ucsd.edu'
    # password='MyDisGenPW'


