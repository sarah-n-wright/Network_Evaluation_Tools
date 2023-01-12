import data_import_tools as dit
from processing_functions import update_nodes, convert_node_ids
import pandas as pd
import csv
from tqdm import tqdm
import requests
from getpass import getpass
import argparse

def convert_genesets(genesets, initial_id, target_id):
        all_nodes = set()
        for set_id in genesets:
            all_nodes = all_nodes.union(genesets[set_id])
        updated_nodes, unmapped1 = update_nodes(all_nodes, id_type=initial_id)
        converted_nodes, unmapped2 = convert_node_ids(updated_nodes.values(),initial_id=initial_id, target_id=target_id)
        node_map = {g: converted_nodes[updated_nodes[g]] for g in updated_nodes if updated_nodes[g] in converted_nodes}
        new_gene_sets = {}
        for set_id in genesets:
            new_gene_sets[set_id] = {node_map[g] for g in genesets[set_id] if g in node_map} 
        print("Unmapped genes:", len(unmapped1 + list(unmapped2)), "/", len(all_nodes))
        return new_gene_sets 
    
def write_gene_sets(genesets, outfile, sep="\t"):
    out_strings = [sep.join([set_id] + list(genesets[set_id])) for set_id in genesets]
    with open(outfile,'w') as tfile:
	    tfile.write('\n'.join(out_strings))
     
def get_disgenet_sets(diseasefile, email, password, outfile):
    # first get a list of diseases that have between 20 and 300 genes
    disease_stats = pd.read_csv(diseasefile, sep="\t")
    disease_stats = disease_stats[disease_stats.diseaseType == 'disease']
    keep_diseases = disease_stats[(disease_stats.NofGenes >= 20) & (disease_stats.NofGenes <= 300)]
    disease_ids = keep_diseases.diseaseId.values
    disease_gene_dfs = []
    for disease in tqdm(disease_ids):
        disease_gene_dfs.append(query_disgenet_disease(email, password, disease))
    gda_df = pd.concat(disease_gene_dfs)
    gda_df.to_csv(outfile, sep="\t")

    
def query_disgenet_disease(email, password, disease_code):
    #For this example we are going to use the python default http library

    #Build a dict with the following format, change the value of the two keys your DisGeNET account credentials, if you don't have an account you can create one here https://www.disgenet.org/signup/ 

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
        gda_response = s.get(api_host+'/gda/disease/'+disease_code, params={'source':'ALL'})
        try:
            result = pd.DataFrame.from_records(gda_response.json())
            result['diseaseId'] = disease_code
            return result
        except ValueError:
            print("ERROR:", gda_response.json())

    if s:
        s.close()


if __name__=="__main__":
    #email = getpass("User Email:")
    parser = argparse.ArgumentParser(description='GetDisGenetSets')
    parser.add_argument('--email', metavar='d', required=True)
    parser.add_argument('--password', metavar="nodeA", required=True)
    args = parser.parse_args()    
    #password = getpass("User Password:")
    get_disgenet_sets('/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/disease_associations.tsv.gz', args.email, args.password, 
                        outfile='/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/updated_gda_2023.tsv')
    #genesets = dit.load_node_sets('Git/Network_Evaluation_Tools/Data/DisGeNET_genesets.txt')
    #new_genesets = convert_genesests(genesets, initial_id='Symbol', target_id='Entrez')
    #write_gene_sets(new_genesets, 'Git/Network_Evaluation_Tools/Data/DisGeNET_genesets_entrez.txt')