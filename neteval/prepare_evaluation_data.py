import neteval.data_import_export_tools as dit
from neteval.gene_mapper import update_nodes, convert_node_ids
import pandas as pd
import csv
from tqdm import tqdm
import requests
from getpass import getpass
import argparse
import re


def check_min_genes_per_network(gene_set, network_file_list, datadir, min_genes=20):
    with open(network_file_list, 'r') as f:
        prefs = f.readlines()
    node_files = [datadir + p.strip()+".nodelist" for p in prefs]
    
    for net_file in node_files:
        net_genes = set(pd.read_csv(net_file).Unique_Nodes.astype(str).values)
        intersect = net_genes.intersection(gene_set)
        if len(intersect) < min_genes:
            print("Failed on", net_file, len(intersect))
            return False  # not all networks have at least min genes in the gene set
    return True


def check_all_genesets_against_network(gene_sets, network_file, min_genes=20):
    net_genes = pd.read_csv(network_file)
    net_genes['Unique_Nodes'] = net_genes['Unique_Nodes'].astype(str)
    remove_sets = []
    for set_id in gene_sets:
        result = check_single_geneset(gene_sets[set_id], net_genes, min_genes)
        if result is False:
            remove_sets.append(set_id)
    for set_id in remove_sets:
        gene_sets.pop(set_id)
    print("Gene sets removed by", network_file, ":", len(remove_sets))
    return gene_sets
    
def check_single_geneset(set_nodes, net_genes, min_genes=20):
    intersection_df = net_genes[net_genes.Unique_Nodes.isin(set_nodes)]
    return len(intersection_df) >= min_genes
    
        

def filter_gene_sets(genesets, network_file_list, datadir, min_genes=20):
    with open(network_file_list, 'r') as f:
        prefs = f.readlines()
    node_files = [datadir + p.strip()+".nodelist" for p in prefs]
    keep_gene_sets = {set_id:genesets[set_id] for set_id in genesets}
    for network in node_files:
        keep_gene_sets = check_all_genesets_against_network(keep_gene_sets, network, min_genes)
    print(len(keep_gene_sets), "/", len(genesets), "retained after filtering.")
    return keep_gene_sets


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


def get_disgenet_sets(diseasefile, email, password, outfile, source, types=["disease"], min_genes=20, max_genes=300):
    # first get a list of diseases that have between 20 and 300 genes
    disease_stats = pd.read_csv(diseasefile, sep="\t")
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
        disease_gene_dfs.append(query_disgenet_disease(email, password, disease, source))
    gda_df = pd.concat(disease_gene_dfs)
    gda_df.to_csv(outfile, sep="\t")

    
def query_disgenet_disease(email, password, disease_code, source):
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
        gda_response = s.get(api_host+'/gda/disease/'+disease_code, params={'source':source})
        try:
            result = pd.DataFrame.from_records(gda_response.json())
            result['diseaseId'] = disease_code
            return result
        except ValueError:
            print("ERROR:", gda_response.json())
    if s:
        s.close()

def process_disgenet_data(datafile, outfile, sep="\t", id_type="geneid", min_genes=20):
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
    
    
def process_gwas_catalog_data(datafile, outfile, pval_th=5e-8, include_intergenic=False, min_genes=20, max_genes=300):
    cols = ['DISEASE/TRAIT','SNP_GENE_IDS', 'INTERGENIC', 'P-VALUE', 'MAPPED_TRAIT']
    if include_intergenic:
        cols = cols + ['UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID', 'UPSTREAM_GENE_DISTANCE', 'DOWNSTREAM_GENE_DISTANCE']
    data = pd.read_csv(datafile, sep="\t", usecols=cols, nrows=1000)
    data = data[data["P-VALUE"] <= pval_th]
    data = data.dropna(subset=['SNP_GENE_IDS'])
    if not include_intergenic:
        data = data[data["INTERGENIC"] == 0]
    trait_counts =  data['MAPPED_TRAIT'].value_counts().to_dict()
    keep_traits = [d for d in trait_counts if ((trait_counts[d] > min_genes) & (trait_counts[d] < max_genes))]
    gwas_sets = {}
    for trait in keep_traits:
        trait_set = process_gwas_genes(data, trait, intergenic=include_intergenic)
        if len(trait_set) > min_genes:
            gwas_sets[trait] = trait_set
    
def process_gwas_genes(data, trait, intergenic=False):
    coding_data = data[((data['MAPPED_TRAIT']==trait) & (data['INTERGENIC']==0))]
    if len(coding_genes) > 0:
        coding_genes = coding_data['SNP_GENE_IDS'].apply(lambda x: clean_gwas_gene_id(x))
        coding_genes = [item for sublist in coding_genes for item in sublist]
        coding_genes = set(coding_genes)
    else:
        coding_genes = set()
    if intergenic:
        intergenic_genes = process_gwas_intergenic(data, trait)
    else:
        intergenic_genes = set()
    return coding_genes.union(intergenic_genes)
        
        
def process_gwas_intergenic(data, trait):
    trait_data = data[((data['MAPPED_TRAIT']==trait) & (data['INTERGENIC']==1))]
    #TODO are there entries with only downstream or upstream only??
    trait_data.dropna(subset=["UPSTREAM_GENE_ID","UPSTREAM_GENE_DISTANCE",'DOWNSTREAM_GENE_DISTANCE','DOWNSTREAM_GENE_ID'])
    if len(trait_data) > 0:
        genes = trait_data.apply(lambda x: clean_gwas_gene_id(x["UPSTREAM_GENE_ID"]) if 
                                    x["UPSTREAM_GENE_DISTANCE"] < x['DOWNSTREAM_GENE_DISTANCE'] else 
                                    clean_gwas_gene_id(x['DOWNSTREAM_GENE_ID']), axis=1)
        genes = [item for sublist in genes for item in sublist]
        genes = set(genes)
    else:
        return set()
    

def clean_gwas_gene_id(geneid):
    #TODO do I want to keep multiples?
    return geneid.split(', ')



if __name__=="__main__":
    #email = getpass("User Email:")
    #parser = argparse.ArgumentParser(description='GetDisGenetSets')
    #parser.add_argument('--email', metavar='d', required=True)
    #parser.add_argument('--password', metavar="nodeA", required=True)
    #args = parser.parse_args()    
    #password = getpass("User Password:")
    #process_disgenet_data("/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/updated_gda_2023_BEFREE.tsv",
    #                        outfile='/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/DisGeNET_genesets_BEFREE_2023.txt')
    #process_disgenet_data("/cellar/users/snwright/Data/Transfer/Disgenet_GWASDB_gda.tsv",
    #                        outfile='/cellar/users/snwright/Data/Transfer/Disgenet_GWASDB_genesets.tsv', id_type="gene_symbol")
    get_disgenet_sets('/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/HPO_diseases.tsv', "snwright@ucsd.edu", "Ballon44!", 
                    outfile='/cellar/users/snwright/Data/Transfer/Disgenet_HPO_gda.tsv', source="HPO",
                    types=["disease", "phenotype"], min_genes=10, max_genes=300)
    process_disgenet_data("/cellar/users/snwright/Data/Transfer/Disgenet_HPO_gda.tsv",
                            outfile='/cellar/users/snwright/Data/Transfer/Disgenet_HPO_genesets.tsv', id_type="gene_symbol", min_genes=10)
    #get_intersecting_genes('Git/Network_Evaluation_Tools/Data/v1_net_prefixes.txt', 
                            #datadir='/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/')
    # process_gwas_catalog_data('/cellar/users/snwright/Data/Network_Analysis/Reference_Data/gwas_cat_Jan13_2023.txt', "", pval_th=5e-8, include_intergenic=False, min_genes=20, max_genes=300)
    
    # datadir='/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/'
    # gitdir='/cellar/users/snwright/Git/Network_Evaluation_Tools/'
    # genesets = dit.load_node_sets(gitdir+ 'Data/DisGeNET_genesets_BEFREE_2023.txt')
    # filtered_genesets = filter_gene_sets(genesets, gitdir+'Data/v1_net_prefixes.txt', min_genes=20,datadir=datadir)
    # write_gene_sets(filtered_genesets, gitdir+'Data/DisGeNET_genesets_BEFREE_2023_v1_filtered_genesets_20.txt')
    
    # genesets = dit.load_node_sets(gitdir+ 'Data/DisGeNET_genesets_BEFREE_2023.txt')
    # filtered_genesets = filter_gene_sets(genesets, gitdir+'Data/v2_net_prefixes.txt', min_genes=20,datadir=datadir)
    # write_gene_sets(filtered_genesets, gitdir+'Data/DisGeNET_genesets_BEFREE_2023_v2_filtered_genesets_20.txt')
    
    # genesets = dit.load_node_sets(gitdir+ 'Data/DisGeNET_genesets_ALL_2023.txt')
    # filtered_genesets = filter_gene_sets(genesets, gitdir+'Data/v1_net_prefixes.txt', min_genes=20,datadir=datadir)
    # write_gene_sets(filtered_genesets, gitdir+'Data/DisGeNET_genesets_ALL_2023_v1_filtered_genesets_20.txt')
    
    # genesets = dit.load_node_sets(gitdir+ 'Data/DisGeNET_genesets_ALL_2023.txt')
    # filtered_genesets = filter_gene_sets(genesets, gitdir+'Data/v2_net_prefixes.txt', min_genes=20,datadir=datadir)
    # write_gene_sets(filtered_genesets, gitdir+'Data/DisGeNET_genesets_ALL_2023_v2_filtered_genesets_20.txt')
    # new_genesets = convert_genesets(genesets, initial_id='Symbol', target_id='Entrez')
    # write_gene_sets(new_genesets, 'Git/Network_Evaluation_Tools/Data/Oncogeneic_Components_genesets_entrez.txt')
