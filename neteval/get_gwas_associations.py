import pandas as pd
from datetime import datetime
import requests
import shutil
import argparse

## TODO for some reason this is not workigng. Might have to separate this kind of stuff out from the main package. 

def download_file(url, out_dir):
    """
    Download a file from a URL and save it to the specified location.

    :param url: str, The URL of the file to download.
    :param out_dir: str, The directory to save the file to.
    """
    # Send a HTTP request to the URL
    current_date = datetime.now().strftime("%b%d_%Y")
    if out_dir[-1] != '/':
        out_dir = out_dir + '/'
    out_path = out_dir + f"gwas_catalog_{current_date}.txt"
    with requests.get(url, stream=True) as response:
        # Raise an exception if the request was unsuccessful
        response.raise_for_status()

        # Open the file at the save location and write the content of the response to it
        with open(out_path, 'wb') as out_file:
            # Write file in chunks in case it's large
            shutil.copyfileobj(response.raw, out_file)

    print(f"File downloaded from {url} to {out_path}")
    return out_path


def clean_gwas_catalog_data(datafile, outfile, pval_th=5e-8, include_intergenic=False):
    cols= ['DATE', 'PUBMEDID', 'DISEASE/TRAIT', 'MAPPED_GENE', 'SNP_GENE_IDS', 'P-VALUE', 'MAPPED_TRAIT', 'MAPPED_TRAIT_URI', 'INTERGENIC']
    if include_intergenic:
        cols = cols + ['UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID', 'UPSTREAM_GENE_DISTANCE', 'DOWNSTREAM_GENE_DISTANCE']
    data = pd.read_csv(datafile, sep="\t", usecols=cols)
    # filter on pval
    data = data[data["P-VALUE"] <= pval_th]
    # filter on gene and trait present
    data = data.dropna(subset=['SNP_GENE_IDS', "MAPPED_TRAIT_URI"])
    # filter out intergenic  
    if not include_intergenic:
        data = data[data["INTERGENIC"] == 0]
    # remove associations with multiple genes
    data = data[~data["SNP_GENE_IDS"].str.contains(",")]
    # remove associations with multiple traits
    data = data[~data["MAPPED_TRAIT_URI"].str.contains(",")]
    # create trait code
    data['TRAIT_CODE'] = data['MAPPED_TRAIT_URI'].apply(lambda x: x.split('/')[-1])
    # write the cleaned file
    data.to_csv(outfile, sep="\t", index=False)


def create_gwas_gene_sets(datafile, outfile, min_genes= 5, max_genes=500, include_intergenic=False, split_date=None):
    data = pd.read_csv(datafile, sep="\t")
    if split_date is not None:
        data['DATE'] = pd.to_datetime(data['DATE'])
    trait_counts =  data['TRAIT_CODE'].value_counts().to_dict()
    keep_traits = [d for d in trait_counts if ((trait_counts[d] > min_genes) & (trait_counts[d] < max_genes))]
    gwas_sets = {}
    for trait in keep_traits:
        if split_date is not None:
            test_set = process_gwas_genes(data[data['DATE'] >= split_date], trait, intergenic=include_intergenic)
            print(test_set)
            if len(test_set) > min_genes:
                gwas_sets[trait] = test_set
                train_set = process_gwas_genes(data[data['DATE'] < split_date], trait, intergenic=include_intergenic)
                gwas_sets[trait+split_date] = train_set
        else:
            trait_set = process_gwas_genes(data, trait, intergenic=include_intergenic)
            if len(trait_set) > min_genes:
                gwas_sets[trait] = trait_set
    if split_date is not None:
        write_gene_sets(gwas_sets, outfile+split_date)
    else:
        write_gene_sets(gwas_sets, outfile)
    
def process_gwas_genes(data, trait, intergenic=False):
    coding_data = data[((data['TRAIT_CODE']==trait) & (data['INTERGENIC']==0))]
    if len(coding_data) > 0:
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


def write_gene_sets(genesets, outfile, sep="\t"):
    out_strings = [sep.join([set_id] + [str(g) for g in genesets[set_id]]) for set_id in genesets]
    with open(outfile,'w') as tfile:
        tfile.write('\n'.join(out_strings))

if __name__=='__main__':
    # add argument parsing
    parser = argparse.ArgumentParser(description='Get and clean latest GWAS Catalog data')
    parser.add_argument('-o', metavar='out_directory', required=True, type=str)
    parser.add_argument('-u', action='store_true')
    parser.add_argument('-m', metavar='min_genes', default=5, type=int)
    parser.add_argument('-M', metavar='max_genes', default=500, type=int)
    parser.add_argument('-D', metavar='split_date', required=False, type=str, default=None)
    parser.add_argument('-G', metavar='gwas_file', required=False, type=str, default=None)
    parser.add_argument('-p', metavar='pval_th', required=False, type=float, default=5e-8)
    
    args = parser.parse_args()
    update_gwas = args.u
    if args.G is None:
        update_gwas = True
    outdir = args.o
    if update_gwas:
        url = 'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative'  # replace with your file URL  # replace with your desired file path
        gwas_file = download_file(url, outdir)
    else:
        gwas_file = args.G
    clean_gwas_catalog_data(gwas_file, gwas_file + '.cleaned', pval_th=args.p, include_intergenic=False)
    create_gwas_gene_sets(gwas_file + '.cleaned', gwas_file+'.genesets', min_genes=args.m, 
                            max_genes=args.M, include_intergenic=False, split_date=args.D)