import argparse
from neteval.processing_functions import *
from neteval.gene_mapper import query_mygene
import pandas as pd
import csv

def parse_species(cols, value):
    if cols == "None":
        species = None
        species_code = None
        return species, species_code
    if ";" in cols:
        cols = cols.split(";")
        for i, col in enumerate(cols):
            if col.isnumeric():
                cols[i] = int(col)
    else:
        if cols.isnumeric():
            cols = int(cols)
                
    if value.isnumeric():
        species_code = int(value)
    else:
        species_code = value
    return cols, species_code

def parse_prefix(pref_sep):
    results = []
    if "[" in pref_sep:
        pairs = pref_sep.split("],")
        for pair in pairs:
            prefix = pair.split("[")[0]
            sep = pair.split("[")[1].split("]")[0]   
            results.append((prefix, sep))
        return results
    elif "," in pref_sep:
        prefs = pref_sep.split(",")
        return prefs
    else:
        return [pref_sep]
        
        



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Process Network Data.')
    parser.add_argument('datafile', metavar='d', help='String with the full file path storing the data')
    parser.add_argument('-A', metavar="nodeA", required=True)
    parser.add_argument('-B', metavar="nodeB", required=True)
    parser.add_argument('--species', default="None")
    parser.add_argument('--species_value', default="None")
    parser.add_argument('--score', default="None")
    parser.add_argument('-i', metavar="input_id", required=True, nargs='*')
    parser.add_argument('-t', metavar="target_id", required=True)
    parser.add_argument('-N', metavar='networkName', required=True)
    parser.add_argument('-o', metavar='outpath', required=True)
    parser.add_argument('--header', default='0')
    parser.add_argument('--testMode', default='0')
    parser.add_argument('--sep', default="\t")
    parser.add_argument('--prefix', default="None")
    args = parser.parse_args()    
    #print(args.A, args.B)
    #print(args.B.isnumeric())
    node_a = int(args.A) if args.A.isnumeric() else args.A
    node_b = int(args.B) if args.B.isnumeric() else None if args.B == "None" else args.B
    #print(type(node_a), type(node_b))

    species, species_code = parse_species(args.species, args.species_value)
    
    score = args.score if (args.score != "None") else None
    score = int(score) if ((score is not None) and (score.isnumeric())) else score
    header = int(args.header) if (args.header != "None") else None
    run_in_test_mode = False if args.testMode == '0' else True
    prefixes = parse_prefix(args.prefix) if (args.prefix != "None") else None
    
    if args.sep == "tab":
        sep = "\t"
    elif args.sep == "space":
        sep = "\s+"
    else:
        sep = args.sep
    
    #print(node_a, node_b, species, score, header)
    #print(type(header))
    print(args)
    if True:
        nd = NetworkData(args.datafile, node_a=node_a, node_b=node_b, species=species, target_id_type=args.t,  identifiers=args.i,
                    score=score, species_code=species_code,header=header, net_name=args.N, test_mode=run_in_test_mode, sep=sep,
                    prefixes=prefixes)
        nd.clean_data()
        original_nodes = nd.get_unique_nodes()
        mygene_fields = {"Symbol": "symbol", "Entrez": 'entrezgene', "Uniprot": "uniprot", "Ensembl": "ensembl.gene",
                    "Refseq":"refseq", "EnsemblProtein":"ensembl.protein"}
        scopes = ",".join([mygene_fields[x] for x in args.i])
        fields = mygene_fields[args.t]
        mapped, unmapped = query_mygene(original_nodes, scopes, fields, retries=10)
        save_file = "/cellar/users/snwright/Data/Network_Analysis/mapping_test.txt"
        with open(save_file, "a") as f:
            f.write("\t".join([args.N, str(len(original_nodes)), str(len(unmapped))])+"\n")

    

    if False:
        nd = NetworkData(args.datafile, node_a=node_a, node_b=node_b, species=species, target_id_type=args.t,  identifiers=args.i,
                    score=score, species_code=species_code,header=header, net_name=args.N, test_mode=run_in_test_mode, sep=sep,
                    prefixes=prefixes)
        nd.clean_data()
        nd.convert_nodes()
        try:
            final_nodes = nd.get_unique_nodes()
        except:
            pass
        nd.remove_duplicates()
        nd.write_network_data(args.o)
        nd.write_stats(args.o)
        print("Processing of", args.N, "completed.")
    
