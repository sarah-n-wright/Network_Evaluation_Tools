import pandas as pd
import argparse
import os


def process_corum_file(corum_file, outdir):
    corum_df = pd.read_csv(corum_file, sep='\t', usecols=['ComplexID', 'ComplexName', 'subunits(Entrez IDs)', 'GO ID', 'GO description'])
    corum_df['NumSubunits'] = corum_df['subunits(Entrez IDs)'].apply(lambda x: len(x.split(';')))
    complex_info = corum_df.loc[:,['ComplexID', 'ComplexName', 'GO ID', 'GO description', 'NumSubunits']]
    complex_info.to_csv(os.path.join(outdir,'corum_complex_info.tsv'), sep='\t')
    all_complexes = []
    for index, row in corum_df.iterrows():
        subunits = row['subunits(Entrez IDs)'].split(';')
        subunits = [int(x) for x in subunits]
        all_complexes.append(pd.DataFrame({'Complex': row['ComplexID'], 'Subunit': subunits}))
    complex_df = pd.concat(all_complexes).reset_index(drop=True)
    complex_df.pivot(index='Subunit', columns='Complex', values='Complex')
    complex_df.to_csv(os.path.join(outdir,'corum_processed.tsv'), sep='\t')
    return corum_df


if __name__=="__main__":
    # argparser = argparse.ArgumentParser(description='Calculate the similarity between two networks')
    # argparser.add_argument('--corum-file', help='Path to the file containing the CORUM network')
    # argparser.add_argument('--netdir', help='Path to the directory containing the networks')    
    # argparser.add_argument('--prefix-file', help='Path to the file containing the prefix of the networks')
    # argparser.add_argument('--outdir', help='Path to the output directory')
    
    # args = argparser.parse_args()
    
    corum_file='/cellar/users/snwright/Data/Network_Analysis/Reference_Data/humanComplexes.txt'
    process_corum_file(corum_file, '')
    
    # if processed corum file os.path.join(outdir,'corum_processed.tsv') exists
    if os.path.exists(os.path.join(args.outdir,'corum_processed.tsv')):
        corum_df = pd.read_csv(os.path.join(args.outdir,'corum_processed.tsv'), sep='\t')
    else:
        corum_df = process_corum_file(args.corum_file, args.outdir)
    
    