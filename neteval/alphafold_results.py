#!/usr/bin/env python
import glob
import os
import json
import pandas as pd

def extract_scores(pair, datadir):
    ptms = []
    iptms = []
    try:
        model_files = glob.glob(f'{datadir}/{pair[0]}_{pair[1]}_scores_rank_*.json')
        for mod in model_files:
            data = json.load(open(mod))
            ptms.append(data['ptm'])
            iptms.append(data['iptm'])
    except FileNotFoundError:
        print(model_files)
        #print(f'Error extracting scores for {pair}, not found')
        return {'ptm':[], 'iptm':[]}
    return {'ptm':ptms, 'iptm':iptms}

def extract_scores_for_pair_file(pair_file, datadir):
    with open(pair_file) as f:
        pairs = f.readlines()
    pairs = [pair.strip() for pair in pairs]
    scores = {}
    for pair in pairs:
        split_pair = pair.split(',')
        scores[pair] = extract_scores(split_pair, datadir)
    return scores

def format_results_to_df(score_dict):
    pair_results = []
    for pair in score_dict:
        df = pd.DataFrame(score_dict[pair])
        df['Model'] = df.index
        df['Pair_str'] = pair
        pair_results.append(df)
    if len(pair_results) > 0:
        pair_df = pd.concat(pair_results)
    pair_df['Confidence'] = 0.8 * pair_df['iptm'] + 0.2 * pair_df['ptm']
    pair_df['NCBI_ID_A'] = pair_df['Pair_str'].apply(lambda x: sorted([int(y) for y in x.split('_')])[0])
    pair_df['NCBI_ID_B'] = pair_df['Pair_str'].apply(lambda x: sorted([int(y) for y in x.split('_')])[1])
    return pair_df

def write_results(pair_df, outdir,pair_file):
    pair_df.to_csv(os.path.join(outdir, os.path.basename(pair_file) +'_alphafold_results.tsv'), sep='\t')    
    print("Results written to: ", os.path.join(outdir, os.path.basename(pair_file) +'_alphafold_results.tsv'))
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Extract AlphaFold scores from JSON files')
    parser.add_argument('pair_file', type=str, help='File containing pairs to extract scores for')
    parser.add_argument('datadir', type=str, help='Directory containing JSON files')
    parser.add_argument('outdir', type=str, help='Directory to write results to')
    args = parser.parse_args()
    score_dict = extract_scores_for_pair_file(args.pair_file, args.datadir)
    pair_df = format_results_to_df(score_dict)
    write_results(pair_df, args.outdir, args.pair_file)