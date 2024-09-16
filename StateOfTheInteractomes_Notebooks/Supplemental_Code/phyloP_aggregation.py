import pandas as pd
import os
import numpy as np
from neteval.query_ensembl import *
from neteval.gene_mapper import *

def load_summary_data(datadir, f):
    data = pd.read_csv(os.path.join(datadir, f), sep='\t', header=None)
    data.columns = ['CHR', 'start', 'stop', 'id', 'x', 'strand', 'mean', 'median', 'stdev', 'mad', 'uniq_bases']
    return data
    
def extract_all_stat(datadir, file_list):
    stats = {"mean":{}, "bases":{}, "median":{}}
    for f in file_list:
        data = load_summary_data(datadir, f)
        data = data[data['uniq_bases'] > 0]
        gene_id = data['id'].unique()
        if len(gene_id) > 1:
            print(gene_id)
        else:
            try:
                # all bases should be treated equally so we calculated weighted values from the statistics of each transcript
                stats["mean"][gene_id[0]] = np.average(data['mean'].astype(float).values, weights=data['uniq_bases'].astype(int).values)
                stats["bases"][gene_id[0]] = np.sum(data['uniq_bases'].astype(int).values)
                stats["median"][gene_id[0]] = weighted_median(data['median'].astype(float).values, data['uniq_bases'].astype(int).values)
            except:
                if len(data) > 0:
                    #print(data.head())
                    print(f)
    return stats

def weighted_median(medians, weights):
    a = []
    for i, med in enumerate(medians):
        a += [med] * weights[i]
    return np.median(a)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Aggregate phyloP across all CDS transcripts of a gene')
    parser.add_argument("datadir", type=str, help='Directory containing the summary.bed files')
    parser.add_argument("filelist", type=str, help='Filename in datadir to the list of all summary.bed files')
    args = parser.parse_args()
    with open(os.path.join(args.datadir, args.filelist), 'r') as f:
        file_list = [x.strip()+'.summary.bed' for x in f.readlines()]
    cons = extract_all_stat(datadir, file_list)
    cons["transcript"] = cons.index
    cons['transcript'] = cons.transcript.apply(lambda x: x.split(':')[-1].split('.')[0])
    converted_df, still_missing = query_mygene(cons.transcript.values, scopes='ensembl.transcript', fields='entrezgene')
    cons_mapped = cons.merge(converted_df.reset_index().dropna(subset=['_id']).loc[:, ('query', 'entrezgene')], left_on='transcript', right_on='query', how='inner')
    cons_mapped.to_csv(os.path.join(datadir, 'all_summarized_scores.txt'), sep='\t')