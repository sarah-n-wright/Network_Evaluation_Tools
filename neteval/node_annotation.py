import numpy as np
import pandas as pd
from collections import defaultdict
import requests
from requests.adapters import HTTPAdapter, Retry
import re
from datetime import datetime
from neteval import gene_mapper
import os
from tqdm import tqdm
from scipy.stats import fisher_exact


CWD = os.path.dirname(os.path.abspath(__file__))


def load_nodes(nodefile):
    nodes = pd.read_csv(nodefile, sep='\t')
    nodes.columns = ['Node']
    return nodes

def perform_permutation_test(nodes, annotation_data, n_permutations=10000, metric='CitationCount', stat=np.mean):
    nodes = [x for x in nodes if x in annotation_data.index]
    all_values = annotation_data[metric].values
    observed_values = annotation_data.loc[nodes][metric].values
    observed_stat = stat(observed_values)
    observed_std = np.std(observed_values)
    permuted_stats = np.zeros(n_permutations)
    for i in tqdm(range(n_permutations)):
        permuted_values = np.random.choice(all_values, len(nodes))
        permuted_stats[i] = stat(permuted_values)
    p_value = np.sum(permuted_stats >= observed_stat) / n_permutations
    if p_value == 0:
        p_value = 1/n_permutations
    return observed_stat, stat(permuted_stats), p_value

def permutation_test_wrapper(stats, anno_df, anno_col):
    nets = []
    observed = []
    permuted = []
    p_values = []
    for net in [x for x in stats.network_names]:
        nets.append(net)
        nodes = load_nodes(stats.node_files[net])
        observed_stat, permuted_stat, p_value = perform_permutation_test(nodes['Node'], anno_df, metric=anno_col, stat=np.median)
        observed.append(observed_stat)
        permuted.append(permuted_stat)
        p_values.append(p_value)
    results = pd.DataFrame({'Network':nets, 'Observed': observed, 'Permuted': permuted, 'P-value': p_values})
    results = results.sort_values(by='Observed', ascending=False)
    results['Significant'] = results['P-value'] < 0.05/len(nets)
    return results

def parse_chrm(chr_str):
    chr_map = {'mitochondria':'MT', 'reserved':'other', "unplaced":'other', 'not on reference assembly':'other'}
    if chr_str in chr_map:
        return chr_map[chr_str]
    try:
        match = re.search(r'^(Y|X|\d+(?:\.\d+)?)', chr_str)
        if match:
            return match.group()
        else:
            print(chr_str)
            return 'other'
    except TypeError:
        return 'other'


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


def load_hgnc(datadir):
    
    hgnc = pd.read_csv(os.path.join(datadir, "HGNC_download_Dec20_2023.txt"), sep="\t")
    hgnc = hgnc.dropna(subset=['NCBI Gene ID(supplied by NCBI)'])
    hgnc['NCBI Gene ID(supplied by NCBI)'] = hgnc['NCBI Gene ID(supplied by NCBI)'].astype(int)
    hgnc.rename(columns={'NCBI Gene ID(supplied by NCBI)':"GeneID"}, inplace=True)
    hgnc.Chromosome = hgnc.Chromosome.apply(parse_chrm)
    hgnc.drop(columns=['Approved symbol', 'Status', 'RefSeq IDs', 'Gene group name', 'NCBI Gene ID', 'UniProt ID(supplied by UniProt)'], inplace=True)
    hgnc.set_index('GeneID', inplace=True)
    hgnc.index.name=None
    return hgnc


def load_ensembl(datadir):
    ensem = pd.read_csv(os.path.join(datadir, "Ensembl_export_Dec20_2023.txt.gz"), sep="\t")
    #ensem.drop(columns=["GO domain"], inplace=True)
    ensem.rename(columns={'NCBI gene (formerly Entrezgene) ID':"GeneID"}, inplace=True)
    ensem.drop_duplicates(inplace=True)
    ensem_no_transcript = ensem.drop(columns=['Transcript stable ID', 'Transcript type', 'Transcript length (including UTRs and CDS)']).drop_duplicates()
    ensem_no_transcript["Gene length"] = ensem_no_transcript['Gene end (bp)'] - ensem_no_transcript['Gene start (bp)']
    #ensem_no_transcript["Protein Coding"] = ensem_no_transcript['Gene type'].apply(lambda x: True if x == "protein_coding" else False)
    ensem_no_transcript.dropna(subset=['GeneID'], inplace=True)
    ensem_no_transcript.GeneID = ensem_no_transcript.GeneID.astype(int)
    ensem_no_transcript.drop(columns=['Gene stable ID', 'Gene start (bp)', 'Gene end (bp)', 'HGNC symbol'], inplace=True)
    deduped_gene_stats = ensem_no_transcript.groupby('GeneID').mean()
    deduped_gene_types = ensem_no_transcript.groupby('GeneID')['Gene type'].agg(lambda x: x.mode()[0])
    out_df = pd.concat([deduped_gene_stats, deduped_gene_types], axis=1)
    out_df.index.name=None
    
    return out_df


def load_citations(datadir):
    cite = pd.read_csv(os.path.join(datadir, "gene_citation_counts_Dec20_2023.txt"), sep="\t", header=None)
    cite.columns = ["GeneID", "CitationCount"]
    cite.set_index("GeneID", inplace=True)
    cite.index.name=None
    return cite


def load_uniprot(datadir):
    uni = pd.read_csv(os.path.join(datadir, "uniprot_data_id_length_mass_2023-12-20.tsv"), sep="\t", 
                    names=["GeneID", "id","aa_length", "mass"], header=0)
    uni.set_index("GeneID", inplace=True)
    uni.index.name=None
    return uni


def get_annot_df(datadir):
    hgnc = load_hgnc(datadir)
    ensem = load_ensembl(datadir)
    cite = load_citations(datadir)
    uni = load_uniprot(datadir)
    return pd.concat([hgnc, ensem, cite, uni], axis=1)


def get_ensembl_annotation_data():
    print('Instructions for downloading Ensembl data:')
    print('1: Go to https://www.ensembl.org/biomart/martview')
    print('2: Select Ensembl Genes')
    print('3: Select Human genes')
    print('4: Go to attributes and select the following:')
    print('Gene stable ID, Transcript stable ID, Gene start (bp) Gene end (bp), Transcript length (including UTRs and CDS), Gene GC content, Gene type, Transcript type, GO domain, NCBI gene (formerly Entrezgene) ID, HGNC symbol')
    print('5: Go to results tab and download results as TSV')

def get_ncbi_citation_data():
    print('Instructions for downloading NCBI citation data:')
    print('1: Go to https://www.ncbi.nlm.nih.gov/public/')
    print('2: Select Gene > DATA > gene2pubmed.gz')
    print('Run the command:')
    print('awk -F\'\\t\' \'NR>1 && $1 == 9606 {count[$2]++} END {for (gene in count) print gene "\\t" count[gene]}\' gene2pubmed > text.txt')
    

def get_uniprot_annotation_data(fields, outdir, index_on=2, page_size=500, max_retries=5, taxid=9606, verbose=False):
    if taxid != 9606:
        raise(NotImplementedError("Only human data is currently supported"))
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=max_retries, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return_fields = ",".join(fields)
    url = 'https://rest.uniprot.org/uniprotkb/search?fields='+return_fields+'&format=tsv&query=%28%28organism_id%3A9606%29%29&size=' + str(page_size)
    results = {}
    i = 0
    for batch, total in get_batch(url, re_next_link, session, verbose=verbose):
        i += 1
        print("batch", str(i))
        for line in batch.text.splitlines()[1:]:
            data = line.split("\t")
            if len(data[index_on]) > 0 :
                if ";" in data[index_on]:
                    geneid = data[index_on].split(";")
                    if geneid[1] != '':
                        continue
                    else:
                        geneid = geneid[0]
                else:
                    geneid = data[index_on]
                try:
                    results[int(geneid)] = {fields[i]: data[i] for i in range(len(fields))}
                except:
                    print(data)
                    raise ValueError
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
        

class ExpressionData:
    def __init__(self, filepath, debug=None, force_remapping=False, datatype='mrna', min_obs=0):
        if debug is not None:
            self.data = debug.data
            self.tissues = debug.tissues
            self.gene_idx = debug.gene_idx
            self.global_medians	= debug.global_medians
            self.min_obs = debug.min_obs
            self.datatype = debug.datatype

        else:
            self.filepath = filepath
            self.min_obs = min_obs
            self.datatype = datatype
            self.data = self._load_data(force_remapping=force_remapping)
            self.data = self._consolidate_transcripts()
            self.gene_idx = self._get_gene_idx()
            self.tissues = self._categorize_tissues()
            self.global_medians = self._global_expression_statistics()
    
    def _load_data(self, force_remapping=False):
        if self.datatype == 'mrna':
            return self._load_mrna(force_remapping=force_remapping)
        elif self.datatype == 'protein':       
            return self._load_protein(force_remapping=force_remapping)

    def _load_mrna(self, force_remapping=False):
        raw_data  = pd.read_csv(self.filepath, sep="\t", skiprows=2)
        raw_data.rename(columns={"Name":"Ensembl_ID", "Description":"Symbol"}, inplace=True)
        if "Ensembl_ID" not in raw_data.columns:
            raw_data = pd.read_csv(self.filepath, sep="\t")
            raw_data.rename(columns={"Name":"Ensembl_ID", "Description":"Symbol"}, inplace=True)
        mapped_data = self._convert_ids(raw_data, force_remapping=force_remapping)
        return mapped_data
    
    def _load_protein(self, force_remapping=False):
        raw_data  = pd.read_csv(self.filepath, sep="\t")
        raw_data.rename(columns={"Gene":"Ensembl_ID", "Gene name":"Symbol"}, inplace=True)
        raw_data = raw_data.loc[raw_data.Reliability.isin(["Enhanced", "Supported","Approved"])]
        tissue_counts = raw_data.groupby("Tissue").Symbol.unique().apply(len).sort_values(ascending=False)
        keep_tissues = tissue_counts[tissue_counts > self.min_obs].index.to_list()
        raw_data = raw_data.loc[raw_data.Tissue.isin(keep_tissues)]
        raw_data = raw_data.loc[raw_data.Level.isin(['Not detected', 'Low', 'Medium', 'High'])]
        raw_data["Level"] = raw_data.Level.map({'Not detected':0, 'Low':1, 'Medium':2, 'High':3})
        mapped_data = self._convert_ids(raw_data, force_remapping=force_remapping)
        tissue_data = mapped_data.pivot_table(index=["Symbol","Ensembl_ID", "Entrez"], columns="Tissue", values="Level", aggfunc='mean').reset_index()
        tissue_data.columns.name = None
        tissue_data.index.name = None
        return tissue_data
    
    def _convert_ids(self, data, force_remapping=False):
        #TODO update to genemapper package
        if ("Entrez" not in data.columns) or (force_remapping):
            all_nodes = list(set([x.split(".")[0] for x in data.Ensembl_ID.to_list()]))
            gmap, unmapped = gene_mapper.update_nodes(nodes=all_nodes, id_type="Ensembl")
            converted_map, missing = gene_mapper.convert_node_ids(list(gmap.values()), "Ensembl", "Entrez")
            mapped_nodes = [gene for gene in [node for node in all_nodes if (node in gmap.keys())] if gene in converted_map.keys()]
            id_map = {node:converted_map[gmap[node]] for node in mapped_nodes}
            data["temp_id"] = data.Ensembl_ID.apply(lambda x: x.split(".")[0])
            data["Entrez"] = data.temp_id.apply(lambda x: id_map[x] if x in id_map.keys() else np.nan)
            data.drop(columns=["temp_id"], inplace=True)
        return data.dropna(subset=["Entrez"])
    
    def _consolidate_transcripts(self):
        dups = self.data[self.data.Entrez.duplicated()].Entrez.to_list()
        if len(dups) == 0:
            return self.data
        else:
            dup_results = []
            if self.datatype == 'protein':
                for gene in dups:
                    gene_data = self.data[self.data.Entrez == gene]
                    consolidated_data = {}
                    for tissue in gene_data.columns[2:-1]:
                        consolidated_data[tissue] = gene_data[tissue].mean()
                    dup_results.append(pd.DataFrame({**consolidated_data, "Entrez":gene}, index=[gene]))
            else:
                for gene in dups:
                    gene_data = self.data[self.data.Entrez == gene]
                    consolidated_data = {}
                    for tissue in gene_data.columns[2:-1]:
                        if min (gene_data[tissue]) > 0.0001:
                            consolidated_data[tissue] = gene_data[tissue].mean()
                        elif max(gene_data[tissue]) < 0.0001:
                            consolidated_data[tissue] = gene_data[tissue].mean()
                        else:
                            consolidated_data[tissue] = max(gene_data[tissue])
                    summary = {}
                    if len(gene_data.Ensembl_ID.unique()) == 1:
                        summary['Ensembl_ID'] = gene_data.Ensembl_ID.unique()
                    else:
                        ids = [x.split(".")[0] for x in gene_data.Ensembl_ID.unique()]
                        if ids[0] == ids[1]:
                            summary['Ensembl_ID'] = ids[0]
                        else:
                            summary['Ensembl_ID'] = "_".join(gene_data.Ensembl_ID.unique())
                    if len(gene_data.Symbol.unique()) == 1:
                        summary['Symbol'] = gene_data.Symbol.unique()
                    else:
                        summary['Symbol'] = "_".join(gene_data.Symbol.unique())
                    dup_results.append(pd.DataFrame({**summary, **consolidated_data, "Entrez":gene}, index=[gene]))
            all_results = pd.concat(dup_results)
            out_data = self.data[~self.data.Entrez.isin(dups)]
            return pd.concat([out_data, all_results]).reset_index(drop=True)
    
    def _categorize_tissues(self):
        if self.datatype == 'mrna':
            tissue_categories = defaultdict(list)
            for tissue in self.data.columns[2:-1]:
                tissue_categories[tissue.split(" - ")[0]].append(tissue)
                tissue_categories['Cardiovascular'] = tissue_categories['Heart'] + tissue_categories['Artery']
                tissue_categories['Digestive'] = tissue_categories['Esophagus'] + tissue_categories['Stomach'] + tissue_categories['Small Intestine'] + tissue_categories['Colon']
        elif self.datatype == 'protein':
            tissue_categories = {'Digestive': ['colon', 'duodenum', 'esophagus', 'gallbladder', 'rectum', 'salivary gland', 'small intestine', 'stomach 1', 'stomach 2', 'oral mucosa'],
                'Liver': ['liver'],
                'Pancreas': ['pancreas'],
                'Appendix' : ['appendix'],
                'Brain': ['caudate', 'cerebellum', 'cerebral cortex', 'hippocampus'],
                'Reproductive (Female)': ['breast', 'cervix', 'endometrium 1', 'endometrium 2', 'fallopian tube', 'ovary', 'placenta', 'vagina'],
                'Reproductive (Male)': ['epididymis', 'prostate', 'seminal vesicle', 'testis'],
                'Respiratory': ['bronchus', 'lung', 'nasopharynx'],
                'Urinary': ['kidney', 'urinary bladder'],
                'Endocrine': ['adrenal gland', 'parathyroid gland', 'thyroid gland'],
                'Muscle': ['skeletal muscle', 'smooth muscle'],
                'Immune': ['bone marrow', 'lymph node', 'spleen', 'tonsil'],
                'Cardiovascular': ['heart muscle'],
                'Skin': ['skin 1', 'skin 2'],
                'Soft tissue': ['soft tissue 1', 'soft tissue 2'],
                'Adipose': ['adipose tissue']}
        else:
            raise ValueError("Datatype must be 'mrna' or 'protein'")
        return tissue_categories

    def _get_gene_idx(self):
        gene_idx = self.data.loc[:, 'Entrez'].to_dict()
    # reverse this dictionary
        gene_idx = {v: k for k, v in gene_idx.items()}
        return gene_idx

    def _subset_tissue_data(self, data, tissue):
        if type(tissue) == str:
            return pd.concat([data.iloc[:, [0,1]], data[self.tissues[tissue]], data.iloc[:, -1]], axis=1)
        else:
            tissue_types = []
            for x in tissue:
                tissue_types += self.tissues[x]
            return pd.concat([data.iloc[:, [0,1]], data[tissue_types], data.iloc[:, -1]], axis=1)

    def get_gene_tissue_exp(self, gene, tissue):
        gene_data = self.get_gene_expression( gene, )
        return self._subset_tissue_data(gene_data, tissue)
    
    def get_gene_global_exp(self, gene):
        gene_data = self.get_gene_expression(gene)
        return gene_data.median()
    
    def get_gene_exp(self, gene):
        return self.data.loc[self.gene_idx[gene], self.data.columns[2:-1]]

    def tissue_expression_median(self, gene , tissue):
        gene_data = self.get_gene_tissue_exp(gene, tissue)
        return gene_data.median()

    def _global_expression_statistics(self):
        #TODO check the axis of this
        statistics_df = self.data.loc[:, self.data.columns[2:-1]].T.describe().T
        statistics_df['median'] = self.data.loc[:, self.data.columns[2:-1]].median(axis=1)
        statistics_df['non_zero'] = self.data.loc[:, self.data.columns[2:-1]].astype(bool).sum(axis=1)
        statistics_df.index = self.data.Entrez
        return statistics_df

    def tissue_expression_statistics(self, tissue):
        tissue_df = self._subset_tissue_data(self.data, tissue)
        # if there are multiple tissues
        if tissue_df.shape[1] > 4:
            # for some reason doing each independently is faster than using the pandas describe function
            statistics_df = tissue_df.loc[:, tissue_df.columns[2:-1]].quantile([0.25, 0.5, 0.75], axis=1).T
            statistics_df['Entrez'] = tissue_df.Entrez
            statistics_df['mean'] = tissue_df.loc[:, tissue_df.columns[2:-1]].mean(axis=1)
            statistics_df['std'] = tissue_df.loc[:, tissue_df.columns[2:-1]].std(axis=1)
            statistics_df['min'] = tissue_df.loc[:, tissue_df.columns[2:-1]].min(axis=1)
            statistics_df['max'] = tissue_df.loc[:, tissue_df.columns[2:-1]].max(axis=1)
            statistics_df['median'] = tissue_df.loc[:, tissue_df.columns[2:-1]].median(axis=1)
            statistics_df['non_zero'] = tissue_df.loc[:, tissue_df.columns[2:-1]].astype(bool).sum(axis=1)
            statistics_df.index = tissue_df.Entrez
            return statistics_df
        else:
            return tissue_df.loc[:, ("Entrez", tissue)].set_index("Entrez", drop=True)

        
def classify_genes(df, min_exp=1, ratio=5, n_group=7, exclude_cols=2):
    tissue_columns = df.columns[exclude_cols:]  # Assuming the first three columns are 'Ensembl_ID', 'Symbol', and 'Gene_Name'
    classifications = {}
    tissue_genes = {tiss:{'Tissue Enriched':[], 'Group Enriched':[], 'Tissue Enhanced':[]} for tiss in tissue_columns}
    for index, row in tqdm(df.iterrows()):
        gene_expression = row[tissue_columns]
        max_expression = gene_expression.max()
        # check if gene is sufficiently expressed in any tissues
        if max_expression <= min_exp:
            classifications[index] = {'Low Expression': []}
            continue
        # Identify the tissue with the maximal expression
        max_tissue = gene_expression.index[np.argmax(gene_expression)]
        mean_expression = gene_expression.mean()
        
        #Tissue Enriched: Genes with an expression level greater than 1 that also have at least 
        #five-fold higher expression levels in a particular tissue compared to all other tissues.
        tissue_enriched = gene_expression[max_tissue] >= ratio * gene_expression.drop(max_tissue)
        if tissue_enriched.all():
            classifications[index] = {'Tissue Enriched': [max_tissue]}
            tissue_genes[max_tissue]['Tissue Enriched'].append(index)
            continue
            
        #Group Enriched: Genes with an expression level greater than 1 that also have at least five-fold higher 
        #expression levels in a group of 2-7 tissues compared to all other tissues, and that are not considered Tissue Enriched.
        top_tissues = list(gene_expression.sort_values(ascending=False).index[0:n_group][::-1])
        group_found = False
        while len(top_tissues) > 1:
            low_tiss = top_tissues[0]
            if gene_expression[low_tiss] < min_exp:
                top_tissues.pop(0)
            else:
                if np.sum(gene_expression[low_tiss] < ratio * gene_expression.drop(low_tiss)) < (len(top_tissues)):
                    group_found = True
                    break
                else:
                    top_tissues.pop(0)
        if group_found:
            classifications[index] = {'Group Enriched': top_tissues}
            for tiss in top_tissues:
                tissue_genes[tiss]['Group Enriched'].append(index)
            continue
            
            
        #Tissue Enhanced: Genes with an expression level greater than 1 that also have at least five-fold higher 
        #expression levels in a particular tissue compared to the average levels in all other tissues, and that 
        #are not considered Tissue Enriched or Group Enriched.
        tissue_enhanced = (gene_expression >= min_exp) & (gene_expression >= [ratio * gene_expression.drop(t).mean() for t in gene_expression.index])
        enhanced_tissues = tissue_enhanced[tissue_enhanced].index.tolist()
        if tissue_enhanced.any():
            classifications[index] = {'Tissue Enhanced': enhanced_tissues}
            #classifications.append(f'Tissue Enhanced in {max_tissue}')
            for tiss in enhanced_tissues:
                tissue_genes[tiss]['Tissue Enhanced'].append(index)
            continue

        classifications[index] = {'Not Classified':[]}
    return classifications, tissue_genes

def get_stats_tissue_enriched(data, nodes, tiss_genes, use_classes = ['Tissue Enriched', 'Group Enriched'], min_tiss_genes=10):
    M = len(data)
    present_nodes = [x for x in nodes if x in data.index.values]
    n = len(present_nodes)
    tissue_specific_genes = get_tissue_genes(tiss_genes, data.index.values, use_classes=use_classes)
    tissue_specific_genes = {t: g for t,g in tissue_specific_genes.items() if len(g) > min_tiss_genes}
    N_vec = pd.Series({t:len(g) for t,g in tissue_specific_genes.items()})
    x_vec = get_shared_genes(tissue_specific_genes, nodes)
    contingency_data = pd.DataFrame({'M': M, 'N': N_vec, 'n': n, 'x': x_vec})
    try:
        contingency_data['P-value'] = contingency_data.apply(lambda y: fisher_exact(np.array([[y.x, y.n - y.x],
                                                                    [y.N - y.x, y.M-(y.n+y.N) + y.x]]))[1], axis=1)
        contingency_data['stat'] = contingency_data.apply(lambda y: fisher_exact(np.array([[y.x, y.n - y.x],
                                                                [y.N - y.x, y.M-(y.n+y.N) + y.x]]))[0], axis=1)
    except ValueError:
        print(contingency_data)
        print(contingency_data.apply(lambda y:np.sum(np.array([[y.x, y.n - y.x],
                                                                [y.N - y.x, y.M-(y.n+y.N) + y.x]])), axis=1))
    return contingency_data

def get_tissue_genes(tiss_genes,background_genes, use_classes=['Tissue Enriched']):
    gene_lists = {}
    for tiss in tiss_genes:
        gene_lists[tiss] = []
        for cl in tiss_genes[tiss]:
            if cl in use_classes:
                gene_lists[tiss] += [g for g in tiss_genes[tiss][cl] if g in background_genes]
    return {t: set(g) for t, g in gene_lists.items()}

def get_shared_genes(gene_lists, nodes):
    out = {}
    for tiss, genes in gene_lists.items():
        out[tiss] = len(set(nodes).intersection(genes))
    return pd.Series(out)
        