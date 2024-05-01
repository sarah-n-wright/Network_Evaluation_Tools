import subprocess
import os
import pandas as pd
import random
from sklearn.metrics import auc
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score,precision_recall_curve
import numpy as np
from tqdm import tqdm
import argparse
import multiprocessing as mp
from neteval.data_import_export_tools import load_edgelist_to_networkx
import networkx as nx
import warnings
import csv
from heapq import heappush, heappop
import pickle
from collections import defaultdict

DIRPATH = os.path.dirname(os.path.realpath(__file__))

def create_network_folds(edgefile, nfolds=10, lcc=True):
    if lcc:
        G = load_edgelist_to_networkx(edgefile, id_type="Entrez", testmode=False, timer=None, delimiter="\t", node_cols=[0,1],
                                keep_attributes=False, verbose=False)
        # get largest connected component
        H = max(nx.connected_components(G), key=len)
        net_df = nx.to_pandas_edgelist(G.subgraph(H))
    else:
        net_df = pd.read_csv(edgefile, sep="\t", usecols=['Entrez_A', 'Entrez_B'])
    net_size = len(net_df)
    shuffled_idx = random.sample([i for i in range(net_size)], net_size)
    folds = [shuffled_idx[i::nfolds] for i in range(nfolds)]
    for i in range(nfolds):
        net_df.iloc[folds[i]].to_csv(edgefile+".fold"+str(i+1)+"_test", sep='\t', index=False)
        net_df.drop(index=folds[i]).to_csv(edgefile+".fold"+str(i+1), sep='\t', index=False, header=False)



def l3_edge_prediction(filepath, name, outpath, execdir):
    """_summary_

    Args:
        filepath (_type_): Current path to edge file
        name (_type_): name to append to output files
        outpath (_type_): Location to save the outputs
        execdir (_type_): Location of the L3.out executable
    """
    subprocess.call([DIRPATH+"/L3_prediction.sh", filepath, name, outpath, execdir])
    

def process_prediction_job(args):
    edgefile, net_name, outpath, execdir = args
    l3_edge_prediction(edgefile, net_name, outpath, execdir)
    #os.remove(edgefile)


def get_node_count(edgefile=None, nodelist=None):
    if edgefile is not None:
        command = ['wc', '-l', edgefile.split("_net.txt")[0]+ ".nodelist"]
    elif nodelist is not None:
        command = ['wc', '-l', nodelist]
    else:
        raise ValueError("Must provide either an edgefile or a nodelist.")
    result = subprocess.run(command, capture_output=True, text=True)
    output = int(result.stdout.strip().split()[0])
    return output - 1

def get_edge_count(edgefile):
    command = ['wc', '-l', edgefile]
    result = subprocess.run(command, capture_output=True, text=True)
    output = int(result.stdout.strip().split()[0])
    return output - 1
    

class RankMPS:
    def __init__(self, filename, name_feature_1,
                name_feature_2):
        self.name_feature_1 = name_feature_1
        self.name_feature_2 = name_feature_2
        self.filename=filename
        self.ranked_items = []
        
    def get_top_k_predictions(self, k=500):
        top_items = self.rank_n(k)
        pred_df = pd.DataFrame(top_items, columns=['Entrez_A', 'Entrez_B', self.name_feature_1, self.name_feature_2])
        return pred_df

    def rank_n(self, predict_n, exclude=None):
        with open(self.filename, 'r', buffering=1) as input_file:
            csv_reader = csv.reader(input_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
            header = next(csv_reader)
            index_feature_1 = -1
            index_feature_2 = -1
            for index, field_name in enumerate(header):
                if field_name == self.name_feature_1:
                    index_feature_1 = index
                if field_name == self.name_feature_2:
                    index_feature_2 = index
            #
            # replace with exceptions...
            if index_feature_1 < 0 or index_feature_2 < 0:
                return None

            for record in csv_reader:
                self.process_record(record, index_feature_1, index_feature_2, 
                        predict_n, exclude=exclude)

            top_items = [(score, item) for score, item in self.ranked_items]
            top_items.sort(reverse=True)  # This sorts by score in descending order
            return [item for score, item in top_items]
        
    def process_record(self, record, index_feature_1, index_feature_2, k, exclude=[]):
        if record[0] == "0":
            # print("the couple already edged in the training graph!!!")
            return 
    #
        u = record[1]
        v = record[2]

        value_feature_1 = float(record[index_feature_1])
        value_feature_2 = float(record[index_feature_2])
        #print(value_feature_1, value_feature_2)

        item = (u, v, value_feature_1, value_feature_2)
        score = (value_feature_1, value_feature_2)
        
        if (int(u), int(v)) not in exclude and (int(v), int(u)) not in exclude:
            # Use a tuple (negative score, item) to keep the smallest scores at the top of the heap
            if len(self.ranked_items) < k:
                heappush(self.ranked_items, (score, item))
            else:
                # If the current score is higher than the smallest score in the heap, replace it
                heappush(self.ranked_items, (score, item))
                _ = heappop(self.ranked_items)



class EdgePredictionResults:
    def __init__(self, edgefile, net_name, pred_method, resultdir, suffix='', execdir=None):
        self.edgefile = edgefile
        self.net_name = net_name
        self.suffix=suffix
        self.resultdir = resultdir
        self.pred_method=pred_method
        if self.pred_method == 'L3':
            assert execdir is not None, "Must provide the location of the L3.out executable."
        self.execdir = execdir
        self.edge_count = get_edge_count(edgefile)
        self.node_count = get_node_count(edgefile=edgefile)
    
    def run_single_prediction(self):
        if self.pred_method == 'L3':
            l3_edge_prediction(self.edgefile, self.net_name +self.suffix, self.resultdir, self.execdir)
            prediction_file = self.resultdir + "L3_predictions_" + self.net_name +self.suffix + ".dat"
        else:
            raise NotImplementedError("Only L3 is currently implemented.")
        return prediction_file
    
    def load_prediction_network_lcc(self):
        G = load_edgelist_to_networkx(self.edgefile, id_type="Entrez", testmode=False, timer=None, delimiter="\t", node_cols=[0,1],
                                keep_attributes=False, verbose=False)
        lcc = max(nx.connected_components(G), key=len)
        return G.subgraph(lcc)
    
    def filter_test_edges(self, test_edges, G, verbose=True):
        n1 = len(test_edges)
        test_edges = test_edges.loc[test_edges.Entrez_A.isin(G.nodes) & test_edges.Entrez_B.isin(G.nodes)]
        n2 = len(test_edges)
        if verbose:
            print((n1-n2), 'of', n1, 'test edges removed.')
        return test_edges

    def filter_gold_standard_test_edges(self, test_edges, G, verbose=True):
        n1 = len(test_edges)
        # check which edges are possible
        test_edges = test_edges.loc[test_edges.Entrez_A.isin(G.nodes) & test_edges.Entrez_B.isin(G.nodes)]
        n2 = len(test_edges)
        if verbose:
            print((n1-n2), 'of', n1, 'test edges removed due to impossibility.')
        # test if edges in G already
        test_edges = test_edges.loc[~test_edges.apply(lambda x: G.has_edge(x.Entrez_A, x.Entrez_B), axis=1)]
        n3 = len(test_edges)
        if verbose:
            print((n2-n3), 'of', n2, 'test edges removed due to already existing edges.')
        present = n2-n3
        possible= n2
        return test_edges, present, possible
    
    def calculate_baseline_AUPRC(self, G, test_edges):
        baseline = len(test_edges) / ((len(G.nodes) * (len(G.nodes) - 1) / 2) - len(G.edges))
        return baseline
    
    def get_top_k_predictions(self, k=500):
        if self.pred_method == 'L3':
            pred_df = pd.read_csv(self.resultdir+'L3_predictions_'+self.net_name+self.suffix+'.dat', sep='\t', header=None, names=['Entrez_A', 'Entrez_B', 'Score'], nrows=k)
        if self.pred_method == 'MPS':
            files = os.listdir(self.resultdir+self.net_name+'/Topological/')
            pred_files = [f for f in files if self.suffix+'.csv' in f]
            if len(pred_files) > 1:
                warnings.warn("Multiple prediction files found. Taking the first one.")
            use_file= pred_files[0]
            mps_pred = RankMPS(self.resultdir+self.net_name+'/Topological/'+use_file, 'sim__jac__MAX__sum','Preferential_Attachment')
            pred_df = mps_pred.get_top_k_predictions(k)
            pred_df['Entrez_A'] = pred_df['Entrez_A'].astype(int)
            pred_df['Entrez_B'] = pred_df['Entrez_B'].astype(int)
        return pred_df
    
    def calculate_k_values(self, G, test_edges, k):
        pk = k
        ptest = len(test_edges)
        p1percent=int(len(G.edges)*0.01)
        print('>>>K', k, 'P@K', pk, 'P@TEST', ptest, 'P@1%', p1percent)
        return {f'p@{k}':pk, 'p@test':ptest, 'p@1%':p1percent}
    
    def calculate_stats_at_k(self, G, test_edges, k):

        pred_df = self.get_top_k_predictions(k)
        if self.pred_method=='MPS':
            pred_df['Score'] = [i for i in range(len(pred_df), 0, -1)]
        test_edges = [e for e in test_edges[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)]        
        total_possible_predictions = ((len(G.nodes) * (len(G.nodes) - 1) / 2) - len(G.edges))
        pAUPRC, baseAUPRC, _ = calculate_pAUPRC(test_edges, total_possible_predictions, pred_df)
        results={}
        results['baseAUPRC'] = baseAUPRC
        results['AUPRC'] = pAUPRC
        #results['AUROC'] = calculate_AUROC(test_edges, pred_df)
        results['P@k'] = calculate_Pk(test_edges, pred_df, k)
        return results
        
    def evalutate_all_metrics(self, k):
        test_edges = load_edge_list(self.edgefile+"_test", edge_tuples=False)
        G = self.load_prediction_network_lcc()
        test_edges = self.filter_test_edges(test_edges, G, verbose=True)
        k_dict = self.calculate_k_values(G, test_edges, k)
        all_results = {}
        for name, k in k_dict.items():
            all_results[name] = self.calculate_stats_at_k(G, test_edges, k)
        results_df = pd.DataFrame.from_dict(all_results, orient='index')
        results_df.reset_index(inplace=True)
        results_df.rename(columns={'index':'partial'}, inplace=True)
        return results_df
    
    def evaluate_gold_standard(self, benchmark, k):
        bench_dict = {'corum':os.path.join(DIRPATH, "../Data/corum.pc_net.txt"),
            'panther': os.path.join(DIRPATH, "../Data/panther.pc_net.txt")}
        test_edges = load_edge_list(bench_dict[benchmark], edge_tuples=False)
        G = self.load_prediction_network_lcc()
        test_edges, already_present, possible = self.filter_gold_standard_test_edges(test_edges, G, verbose=True)
        k_dict = self.calculate_k_values(G, test_edges, k)
        all_results = {}
        for name, k in k_dict.items():
            all_results[name] = self.calculate_stats_at_k(G, test_edges, k)
            all_results[name]['Present'] = already_present
            all_results[name]['Possible'] = possible
        results_df = pd.DataFrame.from_dict(all_results, orient='index')
        results_df.reset_index(inplace=True)
        results_df.rename(columns={'index':'partial'}, inplace=True)
        return results_df


    def run_10_fold_cv(self, pred_method, nfolds=10, verbose=True):
        shuffled_idx = random.sample([i for i in range(self.edge_count)], self.edge_count)
        folds = [shuffled_idx[i::nfolds] for i in range(nfolds)]
        net_data = load_edge_list(self.edgefile)
        if pred_method == 'L3':
            prediction_files = []
            test_files = []
            input_files = []
            for i in tqdm(range(nfolds)):
                net_data.loc[folds[i]].to_csv(self.edgefile+".fold"+str(i+1)+"_test", sep='\t', index=False)
                test_files.append(self.edgefile+".fold"+str(i+1)+"_test")
                net_data.drop(index=folds[i]).to_csv(self.edgefile+".fold"+str(i+1), sep='\t', index=False)
                prediction_files.append(self.resultdir + "L3_predictions_" + self.net_name + ".fold"+str(i+1) + ".dat")
                input_files.append(self.edgefile+".fold"+str(i+1))
            if verbose:    
                print("Running prediction on", mp.cpu_count(), "cores.")
            with mp.Pool(3) as pool:
                pool.map(process_prediction_job, [(self.edgefile+".fold"+str(i+1), self.net_name + ".fold"+str(i+1), self.resultdir, self.execdir) for i in range(nfolds)])

        return prediction_files, test_files, input_files
    

    def evaluate_10_fold_cv_performance(self, prediction_files, test_files, metrics = ["AUROC", "AUPRC", "P@500", "NDCG"], verbose=True):
        assert len(prediction_files) == len(test_files), "There must be the same number of performance files as test files."
        #results = {"fold": [i for i in range(len(prediction_files))], **{metric:[] for metric in metrics}}
        #metric_funcs = {"AUROC":self.calculate_AUROC, "AUPRC":self.calculate_AUPRC, "P@500":self.calculate_P500, "NDCG":self.calculate_NDCG}
        # TODO technically this should be checking nodes as well, because sampling could remove nodes from the network
        nfolds = len(prediction_files)
        results = []
        for i in range(nfolds):
            results.append(process_evaluation_job((i, load_edge_list(test_files[i], edge_tuples=True), load_prediction_file(prediction_files[i]), metrics)))
        #with mp.Pool(mp.cpu_count()) as pool:
        #    results = pool.map(process_evaluation_job, [(i, load_edge_list(test_files[i], edge_tuples=True), load_prediction_file(prediction_files[i]), metrics) for i in range(len(prediction_files))])
        
        # for i in range(len(prediction_files)):
        #     predicted_edges = self.load_prediction_file(prediction_files[i])
        #     test_edges = self.load_test_file(test_files[i])
        #     for metric in metrics:
        #         results[metric].append(metric_funcs[metric](test_edges, predicted_edges))
        #     if verbose:
        #         print("Evaluation for fold", i+1, "complete.")
        results_df = pd.DataFrame(results)
        results_df['network'] = self.net_name
        results_df['baselineAUPRC'] = self.calculate_baseline_AUPRC(folds = nfolds)
        return results_df
    

    def evaluate_gold_standard_performance(self, input_edge_files, prediction_files, standard_file, metrics = ["AUROC", "AUPRC", "P@500", "NDCG"], verbose=True):
        assert len(input_edge_files) == len(prediction_files), "There must be the same number of prediction files as input edge files."
        raw_test_edge_df = load_edge_list(standard_file, edge_tuples=False)
        all_results = []
        for i in range(len(input_edge_files)):
            test_edges = filter_test_edges(raw_test_edge_df, input_edge_files[i])
            result = process_evaluation_job((i, test_edges, load_prediction_file(prediction_files[i]), metrics))
            result['baselineAUPRC'] = self.calculate_baseline_AUPRC(test_edges=test_edges)
            all_results.append(result)
        results_df = pd.DataFrame(all_results)
        results_df['network'] = self.net_name
        return results_df
            

    def calculate_baseline_AUPRC(self, test_edges=None, folds=10):
        if test_edges is None:
            positives = int(self.edge_count / folds)
        else:
            positives = len(test_edges)
        # all possible predicted edges is equal to the pairwise combinations of nodes minus the number of edges in the network (without the held out set)
        all_possible = (self.node_count * (self.node_count - 1) / 2) - (self.edge_count - positives)
        return positives/all_possible
    
    def get_top_k_predictions_not_in_test(self, test_edges, k=100):
        if self.pred_method == 'L3':
            top_k_edges_not_in_test = []
            with open(self.resultdir+'L3_predictions_'+self.net_name+self.suffix+'.dat', 'r') as f:
                for line in f:
                    edge = line.strip().split('\t')
                    if (int(edge[0]), int(edge[1])) not in test_edges and (int(edge[1]), int(edge[0])) not in test_edges:
                        top_k_edges_not_in_test.append((int(edge[0]), int(edge[1])))
                    if len(top_k_edges_not_in_test) > k:
                        break
        elif self.pred_method == 'MPS':
            files = os.listdir(self.resultdir+self.net_name+'/Topological/')
            pred_files = [f for f in files if self.suffix+'.csv' in f]
            if len(pred_files) > 1:
                warnings.warn("Multiple prediction files found. Taking the first one.")
            use_file= pred_files[0]
            mps_pred = RankMPS(self.resultdir+self.net_name+'/Topological/'+use_file, 'sim__jac__MAX__sum','Preferential_Attachment')
            top_k_edges_not_in_test = mps_pred.rank_n(k, exclude=test_edges)
            top_k_edges_not_in_test = [(int(x[0]), int(x[1])) for x in top_k_edges_not_in_test]
        return top_k_edges_not_in_test
    
    def get_k_random_edges_not_in_G(self, G, test_edges,  k=100):
        random_edge_set = []
        while len(random_edge_set) < k:
            edge = (random.choice(list(G.nodes)), random.choice(list(G.nodes)))
            if edge[0] != edge[1] and not G.has_edge(*edge):
                if edge not in test_edges and edge[::-1] not in test_edges:
                    random_edge_set.append(edge)
        return random_edge_set
    
    def load_edge_network_coverage(self, coverage_file):
        if '.pkl' in coverage_file:
            with open(coverage_file, 'rb') as f:
                edge_counts = pickle.load(f)
            return defaultdict(int, edge_counts)
        else:
            raise NotImplementedError("Only pickled edge counts are currently implemented.")

    def evaluate_random_network_coverage(self, random_edges, coverage_dict,  k=100):
        random_coverage = [max(coverage_dict[edge], coverage_dict[edge[::-1]]) for edge in random_edges]
        return random_coverage
    
    def evaluate_network_coverage(self, coverage_file, k=100, permutations=100):
        # get top predictions not in the network
        G = self.load_prediction_network_lcc()
        test_edges = load_edge_list(self.edgefile+"_test", edge_tuples=True)
        top_k_edges_not_in_G = self.get_top_k_predictions_not_in_test(test_edges, k)
        edge_coverage = self.load_edge_network_coverage(coverage_file)
        random_means = []
        random_medians = []
        for _ in range(permutations):
            random_edges = self.get_k_random_edges_not_in_G(G, test_edges, k)
            random_coverage = [max(edge_coverage[edge], edge_coverage[edge[::-1]]) for edge in random_edges]
            random_means.append(np.mean(random_coverage))
            random_medians.append(np.median(random_coverage))
        top_k_coverage = {edge:max(edge_coverage[edge], edge_coverage[edge[::-1]]) for edge in top_k_edges_not_in_G}
        top_k_mean = np.mean([x for x in top_k_coverage.values()])
        top_k_median = np.median([x for x in top_k_coverage.values()])
        unverified_edges = [edge for edge in top_k_edges_not_in_G if top_k_coverage[edge] < 2]
        df = pd.DataFrame({'random_means':random_means, 'random_medians':random_medians, 'top_k_mean':top_k_mean, 'top_k_median':top_k_median, 'k':k,
                           '0-1':len([x for x in top_k_coverage.values() if x <= 1]), '2-5': len([x for x in top_k_coverage.values() if 2 <= x <= 5]), 
                           '6-10': len([x for x in top_k_coverage.values() if 6 <= x <= 10]), '11-20': len([x for x in top_k_coverage.values() if 11 <= x <= 20]),
                           '21+': len([x for x in top_k_coverage.values() if x > 20])})
        return df, unverified_edges

def load_prediction_file(filepath, edge_tuples=False):
    df = pd.read_csv(filepath, header=None, names=['Entrez_A', 'Entrez_B', 'Score'], sep="\t")
    df[['Entrez_A', 'Entrez_B']] = np.sort(df[['Entrez_A', 'Entrez_B']].values, axis=1)  # Sort each edge
    if edge_tuples:
        return set(df.iloc[:, [0,1]].itertuples(index=False, name=None))
    else:
        return df


def load_edge_list(filepath, edge_tuples=False):
    df = pd.read_csv(filepath, sep="\t", usecols=[0,1])
    df.columns = ['Entrez_A', 'Entrez_B']
    try:
        df[['Entrez_A', 'Entrez_B']] = np.sort(df[['Entrez_A', 'Entrez_B']].values, axis=1)
    except KeyError as e:
        print(filepath)
        print(df.head())# Sort each edge
        raise e
    if edge_tuples:
        return set(df.itertuples(index=False, name=None))
    else:
        return df
    
def filter_test_edges(raw_test_edge_df, edgefile):
    base_edges_df = load_edge_list(edgefile, edge_tuples=False)
    nodes = set(base_edges_df.Entrez_A.values).union(set(base_edges_df.Entrez_B.values))
    node_filtered_df = raw_test_edge_df[raw_test_edge_df.Entrez_A.isin(nodes) & raw_test_edge_df.Entrez_B.isin(nodes)]
    node_filtered = set(node_filtered_df.iloc[:, [0,1]].itertuples(index=False, name=None))
    base_edges = set(base_edges_df.iloc[:, [0,1]].itertuples(index=False, name=None))
    edge_filtered = node_filtered - base_edges
    return edge_filtered


    
def process_evaluation_job(args):
    i, test_edges, predicted_edges, metrics = args
    metric_funcs = {"AUROC":calculate_AUROC, "AUPRC":calculate_AUPRC, "P@500":calculate_Pk, "NDCG":calculate_NDCG}
    results = {'fold':i}
    for metric in metrics:
        if metric in metric_funcs:
            results[metric] = metric_funcs[metric](test_edges, predicted_edges)
    return results

def calculate_pAUPRC(test, total_possible, pred_df):
    pred_edges = pred_df[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)
    y_true = [(x in test) or (x[::-1] in test) for x in pred_edges]
    scores = pred_df['Score'].values
    if np.sum(y_true) == 0:
        return 0, 0, 0
    precision, recall, _ = precision_recall_curve(y_true, scores)
    recall = recall * (np.sum(y_true)/len(test))
    baseline = (len(test)/total_possible) * max(recall)
    pAUPRC = auc(recall, precision)
    return pAUPRC, baseline, pAUPRC/baseline
    
            
def calculate_AUROC(test_edges, predicted_edges):
    y_true = [edge in test_edges for edge in predicted_edges[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)]
    scores = predicted_edges['Score'].values
    auroc = roc_auc_score(y_true, scores)
    return auroc

def calculate_AUPRC(test_edges, predicted_edges):
    y_true = [edge in test_edges for edge in predicted_edges[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)]
    scores = predicted_edges['Score'].values
    precision, recall, _ = precision_recall_curve(y_true, scores)
    return auc(recall, precision)

def calculate_Pk(test_edges, predicted_edges, k=500):
    top_k_preds = predicted_edges.nlargest(k, 'Score')[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)
    relevant_at_k = np.sum([(x in test_edges) or (x[::-1] in test_edges) for x in top_k_preds])
    return relevant_at_k / k

def calculate_NDCG(test_edges, predicted_edges):
    y_true = np.array([edge in test_edges for edge in predicted_edges[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)])
    scores = predicted_edges['Score'].values
    return ndcg_at_k(y_true, scores)

def ndcg_at_k(y_true, scores, k=None):
    actual = dcg_at_k(y_true, scores, k)
    best = dcg_at_k(y_true, y_true, k)
    return actual / best

def dcg_at_k( y_true, scores, k=None):
    order = np.argsort(scores)[::-1]
    y_true = np.take(y_true, order[:k])
    gain = 2 ** y_true - 1
    discounts = np.log2(np.arange(len(y_true)) + 2)
    return np.sum(gain / discounts)

def evaluation_wrapper(args, predicted_files, test_files, input_files=None, benchmarks=['self'], verbose=True, outdir=DIRPATH+"../Data/example_outputs/"):
    print('>>>BENCHMARKS:', benchmarks)
    if 'self' in benchmarks:
        out = epr.evaluate_10_fold_cv_performance(predicted_files, test_files)
        out.to_csv(outdir +args.networkprefix + "_L3_results.tsv", sep="\t", index=False)
        if verbose:
            print("Self evaluation complete.")
    if 'corum' in benchmarks:
        assert input_files is not None, "Must provide input files to evaluate against corum."
        corum_file = os.path.join(DIRPATH, "../Data/corum.pc_net.txt")
        corum = epr.evaluate_gold_standard_performance(input_files, predicted_files, corum_file)
        corum.to_csv(outdir +args.networkprefix + "_L3_corum_results.tsv", sep="\t", index=False)
        if verbose: 
            print("Corum evaluation complete.")
    if 'panther' in benchmarks:
        panther_file = os.path.join(DIRPATH, "../Data/panther.pc_net.txt")
        panther = epr.evaluate_gold_standard_performance(input_files, predicted_files, panther_file)
        panther.to_csv(outdir+args.networkprefix + "_L3_panther_results.tsv", sep="\t", index=False)
        if verbose:
            print("Panther evaluation complete.")
        
    


if __name__ == "__main__":
    #edgefile = '/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_fixed/bioplex.v3.293T_net.txt'
    #create_network_folds(edgefile, nfolds=10, lcc=True)
    # Run the edge prediction script
    parser = argparse.ArgumentParser(description='Run edge prediction on a network.')
    parser.add_argument("--networkprefix", type=str, 
		help='Prefix of network to be evaluated. File must be 2-column edge list where each line is a gene interaction separated by a tab delimiter.')
    parser.add_argument("--runwhat", type=str, default="Both", choices=["Prediction", "Evaluation", "Both", "Folds", 'PredictFolds','GoldStandard', 'Coverage'])
    parser.add_argument('--benchmarks',nargs='+', default=['self'], choices=['self', 'corum', 'kegg', 'panther'])
    parser.add_argument('--networksuffix', type=str, default="")
    parser.add_argument('--outdir', type=str, default='/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/L3/L3_outputs_fixed/')
    parser.add_argument('--pred_method', type=str, default='L3' )
    parser.add_argument('--datadir', type=str, default='/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_fixed/')
    parser.add_argument('--execdir', type=str, default="/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/kpisti-L3-ed6b18f/", help='Location of the L3.out executable.')
    args = parser.parse_args()
    
    if args.pred_method == 'MPS':
        assert args.runwhat in ['Evaluation', 'GoldStandard', 'Coverage'], "MPS is only implemented for evaluation and gold standard evaluation."

    print('>>ARGS:', args)
    
    if args.runwhat in ["Folds"]:
        create_network_folds(args.datadir+ args.networkprefix+"_net.txt", nfolds=10, lcc=True)
    if args.runwhat in ["PredictFolds"]:
        epr = EdgePredictionResults(args.datadir+ args.networkprefix+"_net.txt"+args.networksuffix, args.networkprefix, suffix=args.networksuffix, resultdir=args.outdir, pred_method='L3', execdir=args.execdir)
        epr.run_single_prediction()
    if args.runwhat in ["Prediction", "Both"]:
        epr = EdgePredictionResults(args.datadir+ args.networkprefix+"_net.txt", args.networkprefix, suffix=args.networksuffix, resultdir=args.outdir, pred_method='L3', execdir=args.execdir)
        predicted_files, test_files, input_files = epr.run_10_fold_cv("L3", nfolds=10)
    # save the prediciton and test file paths
        pd.DataFrame({"prediction_files":predicted_files, "test_files":test_files, "input_files": input_files}).to_csv(args.outdir +args.networkprefix + "_L3_filepaths.tsv", sep="\t", index=False)
    if args.runwhat in ['Evaluation']:
        print('Creating EdgePredictionResults object')
        epr = EdgePredictionResults(args.datadir+ args.networkprefix+"_net.txt"+args.networksuffix, args.networkprefix, suffix=args.networksuffix, resultdir=args.outdir, pred_method=args.pred_method, execdir=args.execdir)
        print('Calculating statistics')
        results_df = epr.evalutate_all_metrics(k=500)
        print('Writing statistics')
        if args.pred_method=='L3':
            results_df.to_csv(args.outdir +args.networkprefix+args.networksuffix +'_'+args.pred_method + "_results.tsv", sep='\t', index=False)
        elif args.pred_method=='MPS':
            results_df.to_csv(args.outdir +args.networkprefix+'/Predictions/'+args.networkprefix+args.networksuffix +'_'+args.pred_method + "_results.tsv", sep='\t', index=False)
    if args.runwhat in ['GoldStandard']:
        epr = EdgePredictionResults(args.datadir+ args.networkprefix+"_net.txt"+args.networksuffix, args.networkprefix, suffix=args.networksuffix, resultdir=args.outdir, pred_method=args.pred_method, execdir=args.execdir)
        for bench in args.benchmarks:
            print('Evaluating', bench)
            results_df = epr.evaluate_gold_standard(bench, k=500)
            if args.pred_method == 'L3':
                results_df.to_csv(args.outdir +args.networkprefix+args.networksuffix +'_'+args.pred_method + "_"+bench+"_results.tsv", sep='\t', index=False)
            elif args.pred_method == 'MPS':
                results_df.to_csv(args.outdir +args.networkprefix+'/Predictions/'+args.networkprefix+args.networksuffix +'_'+args.pred_method + "_"+bench+"_results.tsv", sep='\t', index=False)
    if args.runwhat in ['Coverage']:
        epr = EdgePredictionResults(args.datadir+ args.networkprefix+"_net.txt"+args.networksuffix, args.networkprefix, suffix=args.networksuffix, resultdir=args.outdir, pred_method=args.pred_method, execdir=args.execdir)
        coverage_file = os.path.join(DIRPATH, '../Data/mapped_edge_counts.pkl')
        results_df, unverified_edges = epr.evaluate_network_coverage(coverage_file, k=100)
        if args.pred_method == 'L3':
            results_df.to_csv(args.outdir +args.networkprefix+args.networksuffix +'_'+args.pred_method + "_coverage_results.tsv", sep='\t', index=False)
            with open(args.outdir +args.networkprefix+args.networksuffix +'_'+args.pred_method + "_unverified_edges.txt", 'w') as f:
                edge_strings = ['\t'.join([str(edge[0]), str(edge[1])]) for edge in unverified_edges]
                f.write('\n'.join(edge_strings) + '\n')
        if args.pred_method == 'MPS':
            results_df.to_csv(args.outdir +args.networkprefix+'/Predictions/'+args.networkprefix+args.networksuffix +'_'+args.pred_method + "_coverage_results.tsv", sep='\t', index=False)
            with open(args.outdir +args.networkprefix+'/Predictions/'+args.networkprefix+args.networksuffix +'_'+args.pred_method + "_unverified_edges.txt", 'w') as f:
                edge_strings = ['\t'.join([str(edge[0]), str(edge[1])]) for edge in unverified_edges]
                f.write('\n'.join(edge_strings) + '\n')
    if args.runwhat in ["Both"]:
        epr = EdgePredictionResults(args.datadir+ args.networkprefix+"_net.txt", args.networkprefix, suffix=args.networksuffix, resultdir=args.outdir, pred_method='L3', execdir=args.execdir)
        if args.runwhat == "Evaluation":
            files_df = pd.read_csv(args.outdir+args.networkprefix + "_L3_filepaths.tsv", sep="\t")
            predicted_files = list(files_df['prediction_files'].values)
            test_files = list(files_df['test_files'].values)
            input_files = list(files_df['input_files'].values)
            evaluation_wrapper(args, predicted_files, test_files, input_files, benchmarks=args.benchmarks, outdir=args.outdir)
        #out = pd.read_csv("/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/" +args.networkprefix + "_L3_results.tsv", sep="\t")
        #out = epr.evaluate_10_fold_cv_performance(predicted_files, test_files)
        #out.to_csv("/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/" +args.networkprefix + "_L3_results.tsv", sep="\t", index=False)
        #corum_file = "/cellar/users/snwright/Data/Network_Analysis/Reference_Data/corum.pc_net.txt"
        #corum = epr.evaluate_gold_standard_performance(input_files, predicted_files, corum_file)
        #corum.to_csv("/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/" +args.networkprefix + "_L3_corum_results.tsv", sep="\t", index=False)