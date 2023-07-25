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

DIRPATH = os.path.dirname(os.path.realpath(__file__))
EXECDIR = "/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/kpisti-L3-ed6b18f/"
DATADIR = "/cellar/users/snwright/Data/Network_Analysis/Processed_Data/v2_2022/"
RESULTDIR = "/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/L3_outputs/"

def l3_edge_prediction(filepath, name, outpath, execdir):
    """_summary_

    Args:
        filepath (_type_): Current path to edge file
        name (_type_): name to append to output files
        outpath (_type_): Location to save the outputs
        execdir (_type_): Location of the L3.out executable
    """
    subprocess.call([DIRPATH+"/../L3_prediction.sh", filepath, name, outpath, execdir])
    
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
    


class EdgePredictionResults:
    def __init__(self, edgefile, net_name):
        self.edgefile = edgefile
        self.net_name = net_name
        self.edge_count = get_edge_count(edgefile)
        self.node_count = get_node_count(edgefile=edgefile)
    
    def run_single_prediction(self, pred_method):
        if pred_method == 'L3':
            l3_edge_prediction(self.edgefile, self.net_name, RESULTDIR, EXECDIR)
            prediction_file = RESULTDIR + "L3_predictions_" + self.net_name + ".dat"
        else:
            raise NotImplementedError("Only L3 is currently implemented.")
        return prediction_file
    
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
                prediction_files.append(RESULTDIR + "L3_predictions_" + self.net_name + "fold"+str(i+1) + ".dat")
                input_files.append(self.edgefile+".fold"+str(i+1))
            if verbose:    
                print("Running prediction on", mp.cpu_count(), "cores.")
            with mp.Pool(5) as pool:
                pool.map(process_prediction_job, [(self.edgefile+".fold"+str(i+1), self.net_name + "fold"+str(i+1), RESULTDIR, EXECDIR) for i in range(nfolds)])

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
    
    
def load_prediction_file(filepath, edge_tuples=False):
    df = pd.read_csv(filepath, header=None, names=['Entrez_A', 'Entrez_B', 'Score'], sep="\t")
    df[['Entrez_A', 'Entrez_B']] = np.sort(df[['Entrez_A', 'Entrez_B']].values, axis=1)  # Sort each edge
    if edge_tuples:
        return set(df.iloc[:, [0,1]].itertuples(index=False, name=None))
    else:
        return df

# def load_test_file(filepath, names=['Entrez_A', 'Entrez_B'], edge_tuples=True):
#     #TODO check the ordering of the nodes?
#     df = pd.read_csv(filepath, sep="\t", usecols=[0,1], names=names)
#     df[['Entrez_A', 'Entrez_B']] = np.sort(df[['Entrez_A', 'Entrez_B']].values, axis=1)
#     if edge_tuples:# Sort each edge
#         return set(df.itertuples(index=False, name=None))
#     else:
#         return df

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
    metric_funcs = {"AUROC":calculate_AUROC, "AUPRC":calculate_AUPRC, "P@500":calculate_P500, "NDCG":calculate_NDCG}
    results = {'fold':i}
    for metric in metrics:
        if metric in metric_funcs:
            results[metric] = metric_funcs[metric](test_edges, predicted_edges)
    return results
    
            
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

def calculate_P500(test_edges, predicted_edges, k=500):
    top_k_preds = predicted_edges.nlargest(k, 'Score')[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)
    relevant_at_k = len(set(top_k_preds) & test_edges)
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


if __name__ == "__main__":
    # Run the edge prediction script
    parser = argparse.ArgumentParser(description='Run edge prediction on a network.')
    parser.add_argument("network_prefix", type=str, 
		help='Prefix of network to be evaluated. File must be 2-column edge list where each line is a gene interaction separated by a tab delimiter.')
    parser.add_argument("run_what", type=str, default="Both", choices=["Prediction", "Evaluation", "Both"])
    args = parser.parse_args()
    # class Args:
    #     def __init__(self):
    #         self.network_prefix = 'bioplex.v3.293T'
    #         self.run_what = "Evaluation"
    # args = Args()
    
    epr = EdgePredictionResults(DATADIR+ args.network_prefix+"_net.txt", args.network_prefix)
    if args.run_what in ["Prediction", "Both"]:
        predicted_files, test_files, input_files = epr.run_10_fold_cv("L3", nfolds=10)
    # save the prediciton and test file paths
        pd.DataFrame({"prediction_files":predicted_files, "test_files":test_files, "input_files": input_files}).to_csv("/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/" +args.network_prefix + "_L3_filepaths.tsv", sep="\t", index=False)
    if args.run_what in ["Both", 'Evaluation']:
        if args.run_what == "Evaluation":
            files_df = pd.read_csv("/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/" +args.network_prefix + "_L3_filepaths.tsv", sep="\t")
            predicted_files = list(files_df['prediction_files'].values)
            test_files = list(files_df['test_files'].values)
            input_files = list(files_df['input_files'].values)
        #out = pd.read_csv("/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/" +args.network_prefix + "_L3_results.tsv", sep="\t")
        out = epr.evaluate_10_fold_cv_performance(predicted_files, test_files)
        out.to_csv("/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/" +args.network_prefix + "_L3_results.tsv", sep="\t", index=False)
        corum_file = "/cellar/users/snwright/Data/Network_Analysis/Reference_Data/corum.pc_net.txt"
        corum = epr.evaluate_gold_standard_performance(input_files, predicted_files, corum_file)
        corum.to_csv("/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/" +args.network_prefix + "_L3_corum_results.tsv", sep="\t", index=False)