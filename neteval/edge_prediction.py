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

class EdgePredictionResults:
    def __init__(self, edgefile, net_name):
        self.edgefile = edgefile
        self.net_name = net_name
        self.edge_count = self.get_edge_count()
        
    def get_edge_count(self):
        command = ['wc', '-l', self.edgefile]
        result = subprocess.run(command, capture_output=True, text=True)
        output = int(result.stdout.strip().split()[0])
        return output - 1
    
    def run_10_fold_cv(self, pred_method, nfolds=10):
        shuffled_idx = random.sample([i for i in range(self.edge_count)], self.edge_count)
        folds = [shuffled_idx[i::nfolds] for i in range(nfolds)]
        net_data = self.load_edge_list(self.edgefile)
        if pred_method == 'L3':
            prediction_files = []
            test_files = []
            for i in tqdm(range(nfolds)):
                net_data.loc[folds[i]].to_csv(self.edgefile+".fold"+str(i+1)+"_test", sep='\t', index=False)
                test_files.append(self.edgefile+".fold"+str(i+1)+"_test")
                net_data.drop(index=folds[i]).to_csv(self.edgefile+".fold"+str(i+1), sep='\t', index=False)
                l3_edge_prediction(self.edgefile+".fold"+str(i+1), self.net_name + "fold"+str(i+1), RESULTDIR, EXECDIR)
                prediction_files.append(RESULTDIR + "L3_predictions_" + self.net_name + "fold"+str(i+1) + ".dat")
                os.remove(self.edgefile+".fold"+str(i+1))
        return prediction_files, test_files
        
    def evaluate_10_fold_cv_performance(self, prediction_files, test_files, metrics = ["AUROC", "AUPRC", "P@500", "NDCG"]):
        assert len(prediction_files) == len(test_files), "There must be the same number of performance files as test files."
        results = {"fold": [i for i in range(len(prediction_files))], **{metric:[] for metric in metrics}}
        metric_funcs = {"AUROC":self.calculate_AUROC, "AUPRC":self.calculate_AUPRC, "P@500":self.calculate_P500, "NDCG":self.calculate_NDCG}
        for i in range(len(prediction_files)):
            predicted_edges = self.load_prediction_file(prediction_files[i])
            test_edges = self.load_test_file(test_files[i])
            for metric in metrics:
                results[metric].append(metric_funcs[metric](test_edges, predicted_edges))
        return pd.DataFrame(results)
                
    def load_prediction_file(self, filepath):
        df = pd.read_csv(filepath, header=None, names=['Entrez_A', 'Entrez_B', 'Score'], sep="\t")
        df[['Entrez_A', 'Entrez_B']] = np.sort(df[['Entrez_A', 'Entrez_B']].values, axis=1)  # Sort each edge
        return df
    
    def load_test_file(self, filepath):
        #TODO check the ordering of the nodes?
        df = pd.read_csv(filepath, sep="\t", usecols=[0,1])
        df[['Entrez_A', 'Entrez_B']] = np.sort(df[['Entrez_A', 'Entrez_B']].values, axis=1)  # Sort each edge
        return set(df.itertuples(index=False, name=None))
    
    def load_edge_list(self, filepath, edge_tuples=False):
        df = pd.read_csv(filepath, sep="\t", usecols=[0,1])
        try:
            df[['Entrez_A', 'Entrez_B']] = np.sort(df[['Entrez_A', 'Entrez_B']].values, axis=1)
        except KeyError as e:
            print(self.edgefile)
            print(df.head())# Sort each edge
            raise e
        if edge_tuples:
            return set(df.iloc[:, [0,1]].itertuples(index=False, name=False))
        else:
            return df
            
    def calculate_AUROC(self, test_edges, predicted_edges):
        y_true = [edge in test_edges for edge in predicted_edges[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)]
        scores = predicted_edges['Score'].values
        auroc = roc_auc_score(y_true, scores)
        return auroc
    
    def calculate_AUPRC(sel, test_edges, predicted_edges):
        y_true = [edge in test_edges for edge in predicted_edges[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)]
        scores = predicted_edges['Score'].values
        precision, recall, _ = precision_recall_curve(y_true, scores)
        return auc(recall, precision)
    
    def calculate_P500(self, test_edges, predicted_edges, k=500):
        top_k_preds = predicted_edges.nlargest(k, 'Score')[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)
        relevant_at_k = len(set(top_k_preds) & test_edges)
        return relevant_at_k / k
    
    def calculate_NDCG(self, test_edges, predicted_edges):
        y_true = np.array([edge in test_edges for edge in predicted_edges[['Entrez_A', 'Entrez_B']].itertuples(index=False, name=None)])
        scores = predicted_edges['Score'].values
        return self.ndcg_at_k(y_true, scores)
    
    def ndcg_at_k(self, y_true, scores, k=None):
        actual = self.dcg_at_k(y_true, scores, k)
        best = self.dcg_at_k(y_true, y_true, k)
        return actual / best

    def dcg_at_k(self, y_true, scores, k=None):
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
    args = parser.parse_args()
    epr = EdgePredictionResults(DATADIR+ args.network_prefix+"_net.txt", args.network_prefix)
    a, b = epr.run_10_fold_cv("L3", nfolds=10)
    out = epr.evaluate_10_fold_cv_performance(a, b)
    out.to_csv("/cellar/users/snwright/Data/Network_Analysis/Edge_Prediction/" +args.network_prefix + "_L3_results.tsv", sep="\t", index=False)