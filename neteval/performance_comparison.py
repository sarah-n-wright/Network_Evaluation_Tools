import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
import numpy as np
from sklearn.linear_model import LinearRegression

class EvaluationResults:
	def __init__(self, eval_dir, prefix_file, metrics=["Performance", "Gain"], genesets=['disgen']):
		with open(prefix_file, 'r') as f:
			self.prefixes = [pref.strip() for pref in f.readlines()]
		# make sure metrics is in list form
		if type(metrics) == str:
			self.metrics = [metrics]
		else:
			self.metrics = metrics
		#make sure genesets is in list form
		if type(genesets) == str:
			self.genesets = [genesets]
		else:	
			self.genesets = genesets
		self.eval_dir = eval_dir
		self.results = {}
		self.network_names = {pref:pref for pref in self.prefixes}
		for metric in self.metrics:
			self.results[metric] = {}
			for geneset in self.genesets:
				self.results[metric][geneset] = self.load_evaluation_files(metric, geneset)
		
		
	def load_evaluation_files(self, use_metric, use_geneset):
		paths = {'Performance': ["Performance/", '.performance.csv'], 
			"Gain": ["Performance_Gain/", ".performance_gain.csv"]}
		files = self.get_file_list(paths[use_metric][0],"."+use_geneset+ paths[use_metric][1])
		results_list = []
		for i, f in enumerate(files):
			try:
				df = pd.read_csv(f)
				df["Network"] = self.prefixes[i]
				results_list.append(df)
			except:
				print("FAILED:", f)
		results = pd.concat(results_list).pivot(columns="Network", index="Unnamed: 0", values="0")
		results.index.name=None
		results.columns.name=None
		return results

	def get_file_list(self,metric, suff):
		file_list = []
		for pref in self.prefixes:
			file_list.append(self.eval_dir+metric+pref+suff)
		return file_list

	def rank_all(self):
		ranks = {}
		for metric in self.metrics:
			for geneset in self.genesets:
				ranks[metric + "-" + geneset] = self.rank_networks(metric, geneset)
		self.rankings = pd.DataFrame(ranks)

	def rank_networks(self, use_metric, use_geneset):
		if use_metric not in self.metrics:
			print("Metric not available")
			return
		if use_geneset not in self.genesets:
			print("Geneset not available")
			return
		results = self.results[use_metric][use_geneset]
		rankings = results.rank(axis=1, ascending=False).mean(axis=0)
		return rankings
	
	def set_network_names(self, name_dict):
		for net in self.network_names:
			self.network_names[net] = name_dict[net]

	def plot_rank_scatter(self, metric, geneset1, geneset2, set1_name=None, set2_name=None, labels=True, savepath=None):
		_ = plt.figure(figsize=(7,7))
		ax = plt.gca()
		col1 = metric + "-" + geneset1
		col2 = metric + "-" + geneset2
		n = len(self.rankings)
		plot_df = self.rankings.loc[:, (col1, col2)]
		plot_df.plot(x=col1, y=col2, kind="scatter", ax=ax, s= 50)
		if labels:
			texts = [plt.text(plot_df[col1][i], plot_df[col2][i], plot_df.index[i], fontsize=7) for i in range(len(plot_df))]
			adjust_text(texts)
			# for i in range(len(plot_df)):
			# 	plt.annotate(plot_df.index[i], (plot_df[col1][i], plot_df[col2][i]), fontsize=7)

		if set1_name is not None:
			plt.xlabel(set1_name+" - "+metric, fontsize=14)
			plt.ylabel(set2_name + " - "+ metric, fontsize=14)
		plt.ylim((1, n))
		plt.xlim((1, n))
		ax.invert_yaxis()
		ax.invert_xaxis()
		ax.spines[['right', 'top']].set_visible(False)
		if savepath is not None:
			plt.savefig(savepath, dpi=600, bbox_inches="tight")

	def size_adjusted_performances(self, sizes):
		sizes = {net: np.random.random() for net in self.rankings.index}
		self.size_adjusted = {}
		rankings = {}
		for metric in self.metrics:
			self.size_adjusted[metric] = {}
			for geneset in self.genesets:
				results = self.results[metric][geneset]
				disease_residuals = {}
				for disease in results.index:
					y = results.loc[disease]
					x = np.array([sizes[net] for net in y.index])
					if len(x.shape) == 1:
						x = x.reshape(-1,1)
					reg = LinearRegression().fit(x, y)
					disease_residuals[disease] = y - reg.predict(x)
				self.size_adjusted[metric][geneset] = pd.DataFrame(disease_residuals)
				rankings[metric + "-" + geneset] = pd.DataFrame(disease_residuals).rank(axis=0, ascending=False).mean(axis=1)
		self.size_adjusted_rankings = pd.DataFrame(rankings)
		

# if __name__=="__main__":
# 	eval_dir = "/cellar/users/snwright/Data/Network_Analysis/Evaluation/"
# 	pref = "/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/v2_net_prefixes.txt"
# 	er = EvaluationResults(eval_dir, pref, metrics=["Performance", "Gain"], genesets=['disgen', "cancer"])
# 	er.rank_all()
# 	er.plot_rank_scatter('Performance', 'disgen', 'cancer', savepath='/cellar/users/snwright/Data/Transfer/performance_scatter.png')
# 	er.size_adjusted_performances({pref:np.random.random() for _ in pref})