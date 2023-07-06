import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
import numpy as np
from sklearn.linear_model import LinearRegression
import seaborn as sns

class EvaluationResults:
	def __init__(self, eval_dir, prefix_file, metrics=["Performance", "Gain"], genesets=['disgen']):
		if type(prefix_file) == str:
			with open(prefix_file, 'r') as f:
				self.prefixes = [pref.strip() for pref in f.readlines()]
		else:
			self.prefixes = prefix_file
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

	def size_adjusted_performances(self, sizes, debug=False):
		if debug:
			sizes = {net:[np.random.randint(1,100)] for net in self.prefixes}
		self.size_adjusted = {}
		self.sizes = sizes
		rankings = {}
		self.fits = {}
		for metric in self.metrics:
			self.size_adjusted[metric] = {}
			self.fits[metric] = {}
			for geneset in self.genesets:
				results = self.results[metric][geneset]
				disease_residuals = {}
				fits = {}
				for disease in results.index:
					y = results.loc[disease]
					x = np.array([sizes[net] for net in y.index])
					if len(x.shape) == 1:
						x = x.reshape(-1,1)
					reg = LinearRegression().fit(x, y)
					fits[disease] = (reg.coef_, reg.intercept_)
					disease_residuals[disease] = y - reg.predict(x)
				self.size_adjusted[metric][geneset] = pd.DataFrame(disease_residuals)
				self.fits[metric][geneset] = fits
				rankings[metric + "-" + geneset] = pd.DataFrame(disease_residuals).rank(axis=0, ascending=False).mean(axis=1)
		self.size_adjusted_rankings = pd.DataFrame(rankings)

	def plot_size_fit(self, metric, geneset, disease):
		# check that size_adjusted has been run
		if metric not in self.size_adjusted:
			print("Metric not available")
			return
		if geneset not in self.size_adjusted[metric]:
			print("Geneset not available")
			return
		if disease not in self.size_adjusted[metric][geneset]:
			print("Disease not available")
			return
		# get data
		raw_results = self.results[metric][geneset].loc[disease]
		if len(self.sizes[self.prefixes[0]]) == 1:
			size_df = pd.DataFrame({"Size":self.sizes})
			size_df['Size'] = size_df["Size"].apply(lambda x: x[0])
		else:
			raise NotImplementedError("Size adjustment plotting only implemented for single size per network")
		raw_results = pd.DataFrame({metric:raw_results}).join(size_df)
		size_results = self.size_adjusted[metric][geneset].loc[:,disease]
		size_results = pd.DataFrame({metric+"_adj":size_results}).join(size_df)
		fit_coeffs = self.fits[metric][geneset][disease]

		self._calculate_fit_line(metric, geneset, disease)
		fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(14,7))
		sns.scatterplot(data=raw_results, x="Size", y=metric, ax=ax1)
		sns.scatterplot(data=size_results, x="Size", y=metric+"_adj", ax=ax2)
		fit_data = self._calculate_fit_line(metric, geneset, disease)
		sns.lineplot(data = fit_data, x="Size", y="Fit", ax=ax1, color="red")

	def _calculate_fit_line(self, metric, geneset, disease):
		try:
			fit_coeffs = self.fits[metric][geneset][disease]
		except KeyError:
			print("Disease not available")
			return
		intercept = fit_coeffs[1]
		coeffs = fit_coeffs[0]
		fit_line_df = pd.DataFrame({'Size':self.sizes}) # this won't work for multiple size data points
		fit_line_df["Size"] = fit_line_df["Size"].apply(lambda x: x[0])
		x = np.array([self.sizes[net] for net in fit_line_df.index])
		fit_line_df['Fit'] = np.multiply(x, coeffs) + intercept
		fit_line_df = fit_line_df.sort_values(by="Size")
		return fit_line_df
		
		

if __name__=="__main__":
	eval_dir = "/cellar/users/snwright/Data/Network_Analysis/Evaluation/"
	pref = "/cellar/users/snwright/Git/Network_Evaluation_Tools/Data/v2_net_prefixes.txt"
	er = EvaluationResults(eval_dir, pref, metrics=["Performance", "Gain"], genesets=['disgen', "cancer"])
	er.rank_all()
	er.plot_rank_scatter('Performance', 'disgen', 'cancer', savepath='/cellar/users/snwright/Data/Transfer/performance_scatter.png')
	er.size_adjusted_performances(sizes=None, debug=True)
	er.plot_size_fit("Performance", "disgen", "trachomatis")
	er.fits.keys()
