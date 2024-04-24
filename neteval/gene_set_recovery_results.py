import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression, HuberRegressor
import seaborn as sns
import scipy.cluster.hierarchy as sch
from sklearn.cluster import SpectralBiclustering
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import os
from collections import defaultdict


def cluster_matrix(mat, **kwargs):
    """Function to cluster a matrix by rows and columns

    Args:
        mat (pandas.DataFrame): Matrix to cluster
        kwargs: Keyword arguments to pass to scipy.cluster.hierarchy.linkage
    Returns:
        pandas.DataFrame: Clustered matrix
    """
    # Clustering on rows
    pairwise_distances_rows = sch.distance.pdist(mat.values)
    linkage_rows = sch.linkage(pairwise_distances_rows, **kwargs)
    cluster_distance_threshold_rows = pairwise_distances_rows.max()/2
    idx_to_cluster_array_rows = sch.fcluster(linkage_rows, cluster_distance_threshold_rows, criterion='distance')
    idx_rows = np.argsort(idx_to_cluster_array_rows)
    
    # Clustering on columns
    pairwise_distances_cols = sch.distance.pdist(mat.T.values)
    linkage_cols = sch.linkage(pairwise_distances_cols, **kwargs)
    cluster_distance_threshold_cols = pairwise_distances_cols.max()/2
    idx_to_cluster_array_cols = sch.fcluster(linkage_cols, cluster_distance_threshold_cols, criterion='distance')
    idx_cols = np.argsort(idx_to_cluster_array_cols)
    
    return mat.iloc[idx_rows, idx_cols]


def bicluster_matrix(mat, **kwargs):
    """Function to bicluster a matrix using spectral biclustering
    args:
        mat (pandas.DataFrame): Matrix to bicluster
        kwargs: Keyword arguments to pass to sklearn.cluster.SpectralBiclustering
    returns:
        pandas.DataFrame: Biclustered matrix
    """
    # Normalize the values to range from 0 to 1
    df_normalized = (mat - mat.min().min()) / (mat.max().max() - mat.min().min())

    # Perform spectral biclustering
    model = SpectralBiclustering(random_state=0, **kwargs)
    model.fit(df_normalized.values)

    # Rearrange the rows and columns of the DataFrame according to the biclusters
    mat_biclustered = mat.iloc[np.argsort(model.row_labels_)]
    mat_biclustered = mat_biclustered.iloc[:, np.argsort(model.column_labels_)]
    return mat_biclustered


def cluster_by_annot(df, annot_types, scale=False):
    """Function to cluster a matrix by column annotation types. First clusters the high-level annotation 
    types, then clusters the columns within each high-level annotation type, 
    then clusters the rows across the whole matrix.

    args:
        df (pandas.DataFrame): Matrix to cluster
        annot_types (pandas.DataFrame): DataFrame with the same index as df and a single column 
            with the high-level annotation type for each column in df
        scale (bool): Whether to scale the values in df before clustering. Default is False.

    returns:
        pandas.DataFrame: Clustered matrix
    """
    # If desired, standardize the values in the matrix
    if scale:
        scaler = StandardScaler()
        df_scaled = pd.DataFrame(scaler.fit_transform(df), index=df.index, columns=df.columns)
    else:
        df_scaled = df.copy()	
    # STEP 1: Cluster the high-level annotation types
    # first we want to order the high level annotation types:
    annot_type_means = df_scaled.groupby(annot_types.iloc[:,0], axis=1).mean()
    linkage_matrix_subtypes = linkage(pdist(annot_type_means.T, 'euclidean'), method='average')
    dendrogram_sub_types = dendrogram(linkage_matrix_subtypes, labels=annot_type_means.columns, no_plot=True)
    ordered_sub_types = dendrogram_sub_types['ivl']
    # STEP 2: Cluster the columns within each high-level annotation type
    #Extract the annotation types in the order they appear in input data
    df_annot_types = annot_types.loc[df.columns]
    # create a list to store the order of the diseases
    type_clusters = []
    # iterate over the high level classifications
    for sub_type in ordered_sub_types:
        # Get the samples in the current high level classification
        samps_in_type = df_annot_types[df_annot_types.iloc[:, 0] == sub_type].index
        df_subset = df_scaled[samps_in_type]
        # Calculate the distance matrix for the subset
        dist_matrix_subset = pdist(df_subset.T, 'cosine') # transpose because we want to cluster columns, not rows
        # Perform hierarchical clustering for the subset
        linkage_matrix_subset = linkage(dist_matrix_subset, method='average')
        # Generate dendrogram for subset
        dendrogram_subset = dendrogram(linkage_matrix_subset, labels=samps_in_type, no_plot=True)
        # Store the cluster order
        type_clusters.extend(dendrogram_subset['ivl'])
    # Create a DataFrame with the clustered high-level annotation types
    df_clustered_types = df_scaled[type_clusters]
    # STEP 3: Cluster the rows
    # Calculate the distance matrix for genes
    dist_matrix_genes = pdist(df_clustered_types, 'euclidean')
    # Perform hierarchical clustering for genes
    linkage_matrix_genes = linkage(dist_matrix_genes, method='average')
    # Generate dendrogram for genes
    dendrogram_genes = dendrogram(linkage_matrix_genes, labels=df_clustered_types.index, no_plot=True)
    # Get the gene order
    gene_clusters = dendrogram_genes['ivl']

    # Final DataFrame with clustered diseases and genes
    df_final = df_clustered_types.loc[gene_clusters]
    # Reorder the columns and rows in the original DataFrame according to the clustering
    row_idx = df_final.index
    col_idx = df_final.columns
    return df.loc[row_idx, col_idx]

class EvaluationResults:
    def __init__(self, eval_dir, prefix_file, metrics=["Performance", "Gain"], genesets=['disgen'], dev=False, verbose=True):
        """Class to store the results of evaluating a set of genesets against a set of networks

        Args:
            eval_dir (str): Path to directory containing evaluation results
            prefix_file (str): Path to file containing prefixes of networks to evaluate
            metrics (list, optional): The metric(s) to include from those generated in the 
                evaluation pipeline. Defaults to ["Performance", "Gain"].
            genesets (list, optional): The name(s) of the source of genesets used in the 
                evaluation pipeline. Defaults to ['disgen'].
        """
        # make sure prefixes is in list form
        self.dev = dev
        self.verbose = verbose
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
        self.rankings = {}
        # initialize the network name attribute
        self.network_names = {pref:pref for pref in self.prefixes}
        # load the results
        if self.dev:
            self.results = self.load_new_evaluation_files()
        else:
            for metric in self.metrics:
                self.results[metric] = {}
                for geneset in self.genesets:
                    self.results[metric][geneset] = self.load_evaluation_files(metric, geneset)
    
    def load_new_evaluation_files(self):
        results = defaultdict(dict)
        for geneset in self.genesets:
            datapath = os.path.join(self.eval_dir, geneset)
            geneset_results = {metric:[] for metric in self.metrics}
            for pref in self.prefixes:
                for metric in self.metrics:
                    filename = os.path.join(datapath, pref+"_" + metric + ".tsv")
                    try:
                        df = pd.read_csv(filename, sep="\t")
                        df['Network'] = pref
                        geneset_results[metric].append(df)
                    except FileNotFoundError:
                        print("Not found:", filename)
            for metric in self.metrics:
                results_df = pd.concat(geneset_results[metric])
                for col in results_df.columns[1:-1]:
                    results[metric + '-' + col][geneset] = results_df.pivot(columns="Network", index="Unnamed: 0", values=col)
                    results[metric + '-' + col][geneset].index.name=None
                    results[metric + '-' + col][geneset].columns.name=None
        return results

    def load_evaluation_files(self, use_metric, use_geneset):
        """Function to load the evaluation files for a given metric and geneset. Uses
        the prefixes attribute to find the files. Displays files that fail to load.

        Args:
            use_metric (str): Metric to load
            use_geneset (str): Geneset to load

        Returns:
            pandas.DataFrame: DataFrame with the evaluation results
        """
        # map metrics to the appropriate file suffixes
        paths = {'Performance': ["Performance/", '.performance.csv'], 
            "Gain": ["Performance_Gain/", ".performance_gain.csv"],
            "AUPRC": ["AUPRCs/", ".auprcs.csv"]}
        # get the list of files
        files = self.get_file_list(paths[use_metric][0],"."+use_geneset+ paths[use_metric][1])
        results_list = []
        # iterate over the files and try to load them, print file names if they fail
        failed_count = 0
        for i, f in enumerate(files):
            try:
                df = pd.read_csv(f)
                df["Network"] = self.prefixes[i]
                results_list.append(df)
            except:
                failed_count += 1
                if self.verbose:
                    print("FAILED:", f)
        if failed_count > 0:
            print("Failed to load", failed_count, "files out of", len(files), "for", use_metric, "and", use_geneset)
        # concatenate the results and pivot to get the correct format
        if use_metric == "AUPRC":
            results = pd.concat(results_list).pivot(columns="Network", index="Unnamed: 0", values="AUPRC")
        else:
            results = pd.concat(results_list).pivot(columns="Network", index="Unnamed: 0", values="0")
        results.index.name=None
        results.columns.name=None
        return results

    def get_file_list(self,metric, suff):
        """Function to get a list of files for a given metric and suffix. Note this does not 
        check if the files exist.

        Args:
            metric (str): The metric to use
            suff (str): The suffix to use

        Returns:
            list: A list of file paths
        """
        file_list = []
        for pref in self.prefixes:
            file_list.append(self.eval_dir+metric+pref+suff)
        return file_list
    
    def get_results_metrics(self):
        return [metric for metric in self.results.keys()]

    def rank_all(self, na_option='keep'):
        """Wrapper function to rank all networks for all metrics and genesets
        
        Returns:
            None
        """
        ranks = {}
        # iterate over all availale metrics
        for metric in self.results:	
            # iterate over all available genesets
            for geneset in self.genesets:
                ranks[metric + "-" + geneset] = self.rank_networks(metric, geneset, na_option)
        # store the rankings in a DataFrame
        self.rankings = pd.DataFrame(ranks)

    def rank_networks(self, use_metric, use_geneset, na_option='keep'):
        """Function to rank networks for a given metric and geneset

        Args:
            use_metric (str): The metric to use
            use_geneset (str): The geneset to use

        Returns:
            pandas.DataFrame: A DataFrame with the rankings
        """
        assert na_option in ['keep', 'top', 'bottom', 'drop', 'center']
        # check if the metric and geneset are available
        if use_metric not in self.results.keys():
            print("Metric not available")
            return
        if use_geneset not in self.genesets:
            print("Geneset not available")
            return
        # get the results for the given metric and geneset
        results = self.results[use_metric][use_geneset]
        if na_option == 'drop':
            results = results.dropna(axis=0, thresh=len(results.columns)-3)
            # rank the networks
            rankings = results.rank(axis=1, ascending=False).mean(axis=0)
        elif na_option == 'center':
            rankings = results.rank(axis=1, ascending=False)
            correction =( results.shape[1] - rankings.isna().sum(axis=1))/2 +0.5
            for col in results.rank(axis=1, ascending=False):
                rankings[col] = rankings[col] - correction
            rankings = rankings.mean(axis=0)
        else:
            rankings = results.rank(axis=1, ascending=False, na_option=na_option).mean(axis=0)
        return rankings
    
    def set_network_names(self, name_dict):
        """Function to set the network names

        Args:	
            name_dict (dict): A dictionary with the network names

        Returns:	
            None
        """
        for net in self.network_names:
            self.network_names[net] = name_dict[net]

    def plot_rank_scatter(self, metric, geneset, metric2=None, geneset2=None, axislabel1=None, axislabel2=None, labels=True, savepath=None):
        """Compares the ranked network performance results of all networks for two sources of genesets for a given performance metric,
        or two performance metrics for a given source of genesets. Creates a scatter plot of the results. 

        Args:
            metric (str): The primary metric to use
            geneset (str): The primary geneset source to use
            metric2 (str, optional): The secondary metric for comparison across a single source of genesets. Defaults to None.
            geneset2 (str, optional): The secondary source of genesets for comparison across a single metric. Defaults to None.
            axislabel1 (str, optional): Axis label for primary metric/geneset soucre. Defaults to None (uses column name).
            axislabel2 (str, optional): Axis label for comparison metric/geneset source. Defaults to None (used column name).
            labels (bool, optional): Should the individual points be labeled with network names?. Defaults to True.
            savepath (str, optional): Full filepath for saving the figure to file. Defaults to None.

        Raises:
            ValueError: If both secondary metric and secondary geneset are specified, or if neither are specified.
        
        Returns:
            None
        """
        # check if both secondary metric and secondary geneset are specified
        if metric2 is not None:
            if geneset2 is not None:
                raise ValueError("Can only specify two metrics or two genesets, not two of each")
            else:
                col1 = metric + "-" + geneset
                col2 = metric2 + "-" + geneset
        # check if neither secondary metric nor secondary geneset are specified
        else:
            if geneset is None:
                raise ValueError("Must specify either two metrics or two genesets to enable comparison")
            else:
                col1 = metric + "-" + geneset
                col2 = metric + "-" + geneset2
        # extract the necessary columns from the rankings DataFrame and plot
        _ = plt.figure(figsize=(4,4))
        ax = plt.gca()
        n = len(self.rankings)
        plot_df = self.rankings.loc[:, (col1, col2)]
        plot_df.plot(x=col1, y=col2, kind="scatter", ax=ax, s= 50)
        # add labels if specified
        if labels:
            _ = [plt.text(plot_df[col1][i], plot_df[col2][i], self.network_names[plot_df.index[i]], fontsize=7) for i in range(len(plot_df))]
        # add labels if specified
        if axislabel1 is not None:
            plt.xlabel(axislabel1, fontsize=14)
        if axislabel2 is not None:
            plt.ylabel(axislabel2, fontsize=14)
        # format the plot and save if specified
        plt.ylim((1, n))
        plt.xlim((1, n))
        ax.invert_yaxis()
        ax.invert_xaxis()
        ax.spines[['right', 'top']].set_visible(False)
        if savepath is not None:
            plt.savefig(savepath, dpi=600, bbox_inches="tight")

    def size_adjusted_performances(self, sizes, debug=False, center=None, na_option='keep'):
        """Regress out network size from performance metrics

        Args:
            sizes (dict): Dictionary of network sizes
            debug (bool, optional): If True, randomly assign size variables. Defaults to False.

        Returns:
            None
        """
        if debug:
            sizes = {net:[np.random.randint(1,100)] for net in self.prefixes}
        size_adjusted = {}
        self.sizes = sizes
        rankings = {}
        self.fits = {}
        # for each metric and geneset, regress out network size
        for metric in self.metrics:
            size_adjusted[metric+'_SzAdj'] = {}
            self.fits[metric] = {}
            for geneset in self.genesets:
                results = self.results[metric][geneset]
                #results = results.fillna(0)
                disease_residuals = {}
                fits = {}
                for disease in results.index:
                    y_original = results.loc[disease]
                    valid_networks = y_original.notna()
                    y = y_original[valid_networks]
                    if y.empty:
                        continue
                    x = np.log10(np.array([sizes[net] for net in y.index]))
                    if len(x.shape) == 1:
                        x = x.reshape(-1,1)
                    # perform a linear fit using Huber regression to account for outliers
                    try:
                        reg = HuberRegressor().fit(x, y)
                    except ValueError:
                        print("HuberRegression failed for " + metric + "-" + geneset + "-" + disease + ", using LinearRegression")
                        reg = LinearRegression().fit(x, y)
                    # calculate the residuals as the adjusted performance metric
                    fits[disease] = (reg.coef_, reg.intercept_)
                    disease_residuals[disease] = y_original.copy()
                    if center is None:
                        disease_residuals[disease][valid_networks] = y - reg.predict(x) #+ reg.intercept_
                    elif isinstance(center, int):
                        disease_residuals[disease][valid_networks] = y - reg.predict(x) + center
                    elif center == 'intercept':
                        disease_residuals[disease][valid_networks] = y - reg.predict(x) + reg.intercept_
                    elif center in ['mean', 'median']:
                        if center == 'mean':
                            center_val = y.mean()
                        else:
                            center_val = y.median()
                        disease_residuals[disease][valid_networks] = y - reg.predict(x) + center_val
                    elif center == 'pivot':
                        disease_residuals[disease][valid_networks] = y - reg.predict(x) + reg.predict(np.min(x[valid_networks]).reshape(1, -1))
                # save the residuals and fits
                size_adjusted[metric+'_SzAdj'][geneset] = pd.DataFrame(disease_residuals).T
                self.fits[metric][geneset] = fits
                # re rank the data
                if na_option == 'drop':
                    results = size_adjusted[metric+'_SzAdj'][geneset].dropna(axis=0, thresh=len(results.columns)-3)
                    # rank the networks
                    ranks = results.rank(axis=1, ascending=False).mean(axis=0)
                elif na_option == 'center':
                    ranks = size_adjusted[metric+'_SzAdj'][geneset].rank(axis=1, ascending=False)
                    correction =(size_adjusted[metric+'_SzAdj'][geneset].shape[1] - ranks.isna().sum(axis=1))/2 +0.5
                    for col in size_adjusted[metric+'_SzAdj'][geneset].rank(axis=1, ascending=False):
                        ranks[col] = ranks[col] - correction
                    ranks = ranks.mean(axis=0)
                else:
                    ranks = size_adjusted[metric+'_SzAdj'][geneset].rank(axis=1, ascending=False, na_option=na_option).mean(axis=0)
                rankings[metric + '-SzAdj'+ "-" + geneset ] = ranks
        # replace or add size adjusted metrics to the results
        for result in size_adjusted: 
            self.results[result] = size_adjusted[result]
        for rank in rankings:
            self.rankings[rank] = rankings[rank]

        #self.size_adjusted_rankings = pd.DataFrame(rankings)

    def plot_size_fit(self, metric, geneset, disease):
        """Plot the fit of network size to performance metric for a specified disease geneset, 
        including the linear fit and the residuals

        Args:
            metric (str): The performance metric to plot
            geneset (str): The geneset source to plot
            disease (str): The disease to plot

        Returns:
            None
        """
        metric = metric
        # check that size_adjusted has been run and that the specified metric, geneset, and disease are available
        if metric not in self.results.keys():
            print("Metric not available")
            return
        if geneset not in self.results[metric]:
            print("Geneset not available")
            return
        if disease not in self.results[metric][geneset].index:
            print("Disease not available")
            return
        # get data
        raw_results = self.results[metric][geneset].loc[disease]
        if type(self.sizes[self.prefixes[0]]) == int:
            size_df= pd.DataFrame({"Size":self.sizes})
        elif (self.sizes[self.prefixes[0]]) == 1:
            size_df = pd.DataFrame({"Size":self.sizes})
            size_df['Size'] = size_df["Size"].apply(lambda x: x[0])
        else:
            raise NotImplementedError("Size adjustment plotting only implemented for single size per network")
        
        # get the results for plotting
        raw_results = pd.DataFrame({metric:raw_results}).join(size_df)
        size_results = self.results[metric+'_SzAdj'][geneset].loc[disease]
        size_results = pd.DataFrame({metric+"_SzAdj":size_results}).join(size_df)

        # calculate the fit line from the fit parameters
        self._calculate_fit_line(metric, geneset, disease)
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(14,7))
        sns.scatterplot(data=raw_results, x="Size", y=metric, ax=ax1)
        sns.scatterplot(data=size_results, x="Size", y=metric+"_SzAdj", ax=ax2)
        fit_data = self._calculate_fit_line(metric, geneset, disease)
        sns.lineplot(data = fit_data, x="Size", y="Fit", ax=ax1, color="red")
        ax1.set_xscale("log")
        ax2.set_xscale("log")

    def _calculate_fit_line(self, metric, geneset, disease):
        """Calculate the fit from the fit parameters and results

        Args:
            metric (str): The performance metric to plot
            geneset (str): The geneset source to plot
            disease (str): The disease to plot

        Returns:
            pandas.DataFrame: The fit line data
        """
        # check that size_adjusted has been run and that the specified metric, geneset, and disease are available
        try:
            fit_coeffs = self.fits[metric][geneset][disease]
        except KeyError:
            print("Disease not available")
            return
        # extract the fit parameters
        intercept = fit_coeffs[1]
        coeffs = fit_coeffs[0]
        # initialize the dataframe
        fit_line_df = pd.DataFrame({'Size':self.sizes}) # this won't work for multiple size data points
        if type(self.sizes[self.prefixes[0]]) == list:
            fit_line_df["Size"] = fit_line_df["Size"].apply(lambda x: x[0])
        # calculate the fit line
        x = np.array([self.sizes[net] for net in fit_line_df.index])
        fit_line_df['Fit'] = np.multiply(np.log10(x), coeffs) + intercept
        fit_line_df = fit_line_df.sort_values(by="Size")
        return fit_line_df

    def plot_clustermap(self, metric, geneset, subset=None, network_subset=None, 
            n_subset=100, savepath=None, display_max=20, display_min=-5, 
            cluster_func=cluster_matrix, col_annot=None, row_annot=None, **cluster_kwargs):
        """Plot a clustermap of the performance metric for each disease and network. Allows for
        subsetting of the genesets and networks to be plotted. Three clustering methods are
        available: `cluster_matrix` (hierarchical clustering of columns and rows), 
        `bicluster_matrix` (spectral biclustering for joint clustering of columns and rows), and
        `cluster_by_annot` (`cluster_matrix` while maintaining and ordering column groupings specified).

        Args:
            metric (str): The performance metric to plot
            geneset (str): The geneset source to plot
            subset (str or list, optional): Method for subsetting the genesets to be plotted. 
                One of "top","varaince", "quartile" to subset based on highest mean score, highest
                variance, or largest interquartile range of scores. A list of genesets to directly 
                specify. Defaults to None, all genesets plotted.
            network_subset (list, optional): Subset of networks to be plotted. Defaults to None.
            n_subset (int, optional): If the subset in ['top', 'variance', 'quartile'], 
                specifies the size of the subset. Defaults to 100.
            savepath (str, optional): Complete file path for saving the figure. Defaults to None.
            display_max (int, optional): Cap on the highest performance value to show in colormap. Defaults to 20.
            display_min (int, optional): Floor on the lowest performance value to show in the colormap. Defaults to -5.
            cluster_func (_type_, optional): Function to use for clustering. One of `cluster_matrix`,
                `bicluster_matrix`, or `cluster_by_annot`. Defaults to cluster_matrix.
            col_annot (pandas.DataFrame, optional): Dataframe specifying colors for annotating columns. Defaults to None.
            row_annot (pandas.DataFrame, optional): Dataframe specifying colors for annotating rows. Defaults to None.
            **cluster_kwargs: Additional keyword arguments to pass to the clustering function.
        Returns:
            None	
        """
        # TODO legend for annotations
        result_df = self.results[metric][geneset]
        # perform subsetting
        if network_subset is not None:
            result_df = result_df.loc[:,network_subset]
        if subset is not None:
            result_df = self._subset_results(result_df, subset=subset, n_subset=n_subset)
        result_df.rename(columns=self.network_names, inplace=True)
        # get ordering of rows and columns from specified function
        clustered_results = cluster_func(result_df.T, **cluster_kwargs)
        # if no annotations, plot heatmap directly
        if (col_annot is None) and (row_annot is None):
            _, ax = plt.subplots(figsize=(10,10))
            sns.heatmap(clustered_results, cmap="RdBu", ax=ax, 
                center=0, vmin=display_min, vmax=display_max)
        # if annotations present, plot using clustermap with no clustering
        else:
            cg = sns.clustermap(clustered_results, row_cluster=False, col_cluster=False,
                cmap='RdBu', center=0, vmin=display_min, vmax=display_max,
                row_colors=row_annot, col_colors=col_annot, figsize=(10,15))
            ax = cg.ax_heatmap
        ax.set_xticks([])
        if savepath is not None:	
            plt.savefig(savepath, dpi=300, bbox_inches='tight')

    def plot_boxplot(self, metric, geneset, subset=None, network_subset=None, 
            n_subset=100, max_display=50, savepath=None, size_adjusted_metric=False):
        """Plot a violin plot of the performance metric for each network. Allows for
        subsetting of the genesets and networks to be plotted. 

        Args:	
            metric (str): The performance metric to plot
            geneset (str): The geneset source to plot
            subset (str or list, optional): Method for subsetting the genesets to be plotted.
                One of "top","varaince", "quartile" to subset based on highest mean score, highest
                variance, or largest interquartile range of scores. A list of genesets to directly
                specify. Defaults to None, all genesets plotted.
            network_subset (list, optional): Subset of networks to be plotted. Defaults to None.
            n_subset (int, optional): If the subset in ['top', 'variance', 'quartile'],
                specifies the size of the subset. Defaults to 100.
            max_display (int, optional): Cap on the highest performance value to show in the plot. Defaults to 50.
            savepath (str, optional): Complete file path for saving the figure. Defaults to None.

        Returns:
            None
        """
        if size_adjusted_metric:
            result_df = self.size_adjusted[metric][geneset]
        else:
            result_df = self.results[metric][geneset]
        # perform subsetting
        if network_subset is not None:
            result_df = result_df.loc[:,network_subset]
        if subset is not None:
            result_df = self._subset_results(result_df, subset=subset, n_subset=n_subset)
        # fromat the dataframe for plotting
        plot_df = result_df.melt(value_vars=result_df.columns, var_name="Network", value_name=metric).sort_values(by=metric, ascending=False)
        plot_df["Network"] = plot_df["Network"].apply(lambda x: self.network_names[x])
        # enforce the max display value
        plot_df[metric] = plot_df[metric].apply(lambda x: min(x, max_display))
        # sort the networks by median score
        plot_order = plot_df.groupby("Network").median().sort_values(by=metric, ascending=False).index.to_list()
        # plot the boxplot
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(2,2*(len(plot_order)/10)))
        sns.boxplot(data=plot_df,y="Network", x=metric, linewidth=0.5,  fliersize=0.5, order=plot_order, color="lightgrey", 
            notch=True, ax=ax)
        # format the figure
        plt.ylabel("")
        plt.xlabel("Performance Score")
        ax.tick_params(axis='x', labelsize=9)
        ax.tick_params(axis='y', pad=0,labelsize=9)
        if savepath is not None:
            plt.savefig(savepath, dpi=300, bbox_inches='tight')

    def _subset_results(self, results, subset, n_subset):
        """_summary_

        Args:
            results (panadas.DataFrame): The results dataframe to subset
            subset (str, or list): The method for subsetting the results. One of "top", "variance", 
                "quartile", or a list of genesets to subset.
            n_subset (int): If the subset is "top", "variance", or "quartile", the number 
                of genesets

        Returns:
            plot_results (pandas.DataFrame): The subsetted results dataframe
        """
        if subset == "top": # subset based on mean score across networks
            subset_sets = results.mean(axis=1).sort_values(ascending=False)[:n_subset].index.to_list()
        elif subset == "variance": # subset based on high variance across networks
            subset_sets = results.var(axis=1).sort_values(ascending=False)[:n_subset].index.to_list()
        elif subset == "quartile": # subset based on high interquartile range across networks
            iqrs = results.quantile([0.25, 0.75], axis=1).T
            iqrs['IQR'] = iqrs[0.75] - iqrs[0.25]
            susbet_sets = iqrs.sort_values(by="IQR", ascending=False)[:n_subset].index.to_list()
        elif type(subset) == list: # subset based on a list of genesets
            subset_sets = subset
        plot_results = results.loc[subset_sets]
        return plot_results

