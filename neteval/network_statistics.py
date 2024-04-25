import pandas as pd
import numpy as np
from collections import defaultdict
import subprocess
import matplotlib.pyplot as plt


def load_network_names(filepath):
    with open(filepath, 'r') as f:
        network_names = f.readlines()
    network_names = [n.strip().split(":") for n in network_names]
    network_dict = {net[0]:net[1] for net in network_names }
    return network_dict


class NetworkStats:
    def __init__(self, prefix_file, datadir, exclude_net=[], node_counts=None, edge_counts=None, debug=False):
        self.datadir = datadir
        self.prefix_file = prefix_file
        if debug:
            self.node_files = debug.node_files
            self.edge_files = debug.edge_files
            self.prefixes = debug.prefixes
            self.node_counts = debug.node_counts
            self.network_node_counts = debug.network_node_counts
            self.network_edge_counts = debug.network_edge_counts
        else:
            if isinstance(prefix_file, str):
                with open(self.prefix_file) as f:
                    network_list = f.readlines()
                self.prefixes = [p.strip() for p in network_list]
            elif isinstance(prefix_file, list):
                self.prefixes = prefix_file
            self.node_files = {net:datadir + net + ".nodelist" for net in self.prefixes if net not in exclude_net}
            self.edge_files = {net:datadir + net + "_net.txt" for net in self.prefixes if net not in exclude_net}
            self.node_counts = defaultdict(int)
            self.network_node_counts = defaultdict(int)
            self.network_edge_counts = defaultdict(int)
            if node_counts is None:
                self._count_nodes()
            if edge_counts is None:
                self._count_edges()
        self.nodes = set(self.node_counts.keys())
        self.network_names = {net:net for net in self.prefixes}
        self.densities = defaultdict(float)
        self._calculate_density()
                
    def _count_nodes(self):
        for pref in self.prefixes:
            node_df = pd.read_csv(self.node_files[pref], sep="\t")
            self.network_node_counts[pref] = len(node_df)
            for node in node_df.Unique_Nodes:
                self.node_counts[node] += 1

    def _count_edges(self):
        for pref in self.prefixes:
            command = ['wc', '-l', self.edge_files[pref]]
            result = subprocess.run(command, capture_output=True, text=True)
            output = int(result.stdout.strip().split()[0])
            self.network_edge_counts[pref] = output
            
    def _calculate_density(self):
        for pref in self.prefixes:
            try:
                self.densities[pref] = self.network_edge_counts[pref]/(self.network_node_counts[pref]*(self.network_node_counts[pref]-1))
            except:
                self.densities[pref] = np.nan
            
    def set_network_names(self, names):
        """Function to set the network names

		Args:	
			names (dict, str): A dictionary with the network names, or a string with the path to a file with the network names

		Returns:	
			None
		"""
        if type(names) == str:
            name_dict = load_network_names(names)
        else:
            name_dict = names
        for net in self.network_names:
            self.network_names[net] = name_dict[net]
    
    def plot_size_scatter(self):
        _ = plt.figure(figsize=(5,5))
        pd.DataFrame({'nodes':self.network_node_counts, 'edges':self.network_edge_counts}).plot.scatter(x='nodes', y='edges')
        plt.yscale("log")
        plt.ylabel("Number of Edges")
        plt.xlabel("Number of Nodes")
        
    def plot_density_scatter(self):
        _ = plt.figure(figsize=(5,5))
        pd.DataFrame({'nodes':self.network_node_counts, 'density':self.densities}).plot.scatter(x='nodes', y='density')
        plt.ylabel("Density")
        plt.xlabel("Number of Nodes")
        plt.yscale('log')
        
    def plot_size_bar(self, colors=None):
        fig = plt.figure(figsize=(7, 4), facecolor='white') # Create matplotlib figure
        if colors is None:
            colors = ['red', 'blue']
        ax = fig.add_subplot(111) # Create matplotlib axes
        ax2 = ax.twinx() # Create another axes that shares the same x-axis as ax.

        width = 0.4

        node_df  = pd.DataFrame({'nodes':self.network_node_counts}).sort_values(by="nodes", ascending=False)
        node_df.plot(kind='bar', color=colors[0], ax=ax, width=width, position=1)
        pd.DataFrame({'edges':self.network_edge_counts}).loc[node_df.index].plot(kind='bar', color=colors[1], ax=ax2, width=width, position=0)

        ax.set_ylabel('Nodes', color=colors[0])
        ax2.set_ylabel('Edges', color=colors[1])
        ax2.set_yscale('log')
        ax.set_xlim(-1, len(node_df)-0.5)
        ax.spines['left'].set_color(colors[0])
        ax2.spines['right'].set_color(colors[1])
        ax2.spines['left'].set_visible(False)
        ax.tick_params(axis='y', colors=colors[0])
        ax2.tick_params(axis='y', colors=colors[1]) 
        plt.show()