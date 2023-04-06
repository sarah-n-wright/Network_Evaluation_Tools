import networkx as nx
import obonet as obo
import warnings
import os

def write_obo_file(graph, filename):
    with open(filename, 'w') as obo_file:
        obo_file.write('format-version: 1.2\n')
        obo_file.write('date: 2023-04-03T09:00:00Z\n')
        obo_file.write('saved-by: Your Name\n\n')

        for node in graph.nodes():
            obo_file.write('[Term]\n')
            obo_file.write('id: {}\n'.format(node))
            obo_file.write('name: {}\n'.format(graph.nodes[node]['name']))
            obo_file.write('namespace: {}\n'.format(graph.nodes[node]['namespace']))
            obo_file.write('def: {}\n'.format(graph.nodes[node]['def']))# replace "namespace" with your desired namespace
            attr = [key for key in graph.nodes[node].keys() if key not in ["name", 'namespace', 'def']]
            for key in attr:
                if type(graph.nodes[node][key]) == list:
                    for value in graph.nodes[node][key]:
                        obo_file.write(key +': {}\n'.format(value))
                else:
                    obo_file.write(key+': {}\n'.format(graph.nodes[node][key]))

            obo_file.write('\n')

    print('OBO file written successfully to {}'.format(filename))

# example usage with a toy graph

if __name__=="__main__":
    G = obo.read_obo('http://purl.obolibrary.org/obo/go.obo')
    u='GO:0043235'
    ancestors = nx.ancestors(G, u)
    descendants = nx.descendants(G, u)
    subgraph_nodes = ancestors.union(descendants).union([u])
    subgraph = G.subgraph(subgraph_nodes)
    dirpath = os.path.dirname(os.path.abspath(__file__))
    write_obo_file(subgraph, dirpath + '/sub_go.obo')
    warnings.warn("You must manually transfer the header information from the original .obo to the new .obo", UserWarning)
    
    